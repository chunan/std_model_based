#include <limits>
#include <climits>
#include <cassert>
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>
#include "ugoc_utility.h"
#include "anguso_arg_parser.h"
#include "feature.h"
#include "hmmlite.h"
#include "query_vite_runner.h"
#include "parm.h"
#include "dtw_util.h"
#include "query_hmm.h"

using std::string;

using StdCommonUtil::QueryProfileList;
using StdCommonUtil::SnippetProfile;

using QueryViteUtil::LogBjOt;
using QueryViteUtil::AccumLogBjOt;
using QueryViteUtil::QueryViteManager;
using QueryViteUtil::QDViteBatchParamSet;
using QueryViteUtil::QueryViteRunner;

using DtwUtil::Path;
using DtwUtil::Parmtype;
using DtwUtil::DtwManager;
using DtwUtil::DtwBatchParamSet;
using DtwUtil::QDDtwBatchParamSet;
using DtwUtil::DtwRunner;
using DtwUtil::FrameDtwRunner;


#define PARTIAL_RESULTS


#define float_inf std::numeric_limits<float>::infinity()


static void usage(const char *program_name) {/*{{{*/
  fprintf(stderr, "Usage: %s [args]\n", program_name);
  fprintf(stderr, " Required args            Description            Default\n");
  fprintf(stderr, "  -load-param  <string>   load from parameter    \n");
  fprintf(stderr, "  -doc-list    <string>   list of doc state file \n");
  fprintf(stderr, "  -query-list  <string>   list of query feature  \n");
  fprintf(stderr, "  -result      <string>   file to store result   \n");
  fprintf(stderr, " Optional args\n");
  fprintf(stderr, "  -doc-dir     <string>   directory of doc       none\n");
  fprintf(stderr, "  -query-dir   <string>   directory of query     none\n");
  fprintf(stderr, "  -nthread     <unsigned> number of thread       1\n");
  fprintf(stderr, "  -ignore      <string>   list of q-d to ignore  none\n");
  fprintf(stderr, "  -debug-ans   <string>   ans indiater @result   none\n");
  fprintf(stderr, "  -dump-path                                     false\n");
}/*}}}*/

void ParseArg(const int argc, const char **argv,/*{{{*/
               C_AngusoArgParser* arg_parser) {
  /* Argument parsing */
  arg_parser->addArgument("load-param", C_AngusoArgParser::A_STRING, true);
  arg_parser->addArgument("doc-list", C_AngusoArgParser::A_STRING, true);
  arg_parser->addArgument("query-list", C_AngusoArgParser::A_STRING, true);
  arg_parser->addArgument("result", C_AngusoArgParser::A_STRING, true);

  arg_parser->addArgument("doc-dir", C_AngusoArgParser::A_STRING, false);
  arg_parser->addArgument("query-dir", C_AngusoArgParser::A_STRING, false);
  arg_parser->addArgument("nthread", C_AngusoArgParser::A_UNSIGNED, false);
  arg_parser->addArgument("ignore", C_AngusoArgParser::A_STRING, false);
  arg_parser->addArgument("debug-ans", C_AngusoArgParser::A_STRING, false);
  arg_parser->addArgument("dump-path", C_AngusoArgParser::A_SWITCH, false);

  arg_parser->setDefaultValue("doc-dir", "");
  arg_parser->setDefaultValue("query-dir", "");
  arg_parser->setDefaultValue("nthread", 1u);
  arg_parser->setDefaultValue("ignore", "");
  arg_parser->setDefaultValue("debug-ans", "");
  arg_parser->setDefaultValue("dump-path", false);

  arg_parser->processArgument(argc, argv);
  if ( arg_parser->fail() ) {
    cerr << arg_parser->getErrorMsg() << endl;
    usage(argv[0]);
    exit(EXIT_FAILURE);
  }

}/*}}}*/



void CopyFeatPtr(const vector<Parm>& parm_vec,/*{{{*/
                 vector<const DenseFeature*>* targ_feat_vec) {
  targ_feat_vec->resize(parm_vec.size());
  for (unsigned i = 0; i < parm_vec.size(); ++i)
    (*targ_feat_vec)[i] = &parm_vec[i].Feat();
}/*}}}*/



struct Printer {/*{{{*/
  void operator()(string& str) {
    cout << str << endl;
  }
};/*}}}*/



class QSDtwBatchParamSet : public DtwBatchParamSet {/*{{{*/
  public:
    QSDtwBatchParamSet(Dispatcher<UPair>* dispatcher,
                       vector<Parm>* query_parm_list,
                       vector<Parm>* doc_parm_list,
                       vector<SnippetProfileList>* candidates,
                       vector<IPair>* vqboundary,
                       vector<vector<vector<float> > >* snippet_dist,
                       vector<vector<vector<Path> > >* paths) {
      dispatcher_ = dispatcher;
      query_parm_list_ = query_parm_list;
      doc_parm_list_ = doc_parm_list;
      candidates_ = candidates;
      snippet_dist_ = snippet_dist;
      paths_ = paths;
      vqboundary_ = vqboundary;

      /* Resize output variables */
      snippet_dist_->resize(candidates_->size());
      if (paths_) paths_->resize(candidates_->size());

      for (unsigned qidx = 0; qidx < candidates_->size(); ++qidx) {
        (*snippet_dist_)[qidx].resize((*candidates_)[qidx].size());
        if (paths_) (*paths_)[qidx].resize((*candidates_)[qidx].size());
      }

    }
    virtual bool InstallDtwRunner(DtwRunner* runner);
  private:
    /* data */
    Dispatcher<UPair>* dispatcher_;
    vector<Parm>* query_parm_list_;
    vector<Parm>* doc_parm_list_;
    vector<SnippetProfileList>* candidates_;
    vector<IPair>* vqboundary_;
    vector<vector<vector<float> > >* snippet_dist_;
    vector<vector<vector<Path> > >* paths_;
};/*}}}*/

bool QSDtwBatchParamSet::InstallDtwRunner(DtwRunner* runner) {/*{{{*/

  UPair* ticket = dispatcher_->GetObjPtr();
  if (ticket) {
    int qidx = ticket->first;
    int sidx = ticket->second;
    const SnippetProfile& pf = (*candidates_)[qidx].GetProfile(sidx);
    int didx = pf.Didx();
    Parm* q_parm = &(*query_parm_list_)[qidx];
    Parm* d_parm = &(*doc_parm_list_)[didx];
    const IPair* qboundary = vqboundary_ ? &(*vqboundary_)[qidx] : NULL;
    const IPair* dboundary = &pf.Boundary();
    vector<float>* snippet_dist = &(*snippet_dist_)[qidx][sidx];
    vector<Path>* paths = paths_ ? &(*paths_)[qidx][sidx] : NULL;

    runner->InitDtw(q_parm, d_parm, snippet_dist, paths, qboundary, dboundary);
  }

  return ticket != NULL;

}/*}}}*/



struct UTriple {/*{{{*/
  UTriple (unsigned q, unsigned s, unsigned p) : qidx(q), sidx(s), pidx(p) {}
  unsigned qidx;
  unsigned sidx;
  unsigned pidx;
};/*}}}*/

class PrfDtwBatchParamSet : public DtwBatchParamSet {/*{{{*/
  public:
    PrfDtwBatchParamSet(Dispatcher<UTriple>* dispatcher,/*{{{*/
                        const vector<Parm>* doc_parm_list,
                        const vector<SnippetProfileList>* candidates,
                        const vector<unsigned>* num_prf) {
      dispatcher_ = dispatcher;
      doc_parm_list_ = doc_parm_list;
      candidates_ = candidates;
      num_prf_ = num_prf;

      /* Resize output variables */
      output_dist_.resize(candidates_->size()); // Number of queries

      for (unsigned qidx = 0; qidx < candidates_->size(); ++qidx) {
        // Number of hypothesized regions
        output_dist_[qidx].resize((*candidates_)[qidx].size());

        for (unsigned sidx = 0; sidx < output_dist_[qidx].size(); ++sidx) {
          // Number of pseudo relevant regions
          output_dist_[qidx][sidx].resize((*num_prf_)[qidx]);
        }
      }

    }/*}}}*/
    virtual bool InstallDtwRunner(DtwRunner* runner);
    void IntegratePrfDist(const vector<float>& prf_weight,
                          vector<SnippetProfileList>* old_snippet_lists);
  private:
    /* data */
    Dispatcher<UTriple>* dispatcher_;
    const vector<Parm>* doc_parm_list_;
    const vector<SnippetProfileList>* candidates_;
    const vector<unsigned>* num_prf_;
    /* local */
    vector<vector<vector<vector<float> > > > output_dist_;
};/*}}}*/

bool PrfDtwBatchParamSet::InstallDtwRunner(DtwRunner* runner) {/*{{{*/

  UTriple* ticket = dispatcher_->GetObjPtr();
  if (ticket) {
    const int qidx = ticket->qidx;
    const int sidx = ticket->sidx;
    const int pidx = ticket->pidx;
    const SnippetProfile& pseu_region = (*candidates_)[qidx].GetProfile(pidx);
    const SnippetProfile& hypo_region = (*candidates_)[qidx].GetProfile(sidx);
    assert(pidx != sidx);
    const Parm* q_parm = &(*doc_parm_list_)[pseu_region.Didx()];
    const Parm* d_parm = &(*doc_parm_list_)[hypo_region.Didx()];
    const IPair* qboundary = &pseu_region.Boundary();
    const IPair* dboundary = &hypo_region.Boundary();
    vector<float>* snippet_dist = &output_dist_[qidx][sidx][pidx];
    vector<Path>* paths = NULL;

    runner->InitDtw(q_parm, d_parm, snippet_dist, paths, qboundary, dboundary);
  }

  return ticket != NULL;
}/*}}}*/

void PrfDtwBatchParamSet::IntegratePrfDist(/*{{{*/
    const vector<float>& prf_weight,
    vector<SnippetProfileList>* old_snippet_lists) {

  for (unsigned qidx = 0; qidx < old_snippet_lists->size(); ++qidx) {

    int n_not_found = 0;

    /* push distance zero if sidx < (*num_prf_)[qidx] */
    for (unsigned pidx = 0; pidx < (*num_prf_)[qidx]; ++pidx)
      output_dist_[qidx][pidx][pidx].push_back(0.0);

    /* for each snippet in qidx list */
    for (unsigned sidx = 0; sidx < (*old_snippet_lists)[qidx].size(); ++sidx) {

      SnippetProfile& snippet = (*old_snippet_lists)[qidx].ProfileRef(sidx);
      vector<vector<float> >& dist_vec = output_dist_[qidx][sidx];

      if (dist_vec[0].size() < (*num_prf_)[qidx]) {
        /* Colloct them in a single vector, namely, output_dist_[qidx][sidx][0] */
        if (dist_vec[0].empty()) dist_vec[0].push_back(-float_inf);
        for (unsigned pidx = 1; pidx < (*num_prf_)[qidx]; ++pidx) {
          if (dist_vec[pidx].empty()) // Cannot find path with PRF(pidx)
            dist_vec[0].push_back(-float_inf);
          else
            dist_vec[0].push_back(dist_vec[pidx][0]);
        }
      }

      // Calculate innerproduct
      snippet.ScoreRef() += inner_product(
          dist_vec[0].begin(), dist_vec[0].end(), prf_weight.begin(), 0.0f);
      if (snippet.ScoreRef() == -float_inf)
        n_not_found++;

    } // end for sidx

    (*old_snippet_lists)[qidx].Sort();
    if (n_not_found > 0) {
      cerr << "Warning: PrfDtwBatchParamSet::IntegratePrfDist():"
        << "qidx = " << qidx << "#not found = " << n_not_found << endl;
      (*old_snippet_lists)[qidx].Resize(-n_not_found);
    }

  } // end for qidx

}/*}}}*/



void GeneratePrfWeight(const float prf_lambda, vector<float>* prf_weight) {/*{{{*/
  float sum = 0.0;
  for (unsigned pidx = 0; pidx < prf_weight->size(); ++pidx) {
    (*prf_weight)[pidx] = pow(0.5, prf_lambda * (pidx + 1));
    sum += (*prf_weight)[pidx];
  }

  /*
  float normalizer = 3.0 / sum;
  cout << "PRF weight = {";
  for (unsigned pidx = 0; pidx < prf_weight->size(); ++pidx) {
    (*prf_weight)[pidx] *= normalizer;
    cout << " " << (*prf_weight)[pidx];
  }
  cout << "}\n";
  */
}/*}}}*/


void InitPrfDispatcher(Dispatcher<UTriple>* prf_dispatcher, /*{{{*/
                       const vector<SnippetProfileList>& candidates,
                       const vector<unsigned>& num_prf) {

  prf_dispatcher->Clear();

  for (unsigned qidx = 0; qidx < candidates.size(); ++qidx) {
    for (unsigned sidx = 0; sidx < candidates[qidx].size(); ++sidx) {
      for (unsigned pidx = 0; pidx < num_prf[qidx]; ++pidx) {
        if (sidx == pidx) continue;
        prf_dispatcher->Push(UTriple(qidx, sidx, pidx));
      }
    }
  }

}/*}}}*/

void InitTrainerDispatcher(Dispatcher<pair<unsigned, int8_t> >* trainer_disp, /*{{{*/
                           const vector<SnippetProfileList>& candidates) {
  for (unsigned qidx = 0; qidx < candidates.size(); ++qidx) {
    trainer_disp->Push(make_pair<unsigned, int8_t>(qidx, 0));
    trainer_disp->Push(make_pair<unsigned, int8_t>(qidx, 1));
  }
}/*}}}*/

/* score[i] >= t
 * score[i + 1] < t
 */
unsigned ScoreSearch(const SnippetProfileList& snippet_list, float t,
                     int l = 0, int h = INT_MAX) {
  int low = max(l, 0);
  int high = min(h, static_cast<int>(snippet_list.size()) - 1);
  assert(low <= high);

  while (low < high) {
    int mid = (low + high + 1) / 2;
    if (snippet_list.GetProfile(mid).Score() >= t) {
      low = mid;
    } else {
      high = mid - 1;
    }
  }

  return low;
}

void SnippetListStat(const SnippetProfileList& snippet_list,
                     float* p_mean, float* p_std,
                     unsigned begin = 0, unsigned end = -1) {

  float& mean = *p_mean;
  float& std = *p_std;

  if (end < 0 || end > snippet_list.size()) end = snippet_list.size();
  if (begin < 0 || begin >= snippet_list.size()) begin = 0;
  int n = end - begin;

  if (n == 1) {
    mean = snippet_list.GetProfile(begin).Score();
    std = 0.0f;
  } else {
    mean = 0.0f;
    std = 0.0f;
    for (unsigned s = begin; s < end; ++s) {
      float score = snippet_list.GetProfile(s).Score();
      mean += score;
      std += score * score;
    }
    mean /= (end - begin);
    std = sqrt(std / n - mean * mean);
  }
}



void Std_viterbi_dtw(C_AngusoArgParser& arg_parser) {/*{{{*/

  Timer timer;
  unsigned tid;
  InstallTimer(DtwUtil::timer, timer);
  InstallTimer(QueryViteUtil::timer, timer);

  /* Get arguments from arg_parser */
  string s_doclist, s_querylist, s_loadparam, s_resultfile, s_ignorelist;
  string s_docdir, s_querydir;
  string s_answerfile;
  unsigned nthread;
  bool doDumpPath;
  arg_parser.getStringArgument("load-param", &s_loadparam);
  arg_parser.getStringArgument("doc-list", &s_doclist);
  arg_parser.getStringArgument("query-list", &s_querylist);
  arg_parser.getStringArgument("result", &s_resultfile);
  arg_parser.getStringArgument("doc-dir", &s_docdir);
  arg_parser.getStringArgument("query-dir", &s_querydir);
  arg_parser.getUnsignedArgument("nthread", &nthread);
  arg_parser.getStringArgument("ignore", &s_ignorelist);
  arg_parser.getStringArgument("debug-ans", &s_answerfile);
  arg_parser.getSwitchArgument("dump-path", &doDumpPath);



  /* Variables */
  vector<string> doc_list;             // Document filenames
  vector<string> query_list;           // Query filenames
  QueryProfileList q_profile_list;     // Query profile list
  AnswerList ans_list;

  HMM_GMM bg_model;                    // HMM
  vector<Gaussian*> gauss_pool;        // Gaussian Pool
  vector<GaussianMixture*> state_pool; // State Pool

  Dispatcher<UPair> dispatcher;

  vector<Parm> query_parm;             // Query parm
  vector<const DenseFeature*> pquery_feat;
  vector<Parm> doc_parm;               // Query parm
  vector<const DenseFeature*> pdoc_feat;
  vector<Labfile> doc_lab;
  QDArray<vector<float> > snippet_like;
  QDArray<vector<Labfile> > query_to_doc_lab;
#ifdef UNCONSTRAINED
  vector<vector<float> > next_trans;
  vector<vector<float> > self_trans;
  vector<LogBjOt> query_logBjOt;
#else
  vector<float> next_trans;
  vector<float> self_trans;
  vector<AccumLogBjOt> query_logBjOt;
#endif
  float trans_weight = 5.0;



  tid = timer.Tic("Parsing Lists");
  StdCommonUtil::ParseList(s_doclist.c_str(), s_docdir, &doc_list, NULL);
  StdCommonUtil::ParseList(s_querylist.c_str(), s_querydir, &query_list,
                           &q_profile_list);
  if (!s_ignorelist.empty())
    ParseIgnore(s_ignorelist.c_str(), doc_list, &q_profile_list);
  if (!s_answerfile.empty())
    ans_list.Init(s_answerfile, q_profile_list, doc_list);
  timer.Toc(tid);

  //for_each(doc_list.begin(), doc_list.end(), Printer());

  tid = timer.Tic("Loading HMMs");
  LoadHMMGMG(s_loadparam, &bg_model, state_pool, gauss_pool);
  bg_model.SyncUsed();
  for(unsigned g = 0; g < gauss_pool.size(); g++)
    if (bg_model.getGisUsed(g)) gauss_pool[g]->InvertCov();
  timer.Toc(tid);


  tid = timer.Tic("Load Parms");
  DtwUtil::LoadParmList(query_list, s_querydir, &query_parm, DtwUtil::FRAME);
  DtwUtil::LoadParmList(doc_list, s_docdir, &doc_parm, DtwUtil::FRAME);
  CopyFeatPtr(query_parm, &pquery_feat);
  CopyFeatPtr(doc_parm, &pdoc_feat);
  timer.Toc(tid);


  tid = timer.Tic("Load doc labels");
  vector<string> doc_rec_list = doc_list;
  for_each(doc_rec_list.begin(), doc_rec_list.end(), ReplaceExt("rec"));
  //for_each(doc_rec_list.begin(), doc_rec_list.end(), Printer());
  QueryViteUtil::LoadDocLab(doc_rec_list, &doc_lab);
  timer.Toc(tid);

  /* Resname doc list */
  KeepBasename(&doc_list);



  /***************************** Viterbi **********************************/

  /* Init QueryViteRunner's Dispatcher */
  InitDispatcher(&dispatcher, q_profile_list, doc_list);
  dispatcher.Verbose();
  dispatcher.SetVerboseInt(dispatcher.size() / 50);


  /* Create ParamSet  */
  QDViteBatchParamSet vite_batch_param(
      &dispatcher,
      &bg_model,
      &pquery_feat,
      &doc_lab,
      &snippet_like,
      &query_to_doc_lab,
      NULL, // vqboundary
      NULL, // vdboundary
      &next_trans,
      &self_trans,
      &query_logBjOt,
      trans_weight);

  vite_batch_param.Init();


  /* Setup QueryViteRunner parameters */
  QueryViteRunner::nsnippet_ = 3;
  QueryViteRunner::ratio_bound_ = 2.0;


  /* Create `nthread' of QueryViteManager */
  vector<QueryViteManager> query_vite_man(nthread);
  vector<QueryViteRunner> query_vite_runner(nthread);
  for (unsigned t = 0; t < nthread; ++t) {
    query_vite_man[t].Install(&vite_batch_param, &query_vite_runner[t]);
  }


  /* Start Viterbi decoding */
  cout << "Viterbi:\n";
  tid = timer.Tic("Vite all");
  CastThreads(query_vite_man);
  timer.Toc(tid);
  cout << "\n";


  /* Push into vite_snippet_lists */
  vector<SnippetProfileList> vite_snippet_lists(query_list.size());
  vector<SnippetProfileList> vite_all_snp_lists(query_list.size());
  for (unsigned qidx = 0; qidx < vite_snippet_lists.size(); ++qidx) {
    /* In each query qidx */
    for (unsigned didx = 0; didx < doc_list.size(); ++didx) {
      /* In each captured snippet */
      for (unsigned s = 0; s < snippet_like(qidx, didx).size(); ++s) {
        Labfile& lab = query_to_doc_lab(qidx, didx)[s];
        IPair bound = IPair(doc_lab[didx].getStartF(lab.getCluster(0)),
                            doc_lab[didx].getEndF(lab.getCluster(-1)));
        vite_snippet_lists[qidx].push_back(
            qidx, didx, s, snippet_like(qidx, didx)[s], bound);
      }
    }
    const unsigned FDtwNRescore = min(500u, vite_snippet_lists[qidx].size());
    //cout << "vite_snippet_lists[" << qidx << "] = \n" << vite_snippet_lists[qidx];
    vite_snippet_lists[qidx].Sort();
    vite_all_snp_lists[qidx] = vite_snippet_lists[qidx];

    if (FDtwNRescore < vite_snippet_lists[qidx].size()) {
      cout << "Q" << qidx << ": Reduce search snippets from "
        << vite_snippet_lists[qidx].size() << " to " << FDtwNRescore
        << " for fDTW " << endl;
    }
    vite_snippet_lists[qidx].Resize(FDtwNRescore);

  }


  vector<const vector<SnippetProfileList>* > v_snp_lists;

#ifdef PARTIAL_RESULTS
  v_snp_lists.push_back(&vite_all_snp_lists);
  DumpResult(s_resultfile + ".vite",
             q_profile_list,
             v_snp_lists,
             doc_list,
             &ans_list);
#endif

  /**************************** Frame-based DTW *******************************/

  /* Output Variables */
  vector<vector<vector<float> > > snippet_dist;
  vector<vector<vector<Path> > > frame_paths;


  /* Init dispatcher */
  InitDispatcher(&dispatcher, vite_snippet_lists);
  dispatcher.Verbose();
  dispatcher.SetVerboseInt(dispatcher.size() / 50);


  /* Install Batch Parameter Set */
  QSDtwBatchParamSet dtw_batch_param(
      &dispatcher,
      &query_parm,
      &doc_parm,
      &vite_snippet_lists,
      NULL, // vqboundary_
      &snippet_dist,
      &frame_paths);

  /* Setup FrameDtwRunner parameter */
  FrameDtwRunner::nsnippet_ = 1;



  /* Multi-thread */
  vector<DtwManager> dtw_man(nthread);
  vector<FrameDtwRunner> frame_runner(nthread);
  for (unsigned t = 0; t < nthread; ++t) {
    dtw_man[t].Install(&dtw_batch_param, &frame_runner[t]);
  }

  cout << "Frame-based DTW:\n";
  tid = timer.Tic("FrameDTW");
  CastThreads(dtw_man);
  timer.Toc(tid);
  cout << "\n";


  vector<SnippetProfileList> dtw_snippet_lists(vite_snippet_lists.size());
  vector<SnippetProfileList> dtw_all_snp_lists(vite_snippet_lists.size());

  /* Change score and boundary according to frame-based DTW */
  for (unsigned qidx = 0; qidx < vite_snippet_lists.size(); ++qidx) {
    /* In each query qidx */
    for (unsigned sidx = 0; sidx < vite_snippet_lists[qidx].size(); ++sidx) {
      if (!snippet_dist[qidx][sidx].empty()) {
        assert(!frame_paths[qidx][sidx].empty());
        const SnippetProfile& vite_snippet = vite_snippet_lists[qidx].GetProfile(sidx);
        dtw_snippet_lists[qidx].push_back(
            vite_snippet.Qidx(),
            vite_snippet.Didx(),
            vite_snippet.NthSnippet(),
            snippet_dist[qidx][sidx][0],
            IPair(frame_paths[qidx][sidx][0].front().second + 1,
                  frame_paths[qidx][sidx][0].back().second));
      }
    }

    const unsigned PrfNRescore = min(500u, dtw_snippet_lists[qidx].size());

    dtw_snippet_lists[qidx].Sort();
    dtw_all_snp_lists[qidx] = dtw_snippet_lists[qidx];

    if (PrfNRescore < dtw_snippet_lists[qidx].size()) {
      cout << "Q" << qidx << ": Reduce search snippets from "
        << dtw_snippet_lists[qidx].size() << " to " << PrfNRescore
        << " for PRF" << endl;
    }
    dtw_snippet_lists[qidx].Resize(PrfNRescore);
  }

#ifdef PARTIAL_RESULTS
  v_snp_lists.push_back(&dtw_all_snp_lists);
  DumpResult(s_resultfile + ".fdtw",
             q_profile_list,
             v_snp_lists,
             doc_list,
             &ans_list);
#endif



  /********************************* PRF **************************************/

  /* Input variables */
  const int fix_num_prf = 7;
  const float prf_lambda = 0.4;
  vector<unsigned> num_prf(doc_list.size(), fix_num_prf);
  vector<float> prf_weight(fix_num_prf);
  GeneratePrfWeight(prf_lambda, &prf_weight);
  Dispatcher<UTriple> prf_dispatcher;


  cout << "PRF weight = " << prf_weight << endl;

  /* Output variables */
  vector<vector<float> > prf_dist;

  /* Init dispatcher */
  InitPrfDispatcher(&prf_dispatcher, dtw_snippet_lists, num_prf);
  prf_dispatcher.Verbose();
  prf_dispatcher.SetVerboseInt(prf_dispatcher.size() / 50);

  /* Install Batch Parameter Set */
  PrfDtwBatchParamSet prf_batch_param(
      &prf_dispatcher,
      &doc_parm,
      &dtw_snippet_lists,
      &num_prf);

  /* Setup FrameDtwRunner parameter */
  FrameDtwRunner::nsnippet_ = 1;

  /* Multi-thread */
  for (unsigned t = 0; t < nthread; ++t) {
    dtw_man[t].Install(&prf_batch_param, &frame_runner[t]);
  }

  cout << "PRF Frame-based DTW:\n";
  tid = timer.Tic("PRF FrameDTW");
  CastThreads(dtw_man);
  timer.Toc(tid);
  cout << "\n";

  PrintCalCount(frame_runner);

  /* PRF snippets */
  vector<SnippetProfileList> prf_snippet_lists = dtw_snippet_lists;

  /* New Scores from PRF */
  prf_batch_param.IntegratePrfDist(prf_weight, &prf_snippet_lists);

#ifdef PARTIAL_RESULTS
  v_snp_lists.push_back(&prf_snippet_lists);
  DumpResult(s_resultfile + ".prf",
             q_profile_list,
             v_snp_lists,
             doc_list,
             &ans_list);
#endif

  /********************************* Query HMMs *******************************/
  /* Input data */
  Dispatcher<pair<unsigned, int8_t> > train_disp;
  InitDispatcher(&dispatcher, prf_snippet_lists);
  InitTrainerDispatcher(&train_disp, prf_snippet_lists);
  train_disp.Verbose();
  train_disp.SetVerboseInt(train_disp.size() / 50);

  /* New statepool and gausspool */
  vector<GaussianMixture*> a_state_pool; // State Pool
  vector<Gaussian*> a_gauss_pool;        // Gaussian Pool
  bg_model.setpStatePool(&a_state_pool, UNUSE);
  bg_model.setpGaussPool(&a_gauss_pool, UNUSE);

  /* Calculate num_prf */
  vector<UPair> irr_range(prf_snippet_lists.size());
  for (unsigned qidx = 0; qidx < prf_snippet_lists.size(); ++qidx) {
    float mean, std;
    SnippetListStat(prf_snippet_lists[qidx], &mean, &std);
    num_prf[qidx] = ScoreSearch(prf_snippet_lists[qidx], mean + 1.9 * std);
    irr_range[qidx].first =
      ScoreSearch(prf_snippet_lists[qidx], mean + 0.25 * std, num_prf[qidx]);
    irr_range[qidx].second =
      ScoreSearch(prf_snippet_lists[qidx], mean - 0.25 * std, irr_range[qidx].first);
    cout << qidx << setprecision(2) << ": (" << mean << ", " << std << "), "
      << num_prf[qidx] << ", "
      << irr_range[qidx].first << " - " << irr_range[qidx].second
      << "(" << irr_range[qidx].second - irr_range[qidx].first + 1
      << ")" << endl;
  }

  /* Calculate train_weight */
  const float hmm_lambda = 0.25;
  vector<float> hmm_weight(50);
  GeneratePrfWeight(hmm_lambda, &hmm_weight);
  cout << "HMM weight = " << hmm_weight << endl;
  vector<vector<float> > train_weight(prf_snippet_lists.size(), hmm_weight);
  vector<vector<float> > targ_loglike;
  vector<vector<float> > anti_loglike;

  QueryHmmParamSet query_hmm_set(
      &train_disp,
      &pquery_feat,
      &pdoc_feat,
      &prf_snippet_lists,
      &num_prf,
      &train_weight,
      &irr_range,
      &doc_lab,
      &bg_model,
      20, // max_iter
      &targ_loglike,
      &anti_loglike);

  vector<QueryHmmTrainer> trainer(nthread);
  vector<QueryHmmDecoder> decoder(nthread);
  vector<QueryHmmManager> query_hmm_man(nthread);
  for (unsigned t = 0; t < nthread; ++t) {
    query_hmm_man[t].Install(&query_hmm_set, &trainer[t], &decoder[t]);
  }

  cout << "Query HMM:\n";
  tid = timer.Tic("QueryHmm");
  CastThreads(query_hmm_man);
  timer.Toc(tid);
  cout << "\n";

  vector<SnippetProfileList> hmm_snippet_lists(prf_snippet_lists);
  for (unsigned qidx = 0; qidx < hmm_snippet_lists.size(); ++qidx) {

    int n_not_found = 0;
    for (unsigned sidx = 0; sidx < hmm_snippet_lists[qidx].size(); ++sidx) {
      hmm_snippet_lists[qidx].ProfileRef(sidx).ScoreRef() =
        targ_loglike[qidx][sidx] - anti_loglike[qidx][sidx];
    }

    hmm_snippet_lists[qidx].Sort();
    if (n_not_found > 0) {
      cerr << qidx << ": not found = " << n_not_found << endl;
      hmm_snippet_lists[qidx].Resize(-n_not_found);
    }

  }

  v_snp_lists.push_back(&hmm_snippet_lists);
  DumpResult(s_resultfile + ".hmm",
             q_profile_list,
             v_snp_lists,
             doc_list,
             &ans_list);
  /********************************** End *************************************/
  /* Print running time information */
  timer.Print();

}/*}}}*/



int main(const int argc, const char **argv) {
  C_AngusoArgParser arg_parser;

  ParseArg(argc, argv, &arg_parser);

  cout << "======== Ready to run STD ==========\n";
  arg_parser.dumpAllArguments(cout);
  cout << "====================================\n";

  Std_viterbi_dtw(arg_parser);

  return 0;
}

