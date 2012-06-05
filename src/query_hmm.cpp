#include "query_hmm.h"

using StdCommonUtil::SnippetProfile;



void HmmPair::Init(const Labfile& lab, const int start, const int end, /*{{{*/
                   HMM_GMM& bg_model, const float del_ratio) {

  // Install target model
  targ_model.setUse(bg_model.getUse());
  targ_model.setpStatePool(bg_model.getpStatePool(USE), USE);
  targ_model.setpStatePool(bg_model.getpStatePool(UNUSE), UNUSE);
  targ_model.setpGaussPool(bg_model.getpGaussPool(USE), USE);
  targ_model.setpGaussPool(bg_model.getpGaussPool(UNUSE), UNUSE);

  // Install anti model
  anti_model.setUse(bg_model.getUse());
  anti_model.setpStatePool(bg_model.getpStatePool(USE), USE);
  anti_model.setpStatePool(bg_model.getpStatePool(UNUSE), UNUSE);
  anti_model.setpGaussPool(bg_model.getpGaussPool(USE), USE);
  anti_model.setpGaussPool(bg_model.getpGaussPool(UNUSE), UNUSE);



  /***************** Setup targ_model *****************/
  targ_model.setNstate(end - start);
  targ_model.ClearTrans();
  targ_model.ClearPi(USE);

  for (int i = start; i < end; ++i) {
    // state number in background model
    int old_sno = lab.getCluster(i);
    // state index in background model
    int old_sid = bg_model.getGMidx(old_sno);
    // state number in target model
    int new_sno = i - start;

    /* Create a state (and Gaussians, handled by NewStateCopy()). */
    int new_sid = NewStateCopy(old_sid, bg_model, NULL, NULL);
    targ_model.setState(new_sno , new_sid);

    /* Set transitions, left, right */
    float self_trans = bg_model.getTrans(old_sno, old_sno);

    if (i != start) { // Not first state
      targ_model.setLeft(new_sno, new_sno - 1);
    } else { // First state
      targ_model.setLeft(new_sno, -1);
    }

    if (i != end - 1) { // Not last state
      targ_model.setRight(new_sno, new_sno + 1);
      targ_model.setTrans(new_sno, new_sno, self_trans);
      targ_model.setTrans(new_sno, new_sno + 1, 1.0 - self_trans);
    } else { // Last state
      targ_model.setRight(new_sno, -1);
      targ_model.setTrans(new_sno, new_sno, 1.0);
    }

  }

  int tot_frames = lab.getDuration(start, end - 1);
  targ_model.setPi(0, 1.0);
  /* Set starting states: Pi != 0 */
  for (int i = start + 1; i < end; ++i) {
    if (lab.getDuration(start, i - 1) < del_ratio * tot_frames)
      targ_model.setPi(i - start, 1.0);
  }
  targ_model.normPi(USE);

  /* Set ending states: Right = -1 */
  for (int i = end -2; i > 0; --i) {
    if (lab.getDuration(i + 1, end - 1) < del_ratio * tot_frames)
      targ_model.setRight(i - start, -1);
  }



  /***************** Setup anti_model *****************/
  anti_model.CopyForThread(targ_model);
  for (int sno = 0; sno < targ_model.getNstate(); ++sno) {
    int targ_sid = targ_model.getGMidx(sno);
    int anti_sid = NewStateCopy(targ_sid, targ_model, NULL, NULL);
    anti_model.setState(sno, anti_sid);
  }


  Gaussian::setVarFloor(0.32);
  Gaussian::setRequiredFrame(15.0);
  GaussianMixture::setRequiredFrame(15.0);
  /*
  cout << "targ_model:\n";
  targ_model.display(stdout);
  cout << "anti_model:\n";
  anti_model.display(stdout);
  */

}/*}}}*/



QueryHmmParamSet::QueryHmmParamSet(/*{{{*/
    Dispatcher<pair<unsigned, int8_t> >* dispatcher,
    const vector<const DenseFeature*>* query_feat_list,
    const vector<const DenseFeature*>* doc_feat_list,
    const vector<SnippetProfileList>* candidate_lists,
    const vector<unsigned>* num_prf,
    const vector<vector<float> >* train_weight,
    const vector<UPair>* irr_range,
    const vector<Labfile>* doc_lab,
    HMM_GMM* bg_model,
    int max_iter,
    vector<vector<float> >* output_like,
    vector<vector<float> >* anti_like) {

  dispatcher_ = dispatcher;
  query_feat_list_ = query_feat_list;
  doc_feat_list_ = doc_feat_list;
  candidate_lists_ = candidate_lists;
  num_prf_ = num_prf;
  train_weight_ = train_weight;
  irr_range_ = irr_range;
  doc_lab_ = doc_lab;
  bg_model_ = bg_model;
  max_iter_ = max_iter;
  output_like_ = output_like;
  anti_like_ = anti_like;


  /* Resize output variables */
  output_like_->resize(candidate_lists_->size()); // Number of queries
  anti_like_->resize(candidate_lists_->size()); // Number of queries

  for (unsigned qidx = 0; qidx < output_like_->size(); ++qidx) {
    // Number of hypothesized regions
    (*output_like_)[qidx].resize((*candidate_lists_)[qidx].size());
    (*anti_like_)[qidx].resize((*candidate_lists_)[qidx].size());
  }

  InitModels();
  InitInstances();

}/*}}}*/

void QueryHmmParamSet::InitModels() {/*{{{*/

  bg_model_->getpStatePool(UNUSE)->assign(
      bg_model_->getpStatePool(USE)->size(), NULL);
  bg_model_->getpGaussPool(UNUSE)->assign(
      bg_model_->getpGaussPool(USE)->size(), NULL);

  hmm_pairs.resize(candidate_lists_->size());

  for (unsigned qidx = 0; qidx < candidate_lists_->size(); ++qidx) {
    vector<int> state_seq;
    vector<float> likelihood_seq;
    const DenseFeature& feat = (*query_feat_list_)[qidx];
    bg_model_->CalLogBgOt<float>(feat.Data().begin(), feat.Data.end(), feat.LF());
    bg_model_->CalLogBjOtPxs(feat.LT());
    bg_model_->CalLogDelta(state_seq, &likelihood_seq, end_f);
    // FIXME: we need segment boundary right here!

  }

#ifdef NODEF
  for (unsigned qidx = 0; qidx < candidate_lists_->size(); ++qidx) {
    int s_frame = (*candidate_lists_)[qidx].Front().Boundary().first;
    int e_frame = (*candidate_lists_)[qidx].Front().Boundary().second;
    const Labfile& lab = (*doc_lab_)[(*candidate_lists_)[qidx].Front().Didx()];
    int first_lab, last_lab;
    for (first_lab = 0; first_lab < lab.getNumLab(); ++first_lab) {
      if (s_frame <= lab.getEndF(first_lab)) break;
    }
    for (last_lab = first_lab; last_lab < lab.getNumLab(); ++last_lab) {
      if (e_frame <= lab.getEndF(last_lab)) break;
    }

    cout << "lab: (" << first_lab << ", " << last_lab << ")" << endl;
    /* lab.getCluster(first_lab ~ last_lab) are the states we want */
    HmmPair& pair = hmm_pairs[qidx];
    pair.Init(lab, first_lab, last_lab + 1, *bg_model_, 0.2);
  }
#endif
}/*}}}*/

void QueryHmmParamSet::InitInstances() {/*{{{*/
  /* Init instances_ */
  instances.resize(candidate_lists_->size());
  for (unsigned qidx = 0; qidx < instances.size(); ++qidx) {

    instances[qidx].resize(2);

    /* target model's instances */
    float query_train_weight = 0.0;
    for (unsigned p = 0; p < (*num_prf_)[qidx]; ++p) {
      const SnippetProfile& snippet = (*candidate_lists_)[qidx].GetProfile(p);
      instances[qidx][0].push_back(
          Instance((*doc_feat_list_)[snippet.Didx()],
                   snippet.Boundary(),
                   (*train_weight_)[qidx][p]));
      query_train_weight += instances[qidx][0].back().weight_;
    }
    query_train_weight = 1.0;
    instances[qidx][0].push_back(
        Instance((*query_feat_list_)[qidx],
                 IPair(0, (*query_feat_list_)[qidx]->LT() - 1),
                 query_train_weight));

    cout << "Query[" << qidx << "] weight = " << query_train_weight << endl;
    /* anti model's  instances */
    for (unsigned p = (*irr_range_)[qidx].first;
         p <= (*irr_range_)[qidx].second; ++p) {
      const SnippetProfile& snippet =
        (*candidate_lists_)[qidx].GetProfile(p);
      instances[qidx][1].push_back(
          Instance((*doc_feat_list_)[snippet.Didx()],
                   snippet.Boundary(), 1.0));
    }
  }
}/*}}}*/

bool QueryHmmParamSet::InstallTrainerDecoder(QueryHmmTrainer* trainer, /*{{{*/
                                             QueryHmmDecoder* decoder) {

  pair<unsigned, int8_t>* ticket = dispatcher_->GetObjPtr();

  if (ticket) {
    int qidx = ticket->first;
    if (ticket->second == 0) { // target model
      trainer->InitTrainer(&hmm_pairs[qidx].targ_model,
                           &instances[qidx][0]);
      decoder->InitDecoder(&hmm_pairs[qidx].targ_model,
                           &(*candidate_lists_)[qidx],
                           doc_feat_list_,
                           &(*output_like_)[qidx]);

    } else if (ticket->second == 1) { // anti model
      trainer->InitTrainer(&hmm_pairs[qidx].anti_model,
                           &instances[qidx][1]);
      decoder->InitDecoder(&hmm_pairs[qidx].anti_model,
                           &(*candidate_lists_)[qidx],
                           doc_feat_list_,
                           &(*anti_like_)[qidx]);
    } else {
      assert(false);
    }
  }

  return ticket != NULL;
}/*}}}*/


void QueryHmmTrainer::Train(int max_iter) {/*{{{*/

  float dLogP = 10.0, oldLogP = LZERO, logP;
  const UpdateType udtype = UpdateAll;

  /* Calculate total number of frames */
  int nframe = 0;
  for (unsigned nprf = 0; nprf < instances_->size(); ++nprf) {
    nframe += (*instances_)[nprf].length();
  }


  /* Iterate EM */
  for (int iter = 0; iter < max_iter && (dLogP > 0.03 || dLogP < 0.0); ++iter) {

    /* E-step */
    model_->EMInitIter();
    logP = 0.0;
    for (unsigned nprf = 0; nprf < instances_->size(); ++nprf) {
      const Instance& ins = (*instances_)[nprf];
      logP += model_->EMObs<float>(
          ins.feat_->Data().begin() + ins.bound_.first,
          ins.feat_->Data().begin() + ins.bound_.second,
          ins.feat_->LF(), ins.weight_, udtype);
    }

    /* M-step */
    model_->EMUpdate(NULL, 0.15, udtype);
    logP /= nframe;

    dLogP = logP - oldLogP;
    oldLogP = logP;
  }

}/*}}}*/

template<typename T1, typename T2>
ostream& operator<<(ostream& os, const pair<T1, T2>& p) {
  os << "(" << p.first << ", " << p.second << ")";
  return os;
}


void QueryHmmDecoder::Decode() {/*{{{*/

  for (unsigned sidx = 0; sidx < candidates_->size(); ++sidx) {

    const SnippetProfile& snippet = candidates_->GetProfile(sidx);
    const IPair& bound = snippet.Boundary();
    const DenseFeature* feat = (*doc_feat_list_)[snippet.Didx()];
    int nframe = bound.second - bound.first + 1;

    model_->CalLogBgOt<float>(
        feat->Data().begin() + bound.first,
        feat->Data().begin() + bound.second + 1,
        feat->LF());
    model_->CalLogBjOtPxs(nframe);
    model_->CalLogAlpha(nframe);
    model_->CalLogPrO(nframe);
    (*scores_)[sidx] = model_->getPrO() / nframe;
  }
}/*}}}*/


void* QueryHmmManager::Run() {/*{{{*/

  /* Training */
  while (paramset_->InstallTrainerDecoder(trainer_, decoder_)) {
    trainer_->Train(paramset_->max_iter());
    decoder_->Decode();
  }

  return NULL;
}/*}}}*/

