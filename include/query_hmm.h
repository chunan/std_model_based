#ifndef __QUERY_HMM_H__
#define __QUERY_HMM_H__

#include "hmmlite.h"
#include "feature.h"
#include "std_common.h"
#include "thread_util.h"

using StdCommonUtil::SnippetProfileList;
using StdCommonUtil::IPair;

struct HmmPair {/*{{{*/
  void Install(int nstate, int use, HMM_GMM& bg_model);
  void Init(const Labfile& lab, const int start, const int end,
            HMM_GMM& bg_model, const float del_ratio);

  HMM_GMM targ_model;
  HMM_GMM anti_model;
};/*}}}*/

// TODO:
//  v1. Write QueryHmmTrainer::Train(int max_iter)
//  v2. Write QueryHmmViterbi::Decode()
//  v3. Write QueryHmmParamSet::InstallTrainer()
//  v4. Write QueryHmmParamSet::InstallDecoder()


struct Instance {/*{{{*/
  public:
    Instance(const DenseFeature* feat,
             const IPair bound,
             const float weight) {
      feat_ = feat;
      bound_ = bound;
      weight_ = weight;
    }
    int length() const { return bound_.second - bound_.first; }
    const DenseFeature* feat_;
    IPair bound_;
    float weight_;
};/*}}}*/


class QueryHmmTrainer {
  public:
    void Train(int max_iter);
    void InitTrainer(HMM_GMM* model,
                     const vector<Instance>* instances) {
      instances_ = instances;
      model_ = model;
    }
  private:
    const vector<Instance>* instances_;
    HMM_GMM* model_;
};

class QueryHmmDecoder {
  public:
    void Decode();
    void InitDecoder(HMM_GMM* model,
                     const SnippetProfileList* candidates,
                     const vector<const DenseFeature*>* doc_feat_list,
                     vector<float>* scores) {

      model_ = model;
      candidates_ = candidates;
      scores_ = scores;
      doc_feat_list_ = doc_feat_list;
    }

  private:
    HMM_GMM* model_;
    const SnippetProfileList* candidates_;
    const vector<const DenseFeature*>* doc_feat_list_;
    vector<float>* scores_;
};


class QueryHmmParamSet {/*{{{*/
  public:
    QueryHmmParamSet(/*{{{*/
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
        vector<vector<float> >* anti_like);/*}}}*/
    bool InstallTrainerDecoder(QueryHmmTrainer* trainer,
                               QueryHmmDecoder* decoder);
    int max_iter() const { return max_iter_; }
    //bool InstallDecoder(QueryHmmDecoder* decoder);
  private:
    void InitModels();
    void InitInstances();

    /* Input data */
    Dispatcher<pair<unsigned, int8_t> >* dispatcher_;
    const vector<const DenseFeature*>* query_feat_list_;
    const vector<const DenseFeature*>* doc_feat_list_;
    const vector<SnippetProfileList>* candidate_lists_;
    const vector<unsigned>* num_prf_;
    const vector<vector<float> >* train_weight_;
    const vector<UPair>* irr_range_;
    const vector<Labfile>* doc_lab_;
    HMM_GMM* bg_model_;
    int max_iter_;

    /* Output data */
    vector<vector<float> >* output_like_;
    vector<vector<float> >* anti_like_;

    /* HMMs and training samples */
    vector<HmmPair> hmm_pairs;
    vector<vector<vector<Instance> > > instances;

};/*}}}*/


class QueryHmmManager : public ThreadRunner {/*{{{*/
  public:
    void Install(QueryHmmParamSet* paramset,
                 QueryHmmTrainer* trainer,
                 QueryHmmDecoder* decoder) {
      paramset_ = paramset;
      trainer_ = trainer;
      decoder_ = decoder;
    }
    virtual void* Run();
  private:
    QueryHmmParamSet* paramset_;
    QueryHmmTrainer* trainer_;
    QueryHmmDecoder* decoder_;
};/*}}}*/

#endif
