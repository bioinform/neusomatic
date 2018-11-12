#ifndef NEU_SOMATIC_MSA_UTIL 
#define NEU_SOMATIC_MSA_UTIL 

#include <vector>
#include <string>

#include "change_coordinates.hpp"

namespace neusomatic{

class Col{
  private:
    using Idx = unsigned;
    using Base = int;
    std::vector<Base> bases_;
    std::vector<int> bquals_;

  public:
    static const int ALPHABET_SIZE = 6; // a, t, c, g, gap and missing_char
    static const int TAG_SIZE = 5; // number of tag used. 
    Col() = delete;
    explicit Col(size_t nbases): bases_(nbases), bquals_(nbases), base_freq_(ALPHABET_SIZE), //base_rids_(ALPHABET_SIZE),
                                 bqual_mean(ALPHABET_SIZE), mqual_mean(ALPHABET_SIZE), strand_mean(ALPHABET_SIZE), lsc_mean(ALPHABET_SIZE),
                                 rsc_mean(ALPHABET_SIZE), tag_mean(ALPHABET_SIZE, std::vector<float>(TAG_SIZE)) {}

    std::vector<int> base_freq_;
    std::vector<float> bqual_mean;
    std::vector<float> mqual_mean;
    std::vector<float> strand_mean;
    std::vector<int> lsc_mean;
    std::vector<int> rsc_mean;
    std::vector<std::vector<float>> tag_mean;
    decltype(auto) bases() const {return (bases_);}
    decltype(auto) bquals() const {return (bquals_);}

    void emplace(const Idx& id, const Base& b) {
      bases_[id] = b;
    }

    void emplace_qual(const Idx& id, const int& b) {
      bquals_[id] = b;
    }

    decltype(auto) at(Idx i) const {
      return bases_.at(i);
    }

    decltype(auto) size() const {
      return bases_.size();
    }

    decltype(auto) operator[](Idx i) {
      return (bases_[i]);
    }
};

template<typename RefGap, typename Base>
class CondensedArray{
public:
  using Idx = unsigned;
  static const unsigned char missing_chr_ = '~';

  static int DnaCharToDnaCode(const char& dna) {
    switch(dna) {
      case 'A':
      case 'a':
        return 0;
      case 'C':
      case 'c':
        return 1;
      case 'G':
      case 'g':
        return 2;
      case 'T':
      case 't':
        return 3;
      case '-':
      //case 'N': // do not want N to count
        return 4;
      default:
        return 5;
    }
  }

public:
  using Val = Base;
  using MatrixIdx =  Idx;
  using Row = std::unordered_map<Idx, Val>;
  using ColSpace  = std::vector<Col>;
  using RowSpace = std::vector<Row>;
  using TId = unsigned;
  using SnvTuple = std::tuple<double, std::vector<Idx>>;

  decltype(auto) GetColSpace() const {
    return (cspace_);
  }

  decltype(auto) GetColSpaceMQual() const {
    return (cspace_);
  }
  decltype(auto) GetColSpaceStrand() const {
    return (cspace_);
  }
  decltype(auto) GetColSpaceLSC() const {
    return (cspace_);
  }
  decltype(auto) GetColSpaceRSC() const {
    return (cspace_);
  }
  decltype(auto) GetColSpaceTag() const {
    return (cspace_);
  }


  size_t ncol() const {
    return cspace_.size();
  }

  size_t nrow() const {
    return nrow_;
  }

  template<typename GInv>
  explicit CondensedArray(const std::vector<std::string>& msa, const int& total_cov, const GInv& ginv, const RefGap& refgap, const int num_thread = 4):
      nrow_(msa.size()), cspace_(msa[0].size(), Col(msa.size())), bound_(ginv.left(), ginv.right()), //col_major_bases_(msa[0].size()), 
     total_cov_(total_cov)
  {
    _CheckInput(msa); 
    #pragma omp parallel for schedule(dynamic, 256) num_threads(num_thread)
    for (size_t i = 0; i < msa.size(); ++i) {
      auto dna5qseq = _StringToDnaInt(msa[i]);
      this->row_push(dna5qseq.begin(), dna5qseq.end(), i);
    }
  }

  template<typename GInv>
  explicit CondensedArray(const std::vector<std::string>& msa, const std::vector<std::string>& bqual, 
                    const std::vector<int>& mqual, const std::vector<int>& strand, 
                    const std::vector<int>& lsc, const std::vector<int>& rsc,
                    const std::vector<std::vector<int>>& tags, 
                    const int& total_cov, const GInv& ginv, const RefGap& refgap, const int num_thread = 4): 
            nrow_(msa.size()), 
            total_cov_(total_cov),
            cspace_(msa[0].size(), Col(msa.size())),
            bound_(ginv.left(), ginv.right()), 
            mquals_(nrow_),
            strands_(nrow_), 
            lsc_(nrow_),
            rsc_(nrow_),
            tags_(nrow_, std::vector<int>(5))
  {

    _CheckInput(msa); 

    #pragma omp parallel num_threads(num_thread)
    {
      //if (omp_get_thread_num() == 0) {
        //std::cerr << "using " << omp_get_num_threads() << " threads\n";
      //}

    #pragma omp for schedule(dynamic, 256) 
    for (size_t i = 0; i < msa.size(); ++i) {
      auto dna5qseq = _StringToDnaInt(msa[i]);
      this->row_push(dna5qseq.begin(), dna5qseq.end(), i);
      this->row_push_bqual(bqual[i].begin(), bqual[i].end(), i);
      this->row_push_lsc(lsc[i], i);
      this->row_push_rsc(rsc[i], i);
      this->row_push_mqual(mqual[i], i);
      this->row_push_strand(strand[i], i);
      this->row_push_tag(tags[i], i);
    }

    } //end parallel
  }


  template<typename BaseItr>
  void row_push(BaseItr b, BaseItr e, const Idx rid) {
    BaseItr it = b;
    size_t cid = 0;
    for (; it != e; ++it, ++cid) {
      Push_(rid, cid, *it);
    }
  }

  template<typename BaseItr>
  void row_push_bqual(BaseItr b, BaseItr e, const Idx rid) {
    BaseItr it = b;
    size_t cid = 0;
    for (; it != e; ++it, ++cid) {
      PushBQual_(rid, cid, *it);
    }
  }

  void row_push_mqual(int mqual, const Idx rid) {
    mquals_[rid] = mqual;
  }

  void row_push_strand(int strand, const Idx rid) {
    strands_[rid] = strand;
  }

  void row_push_lsc(int pos, const Idx rid) {
    lsc_[rid] = pos;
  }

  void row_push_rsc(int pos, const Idx rid) {
    rsc_[rid] = pos;
  }

  void row_push_tag(std::vector<int> tags, const Idx rid) {
    for (size_t i = 0; i < tags.size(); ++i) {
      tags_[rid][i] = tags[i];
    } 
  }

  template<class BaseItr>
  void col_push(BaseItr b, BaseItr e, const Idx cid) {
    BaseItr it = b;
    size_t i = 0;
    for (; it != e; ++it, ++i) {
      Push_(i, cid, *it);
    }
  }

  void Init(const int num_thread) {
    for (size_t i = 0; i < ncol(); ++i) {
      auto& col = cspace_[i];
      for (size_t j = 0; j < nrow(); ++j) {
        if (col.bases()[j] == missing_chr_) continue;
        col.base_freq_[col.bases()[j]]++;
      }
    }
  }

  void InitWithAlnMetaData(const int num_thread) {
    for (size_t i = 0; i < ncol(); ++i) {
      //column-wise 
      auto& col = cspace_[i];

      for (auto b = 0; b < Col::ALPHABET_SIZE; ++ b) {
        col.bqual_mean[b]=0;
        col.mqual_mean[b]=0;
        col.strand_mean[b]=0;
        col.lsc_mean[b]=0;
        col.rsc_mean[b]=0;
        for (size_t ii = 0; ii < Col::TAG_SIZE; ++ii) {
          col.tag_mean[b][ii]=0;
        }
      }

      //element-wise
      for (size_t j = 0; j < nrow(); ++j) {
        col.base_freq_[col.bases()[j]]++;
        col.bqual_mean[col.bases()[j]]+=float(col.bquals()[j]-33);
        col.mqual_mean[col.bases()[j]]+=mquals_[j];
        col.strand_mean[col.bases()[j]]+=strands_[j];
        col.lsc_mean[col.bases()[j]] += lsc_[j] == i ? 1 : 0;
        col.rsc_mean[col.bases()[j]] += rsc_[j] == i ? 1 : 0;
        for (size_t ii = 0; ii < Col::TAG_SIZE; ++ii) {
          col.tag_mean[col.bases()[j]][ii] += tags_[j][ii];
        }
      }

      for (auto s = 0; s < Col::ALPHABET_SIZE; ++ s) {
        if (col.base_freq_[s] == 0) continue;
        col.bqual_mean[s]/=col.base_freq_[s];
        col.mqual_mean[s]/=col.base_freq_[s];
        col.strand_mean[s]/=col.base_freq_[s];
        col.strand_mean[s]*=100.0;
        for (size_t ii = 0; ii < Col::TAG_SIZE; ++ii) {
          col.tag_mean[s][ii]/=col.base_freq_[s];
        }
      }
    }
  }


  decltype(auto) total_cov() const {
    return (total_cov_);
  }

  decltype(auto) cspace() const {
    return (cspace_);
  }

  decltype(auto) bound() const {
    return (bound_);
  }


private:

  // used only when use seqan interface
  // now obsolete
  const std::vector<TId> aids_; 

  size_t nrow_;
  const int total_cov_;
  ColSpace cspace_;
  const std::pair<Idx, Idx> bound_;
  std::vector<int> mquals_;
  std::vector<int> strands_;
  std::vector<int> lsc_;
  std::vector<int> rsc_;
  std::vector<std::vector<int>> tags_;

  void _CheckInput(const std::vector<std::string>& msa) {
    unsigned ncol = 0;
    if (msa.empty()) {
      throw std::runtime_error("empty msa");
    }
    for(auto const& row : msa) {
      if (ncol == 0) ncol = row.length(); 
      else {
        if (ncol != row.length()) {
          throw std::runtime_error("input msa has unequal row lengthes");
        }
      }
    }
  }

  std::vector<Base> _StringToDnaChar(const std::string &s) {
    std::vector<Base> dna5qseq(s.size());
    for (size_t j = 0; j < s.size(); ++j) {
      switch(s[j]) {
        case 'A':
        case 'a':
          dna5qseq[j] = 'A';
          break;
        case 'C':
        case 'c':
          dna5qseq[j] = 'C';
          break;
        case 'T':
        case 't':
          dna5qseq[j] = 'T';
          break;
        case 'G':
        case 'g':
          dna5qseq[j] = 'G';
          break;
        case '-':
          dna5qseq[j] = '-'; 
          break;
        default:  
          dna5qseq[j] = missing_chr_; 
          break;
      }
    }
    return dna5qseq;
  }

  std::vector<Base> _StringToDnaInt(const std::string &s) {
    std::vector<Base> dna5qseq(s.size());
    for (size_t j = 0; j < s.size(); ++j) {
      dna5qseq[j] = DnaCharToDnaCode(s[j]);
    }
    return dna5qseq;
  }

  void PushBQual_(Idx row, Idx col, Val val) {
    cspace_[col].emplace_qual(row, val);
  }


  void Push_(Idx row, Idx col, Val val) {
    cspace_[col].emplace(row, val);
  }

};

template<typename RefGap, typename GInv>
decltype(auto) CreateCondensedArray(const std::vector<std::string>& msa, const int total_cov, const GInv& ginv, const RefGap& refgap, const int num_threads) {
  using CondensedArray = neusomatic::CondensedArray<RefGap, int>;
  CondensedArray condensed_array(msa, total_cov, ginv, refgap, num_threads);
  condensed_array.Init(num_threads);
  return condensed_array;
}


template<typename RefGap, typename GInv>
decltype(auto) CreateCondensedArray(const std::vector<std::string>& msa, const std::vector<std::string>& bqual,  const std::vector<int>& mqual, const std::vector<int>& strand,
                              const std::vector<int>& lscs, const std::vector<int>& rscs,
                              const std::vector<std::vector<int>>& tag,
                              const int total_cov, const GInv& ginv, const RefGap& refgap, const int num_threads) {
  using CondensedArray = neusomatic::CondensedArray<RefGap, int>; 
  CondensedArray condensed_array(msa, bqual, mqual, strand, lscs, rscs, tag, total_cov, ginv, refgap, num_threads);
  condensed_array.InitWithAlnMetaData(num_threads);
  return condensed_array;
}
}// end neusomatic

#endif
