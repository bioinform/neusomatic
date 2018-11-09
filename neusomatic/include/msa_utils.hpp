#ifndef NEU_SOMATIC_MSA_UTIL 
#define NEU_SOMATIC_MSA_UTIL 

#include <vector>
#include <string>

#include "change_coordinates.hpp"

namespace neusomatic{

template<typename Base>
inline bool IsGap(const Base& b) {
  if (b == 'N' || b == '-' || b == 4) {
    return true;
  }
  return false;
}

template<class Itr>
auto SortedIndices(const Itr b, const Itr e) {
  std::vector<size_t> idx(std::distance(b,e));
  std::iota(idx.begin(), idx.end(), 0);

  sort(idx.begin(), idx.end(),
       [&b](size_t i1, size_t i2) {return *(b+i1) < *(b+i2);});
  return idx;
}

class Col{
  private:
    using Idx = unsigned;
    using Base = int;
    float weight_;
    std::vector<Base> bases_;
    std::vector<int> bquals_;

  public:
    static const int alphabet_size_ = 6; // a, t, c, g, gap and missing_char
    Col() = delete;
    explicit Col(size_t nbases): weight_(0.0), bases_(nbases), bquals_(nbases), base_freq_(alphabet_size_), base_rids_(alphabet_size_),
                                 bqual_mean(alphabet_size_), mqual_mean(alphabet_size_), strand_mean(alphabet_size_), lsc_mean(alphabet_size_),
                                 rsc_mean(alphabet_size_), tag_mean(alphabet_size_) {}

    //std::map<Base, int> base_freq_;
    std::vector<std::vector<Idx>> base_rids_;
    std::vector<Base> base_freq_;
    std::vector<float> bqual_mean;
    std::vector<float> mqual_mean;
    std::vector<float> strand_mean;
    std::vector<int> lsc_mean;
    std::vector<int> rsc_mean;
    std::vector<float> tag_mean;
    //decltype(auto) weight() { return (weight_); }
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

    decltype(auto) weight() const {
      return (weight_);
    }

    decltype(auto) weight(){
      return (weight_);
    }

    //friend std::ostream& operator<<(std::ostream& dest, Col const& col) {
      //for (auto const& b_c: col.base_freq_) {
        //dest<<"("<<b_c.first<<"): "<<b_c.second<<std::endl;
      //}
      //return dest;
    //}
};

template<typename RefGap, typename Base>
class CondensedArray{
public:
  using Idx = unsigned;
  //using DnaBase = Base;
  //static_assert(std::is_integral<Idx>::value, "integer required");
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
      case 'N':
        return 4;
      default:
        return 5;
    }
  }

public:
  class LinkedCol{
  public:
    using Dimer = std::pair<Base, Base>;
    explicit LinkedCol(const Col& left, const Col& right, const int& l, const int& r):
        lpos_(l), rpos_(r){
      if (left.size() != right.size()) {
        throw std::runtime_error("Link two columns with unequal lengths");
      }
      size_t ncol = left.size();
      for (size_t i = 0; i < ncol; i++) {
        dimer_count_[std::make_pair(left.at(i), right.at(i))] ++ ;
        dimer_rid_[std::make_pair(left.at(i), right.at(i))].push_back(i);
      }
    }

    friend std::ostream& operator<< (std::ostream& dest, LinkedCol const& obj) {
      dest<<"cols: "<<obj.lpos_ + 1<<","<<obj.rpos_ + 1<<std::endl;
      for (auto const& di_c: obj.dimer_count_) {
        dest<<"("<<di_c.first.first<<","<<di_c.first.second<<"): "<<di_c.second<<std::endl;
      }
      return dest;
    }

    int lpos_;
    int rpos_;
    LinkedCol() = default;
    std::map<Dimer, int> dimer_count_; 
    std::map<Dimer, std::vector<Idx>> dimer_rid_;
  };

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
  decltype(auto) GetColSpaceBQual() const {
    return (cspace_bqual_);
  }
  decltype(auto) GetColSpaceMQual() const {
    return (cspace_mqual_);
  }
  decltype(auto) GetColSpaceStrand() const {
    return (cspace_strand_);
  }
  decltype(auto) GetColSpaceLSC() const {
    return (cspace_lsc_);
  }
  decltype(auto) GetColSpaceRSC() const {
    return (cspace_rsc_);
  }
  decltype(auto) GetColSpaceTag() const {
    return (cspace_tag_);
  }


  decltype(auto) GetCC() const {
    return (cc_);
  }
  size_t ncol() const {
    return cspace_.size();
  }

  size_t nrow() const {
    return nrow_;
  }

  template<typename GInv>
  explicit CondensedArray(const std::vector<std::string>& msa, const int& total_cov, const GInv& ginv, const RefGap& refgap, const int num_thread = 4):
      nrow_(msa.size()), cspace_(msa[0].size(), Col(msa.size())), bound_(ginv.left(), ginv.right()), col_major_bases_(msa[0].size()), ref_gaps_(refgap), 
      cc_(ref_gaps_), total_cov_(total_cov)
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
                    const std::vector<int>& mqual, const std::vector<bool>& strand, 
                    const std::vector<std::string>& lsc, const std::vector<std::string>& rsc,
                    const std::vector<std::vector<int>>& tag, 
                    const int& total_cov, const GInv& ginv, const RefGap& refgap, const int num_thread = 4): 
      nrow_(msa.size()), cspace_(msa[0].size(), Col(msa.size())), cspace_bqual_(msa[0].size(), Col(msa.size())), cspace_mqual_(msa[0].size(), Col(msa.size())),
            cspace_strand_(msa[0].size(), Col(msa.size())), cspace_lsc_(msa[0].size(), Col(msa.size())),
            cspace_rsc_(msa[0].size(), Col(msa.size())),
            cspace_tag_(5,ColSpace(msa[0].size(), Col(msa.size()))),
            bound_(ginv.left(), ginv.right()), col_major_bases_(msa[0].size()), ref_gaps_(refgap), cc_(ref_gaps_), total_cov_(total_cov)
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
      this->row_push_lsc(lsc[i].begin(), lsc[i].end(), i);
      this->row_push_rsc(rsc[i].begin(), rsc[i].end(), i);
      this->row_push_mqual(mqual[i], bqual[i].size(), i);
      this->row_push_strand(strand[i], bqual[i].size(), i);
      this->row_push_tag(tag[i], bqual[i].size(), i);
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

  void row_push_mqual(int mqual, const int ncol, const Idx rid) {
    for (size_t cid = 0; cid < ncol; cid++) {
      PushMQual_(rid, cid, (unsigned char) mqual);
    }
  }

  void row_push_strand(bool strand, const int ncol, const Idx rid) {
    for (size_t cid = 0; cid < ncol; cid++) {
      PushStrand_(rid, cid, (unsigned char) strand);
    }
  }

  template<typename BaseItr>
  void row_push_lsc(BaseItr b, BaseItr e, const Idx rid) {
    BaseItr it = b;
    size_t cid = 0;
    for (; it != e; ++it, ++cid) {
      PushLSC_(rid, cid, *it);
    }
  }

  template<typename BaseItr>
  void row_push_rsc(BaseItr b, BaseItr e, const Idx rid) {
    BaseItr it = b;
    size_t cid = 0;
    for (; it != e; ++it, ++cid) {
      PushRSC_(rid, cid, *it);
    }
  }

  void row_push_tag(std::vector<int> tag, const int ncol, const Idx rid) {
    for (size_t cid = 0; cid < ncol; cid++) {
      PushTag_(rid, cid, tag);
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

  void Init() {
    for (size_t i = 0; i < ncol(); ++i) {
      auto& col = cspace_[i];
      for (size_t j = 0; j < nrow(); ++j) {
        if (col.bases()[j] == missing_chr_) continue;
        col.base_freq_[col.bases()[j]]++;
        col.base_rids_[col.bases()[j]].push_back(j);
      }
      col_major_bases_[i] = ColMajority_(i);
    }
  }

  void InitWithAlnMetaData() {
    for (size_t i = 0; i < ncol(); ++i) {
      //column-wise 
      auto& col = cspace_[i];
      auto& col_bqual = cspace_bqual_[i];
      auto& col_mqual = cspace_mqual_[i];
      auto& col_strand = cspace_strand_[i];
      auto& col_lsc = cspace_lsc_[i];
      auto& col_rsc = cspace_rsc_[i];

      col_major_bases_[i] = ColMajority_(i);// may not be needed

      //std::list<Base> bases_list = { 'A', 'C', 'G', 'T', '~', '-' };
      for (auto b = 0; b < Col::alphabet_size_; ++ b) {
        col_bqual.bqual_mean[b]=0;
        col_mqual.mqual_mean[b]=0;
        col_strand.strand_mean[b]=0;
        col_lsc.lsc_mean[b]=0;
        col_rsc.rsc_mean[b]=0;
        for (size_t ii = 0; ii < 5; ++ii) {
          cspace_tag_[ii][i].tag_mean[b]=0;
        }
      }

      //element-wise
      for (size_t j = 0; j < nrow(); ++j) {
        //if (calculate_entropy && col.bases()[j] == missing_chr_) continue;
        col.base_freq_[col.bases()[j]]++;
        col.base_rids_[col.bases()[j]].push_back(j);
        col.bqual_mean[col.bases()[j]]+=float(int(col_bqual.bases()[j])-33);
        col.mqual_mean[col.bases()[j]]+=float(int(col_mqual.bases()[j])/70.0);
        col.strand_mean[col.bases()[j]]+=float(int(col_strand.bases()[j]));
        col_lsc.lsc_mean[col.bases()[j]]+=int(col_lsc.bases()[j]=='1');
        col_rsc.rsc_mean[col.bases()[j]]+=int(col_rsc.bases()[j]=='1');
        for (size_t ii = 0; ii < 5; ++ii) {
          cspace_tag_[ii][i].tag_mean[col.bases()[j]]+=float(int(cspace_tag_[ii][i].bases()[j]))/100.0;
        }
      }

      //for(auto it = col.base_freq_.cbegin(); it != col.base_freq_.cend(); ++it) {
      for (auto s = 0; s < Col::alphabet_size_; ++ s) {
        //col.bqual_mean[it->first]/=col.base_freq_[it->first];
        col.bqual_mean[s]/=col.base_freq_[s];
        col.mqual_mean[s]/=col.base_freq_[s];
        col.mqual_mean[s]*=70.0;
        col.strand_mean[s]/=col.base_freq_[s];
        col.strand_mean[s]*=100.0;
        for (size_t ii = 0; ii < 5; ++ii) {
          cspace_tag_[ii][i].tag_mean[s]/=col.base_freq_[s];
          cspace_tag_[ii][i].tag_mean[s]*=100.0;
        }
      }
    }
  }

  decltype(auto) sorted_entropy_cidx_pairs() const {
    return (sorted_entropy_cidx_pairs_);
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

  decltype(auto) change_coord() const {
    return (cc_);
  }

  decltype(auto) ref_normalized_entropy() const {
    return (ref_nomalized_entropy_);
  }

private:

  // used only when use seqan interface
  // now obsolete
  const std::vector<TId> aids_; 

  size_t nrow_;
  ColSpace cspace_;
  ColSpace cspace_bqual_;
  ColSpace cspace_mqual_;
  ColSpace cspace_strand_;
  ColSpace cspace_lsc_;
  ColSpace cspace_rsc_;
  std::vector<ColSpace> cspace_tag_;
  const std::pair<Idx, Idx> bound_;
  std::vector<Base> col_major_bases_;
  std::vector<float> ref_nomalized_entropy_;
  std::vector<std::pair<float, Idx>> sorted_entropy_cidx_pairs_;
  const RefGap ref_gaps_;
  const neusomatic::ChangeCoordinates<RefGap> cc_;
  const int total_cov_;
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
      dna5qseq[j] = DnaCharToDnaCode(j);
    }
    return dna5qseq;
  }

  void PushBQual_(Idx row, Idx col, Val val) {
    cspace_bqual_[col].emplace(row, val);
  }

  void PushMQual_(Idx row, Idx col, Val val) {
    cspace_mqual_[col].emplace(row, val);
  }

  void PushStrand_(Idx row, Idx col, Val val) {
    cspace_strand_[col].emplace(row, val);
  }

  void PushLSC_(Idx row, Idx col, Val val) {
    cspace_lsc_[col].emplace(row, val);
  }

  void PushRSC_(Idx row, Idx col, Val val) {
    cspace_rsc_[col].emplace(row, val);
  }

  void PushTag_(Idx row, Idx col, std::vector<int> val) {
    cspace_tag_[0][col].emplace(row, (unsigned char) val[0]);
    cspace_tag_[1][col].emplace(row, (unsigned char) val[1]);
    cspace_tag_[2][col].emplace(row, (unsigned char) val[2]);
    cspace_tag_[3][col].emplace(row, (unsigned char) val[3]);
    cspace_tag_[4][col].emplace(row, (unsigned char) val[4]);
  }



  void Push_(Idx row, Idx col, Val val) {
    cspace_[col].emplace(row, val);
  }

  Base ColMajority_(Idx ci) const {
    Base result = 4;
    int m = 0;
    auto const& col = cspace_[ci];
    //for (auto const& b_f: col.base_freq_) {
    for (int i = 0; i < (int) col.base_freq_.size(); ++i) {
      if (col.base_freq_[i] > m) {
        result = i;
        m = col.base_freq_[i];
      }
    }
    return result;
  }

  void SortCols_() { 
    sorted_entropy_cidx_pairs_.clear();
    sorted_entropy_cidx_pairs_.resize(ncol());
    for (Idx i = 0; i < ncol(); ++i) {
      sorted_entropy_cidx_pairs_[i] = std::make_pair(cspace_[i].weight(), i); 
    }
    //sort in decreasing order
    std::sort(sorted_entropy_cidx_pairs_.begin(), sorted_entropy_cidx_pairs_.end(),
        [](const std::pair<float, Idx> &left,
          const std::pair<float, Idx> &right) {return left.first > right.first;});
  }

};

template<typename RefGap, typename GInv>
decltype(auto) CreateCondensedArray(const std::vector<std::string>& msa, const int total_cov, const GInv& ginv, const RefGap& refgap, const int num_threads) {
  using CondensedArray = neusomatic::CondensedArray<RefGap, int>;
  CondensedArray condensed_array(msa, total_cov, ginv, refgap, num_threads);
  condensed_array.Init();
  return condensed_array;
}


template<typename RefGap, typename GInv>
decltype(auto) CreateCondensedArray(const std::vector<std::string>& msa, const std::vector<std::string>& bqual,  const std::vector<int>& mqual, const std::vector<bool>& strand,
                              const std::vector<std::string>& lscs, const std::vector<std::string>& rscs,
                              const std::vector<std::vector<int>>& tag,
                              const int total_cov, const GInv& ginv, const RefGap& refgap, const int num_threads) {
  using CondensedArray = neusomatic::CondensedArray<RefGap, int>; 
  CondensedArray condensed_array(msa, bqual, mqual, strand, lscs, rscs, tag, total_cov, ginv, refgap, num_threads);
  condensed_array.InitWithAlnMetaData();
  return condensed_array;
}
}// end neusomatic

#endif
