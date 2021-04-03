#ifndef NEU_SOMATIC_MSA_UTIL 
#define NEU_SOMATIC_MSA_UTIL 

#include <vector>
#include <string>

#include "change_coordinates.hpp"
#include "msa.hpp"

namespace neusomatic{

class Col{
  private:
    using Idx = unsigned;
    using Base = int;
    std::vector<Base> bases_;
    std::vector<int> bquals_;

  public:
    static const int ALPHABET_SIZE = 6; // a, t, c, g, gap and missing_char
    static const int TAG_SIZE = 5;
    explicit Col(size_t nbases): bases_(nbases), bquals_(nbases), base_freq_(ALPHABET_SIZE), //base_rids_(ALPHABET_SIZE),
                                 bqual_mean(ALPHABET_SIZE), mqual_mean(ALPHABET_SIZE), strand_mean(ALPHABET_SIZE), lsc_mean(ALPHABET_SIZE),
                                 rsc_mean(ALPHABET_SIZE), tag_mean(ALPHABET_SIZE, std::vector<float>(TAG_SIZE)) {}

    Col(): base_freq_(ALPHABET_SIZE), 
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

template<typename Base>
class CondensedArray{
public:
  using Idx = unsigned;
  static const unsigned char missing_chr_ = '~';
  static const int DEL_CODE = 4;
  static const int GAP_CODE = 4;
  static const int MISSING_CODE = 5;

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
  using ColSpace  = std::vector<Col>;
  using TId = unsigned;

  decltype(auto) GetColSpace() const {
    return (cspace_);
  }

  decltype(auto) GetBaseFreq(int i) const {
    return (cspace_[i].base_freq_);
  }

  decltype(auto) GetBQMean(int i) const {
    return (cspace_[i].bqual_mean);
  }

  decltype(auto) GetMQMean(int i) const {
    return (cspace_[i].mqual_mean);
  }

  decltype(auto) GetStrandMean(int i) const {
    return (cspace_[i].strand_mean);
  }

  decltype(auto) GetLSCMean(int i) const {
    return (cspace_[i].lsc_mean);
  }

  decltype(auto) GetRSCMean(int i) const {
    return (cspace_[i].rsc_mean);
  }

  decltype(auto) GetTagMean(int i) const {
    return (cspace_[i].tag_mean);
  }

  size_t ncol() const {
    return cspace_.size();
  }

  size_t nrow() const {
    return nrow_;
  }

  decltype(auto) GetGappedRef() const {
    return (gapped_ref_str_);
  }

  CondensedArray() : nrow_(0)  {} 

  template<typename GInv>
  explicit CondensedArray(const std::vector<std::string>& msa, const int& total_cov, const GInv& ginv, const std::string& refgap) :
      nrow_(msa.size()), 
      cspace_(msa[0].size(), 
      Col(msa.size())),
      gapped_ref_str_ (refgap)
  {
    _CheckInput(msa); 
    for (size_t i = 0; i < msa.size(); ++i) {
      auto dna5qseq = _StringToDnaInt(msa[i]);
      this->row_push(dna5qseq.begin(), dna5qseq.end(), i);
    }
  }

  template<typename MSA>
  explicit CondensedArray(const MSA& msa, const bool calculate_qual) :
    nrow_(msa.size()), 
    cspace_(msa.ncol()),
    gapped_ref_str_(msa.gapped_ref_str())
  { 

    // auto start_3p0 = std::chrono::system_clock::now();


    for (int pp = 0; pp < msa.ncol(); ++pp) {
      cspace_[pp].base_freq_[MISSING_CODE] = nrow_;
    }

    for (size_t i = 0; i < nrow_; ++i) {
      // auto start_3p0p0 = std::chrono::system_clock::now();

      auto const& r = msa.bam_records()[i];
      auto const& cigar = r.GetCigar();        
      const auto seq_qual = msa.GetGappedSeqAndQual(r);
      const auto& seq = seq_qual.first;

      // -1,+1 for INS at the begin or end of a read
      auto s = msa.GapPosition(std::max(r.Position() - 1 , msa.ref_gaps().left()) - msa.ref_gaps().left());
      auto e = msa.GapPosition(std::min(r.PositionEnd() + 1, msa.ref_gaps().right()) - msa.ref_gaps().left() - 1);

      // for (int pp = 0; pp < s; ++pp) {
      //   ++cspace_[pp].base_freq_[MISSING_CODE];
      // }


      // for (int pp = e + 1; pp < msa.ncol(); ++pp){
      //   ++cspace_[pp].base_freq_[MISSING_CODE];
      // }

      // auto end_3p0p0 = std::chrono::system_clock::now();
      // std::chrono::duration<double> elapsed_seconds_3p0p0 = end_3p0p0-start_3p0p0;
      // std::time_t end_time_3p0p0 = std::chrono::system_clock::to_time_t(end_3p0p0);
      // std::cout << "elapsed_3p0p0 time: " << elapsed_seconds_3p0p0.count() << "s\n";      
      // auto start_3p0p1 = std::chrono::system_clock::now();

      for (int pp = s; pp <= e; ++pp) {
        --cspace_[pp].base_freq_[MISSING_CODE];
        ++cspace_[pp].base_freq_[DnaCharToDnaCode(seq[pp])];
      }


      // for (int pp = s; pp <= e; ++pp) {
      //   std::map<int, int> alt_counts;
      //   int dp = 0;
      //   auto ref_base = ref[pp];
      //   Compute_ALTs(alt_counts, cspace_[pp].base_freq_, dp, ref_base, report_all_alleles, del_min_af, snp_min_af, ins_min_af);


      // auto end_3p0p1 = std::chrono::system_clock::now();
      // std::chrono::duration<double> elapsed_seconds_3p0p1 = end_3p0p1-start_3p0p1;
      // std::time_t end_time_3p0p1 = std::chrono::system_clock::to_time_t(end_3p0p1);
      // std::cout << "elapsed_3p0p1 time: " << elapsed_seconds_3p0p1.count() << "s\n";      
      // auto start_3p0p2 = std::chrono::system_clock::now();

      if (calculate_qual) {
        // auto start_3p0p2p0 = std::chrono::system_clock::now();

        const auto& qual = seq_qual.second;
        int strand = (int) !r.ReverseFlag();

        // auto end_3p0p2p0 = std::chrono::system_clock::now();
        // std::chrono::duration<double> elapsed_seconds_3p0p2p0 = end_3p0p2p0-start_3p0p2p0;
        // std::time_t end_time_3p0p2p0 = std::chrono::system_clock::to_time_t(end_3p0p2p0);
        // std::cout << "elapsed_3p0p2p0 time: " << elapsed_seconds_3p0p2p0.count() << "s\n";      
        // auto start_3p0p2p1 = std::chrono::system_clock::now();

        const auto clips = msa.GetClipping(r);
        const auto tags = msa.GetTags(r, 5);

        // auto end_3p0p2p1 = std::chrono::system_clock::now();
        // std::chrono::duration<double> elapsed_seconds_3p0p2p1 = end_3p0p2p1-start_3p0p2p1;
        // std::time_t end_time_3p0p2p1 = std::chrono::system_clock::to_time_t(end_3p0p2p1);
        // std::cout << "elapsed_3p0p2p1 time: " << elapsed_seconds_3p0p2p1.count() << "s\n";      
        // auto start_3p0p2p2 = std::chrono::system_clock::now();

        if (clips.first != -1) cspace_[clips.first].lsc_mean[DnaCharToDnaCode(seq[clips.first])] ++;
        if (clips.second != -1) cspace_[clips.second].rsc_mean[DnaCharToDnaCode(seq[clips.second])] ++;

        // auto end_3p0p2p2 = std::chrono::system_clock::now();
        // std::chrono::duration<double> elapsed_seconds_3p0p2p2 = end_3p0p2p2-start_3p0p2p2;
        // std::time_t end_time_3p0p2p2 = std::chrono::system_clock::to_time_t(end_3p0p2p2);
        // std::cout << "elapsed_3p0p2p2 time: " << elapsed_seconds_3p0p2p2.count() << "s\n";      
        // auto start_3p0p2p3 = std::chrono::system_clock::now();

        for (int pp = s; pp <= e; ++pp) {
          cspace_[pp].bqual_mean[DnaCharToDnaCode(seq[pp])] += float(qual[pp] - 33) ;
          cspace_[pp].strand_mean[DnaCharToDnaCode(seq[pp])] += strand;
          cspace_[pp].mqual_mean[ DnaCharToDnaCode(seq[pp]) ] += r.MapQuality();
          for (size_t ii = 0; ii < Col::TAG_SIZE; ++ii) {
            cspace_[pp].tag_mean[ DnaCharToDnaCode(seq[pp]) ][ii] += tags[ii];
          }
        }

        // auto end_3p0p2p3 = std::chrono::system_clock::now();
        // std::chrono::duration<double> elapsed_seconds_3p0p2p3 = end_3p0p2p3-start_3p0p2p3;
        // std::time_t end_time_3p0p2p3 = std::chrono::system_clock::to_time_t(end_3p0p2p3);
        // std::cout << "elapsed_3p0p2p3 time: " << elapsed_seconds_3p0p2p3.count() << "s\n";      

      }

      // auto end_3p0p2 = std::chrono::system_clock::now();
      // std::chrono::duration<double> elapsed_seconds_3p0p2 = end_3p0p2-start_3p0p2;
      // std::time_t end_time_3p0p2 = std::chrono::system_clock::to_time_t(end_3p0p2);
      // std::cout << "elapsed_3p0p2 time: " << elapsed_seconds_3p0p2.count() << "s\n";      


    }//end for

    // auto end_3p0 = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_3p0 = end_3p0-start_3p0;
    // std::time_t end_time_3p0 = std::chrono::system_clock::to_time_t(end_3p0);
    // std::cout << "elapsed_3p0 time: " << elapsed_seconds_3p0.count() << "s\n";      
    // auto start_3p1 = std::chrono::system_clock::now();

    if (calculate_qual) {
      for (size_t ii = 0; ii < msa.ncol(); ++ii) {
        auto & col = cspace_[ii];
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

    // auto end_3p1 = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_3p1 = end_3p1-start_3p1;
    // std::time_t end_time_3p1 = std::chrono::system_clock::to_time_t(end_3p1);
    // std::cout << "elapsed_3p1 time: " << elapsed_seconds_3p1.count() << "s\n";      

  }


  template<typename GInv>
  explicit CondensedArray(const std::vector<std::string>& msa, const std::vector<std::string>& bqual, 
                    const std::vector<int>& mqual, const std::vector<int>& strand, 
                    const std::vector<int>& lsc, const std::vector<int>& rsc,
                    const std::vector<std::vector<int>>& tags, 
                    const int& total_cov, const GInv& ginv, const std::string& refgap): 
            nrow_(msa.size()), 
            cspace_(msa[0].size(), Col(msa.size())),
            gapped_ref_str_(refgap),
            mquals_(nrow_),
            strands_(nrow_), 
            lsc_(nrow_),
            rsc_(nrow_),
            tags_(nrow_, std::vector<int>(5))
  {

    _CheckInput(msa); 

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

  void Init() {
    for (size_t i = 0; i < ncol(); ++i) {
      auto& col = cspace_[i];
      for (size_t j = 0; j < nrow(); ++j) {
        col.base_freq_[col.bases()[j]]++;
      }
    }
  }

  void InitWithAlnMetaData() {
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
    return (nrow_);
  }

  decltype(auto) cspace() const {
    return (cspace_);
  }
  template<typename T1>
  friend std::ostream& operator<<(std::ostream&, const CondensedArray<T1>&);

private:

  size_t nrow_;
  ColSpace cspace_;
  std::string gapped_ref_str_;
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

template<typename Base>
std::ostream& operator<<(std::ostream& os, const CondensedArray<Base>& ca) {
  for (const auto& col: ca.cspace_) {
    std::string col_s;
    for (const auto& b : col.base_freq_) {
      col_s += std::to_string(b) + ":";
    }
    col_s += "\n";
    os << col_s;
  } 
  return os;
}

template<typename RefGap, typename GInv>
decltype(auto) CreateCondensedArray(const std::vector<std::string>& msa, const int total_cov, const GInv& ginv, const RefGap& refgap) {
  using CondensedArray = neusomatic::CondensedArray<int>;
  CondensedArray condensed_array(msa, total_cov, ginv, refgap);
  condensed_array.Init();
  return condensed_array;
}


template<typename RefGap, typename GInv>
decltype(auto) CreateCondensedArray(const std::vector<std::string>& msa, const std::vector<std::string>& bqual,  const std::vector<int>& mqual, const std::vector<int>& strand,
                              const std::vector<int>& lscs, const std::vector<int>& rscs,
                              const std::vector<std::vector<int>>& tag,
                              const int total_cov, const GInv& ginv, const RefGap& refgap) {
  using CondensedArray = neusomatic::CondensedArray<int>; 
  CondensedArray condensed_array(msa, bqual, mqual, strand, lscs, rscs, tag, total_cov, ginv, refgap);
  condensed_array.InitWithAlnMetaData();
  return condensed_array;
}


std::string add_qual_col(auto  & data_array, bool is_int=false){
  auto sep = ":";
  int order [5] = { 4, 0, 1, 2, 3 }; 
  std::string ret = "";
  for ( int n=0 ; n<5 ; ++n )
  {
    if (is_int){
      ret += std::to_string(data_array[order[n]]);
    }else{
      ret += std::to_string(int(round(data_array[order[n]])));
    }
    if (n < 4){
      ret+=":";
    }
  }
  return ret;
}


void add_qual_col(auto & count_writer, auto  & data_array, bool is_int, const std::string& end_str ){
  auto sep = ":";
  int order [5] = { 4, 0, 1, 2, 3 };
  if (is_int){
    count_writer << (data_array[order[0]])<<":"<<(data_array[order[1]])<<":"<<(data_array[order[2]])<<":"<<(data_array[order[3]])<<":"<<(data_array[order[4]])<<end_str;
  }else{
    count_writer << int(round(data_array[order[0]]))<<":"<<int(round(data_array[order[1]]))<<":"<<int(round(data_array[order[2]]))<<":"<<int(round(data_array[order[3]]))<<":"<<int(round(data_array[order[4]]))<<end_str;
  }
}


std::string add_tag_col(auto  & data_array, bool is_int=false, int idx=0){
  auto sep = ":";
  int order [5] = { 4, 0, 1, 2, 3 }; 
  std::string ret = "";
  for ( int n=0 ; n<5 ; ++n )
  {
    if (is_int){
      ret += std::to_string(data_array[order[n]][idx]);
    }else{
      ret += std::to_string(int(round(data_array[order[n]][idx])));
    }
    if (n < 4){
      ret+=":";
    }
  }
  return ret;
}


void add_tag_col(auto & count_writer, auto  & data_array, bool is_int, int idx, const std::string& end_str){
    int order [5] = { 4, 0, 1, 2, 3 };
  if (is_int){
    count_writer << (data_array[order[0]][idx])<<":"<<(data_array[order[1]][idx])<<":"<<(data_array[order[2]][idx])<<":"<<(data_array[order[3]][idx])<<":"<<(data_array[order[4]][idx])<<end_str;
  }else{
    count_writer << int(round(data_array[order[0]][idx]))<<":"<<int(round(data_array[order[1]][idx]))<<":"<<int(round(data_array[order[2]][idx]))<<":"<<int(round(data_array[order[3]][idx]))<<":"<<int(round(data_array[order[4]][idx]))<<end_str;
  }
}


std::string add_qual_col(auto  & data_array, bool is_int, const std::string& end_str ){
  auto sep = ":";
  int order [5] = { 4, 0, 1, 2, 3 };
  std::string ret = "";
  if (is_int){
    // std::stringstream result;
    // std::copy(data_array.begin(), data_array.end(), std::ostream_iterator<int>(result, ":"));
    // cout << result<<std::endl;

    ret = std::to_string((data_array[order[0]]))+sep+std::to_string((data_array[order[1]]))+sep+std::to_string((data_array[order[2]]))+sep+std::to_string((data_array[order[3]]))+sep+std::to_string((data_array[order[4]]))+end_str;
  }else{
    ret = std::to_string(int(round(data_array[order[0]])))+sep+std::to_string(int(round(data_array[order[1]])))+sep+std::to_string(int(round(data_array[order[2]])))+sep+std::to_string(int(round(data_array[order[3]])))+sep+std::to_string(int(round(data_array[order[4]])))+end_str;
  }
  return ret;
}

std::string add_tag_col(auto  & data_array, bool is_int, int idx, const std::string& end_str){
  auto sep = ":";
  int order [5] = { 4, 0, 1, 2, 3 };
  std::string ret = "";
  if (is_int){
    ret = std::to_string((data_array[order[0]][idx]))+sep+std::to_string((data_array[order[1]][idx]))+sep+std::to_string((data_array[order[2]][idx]))+sep+std::to_string((data_array[order[3]][idx]))+sep+std::to_string((data_array[order[4]][idx]))+end_str;
  }else{
    ret = std::to_string(int(round(data_array[order[0]][idx])))+sep+std::to_string(int(round(data_array[order[1]][idx])))+sep+std::to_string(int(round(data_array[order[2]][idx])))+sep+std::to_string(int(round(data_array[order[3]][idx])))+sep+std::to_string(int(round(data_array[order[4]][idx])))+end_str;
  }
  return ret;
}


// void Compute_ALTs(auto & alt_counts, auto & base_freq_, auto  & dp, auto ref_base, bool report_all_alleles){
//   ref_base = std::toupper(ref_base);
//   auto ref_code = DnaCharToDnaCode(ref_base);

//   if (ref_base == 'N') {
//     ref_base = '-';
//   }

//   base_freq_.erase(base_freq_.begin() + 5);

//   std::vector<int> pileup_counts(base_freq_.size());
//   int total_count=0;
//   for (int base = 0; base < (int) base_freq_.size(); ++base) {
//     pileup_counts[base] = base_freq_[base];
//     total_count+=base_freq_[base];
//   }
  
//   if (total_count==0) {
//     return;
//   }
//   auto ref_count = base_freq_[ref_code];
//   auto var_code = ref_code; 
//   int var_count = 0;
//   dp = ref_count;
//   if (report_all_alleles and ref_base != '-'){
//     for (int row = 0;  row < base_freq_.size(); ++row) {
//       auto alt_cnt = base_freq_[row];
//       if (( row != ref_code) and (alt_cnt > 0)){
//         auto af = alt_cnt/float(alt_cnt+ref_count);
//         if ((alt_cnt >= ref_count) or ((row == 4 and  af > del_min_af ) or
//                                         (row != 4 and ref_base != '-' and af > snp_min_af ) or
//                                         (ref_base =='-' and af > ins_min_af))){
//           alt_counts.insert(std::pair<int, int>(row, alt_cnt));
//           dp += alt_cnt;
//         }
//       }
//     }
//   }else{
//     int major = -1;
//     int major_count = 0;
//     int minor = -1;
//     int minor_count = 0;
//     int minor2 = -1;
//     int minor2_count = 0;

//     for (int row = 0;  row < base_freq_.size(); ++row) {
//       if (base_freq_[row] > major_count) {
//         minor2 = minor;
//         minor2_count = minor_count;
//         minor_count = major_count;
//         minor = major;
//         major_count = base_freq_[row];
//         major = row;
//       } else if (base_freq_[row] > minor_count) {
//         minor2 = minor;
//         minor2_count = minor_count;
//         minor_count = base_freq_[row];
//         minor = row;
//       } else if (base_freq_[row] > minor2_count) {
//         minor2_count = base_freq_[row];
//         minor2 = row;
//       }
//     }

//     if (minor != -1 and major != -1){
//       if (minor2 != -1 and ref_code == major and minor == 4 and ref_code != 4 ){
//         if (minor2_count>0.5*minor_count){
//           minor = minor2;
//           minor_count = minor2_count;
//         }
//       }
//     }
//     auto af = minor_count/float(major_count+minor_count);
//     if (major != ref_code){
//       var_code = major;
//       var_count = major_count;
//     } else if (minor != ref_code and ( (minor == 4 and  af > del_min_af ) or
//                                     (minor != 4 and ref_base != '-' and af > snp_min_af ) or
//                                     (ref_base =='-' and af > ins_min_af))){
//       var_code = minor;
//       var_count = minor_count;
//     }
//     if (var_count > 0) { 
//       alt_counts.insert(std::pair<int, int>(var_code,var_count));
//       dp += var_count;
//     }
//   }
// }


}// end neusomatic

#endif
