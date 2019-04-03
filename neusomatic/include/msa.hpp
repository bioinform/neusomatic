#ifndef ALLAMERICAN_MSA_HPP
#define ALLAMERICAN_MSA_HPP

#include <vector>
#include <string>
#include <limits>

#include "bam_variant.hpp"
namespace neusomatic{


template<typename Base, typename GInv>
class GappedSeq {
  /*
   * This data struct reprents an list of alternative Base 
   * and Gap. The Base always starts first.
   */
  std::vector<Base> bases_;
  std::vector<int> gaps_;
  const GInv ginv_;
  std::string seq_str_;

public:
  GappedSeq(unsigned rlen): bases_(rlen), gaps_(rlen) {}
  
  GappedSeq(const std::string& rseq, const GInv& ginv): bases_(rseq.begin(), rseq.end()), gaps_(rseq.length()), ginv_(ginv) {
    /*
     * Create a GappedSeq from a ref string (no gaps) and its corresponding genomic interval
     * */
    if (rseq.size() + ginv.left()  != ginv.right()) {
      throw std::runtime_error("seq length does not match the interval length");
    }
  }
  void SetGap(const size_t i, const int32_t len) {
    if (len > gaps_[i]) {
      gaps_[i] = len; 
    }
  }

  decltype(auto) gaps() {
    return  (gaps_);
  }

  decltype(auto) ginv() const {return (ginv_);}
  int32_t left() const {return ginv_.left();}
  int32_t right() const {return ginv_.right();}
  std::vector<Base> bases() const {return bases_;}
  std::vector<int> gaps() const {return gaps_;}

  std::string to_string() const {
    std::string result;
    for (size_t i = 0; i < bases_.size(); ++i) {
      if (i != 0) {
        for (int32_t j = 0; j < gaps_[i]; ++j) {
          result += '-';
        }
      }
      result += bases_[i];
    }
    return result;
  }

  int length() const {
    int res = 0;
    for (size_t i = 0; i < bases_.size(); ++i) {
      if (i != 0) {
         res += gaps_[i];
      }
    }
    res += bases_.size();
    return res;
  }
};

template<typename Base, typename GInv>
inline std::ostream& operator<<(std::ostream& os, const GappedSeq<Base, GInv> gs) {
  os << gs.ginv() <<":\t";
  for (size_t i = 0; i < gs.bases().size(); ++i) {
    if (i != 0) {
      os << std::string(gs.gaps()[i], '-');
    }
    os << gs.bases()[i]; 
  }
  return os;
}

template<typename Base, typename Variant>
class ReadSeqGap {
  /*
   * Represent a gapped read in the msa. 
   * Alternate Base with gaps. The Base always comes first. 
   */
  std::vector<Base> bases_; 
  std::vector<std::string>  gapstrs_; // the first element in this vector is not in use. 
  std::string bases_str_;
  static const unsigned char missing_chr_ = '~';
  static const unsigned char gapchar_ = '-';
public:

  ReadSeqGap(const GappedSeq<char, typename Variant::TInv>& refgap, 
             const std::vector<Variant>& vars, 
             const int record_begin, 
             const int record_end) {

    // starts from the reference bases and then fill in mutated bases and gaps.
    
    int msa_len=0;
    for (size_t i = 0; i < refgap.gaps().size(); ++ i) {
      const auto& gap = refgap.gaps()[i];
      gapstrs_[i] = std::string(gap, gapchar_);
      msa_len += 1+gap;
    }
    for (auto const& v : vars) {

      if (v.Type() == neusomatic::bio::VarType::SUB) {
        if (!neusomatic::bio::IsOverlapped(v.ginv(), refgap.ginv())) continue; 
        int left = std::max(v.left(), refgap.left());
        int right = std::min(v.right(), refgap.right());
        for (size_t i = left; i < right; ++i) {
          bases_[i - refgap.left()] = v.allele()[i - v.left()];
        }
      } else if (v.Type() == neusomatic::bio::VarType::DEL) {
        if (!neusomatic::bio::IsOverlapped(v.ginv(), refgap.ginv())) continue; 
        int left = std::max(v.left(), refgap.left());
        int right = std::min(v.right(), refgap.right());
        for (size_t i = left; i < right; ++i) {
          bases_[i - refgap.left()] = gapchar_;
        }
      } else if (v.Type() == neusomatic::bio::VarType::INS) {
        if (!neusomatic::bio::IsContainedIn(v.ginv(), refgap.ginv())) continue; 
        for (size_t i = 0; i < v.allele().length(); ++i) {
          gapstrs_[v.left() - refgap.left()][i] = v.allele()[i];
        }
      }
    }
    for (size_t i = refgap.left(); i < record_begin; ++i) {
      bases_[i - refgap.left()] = missing_chr_;
    }

    for (size_t i = record_end; i < refgap.right(); ++i) {
      bases_[i - refgap.left()] = missing_chr_;
    }
  };

  std::string to_string() const {
    std::string result;
    for (size_t i = 0; i < bases_.size(); ++i) {
      if (i != 0) {
        result += bases_[i] == missing_chr_ ? std::string(gapstrs_[i].size(), missing_chr_) : gapstrs_[i];
      } 
      result += bases_[i];
    }
    if (result[0] == missing_chr_) {
      int i = 0;
      while (result[++i] == missing_chr_); 
      for (; result[i] == gapchar_; ++i) result[i] = missing_chr_;
    }
    return result;
  }

};

template<typename BamRecord, typename GInv>
class MSABuilder {
  const std::vector<BamRecord>& bam_records_;
  //std::vector<size_t> rids_;
  GappedSeq<char, GInv> ref_gaps_; 
  std::vector<std::string>  gapstrs_; // the first element in this vector is not in use. 
  int ncol_;
  std::string gapped_ref_str_;
  std::vector<int> gap_positions_;
  static const unsigned char missing_chr_ = '~';
  static const unsigned char gapchar_ = '-';

public:
  MSABuilder() = delete; 
  using ContigGaps = GappedSeq<char,GInv>; 

  decltype(auto) bam_records() const {
    return (bam_records_);
  }

  decltype(auto) ref_gaps() const {
    return (ref_gaps_);
  }

  size_t size() const {
    return bam_records_.size();
  }

  const int ncol() const {
    return ncol_;
  }

  const std::string gapped_ref_str() const {
    return gapped_ref_str_;
  }

  const int GapPosition(const int p) const {
    return gap_positions_[p];
  }

  static void BuildRefGap(const std::vector<BamRecord>& bams, GappedSeq<char, GInv>& refgaps) {
    // This loop setup the reference in gapped space. 
    for (auto const& r : bams) {
      auto vars = neusomatic::GetIndels<BamRecord, int>(r);
      for (auto const& v : vars) {
        if (neusomatic::bio::IsContainedIn(v.ginv(), refgaps.ginv())) {
          if (v.Type() == neusomatic::bio::VarType::INS) {
            int idx = v.left() - refgaps.ginv().left();
            refgaps.SetGap(idx, v.allele().length());
          }
        }
      }
    }
  }


  auto GetGappedSeqAndQual(const BamRecord& r) const {

    std::string msa_bases(ncol(), missing_chr_);
    std::string msa_bquals(ncol(), 33);

    auto const& cigar = r.GetCigar(); 
    auto const& seq = r.Sequence();
    auto const& qual = r.Qualities();
    int32_t ref_pos = r.Position();
    int32_t read_pos = r.AlignmentPosition();
    auto gapl = GapPosition(std::max(r.Position() , ref_gaps_.left()) - ref_gaps_.left());
    auto gapr = GapPosition(std::min(r.PositionEnd(), ref_gaps_.right()) - 1 - ref_gaps_.left());
    if (gapl > gapr) {
      throw std::runtime_error("Invalid read " + r.Qname() + "\n");
    }
    std::fill(msa_bases.begin() + gapl, msa_bases.begin() + gapr ,'-');

    if (cigar.begin()->Type() == 'H') {
      read_pos -= cigar.begin()->Length();
    }
    for (auto c = cigar.begin(); c != cigar.end(); ++c) {
      if (c->Type() == 'M') {
        const int l = std::max(ref_pos, ref_gaps_.left());
        const int r = std::min(int32_t(ref_pos + c->Length()), ref_gaps_.right());
        if (l < r) {
          const int s = read_pos + l - ref_pos;       
          auto ss = seq.substr(s, l - r);
          auto qq = qual.substr(s, l -r); 
          
          for (int pp = l; pp < r; ++pp) {
            int gapp = GapPosition(pp - ref_gaps_.left());
            msa_bases[gapp] = ss[pp - l]; 
            msa_bquals[gapp] = qq[pp - l];
            if (pp != r - 1) {
              for (auto gg = GapPosition(pp - ref_gaps_.left()) + 1 ; gg < GapPosition(pp - ref_gaps_.left() + 1); gg++){
                //msa_bases[gg] = gapchar_;
                msa_bquals[gg] = qq[pp - l];
              }
            }
          }
        }

        ref_pos += c->Length();
        read_pos += c->Length();
      }
      if (c->Type() == 'I') {
        if (ref_pos > ref_gaps_.left() && ref_pos < ref_gaps_.right()) {
          const auto gapp = GapPosition(ref_pos - 1 - ref_gaps_.left());
          const auto gapp_end = GapPosition(ref_pos - ref_gaps_.left());
          auto s = seq.substr(read_pos, c->Length());
          auto q = qual.substr(read_pos, c->Length()); 
          int pp = 0;
          for (; pp < s.size(); ++pp) {
            msa_bases[gapp + pp + 1] = s[pp];
            msa_bquals[gapp + pp + 1] = q[pp];
          }
          for (; gapp + pp + 1 < gapp_end ; ++pp) {
            msa_bquals[gapp + pp + 1] = q.back();
          }
        }
        read_pos += c->Length();
      }
      if (c->Type() == 'D') {
        const int l = std::max(ref_pos, ref_gaps_.left());
        const int r = std::min(int32_t(ref_pos + c->Length()), ref_gaps_.right());
        if (l < r) {
          const unsigned char q = qual[read_pos];
          for (int pp = l; pp < r; ++pp) {
            int gapp = GapPosition(pp - ref_gaps_.left());
            msa_bases[gapp] = gapchar_; 
            msa_bquals[gapp] = q;
            for (auto gg = GapPosition(pp - ref_gaps_.left() - 1) + 1 ; gg < GapPosition(pp - ref_gaps_.left()); gg++){
              //msa_bases[gg] = gapchar_;
              msa_bquals[gg] = q;
            }
          }
          if (r < ref_gaps_.right()) {
            for (auto gg = GapPosition(r - ref_gaps_.left() - 1) + 1 ; gg < GapPosition(r - ref_gaps_.left()); gg++){
              msa_bquals[gg] = q;
            }
          }
        }
        ref_pos += c->Length();
      }
    }
    
    return std::make_pair(msa_bases, msa_bquals);
  }
  
  MSABuilder(const GInv& ginv, const std::vector<BamRecord>& bams, const std::string& refstr):  
    bam_records_(bams), ref_gaps_(refstr, ginv) {
    BuildRefGap(bams, ref_gaps_);      
    ncol_ = ref_gaps_.length();
    gapped_ref_str_ = ref_gaps_.to_string();
    gapstrs_.resize(ref_gaps_.bases().size());
    for (size_t i = 0; i < ref_gaps_.gaps().size(); ++ i) {
      const auto& gap = ref_gaps_.gaps()[i];
      gapstrs_[i] = std::string(gap, gapchar_);
    }

    gap_positions_.resize(ref_gaps_.bases().size());
    int c=0;
    for (size_t i = 0; i < ref_gaps_.bases().size(); ++i) {
      if (i !=0) {
        c += ref_gaps_.gaps()[i];
      } 
      gap_positions_[i] = c;
      c += 1;
    }
  }

  std::vector<std::string> GetMSA() const { 
    std::vector<std::string> result(bam_records_.size());
    for (size_t i = 0; i < bam_records_.size(); ++i) {
      auto const& r = bam_records_[i];
      auto vars = neusomatic::GetVars<BamRecord, int>(r);
      ReadSeqGap<char, neusomatic::bio::Variant<std::string, int>> row(ref_gaps_, vars, r.Position(), r.PositionEnd());
      result[i] = row.to_string();
    }
    return result;
  } 

  auto GetClipping(const BamRecord& r) const {
    int lsc = -1, rsc = -1;
    auto const& cigar = r.GetCigar();        
    auto front_cigar=cigar.front().Type();
    auto end_cigar=cigar.back().Type();
    int pos_lsc = ((front_cigar=='H')||(front_cigar=='S')) ? r.Position() : -1;
    int pos_rsc = ((end_cigar=='H')||(end_cigar=='S')) ? r.PositionEnd()-1 : -1;

    if (ref_gaps_.left() <= pos_lsc && pos_lsc<ref_gaps_.right()){
      lsc = GapPosition(pos_lsc - ref_gaps_.left());
    }
    if (ref_gaps_.left() <= pos_rsc && pos_rsc<ref_gaps_.right()) {
      rsc = GapPosition(pos_rsc - ref_gaps_.left());
    }   
    return std::make_pair(lsc, rsc);
  }

  auto GetTags(const BamRecord& r, int tag_size) const {
    std::vector<int> tags(tag_size);
    const auto length= r.Sequence().size();
    int32_t nm=0;
    uint8_t* p1 = bam_aux_get(r.shared_pointer().get(), "NM");
    if (p1){
      nm = bam_aux2i(p1);
    }
    int32_t as=length;
    uint8_t* p2 = bam_aux_get(r.shared_pointer().get(), "AS");
    if (p2){
      as = bam_aux2i(p2);
    }
    int32_t xs=0;
    uint8_t* p3 = bam_aux_get(r.shared_pointer().get(), "XS");
    if (p3){
      xs = bam_aux2i(p3);
    }
    tags[0]=std::max(std::min(100 - (int) std::rint(((float) nm)/(0.001+length)*100),100),0);
    tags[1]=std::max(std::min((int) std::rint(((float) as)/(0.001+length)*100),100),0);
    tags[2]=std::max(std::min(100 - (int) std::rint(((float) xs)/(0.001+length)*100),100),0); 
    tags[3]=r.ProperPair() ? 100 : 0; //bamrecord is mapped in proper pair
    tags[4]=(r.NumClip()>0) ? 100: 0; //bamrecord has soft/hard clips
    return tags;
  }

  auto GetMSAwithQual() const { 


    std::vector<std::string> result(bam_records_.size());
    std::vector<std::string> bquals(bam_records_.size());
    std::vector<int> lscs(bam_records_.size(), -1);
    std::vector<int> rscs(bam_records_.size(), -1);
    std::vector<int> mquals(bam_records_.size());
    std::vector<int> strands(bam_records_.size());
    std::vector<std::vector<int>> tags(bam_records_.size(),std::vector<int>(5,0));
    for (size_t i = 0; i < bam_records_.size(); ++i) {
      auto const& r = bam_records_[i];
      auto const& cigar = r.GetCigar();        
      const auto seq_qual = GetGappedSeqAndQual(r);
      result[i] = seq_qual.first;
      bquals[i] = seq_qual.second;
      
      const auto clips = GetClipping(r);
      lscs[i] = clips.first; 
      rscs[i] = clips.second;
      mquals[i] = r.MapQuality();
      strands[i] = (int) !r.ReverseFlag();
      const auto tag = GetTags(r, 5);
      tags[i] = tag;
    }
    return std::tuple< std::vector<std::string>, std::vector<std::string>, std::vector<int>, 
                       std::vector<int>, std::vector<int>, std::vector<int>,
                       std::vector<std::vector<int>> > (result,bquals,mquals,strands,lscs,rscs,tags);
  } 


};

} // namespace neusomatic
#endif /* ALLAMERICAN_HPP */
