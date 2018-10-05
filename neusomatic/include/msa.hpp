#ifndef ALLAMERICAN_MSA_HPP
#define ALLAMERICAN_MSA_HPP

#include <vector>
#include <string>
#include <limits>

#include "bam_variant.hpp"
namespace neusomatic{

class Gap {
  /*
   * Represent a gap in the msa.
   */
  int32_t pos_; 
  int32_t len_;

public:
  Gap(int32_t p, int32_t l): pos_(p), len_(l) {};     
  Gap(): pos_(-1), len_(0) {};
  bool IsNull() const { return pos_ == -1 && len_ == 0; }
  void Set(const int32_t p, const int32_t l) {pos_ = p; len_ = l;}
  int32_t len() const {return len_;}
  int32_t pos() const {return pos_;}
};

inline std::ostream& operator<<(std::ostream& os, const Gap& gap)
{
  os << gap.pos() + 1 <<": "<<gap.len();
  return os;
}

template<typename Base, typename Gap, typename GInv>
class GappedSeq {
  /*
   * This data struct reprents an list of alternative Base 
   * and Gap. The Base always starts first.
   */
  std::vector<Base> bases_;
  std::vector<Gap> gaps_;
  const GInv ginv_;
  std::string seq_str_;

public:
  GappedSeq(unsigned rlen): bases_(rlen), gaps_(rlen) {}
  
  GappedSeq(const std::string& rseq, const GInv& ginv): bases_(rseq.begin(), rseq.end()), gaps_(rseq.length()), ginv_(ginv) {}

  void SetGap(const size_t i, const int32_t pos, const int32_t len) {
    if (gaps_[i].IsNull()) {
      gaps_[i].Set(pos, len);
    } else {
      if (len > gaps_[i].len()) gaps_[i].Set(pos, len); 
    }
  }

  decltype(auto) gaps() {
    return  (gaps_);
  }

  decltype(auto) ginv() const {return (ginv_);}
  int32_t left() const {return ginv_.left();}
  int32_t right() const {return ginv_.right();}
  decltype(auto) bases() const {return (bases_);}
  decltype(auto) gaps() const {return (gaps_);}

  std::string to_string() const {
    std::string result;
    for (size_t i = 0; i < bases_.size(); ++i) {
      if (i != 0) {
        for (int32_t j = 0; j < gaps_[i].len(); ++j) {
          result += '-';
        }
      }
      result += bases_[i];
    }
    return result;
  }
};

template<typename Base, typename Gap, typename GInv>
inline std::ostream& operator<<(std::ostream& os, const GappedSeq<Base, Gap, GInv> gs) {
  os << gs.ginv() <<":\t";
  for (size_t i = 0; i < gs.bases().size(); ++i) {
    if (i != 0) {
      os << std::string(gs.gaps()[i].len(), '-');
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
  std::vector<unsigned char> bquals_;

  std::vector<std::string>  gapstrs_; // the first element in this vector is not in use. 
  std::string bases_str_;
  const unsigned char missing_chr_;
  const unsigned char gapchar_;
public:

  ReadSeqGap(const GappedSeq<char, Gap, typename Variant::TInv>& refgap, const std::vector<Variant>& vars, 
             const int record_begin, const int record_end, const unsigned char gapchar = '-', const unsigned char missing_chr = '~'):
              ReadSeqGap(refgap, vars, std::vector<unsigned char>(1,0), false, record_begin, record_end, gapchar, missing_chr){}


  ReadSeqGap(const GappedSeq<char, Gap, typename Variant::TInv>& refgap, const std::vector<Variant>& vars, 
             const std::vector<unsigned char>& bquals, bool calculate_qual_stat,
             const int record_begin, const int record_end, const unsigned char gapchar = '-', const unsigned char missing_chr = '~'):
              bases_(refgap.bases()), gapstrs_(bases_.size()), gapchar_(gapchar), missing_chr_(missing_chr) {
    // starts from the reference bases and then fill in mutated bases and gaps.
    
    int msa_len=0;
    for (size_t i = 0; i < refgap.gaps().size(); ++ i) {
      const auto& gap = refgap.gaps()[i];
      gapstrs_[i] = std::string(gap.len(), gapchar_);
      msa_len += 1+gap.len();
    }
    for (auto const& v : vars) {

      if (v.Type() == neusomatic::bio::VarType::SUB) {
        if (!neusomatic::bio::IsContainedIn(v.ginv(), refgap.ginv())) continue; 
        for (size_t i = 0; i < v.length(); ++i) {
          bases_[v.left() - refgap.left() + i] = v.allele()[i];
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
        auto const& gap = refgap.gaps()[v.left() - refgap.left()]; 
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
    if (calculate_qual_stat){
      bquals_.assign(msa_len,33); 
      int j=0;
      int c=0;
      for (size_t i = 0; i < bases_.size(); ++i) {
        if (record_begin > (refgap.left()+i) ){
          c+=(gapstrs_[i].size()+1);
          continue;
        }
        if (record_end <= (refgap.left()+i) ){
          break;
        }
        if (i != 0) {
          for (size_t k = 0; k < gapstrs_[i].size(); ++k) {
            if (gapstrs_[i][k]!=gapchar_){
              assert((33<=bquals[j])&(bquals[j]<=74));
              bquals_[c+k]=bquals[j++];
            }else{
              if (j>0){
                bquals_[c+k]=bquals[j-1];
              }else{
                bquals_[c+k]=bquals[0];
              }
            }
          }
        }
        c+=gapstrs_[i].size();
        if ((bases_[i]!=missing_chr_) & (bases_[i]!=gapchar_)) {
          assert((33<=bquals[j])&(bquals[j]<=74));
          bquals_[c]=bquals[j++];
        }else if ((bases_[i]==gapchar_)) {
            if (j>0){
              bquals_[c]=bquals[j-1];
            }else{
              bquals_[c]=bquals[0];
            }
        }
        c+=1;
      }
    }
  };

  std::string get_bqual() const {
    std::string bquals_str(bquals_.begin(), bquals_.end());
    return bquals_str;
  }
  std::string get_lsc(const GappedSeq<char, Gap, typename Variant::TInv>& refgap, int pos_lsc) const {
    std::vector<unsigned char> lsc(bquals_.size(),'0');
    if (refgap.left() <= pos_lsc && pos_lsc<refgap.right()){
      int c=0;
      for (size_t i = 0; i < bases_.size(); ++i) {
        if (i !=0) {
          c += gapstrs_[i].size();
        } 
        if ((refgap.left()+i)==pos_lsc){
          lsc[c]='1';  
          break;
        }
        c += 1;
      }
    }
    std::string lsc_str(lsc.begin(), lsc.end());
    return lsc_str;
  }
  std::string get_rsc(const GappedSeq<char, Gap, typename Variant::TInv>& refgap, int pos_rsc) const {
    std::vector<unsigned char> rsc(bquals_.size(),'0');
    if (refgap.left() <= pos_rsc && pos_rsc<refgap.right()){
      int c=0;
      for (size_t i = 0; i < bases_.size(); ++i) {
        if (i !=0) {
          c += gapstrs_[i].size();
        } 
        if ((refgap.left()+i)==pos_rsc){
          rsc[c]='1';  
        }
        c += 1;
      }
    }
    std::string rsc_str(rsc.begin(), rsc.end());
    return rsc_str;
  }

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

template<typename BamRecord, typename Variant>
class MSABuilder {
  using GInv = typename Variant::TInv;
  GInv ginv_;
  const std::vector<BamRecord>& bam_records_;
  //std::vector<size_t> rids_;
  GappedSeq<char, Gap, GInv> ref_gaps_; 

public:
  MSABuilder() = default; 
  using ContigGaps = GappedSeq<char, Gap, GInv>; 

  size_t size() {
    return bam_records_.size();
  }

  MSABuilder(const GInv& ginv, const std::vector<BamRecord>& bams, const std::string& refstr):  
      ginv_(ginv), bam_records_(bams), ref_gaps_(refstr, ginv) {
        
    // This loop setup the reference in gapped space. 
    for (auto const& r : bam_records_) {
      auto vars = neusomatic::GetIndels<BamRecord, int>(r);
      for (auto const& v : vars) {
        if (neusomatic::bio::IsContainedIn(v.ginv(), ginv_)) {
          if (v.Type() == neusomatic::bio::VarType::INS) {
            int idx = v.left() - ginv_.left();
            ref_gaps_.SetGap(idx, v.left(), v.allele().length());
          }
        }
      }
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

  auto GetMSAwithQual() const { 
    std::vector<std::string> result(bam_records_.size());
    std::vector<std::string> bquals(bam_records_.size());
    std::vector<std::string> lscs(bam_records_.size());
    std::vector<std::string> rscs(bam_records_.size());
    std::vector<int> mquals(bam_records_.size());
    std::vector<bool> strands(bam_records_.size());
    std::vector<std::vector<int>> tags(bam_records_.size(),std::vector<int>(5,0));
    for (size_t i = 0; i < bam_records_.size(); ++i) {
      auto const& r = bam_records_[i];
      auto const& cigar = r.GetCigar();        
      auto vars = neusomatic::GetVars<BamRecord, int>(r);
      auto read_bqual_string = r.Qualities();
      auto AlignmentEndPosition = r.AlignmentEndPosition();
      auto AlignmentPosition = r.AlignmentPosition();
      if (cigar.front().Type()=='H'){
        AlignmentPosition-=cigar.front().Length();
      }
      if (cigar.back().Type()=='H'){
        AlignmentEndPosition+=cigar.back().Length();
      }

      int start_qual_pos=AlignmentPosition;
      int end_qual_pos=AlignmentEndPosition;
      if (r.Position()<ref_gaps_.left()){
        if (r.PositionEnd()<ref_gaps_.right()) {
          auto const& rev_cigar = r.GetReverseCigar();  
          int e=0;
          int p=r.PositionEnd();
          for (auto c = rev_cigar.begin(); c != rev_cigar.end(); ++c) {
            switch (c->Type()) {
              case 'M':
                e+=c->Length();
                p-=c->Length();
                break;
              case 'I':
                e+=c->Length();
                break;
              case 'D':
                p-=c->Length();
                break;
              default:
                break;
            }
            if (p <= ref_gaps_.left()){
              if (!(c->Type()=='D')){
                e-=(ref_gaps_.left()-p);                
              }
              break;
            }
          }
          start_qual_pos=AlignmentEndPosition-e;
        }else{
          int e=0;
          int p=r.Position();
          for (auto c = cigar.begin(); c != cigar.end(); ++c) {
            switch (c->Type()) {
              case 'M':
                e+=c->Length();
                p+=c->Length();
                break;
              case 'I':
                e+=c->Length();
                break;
              case 'D':
                p+=c->Length();
                break;
              default:
                break;
            }
            if (p > ref_gaps_.left()){
              if (!(c->Type()=='D')){
                e-=(p-ref_gaps_.left());
              }
              break;
            }
          }
          start_qual_pos=AlignmentPosition+e;
          start_qual_pos=e;
        }
      }
      std::vector<unsigned char> read_bquals_(AlignmentEndPosition-start_qual_pos,74);
      int c=0;
      if (! read_bqual_string.empty()){
        if (read_bqual_string[0]!=' '){
          for (auto j = start_qual_pos; j < AlignmentEndPosition; ++j) {
            read_bquals_[c++]=read_bqual_string[j];
          }
        }
      }
      std::string read_bquals_str(read_bquals_.begin(), read_bquals_.end());   
      ReadSeqGap<char, neusomatic::bio::Variant<std::string, int>> row(ref_gaps_, vars, read_bquals_, true, r.Position(), r.PositionEnd());
      result[i] = row.to_string();
      bquals[i] = row.get_bqual();

      auto front_cigar=cigar.front().Type();
      auto end_cigar=cigar.back().Type();
      int pos_lsc = ((front_cigar=='H')||(front_cigar=='S')) ? r.Position() : -1;
      int pos_rsc = ((end_cigar=='H')||(end_cigar=='S')) ? r.PositionEnd()-1 : -1;

      lscs[i] = row.get_lsc(ref_gaps_, pos_lsc);
      rscs[i] = row.get_rsc(ref_gaps_, pos_rsc);
      mquals[i] = r.MapQuality();
      strands[i] = !r.ReverseFlag();

      auto length=r.Sequence().size();
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
      tags[i][0]=std::max(std::min(100 - (int) std::rint(((float) nm)/(0.001+length)*100),100),0);
      tags[i][1]=std::max(std::min((int) std::rint(((float) as)/(0.001+length)*100),100),0);
      tags[i][2]=std::max(std::min(100 - (int) std::rint(((float) xs)/(0.001+length)*100),100),0); 
      tags[i][3]=r.ProperPair() ? 100 : 0; //BamRecord is mapped in proper pair
      tags[i][4]=(r.NumClip()>0) ? 100: 0; //BamRecord has soft/hard clips

    }
    return std::tuple< std::vector<std::string>, std::vector<std::string>, std::vector<int>, 
                       std::vector<bool>, std::vector<std::string>, std::vector<std::string>,
                       std::vector<std::vector<int>> > (result,bquals,mquals,strands,lscs,rscs,tags);
  } 


  decltype(auto) ref_gaps() const {
    return (ref_gaps_);
  };
           
};

} // namespace neusomatic
#endif /* ALLAMERICAN_HPP */
