#ifndef VARIANT_HPP
#define VARIANT_HPP

#include<map>
#include<memory>
#include "SeqanUtils.h"
#include "Interval.hpp"
//#include<cpputil/bio/RefSeqs.h>

/*
 * Standard: any positions(or interval) of a vcf variant should 
 * refer to reference positions.
 */

namespace neusomatic { namespace bio {


enum class VarType{
  DEL = 0,
  INS,
  SUB   
};

inline std::ostream& operator<<(std::ostream &os, const VarType& var)
{
  switch(var)
  {
    case VarType::DEL:
      os<<"D";
      break;
    case VarType::INS:
      os<<"I";
      break;
    case VarType::SUB:
      os<<"S";
      break;
  }
  return os;
} 


/////////////////////
////////////////////



template<typename T1, typename T2=std::string>
class Variant{
public:
  typedef T1 TSeq; 
  typedef T2 TContig;
  typedef neusomatic::bio::GenomicInterval<TContig> TInv;


  /*
   * Deletion type variant ctor
   * 6 paramters
   */
  //Variant(): del_len_(0), ins_len_(0){
    //std::cerr<<"in default ctor ginv is: "<<ginv_<<std::endl;
  //}

  /*
   * variant ctor
   * 8 parameters
   */
  Variant() = default;
  Variant(TContig ctg,
      StrandType strand,
      int begin,
      unsigned int len,
      TSeq seq)
  {
    //createVariant_(ctg, strand, begin, len, vt, dl, il, seq);
    ginv_ = TInv(ctg, strand, begin, begin + len);
    allele_ = seq;
  }

  Variant(TContig ctg,
      int begin,
      unsigned int len,
      TSeq seq) : Variant(ctg, StrandType::UNKNOWN, begin, len, seq)
  {}
/*
 * bunch of getters and setters
 */ 

  bool Valid() const { return ginv_.Valid();}
  decltype(auto) id(){return (id_);}
  
  decltype(auto) left() const {return (ginv_.left());}

  decltype(auto) right() const {return (ginv_.right());}
  
  decltype(auto) contig() const {return (ginv_.contig());}

  decltype(auto) strand() const {return (ginv_.strand());}

  size_t length() const {return right() - left();}

  decltype(auto) ginv() const {return (ginv_);}

  decltype(auto) allele() const {return (allele_);}

  VarType Type() const {
    if (!Valid()) {
      throw std::runtime_error("Uninitilized Variant!");
    }
    if (ginv_.length() > allele_.length()) {
      return VarType::DEL; 
    } else if (ginv_.length() < allele_.length()) {
      return VarType::INS;
    } else {
      return VarType::SUB;
    }
  }

//  bool IsSnp() const {
//    return length() == seqan::length(allele_);
//  }


/*
 *
 */
private:
  int id_ = INT_MIN;
  TInv ginv_;

  TSeq allele_;
  //bool final_;
};

template<typename TSeq, typename TContig>
inline bool operator==(const Variant<TSeq, TContig>& lhs, const Variant<TSeq, TContig>& rhs)
{
  if (lhs.ginv() != rhs.ginv()) return false;
  if (lhs.allele() != rhs.allele()) return false;
  return true;
}

template<typename TSeq, typename TContig>
inline bool operator!=(const Variant<TSeq, TContig>& lhs, const Variant<TSeq, TContig>& rhs)
{
  return !(lhs == rhs);
}

template<typename TSeq, typename TContig>
inline bool operator<(const Variant<TSeq, TContig>& lhs, const Variant<TSeq, TContig>& rhs)
{
  return lhs.ginv() < rhs.ginv();
}


template<typename T1, typename T2>
inline std::ostream& operator<<(std::ostream& os, const Variant<T1,T2>& var)
{
  os << var.ginv() 
     <<":"<<neusomatic::bio::to_string(var.allele());
  return os;
}

/*
 *
 *  Class: VariantGroup
 *
 *
 */

template<typename TVariant>
class VariantGroup{
/*
 * TGaps should be in compliance with seqan Gaps type
 */
public:
  using Variant = TVariant;
  typedef typename TVariant::TSeq TSeq; 
  typedef typename TVariant::TContig TContig;
  typedef typename TVariant::TInv TInv;
//TContig chromosome
  VariantGroup(TContig ctg, StrandType strand, const neusomatic::bio::RefSeqs& ref_seqs):
      ginv_(ctg, strand), ref_seqs_(&ref_seqs){}

  VariantGroup(TContig ctg, const neusomatic::bio::RefSeqs& ref_seqs):
    ginv_(ctg, StrandType::UNKNOWN), ref_seqs_(&ref_seqs){}

  typename TVariant::TSeq GetRefAllele() const
  {

    const auto& ref_seq = (*ref_seqs_)[contig()];

    typename TVariant::TSeq ref_allele;
    int vgleft = ginv_.left();
    int vgright = ginv_.right();

    //if (IsSnp()) ++vgleft;

    for(int i = vgleft; i < vgright; ++i){
      ref_allele += ref_seq[i];
    }
    return ref_allele;
  }

  decltype(auto) GetAltAlleles() const {
    std::map<size_t, typename TVariant::TSeq> var_seq_map;
    for(auto itr = cid_vid_.begin(); itr != cid_vid_.end(); ++itr) {
      auto search  = var_seq_map.find(itr->second);
      if (search == var_seq_map.end()) {
        var_seq_map.emplace(itr->second, GetAltAlleles_(itr->second));
      } else {
        //std::cerr<<"varaints have existed in other consensus. skipping..."<<std::endl;
      }
    }
    return var_seq_map;
  }


  //bool IsSnp() const;

  int left() const { return ginv_.left();}
  int right() const { return ginv_.right();}
  size_t length() const {return ginv_.length();}

  decltype(auto) variants() const { return (variants_);}

  void AddVariant(const TVariant& new_var, const size_t sid)
  {

    //first merge variant if need
    TVariant var_to_add(new_var);
    auto search = cid_vid_.find(sid);
    if (search != cid_vid_.end()) {
      var_to_add = MergeVariants_(variants_[search->second], new_var);
      if (!var_to_add.Valid()) return;
    }

    //then add the variant
    size_t vid;
    const auto& found = std::find(variants_.begin(), variants_.end(), var_to_add);
    if (found != variants_.end()) { //existing variants.
      vid = std::distance(variants_.begin(), found);
      cid_vid_[sid] = vid;
    } else {  // new variant
      this->ginv_.left() = std::min(this->left(), new_var.left());
#ifdef DEBUG
      if (new_var.contig() != this->contig()) {
        std::cerr<<"Warning: a variant from different contig cannot be add to the VariantGroup!\n";
        return;
      }
#endif
      this->ginv_.right() = std::max(this->right(), new_var.right());

      vid = variants_.size();
      variants_.push_back(var_to_add);
      variants_.back().id() = vid;
      cid_vid_[sid] = vid;
    }
  }


  bool empty() const {
    return variants_.empty();
  }

  size_t size() const {
    return variants_.size();
  }

  decltype(auto) ginv() const {return (ginv_);}


  decltype(auto) contig() const {return (ginv_.contig());}
  decltype(auto) vid_cid() const {return (vid_cid_);}

  decltype(auto) Init() {
    for (auto it = cid_vid_.cbegin(); it != cid_vid_.cend(); ++it) {
      vid_cid_[it->second].push_back(it->first);
    }
  }

private:

  //unique set. But it may contain deprecated variants due to MergeVaraint()
  std::vector<TVariant> variants_;
  std::map<size_t, size_t> cid_vid_; //map consensus index to variant index 
  std::map<size_t, std::vector<size_t>> vid_cid_; // map var index to list of consensus ids 

  TInv ginv_;

  neusomatic::bio::RefSeqs const* ref_seqs_;

  decltype(auto) MergeVariants_(const TVariant& left, const TVariant& right) const {
    typename TVariant::TSeq seq;
#ifdef DEBUG
    if (left.contig() != left.contig() || right < left) {
      std::cerr<<"Warning: variant "<<left<<" and "<<right<<" cannot be merged in VariantGroup!\n";
      return TVariant();
    }
#endif

    seqan::append(seq, left.allele());
    int cursor = left.right();
    for (; cursor < right.left(); ++cursor) {
      seqan::appendValue(seq, (*ref_seqs_)[contig()][cursor]);
    }
    seqan::append(seq, right.allele());
    TVariant ret(contig(), left.left(),  right.right() - left.left(), seq);
    //std::cerr<< "merge seq: "<<seq<<std::endl;
    return ret;
  }

  decltype(auto) GetAltAlleles_(size_t vid) const{

    typename TVariant::TSeq alt_allele;
    const auto& ref_seq = (*ref_seqs_)[contig()];

    int curr = this->left();

    auto const& variant = variants_[vid];
    for (; curr < variant.left(); ++curr) {
      alt_allele += ref_seq[curr];
    }
    // add variant representation

    alt_allele += variant.allele();

    curr = variant.right();

    for (; curr < this->right(); ++curr) {
        alt_allele += ref_seq[curr];
    }
    return alt_allele;
  }
};

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const VariantGroup<T>& vg)
{
  int i=1;
  os << "Variant Group: " <<vg.ginv()<<" has "<<vg.size()<<" variants.";
  for(const auto& var : vg.variants()){
    os <<"\n";
    os << "variant "<<std::to_string(i)<<": "<<var <<"\n";  
    ++i;
  }
  return os;
}


}//namespace bio
}//namespace neusomatic
#endif /* VARIANT_HPP */
