#ifndef LIB_INCLUDE_FRAGSTORESEQAN_HPP_
#define LIB_INCLUDE_FRAGSTORESEQAN_HPP_

#include <vector>
#include <map>
#include <unordered_map>

#include <SeqLib/BamReader.h>
#include "Options.h"

namespace neusomatic{

template<typename GInv, typename BamRecord, typename BamReader>
class CaptureLayout {
  size_t region_idx_;
  const neusomatic::Options & opts_; 
  BamReader bam_reader_;
  std::vector<GInv> ginvs_;
  //std::map <std::string, std::vector<GInv>> contig_ginv_;
  std::map <GInv, std::vector<BamRecord>> ginv_records_;

public:
  template<typename GInvString>
  CaptureLayout(const std::string& bam_path, const std::vector<GInvString>& regions, const neusomatic::Options& o): 
      region_idx_(0), opts_(o) {
      bam_reader_.Open(bam_path);
      LoadRegions_(regions);
  }

  decltype(auto) ginv_records() const {
    return (ginv_records_);
  }

  decltype(auto) Reader() const {
    return (bam_reader_);
  }

  std::map<GInv, std::vector<BamRecord>> GetRegion(const int region_idx, const bool full_contained) {
    if (region_idx < 0 || region_idx >= ginvs_.size()) {
      throw std::runtime_error("region idx not exist");
    }
    std::map <GInv, std::vector<BamRecord>> ginv_records;
    const GInv& curr_ginv = ginvs_[region_idx];
    SeqLib::GenomicRegion gr(curr_ginv.contig(), curr_ginv.left(), curr_ginv.right());
    bam_reader_.SetRegion(gr); 
    std::vector<BamRecord> records;
    BamRecord rec;
    while (bam_reader_.GetNextRecord(rec)) {
      if (full_contained) {
        if (rec.Position() <= curr_ginv.left() && rec.PositionEnd() + 1 >= curr_ginv.right() && rec.MapQuality()>=opts_.min_mapq() && (!rec.SecondaryFlag()) ) {
          records.push_back(rec);
        }
      } else {
        if (rec.Position() < curr_ginv.right() && curr_ginv.left() <= rec.PositionEnd() && rec.MapQuality() >= opts_.min_mapq() && (!rec.SecondaryFlag())) { // overlapped
          records.push_back(rec);
        }
      }
    }
    ginv_records[curr_ginv] = records;
    return ginv_records;
  }

  size_t NumRegion() const {
    return ginvs_.size();
  }

  bool NextRegion(const bool full_contained) {
    if (region_idx_  == ginvs_.size()) return false;
    ginv_records_.clear();

    const GInv& curr_ginv = ginvs_[region_idx_];
    SeqLib::GenomicRegion gr(curr_ginv.contig(), curr_ginv.left(), curr_ginv.right());
    bam_reader_.SetRegion(gr); 

    std::vector<BamRecord> records;
    BamRecord rec;
    while (bam_reader_.GetNextRecord(rec)) {
      if (full_contained) {
        if (rec.Position() <= curr_ginv.left() && rec.PositionEnd() + 1 >= curr_ginv.right() && rec.MapQuality()>=opts_.min_mapq() && (!rec.SecondaryFlag())) {
          records.push_back(rec);
        }
      } else {
        if (rec.Position() < curr_ginv.right() && curr_ginv.left() <= rec.PositionEnd() && rec.MapQuality()>=opts_.min_mapq() && (!rec.SecondaryFlag())) { // overlapped
          records.push_back(rec);
        }
      }
    }
    ginv_records_[curr_ginv] = records;

    region_idx_ ++;

    return true;
  }

private:
  template<typename GInvString>
  void LoadRegions_(const std::vector<GInvString>& regions) {
    ginvs_.resize(regions.size());    
    for (size_t i = 0; i < regions.size(); ++i) {
      ginvs_[i] = GInv(bam_reader_.Header().Name2ID(regions[i].contig()), regions[i].left(), regions[i].right());
    }
  }

};

}//namespace neusomatic

#endif /* LIB_INCLUDE_FRAGSTORESEQAN_HPP_ */
