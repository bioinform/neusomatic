#ifndef LIB_INCLUDE_BEDIO_HPP
#define LIB_INCLUDE_BEDIO_HPP

#include <seqan/bed_io.h>

#include <string>
#include <vector>

#include "Interval.hpp"

namespace neusomatic {

template <typename Interval>
class BedIO {
public:
  BedIO(const std::string& filepath) : filepath_(filepath) {}

  static void WriteToBed3(const std::vector<Interval>& invs, seqan::BedFileOut& out) {
    for (const auto& inv : invs) {
      seqan::BedRecord<seqan::Bed3> record;
      record.beginPos = inv.left();
      record.endPos = inv.right();
      record.ref = inv.contig();
      seqan::writeRecord(out, record);
    }
  }

  std::vector<Interval> ReadBed3() {
    seqan::BedFileIn bed_in_(filepath_.c_str());
    if (seqan::atEnd(bed_in_)) {
      throw std::runtime_error("It seems to be a corrupted or wrong bed file: " +
                               filepath_);
    }

    std::vector<Interval> result;
    while (!seqan::atEnd(bed_in_)) {
      seqan::readRecord(record_, bed_in_);
      // always left <= right
      int l = std::min(record_.beginPos, record_.endPos);
      int r = std::max(record_.endPos, record_.endPos);
      auto refcontig = seqan::toCString(record_.ref);
      Interval inv(refcontig, l, r);
      result.push_back(inv);
    }
    return result;
  }

  std::vector<Interval> ReadBed3_windowed(int window_size) {
    seqan::BedFileIn bed_in_(filepath_.c_str());
    if (seqan::atEnd(bed_in_)) {
      throw std::runtime_error("It seems to be a corrupted or wrong bed file: " +
                               filepath_);
    }
    std::vector<Interval> result;
    while (!seqan::atEnd(bed_in_)) {
      seqan::readRecord(record_, bed_in_);
      // always left <= right
      int l = std::min(record_.beginPos, record_.endPos);
      int r = std::max(record_.endPos, record_.endPos) + 1;
      int span = r - l;
      auto refcontig = seqan::toCString(record_.ref);
      for (int i = 0; i <= span / window_size; i++) {
        int l_2 = std::max(l + i * window_size - 1, l);
        int r_2 = std::min(r, l + (i + 1) * window_size);
        if ((r - r_2) < window_size) {
          r_2 = r;
        }
        Interval inv(refcontig, l_2, r_2);
        result.push_back(inv);
        if (r_2 == r) {
          break;
        }
      }
    }
    return result;
  }

  std::vector<Interval> ReadAndReduceBedInvs() {
    std::vector<Interval> result;
    auto splited_invs = ReadAndSplitBedInvs();
    for (auto it = splited_invs.cbegin(); it != splited_invs.cend(); ++it) {
      neusomatic::bio::IRanges<Interval> iranges(it->second);
      neusomatic::bio::IRanges<Interval> reduced_iranges = iranges.reduce();

      for (auto inv = reduced_iranges.cbegin(); inv != reduced_iranges.cend(); ++inv) {
        Interval ginv(seqan::toCString(it->first), inv->left(), inv->right());
        result.push_back(ginv);
      }
    }
    return result;
  }

  std::vector<Interval> ReadAndDisjointBedInvs() {
    std::vector<Interval> result;
    auto splited_invs = ReadAndSplitBedInvs();
    for (auto it = splited_invs.cbegin(); it != splited_invs.cend(); ++it) {
      neusomatic::bio::IRanges<Interval> iranges(it->second);
      neusomatic::bio::IRanges<Interval> reduced_iranges = iranges.disjoint();

      for (auto inv = reduced_iranges.cbegin(); inv != reduced_iranges.cend(); ++inv) {
        Interval ginv(seqan::toCString(it->first), inv->left(), inv->right());
        result.push_back(ginv);
      }
    }
    return result;
  }

  decltype(auto) ReadAndSplitBedInvs() {
    std::vector<Interval> bed_invs = ReadBed3();
    return neusomatic::bio::SplitByContig(bed_invs);
  }

private:
  const std::string& filepath_;
  // seqan::BedFileIn bed_in_;
  seqan::BedRecord<seqan::Bed3> record_;  // place holder
};

}  // namespace neusomatic

#endif /* LIB_INCLUDE_BEDIO_HPP */
