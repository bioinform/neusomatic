#ifndef CHANGE_COORDINATES_HPP
#define CHANGE_COORDINATES_HPP

#include "SeqanUtils.h"

namespace neusomatic {

class ChangeCoordinates {
  std::vector<int> ungap_pos_;
  std::vector<int> gap_pos_; // reduce to the left-most position
  std::string ref_seq_;
  static const unsigned char gap_chr_ = '-';
public:
  ChangeCoordinates(const std::string& refgap) {
    int i = 0;
    int ungap_cursor = 0;
    auto num_gaps = std::distance(begin(refgap), end(refgap));
    ungap_pos_.resize(num_gaps);
    gap_pos_.reserve(num_gaps);

    for (auto it = begin(refgap); it != end(refgap); ++it, ++i) {
      if (*it == gap_chr_) {
        ungap_pos_[i] = ungap_cursor; 
      } else {
        gap_pos_.push_back(i);
        ungap_pos_[i] = ungap_cursor; 
        ref_seq_ += neusomatic::bio::Base2String(*it, "N");
        ungap_cursor++;
      }
    }
  }

  int UngapPos(const int gap_pos) const {
    return ungap_pos_[gap_pos];  
  }

  int GapPos(const int ref_pos) const {
    return gap_pos_[ref_pos];
  }

  const std::string& RefSeq() const {
    return ref_seq_; 
  }
};

}// namespace neusomatic



#endif /* ALlAMERICAN_CHANGE_COORDINATES_HPP */
