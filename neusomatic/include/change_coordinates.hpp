#ifndef ALlAMERICAN_CHANGE_COORDINATES_HPP
#define ALlAMERICAN_CHANGE_COORDINATES_HPP

#include "SeqanUtils.h"

namespace neusomatic {

template<typename RefGap>
class ChangeCoordinates {
  const RefGap& ref_gaps_;
  std::vector<int> ref_pos_;
  std::vector<int> gap_pos_; // reduce to the left-most position
  std::string ref_seq_;
  static const unsigned char gap_chr_ = '-';
public:
  ChangeCoordinates(const RefGap& refgap): ref_gaps_(refgap) {
    int i = 0;
    int ref_cursor = 0;
    auto num_gaps = std::distance(begin(refgap), end(refgap));
    ref_pos_.resize(num_gaps);
    gap_pos_.reserve(num_gaps);

    for (auto it = begin(refgap); it != end(refgap); ++it, ++i) {
      if (*it == gap_chr_) {
        ref_pos_[i] = ref_cursor; 
      } else {
        gap_pos_.push_back(i);
        ref_pos_[i] = ref_cursor; 
        ref_seq_ += neusomatic::bio::Base2String(*it, "N");
        ref_cursor++;
      }
    }
  }

  int RefPos(const int gap_pos) const {
    return ref_pos_[gap_pos];  
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
