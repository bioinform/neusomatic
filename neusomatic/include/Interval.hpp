#ifndef INTERVAL_HPP
#define INTERVAL_HPP

#include <seqan/misc/interval_tree.h>
#include <seqan/seq_io.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

namespace neusomatic {
namespace bio {

enum class StrandType { PLUS = 0, MINUS, UNKNOWN };

inline std::ostream& operator<<(std::ostream& os, const StrandType& st) {
  switch (st) {
  case StrandType::PLUS:
    os << "+";
    break;
  case StrandType::MINUS:
    os << "-";
    break;
  case StrandType::UNKNOWN:
    os << "?";
    break;
  }
  return os;
}

//
// Utility functions:
//
template <typename Interval, bool half_open>
inline decltype(auto) coverage(const std::vector<Interval>& invs)
    -> std::vector<typename Interval::Depth> {
  typedef typename Interval::Depth CovType;

  const auto& lmin = std::min_element(
      invs.begin(), invs.end(),
      [&](const auto& lhs, const auto& rhs) { return lhs.left() < rhs.left(); });
  const auto& rmax = std::max_element(
      invs.begin(), invs.end(),
      [&](const auto& lhs, const auto& rhs) { return lhs.right() < rhs.right(); });
  int left = lmin->left();
  int right = rmax->right();
  assert(left < right);

  std::vector<CovType> coverage;

  if (half_open)
    coverage.resize(right - left, 0);
  else
    coverage.resize(right - left + 1, 0);

  for (auto s = invs.cbegin(); s != invs.cend(); ++s) {
    int end = 0;

    if (half_open)
      end = s->right();
    else
      end = s->right() + 1;

    for (int p = s->left(); p < end; ++p) {
      coverage[p - left] += s->depth();
    }
  }

  return coverage;
}

struct BasicInterval {
  /*
   * half open interval encoding.
   * e.g. [left, right) -> coverage
   */
  using Depth = int;
  BasicInterval()
      : i1(std::numeric_limits<int>::max()), i2(std::numeric_limits<int>::min()) {}

  BasicInterval(int l, int r) : i1(l), i2(r) {}

  virtual ~BasicInterval() {}  // virtual distructor
  virtual const int left() const { return i1; }
  virtual const int right() const { return i2; }
  virtual int& left() { return i1; }
  virtual int& right() { return i2; }
  virtual int depth() const { return 1; }
  virtual size_t length() const {
    if (right() < left())
      return 0;
    else
      return right() - left();
  }
  virtual bool Valid() const {
    return (i1 != std::numeric_limits<int>::max() &&
            i2 != std::numeric_limits<int>::min());
  }

  int i1;  // compatible with seqan
  int i2;
};

inline std::ostream& operator<<(std::ostream& os, const BasicInterval& bi) {
  os << "interval [" << bi.left() + 1 << "," << bi.right() + 1 << ")" << std::endl;
  return os;
}

inline bool operator<(const BasicInterval& lhs, const BasicInterval& rhs) {
  if (lhs.left() < rhs.left()) {
    return true;
  } else if (lhs.left() > rhs.left()) {
    return false;
  } else {
    if (lhs.right() < rhs.right())
      return true;
    else
      return false;
  }
}

template <typename Contig, typename Haplo = seqan::Dna5String>
class GenomicInterval : public BasicInterval {
  /*
   * In compliance with seqan: 0-based, half open interval [a, b).
   * Contig is the
   */
public:
  using TContig = Contig;

  GenomicInterval() : contig_(), strand_(StrandType::UNKNOWN) {}
  GenomicInterval(Contig ctg, StrandType strand) : contig_(ctg), strand_(strand) {}

  GenomicInterval(const Contig& ctg, int l, int r)
      : BasicInterval(l, r), contig_(ctg), strand_(StrandType::UNKNOWN) {}

  GenomicInterval(const Contig& ctg, const StrandType& s, int l, int r)
      : BasicInterval(l, r), contig_(ctg), strand_(s) {}

  GenomicInterval(const Contig& ctg, const StrandType& s, int l, int r, const Haplo& hap)
      : BasicInterval(l, r), contig_(ctg), strand_(s), haplotype_(hap) {}

  // int lenght() const {return right() - left();}

  decltype(auto) contig() const { return contig_; }
  decltype(auto) contig() { return (contig_); }

  decltype(auto) strand() const { return strand_; }

  decltype(auto) haplotype() const { return haplotype_; }
  decltype(auto) haplotype() { return haplotype_; }

  std::string string() const {
    std::string result;
    result += this->contig();
    result += ":";
    result += std::to_string(this->left());
    result += "-";
    result += std::to_string(this->right());
    return result;
  }

private:
  Contig contig_;
  StrandType strand_;
  Haplo haplotype_;
};

template <typename Contig, typename Haplo = seqan::Dna5String>
inline bool operator==(const GenomicInterval<Contig, Haplo>& lhs,
                       const GenomicInterval<Contig, Haplo>& rhs) {
  if (lhs.contig() != rhs.contig())
    return false;
  if (lhs.strand() != rhs.strand())
    return false;
  if (lhs.left() != rhs.left())
    return false;
  if (lhs.right() != rhs.right())
    return false;
  if (lhs.haplotype() != rhs.haplotype())
    return false;
  return true;
}

template <typename Contig, typename Haplo = seqan::Dna5String>
inline bool operator!=(const GenomicInterval<Contig, Haplo>& lhs,
                       const GenomicInterval<Contig, Haplo>& rhs) {
  return !(lhs == rhs);
}

template <typename Contig, typename Haplo = seqan::Dna5String>
inline std::ostream& operator<<(std::ostream& os,
                                const GenomicInterval<Contig, Haplo>& gInv) {
  // print based on 1-base coordinates;
  os << gInv.contig() << ":" << gInv.strand() << " [" << gInv.left() + 1 << ","
     << gInv.right() + 1 << ") " << gInv.haplotype();
  return os;
}

template <typename Contig, typename Haplo = seqan::Dna5String>
inline bool operator<(const GenomicInterval<Contig, Haplo>& lhs,
                      const GenomicInterval<Contig, Haplo>& rhs) {
  // sort contig by alphabetic order
  if (lhs.contig() < rhs.contig()) {
    return true;
  } else if (lhs.contig() > rhs.contig()) {
    return false;
  } else {
    if (lhs.left() < rhs.left()) {
      return true;
    } else if (lhs.left() > rhs.left()) {
      return false;
    } else {
      if (lhs.right() < rhs.right())
        return true;
      else
        return false;
    }
  }
}

template <typename Interval = BasicInterval, bool half_open = true>
class IRanges {
  /*
   * Use BasicInterval as default interval representation inside class
   */

private:
  std::vector<Interval> invs_;
  int left_;
  int right_;
  typename Interval::TContig contig_;

  void Setup_() {
    if (invs_.empty()) {
      throw std::runtime_error("IRangs cannot build empty interval set");
    }
    contig_ = invs_[0].contig();
    for (auto const& inv : invs_) {
      if (inv.contig() != contig_) {
        throw std::runtime_error("IRangs intervals cannot come from different contigs");
      }
    }
    const auto& lmin = std::min_element(
        invs_.begin(), invs_.end(),
        [&](const auto& lhs, const auto& rhs) { return lhs.left() < rhs.left(); });

    const auto& rmax = std::max_element(
        invs_.begin(), invs_.end(),
        [&](const auto& lhs, const auto& rhs) { return lhs.right() < rhs.right(); });

    left_ = lmin->left();
    right_ = rmax->right();
  }

public:
  typedef typename std::vector<Interval>::const_iterator const_iterator;
  typedef typename std::vector<Interval>::size_type size_type;

  template <typename IntVec = std::initializer_list<int>>
  IRanges(IntVec starts, IntVec ends) {
    assert(starts.size() == ends.size());
    auto e = ends.cbegin();
    for (auto s = starts.cbegin(); s != starts.cend(); ++s) {
      int end_pos = *e;
      if (!half_open)
        ++end_pos;  // if input is a close interval, convert to half open interval
      Interval inv(*s, end_pos);
      invs_.push_back(move(inv));
      ++e;
    }
    Setup_();
  }

  IRanges(std::vector<Interval> invs) : invs_(invs) {
    if (!half_open) {
      for (Interval& inv : invs_) {
        inv.right()++;
      }
    }
    Setup_();
  }

  size_type size() { return invs_.size(); }

  const_iterator cbegin() { return invs_.cbegin(); }
  const_iterator cend() { return invs_.cend(); }

  std::vector<Interval> reduce() const {
    /*
     *  Merge redundant	ranges,	and return the minimum
     *  non-overlapping ranges covering all the input ranges.
     */
    std::vector<typename Interval::Depth> cov = coverage<Interval, true>(invs_);
    std::vector<Interval> result;

    bool prev_base_is_covered = false;
    bool open = false;

    for (auto it = cov.cbegin(); it != cov.cend(); ++it) {
      if (*it > 0) {
        if (prev_base_is_covered) {
        } else {
          Interval inv;
          inv.contig() = contig_;
          inv.left() = left_ + std::distance(cov.cbegin(), it);
          result.push_back(inv);
          open = true;
        }
        prev_base_is_covered = true;
      } else {
        if (prev_base_is_covered) {
          result.back().right() = left_ + std::distance(cov.cbegin(), it);
          open = false;
        } else {
        }
        prev_base_is_covered = false;
      }
    }
    if (open) {
      result.back().right() = left_ + cov.size();
    }

    if (!half_open) {
      for (auto& res : result) res.right() = res.right() - 1;
    }

    return result;
  }

  std::vector<Interval> disjoint() const {
    /*
     * return non-overlapping intervals
     */
    std::vector<typename Interval::Depth> cov = coverage<Interval, true>(invs_);
    std::vector<int> bars;
    std::vector<Interval> result;
    for (const auto& inv : invs_) {
      bars.push_back(inv.left());
      bars.push_back(inv.right());
    }
    sort(bars.begin(), bars.end());
    const auto& last = unique(bars.begin(), bars.end());
    bars.erase(last, bars.end());

    bool is_left = true;
    for (auto it = bars.cbegin(); it != bars.cend(); ++it) {
      if (is_left) {
        Interval sub_inv;
        sub_inv.contig() = contig_;
        sub_inv.left() = *it;
        result.push_back(sub_inv);
        is_left = false;
      } else {
        if (*it == result.back().left()) {
          result.pop_back();
        } else {
          result.back().right() = *it;
          if (cov[*it - left_] > 0)
            --it;
        }
        is_left = true;
      }
    }
    if (!is_left) {
      result.pop_back();
    }

    if (!half_open) {
      for (auto& res : result) res.right() = res.right() - 1;
    }
    return result;
  }
};

/*
 * Below are utility functions:
 * Those functions must take half closed interval [a,b)
 * as input arguments
 */

template <typename TInterval>
inline bool IsOverlapped(const TInterval& x, const TInterval& y) {
  // intervally this function use closed interval
  int yi2 = y.i2 - 1;
  int xi2 = x.i2 - 1;
  if (x.i1 <= yi2 && y.i1 <= xi2)
    return true;
  else
    return false;
}

template <typename TInterval>
inline bool IsContainedIn(const TInterval& x, const TInterval& y) {
  if (x.contig() != y.contig())
    return false;
  return (x.i1 >= y.i1 && x.i2 <= y.i2 && x.i1 < y.i2 && x.i2 > y.i1);
}

template <typename TInterval>
inline bool IsNearTheStartOf(const TInterval& x, const TInterval& y, int wiggleroom) {
  // return true if x is near the start of y
  if (x.i1 < y.i1 + wiggleroom)
    return true;
  else
    return false;
}

template <typename TInterval>
inline bool IsNearTheEndOf(const TInterval& x, const TInterval& y, int wiggleroom) {
  // return true if x is near the end of y
  if (x.i2 > y.i2 - wiggleroom)
    return true;
  else
    return false;
}

template <typename TInterval>
inline TInterval Intersect(const TInterval& x, const TInterval& y) {
  // return a 0 length interval if non-overlap
  TInterval ret;
  ret.i1 = std::max(x.i1, y.i1);
  ret.i2 = std::min(x.i2, y.i2);

  if (ret.i2 <= ret.i1) {
    return TInterval(ret.i1, ret.i1, x.cargo);
  } else {
    ret.cargo = x.cargo;
    return ret;
  }
}

template <typename Contig>
inline neusomatic::bio::GenomicInterval<Contig> Intersect(
    const neusomatic::bio::GenomicInterval<Contig>& x,
    const neusomatic::bio::GenomicInterval<Contig>& y) {
  // return a 0 length interval if non-overlap
  neusomatic::bio::GenomicInterval<Contig> ret;
  if (x.contig() != y.contig())
    return ret;
  ret.left() = std::max(x.i1, y.i1);
  ret.right() = std::min(x.i2, y.i2);
  if (ret.i2 <= ret.i1)
    ret.i2 = ret.i1;
  return ret;
}

template <typename TIntervalTree, typename TInterval>
inline std::vector<TInterval> FindIntervals(const TIntervalTree& inv_tree,
                                            const std::vector<TInterval>& invs) {
  std::vector<TInterval> results;

  for (auto i = invs.cbegin(); i != invs.cend(); ++i) {
    seqan::String<TInterval> overlapped_invs;
    seqan::findIntervals(overlapped_invs, inv_tree, i->left(), i->right());

    for (unsigned i = 0; i < seqan::length(overlapped_invs); ++i) {
      results.push_back(overlapped_invs[i]);
    }
  }

  return results;
}

template <typename TInterval>
inline decltype(auto) SplitByContig(const std::vector<TInterval> invs) {
  using TContig = typename TInterval::TContig;
  std::map<TContig, std::vector<TInterval>> result;
  for (auto it = invs.cbegin(); it != invs.cend(); ++it) {
    auto search = result.find(it->contig());
    if (search == result.end()) {
      std::vector<TInterval> invs_for_one_ctg;
      invs_for_one_ctg.push_back(*it);
      result.emplace(it->contig(), invs_for_one_ctg);
    } else {
      result[it->contig()].push_back(*it);
    }
  }
  return result;
}

template <typename TInterval>
inline auto Length(const TInterval& inv) {
  if (inv.i1 > inv.i2)
    return 0;
  return inv.i2 - inv.i1;
}

}  // namespace bio
}  // namespace neusomatic

#endif /* INTERVAL_HPP */
