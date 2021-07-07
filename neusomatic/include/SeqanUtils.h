#ifndef SEQAN_UTILS_HPP
#define SEQAN_UTILS_HPP

#include <seqan/find.h>
#include <seqan/misc/interval_tree.h>
#include <seqan/seq_io.h>

#include <vector>

namespace neusomatic {
namespace bio {

template <typename Base>
std::string Base2String(const Base& base, const std::string& gap_symbol) {
  // std::cout<<static_cast<unsigned char>(base)<<" :static cast: ";
  // std::cout<<(unsigned)seqan::ordValue(base)<<std::endl;

  switch ((unsigned)seqan::ordValue(base)) {
  case 0:
    return "A";
  case 1:
    return "C";
  case 2:
    return "G";
  case 3:
    return "T";
  default:
    return gap_symbol;
  }
  // should never reach here
  return "";
}

template <>
std::string Base2String<char>(const char& base, const std::string& gap_symbol) {
  switch (base) {
  case 'a':
  case 'A':
  case 0:
    return "A";

  case 'c':
  case 'C':
  case 1:
    return "C";

  case 'g':
  case 'G':
  case 2:
    return "G";
  case 't':
  case 'T':
  case 3:
    return "T";
  default:
    return gap_symbol;
  }
}

template <typename Itr>
std::string to_string(Itr b, Itr e, const std::string& gap_symbol = "N") {
  std::string result;
  auto it = b;
  for (; it < e; ++it) {
    result += Base2String(*it, gap_symbol);
  }
  return result;
}

inline std::string to_string(const seqan::Dna5String& seq,
                             const std::string& gap_symbol = "N") {
  std::string result;
  for (unsigned int i = 0; i < seqan::length(seq); ++i) {
    result += Base2String(seq[i], gap_symbol);
  }
  return result;
}

template <typename TInterval>
inline decltype(auto) ToSeqanIntervals(const std::vector<TInterval>& invs) {
  typedef seqan::IntervalAndCargo<int, TInterval> TIntervalAndCargo;
  std::vector<TIntervalAndCargo> seqan_intervals;
  for (auto it = invs.cbegin(); it != invs.cend(); ++it) {
    TIntervalAndCargo i(it->left(), it->right(), *it);
    seqan_intervals.push_back(i);
  }
  return seqan_intervals;
}

inline std::vector<seqan::CharString> SeqanSplit(seqan::CharString target,
                                                 std::string splitter) {
  /*
    Return a vec of size 1 containing the intact target if splitter is not found in the
    target.
  */
  if (seqan::length(target) == 0 || seqan::length(splitter) == 0) {
    throw std::runtime_error("Invalid string and pattern to split.");
  }
  using TPos = typename seqan::Position<seqan::CharString>::Type;

  std::vector<seqan::CharString> result;
  seqan::Finder<seqan::CharString> finder(target);
  seqan::Pattern<seqan::CharString, seqan::Horspool> pattern(splitter);
  std::vector<std::pair<TPos, TPos>> searches;
  while (seqan::find(finder, pattern)) {
    searches.emplace_back(seqan::beginPosition(finder), seqan::endPosition(finder));
  }

  TPos lastp = 0;
  for (auto const& p : searches) {
    seqan::CharString unit = seqan::infix(target, lastp, p.first);
    if (!seqan::empty(unit))
      result.push_back(unit);
    lastp = p.second;
  }
  if (lastp < seqan::length(target) - 1) {
    seqan::CharString unit = seqan::infix(target, lastp, seqan::length(target));
    result.push_back(unit);
  }
  return result;
}

}  // namespace bio
}  // namespace neusomatic

#endif /* SEQAN_UTILS_HPP */
