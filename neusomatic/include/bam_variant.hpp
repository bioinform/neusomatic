#ifndef BAM_VARIANT_HPP
#define BAM_VARIANT_HPP

// For an alignment record, get the variants based on the MD tags. 

#include <vector>
#include <string>
#include "variant.hpp"

namespace neusomatic{

inline bool IsNuc( char c )
{
   static const char alpha[] = "ACGTRYKMSWBDHVN";
   return ( std::strchr( alpha, c ) != NULL );
}

inline bool IsNumber(const std::string& s)
{
  return !s.empty() && std::find_if(s.begin(), s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
} 

enum class MDType { MATCH = 0, SUB, DEL };
inline std::ostream& operator<<(std::ostream & os, MDType & mdt) {
  switch(mdt) {
    case MDType::MATCH:
      os << 'M';
      break;
    case MDType::SUB:
      os << 'X';
      break;
    case MDType::DEL:
      os << 'D';
      break;
    default:
      os << "ERROR";
      break;
  }
  return os;
}

class MDTags {
private:
  const std::string tag_string_;
  int16_t curr_pos_;
  size_t size;
public:
  MDTags(const std::string& tags) : tag_string_(tags), curr_pos_(0), size(tags.length()) {
  }

   bool Next(std::pair<MDType, std::string>& result) {
    if (curr_pos_ == size) {return false;}
    if (std::isdigit(tag_string_[curr_pos_])) { // match
      if (tag_string_[curr_pos_] == '0') {
        curr_pos_ ++;
      }
      else {
        int len = 0;
        while (curr_pos_ + len < size && std::isdigit(tag_string_[curr_pos_ + len])) {
          ++len;
        } 
        result = std::make_pair(MDType::MATCH, tag_string_.substr(curr_pos_, len));
        curr_pos_ += len ;
        return true;
      }
    } 

    if (tag_string_[curr_pos_] == '^') { // deletion
      int len = 0;
      ++curr_pos_;
      while (curr_pos_ + len < size && IsNuc(tag_string_[curr_pos_ + len])) {
        ++len;
      } 
      result = std::make_pair(MDType::DEL, tag_string_.substr(curr_pos_, len));
      curr_pos_ += len;
      return true;
    }

    if (IsNuc(tag_string_[curr_pos_])) { // substitution
      int len = 0;
      std::string s;
      while (curr_pos_ + len < size && (IsNuc(tag_string_[curr_pos_ + len]) || tag_string_[curr_pos_ + len] == '0')) {
        if (IsNuc(tag_string_[curr_pos_ + len])) {
          s += tag_string_.substr(curr_pos_ + len, 1); 
        }
        ++len;
      }
      result = std::make_pair(MDType::SUB, s);
      curr_pos_ += len;
      return true;
    }
    // should not reach here
    return false;
  }
};

template<typename BamRecord, typename Contig>
std::vector<neusomatic::bio::Variant<std::string, Contig>> GetIndels(const BamRecord& bam) {
  std::vector<neusomatic::bio::Variant<std::string, Contig>> result;
  if (bam.Position() < 0) return result;
  auto const& cigar = bam.GetCigar(); 
  auto const& seq = bam.Sequence();
  int32_t ref_pos = bam.Position();
  int32_t read_pos = bam.AlignmentPosition();
  for (auto c = cigar.begin(); c != cigar.end(); ++c) {
    std::string s;
    switch (c->Type()) {
      case 'H':
        if (c == cigar.begin()){
          read_pos -= c->Length();
        }
        break;
      case 'M':
        ref_pos += c->Length();
        read_pos += c->Length();
        break;
      case 'I':
        s = bam.Sequence().substr(read_pos, c->Length());
        result.emplace_back(bam.ChrID(), ref_pos, 0, s);
        read_pos += c->Length();
        break;
      case 'D':
        result.emplace_back(bam.ChrID(), ref_pos, c->Length(), s);
        ref_pos += c->Length();
        break;
      default:
        break;
    }
  }

  return result;
}

template<typename BamRecord, typename Contig>
std::vector<neusomatic::bio::Variant<std::string, Contig>> GetSNVs(const BamRecord& bam) {
  std::vector<neusomatic::bio::Variant<std::string, Contig>> result;
  std::string mdstr;
  bool status = bam.GetZTag("MD", mdstr);
  if (!status) {
    std::cerr<<"Warning: no MD tag!\n";
    return result;
  }
  if (IsNumber(mdstr)) {
    return result;
  }
  std::transform(mdstr.begin(), mdstr.end(),mdstr.begin(), ::toupper);
  if (bam.Position() < 0) return result;
  auto const& cigar = bam.GetCigar(); 
  auto const& seq = bam.Sequence();
  int32_t ref_pos = bam.Position();
  int32_t read_pos = bam.AlignmentPosition();

  MDTags md(mdstr);
  std::pair<MDType, std::string> unit;

  std::vector<neusomatic::bio::Variant<std::string, Contig>> snvs_MDTags;
  int sum_M_len_before_D = 0;
  int match_len = 0;
  int ref_consumer = 0; 

  for (auto c = cigar.begin(); c != cigar.end(); ++c) {
    match_len = 0;
    ref_consumer = 0; 

    switch (c->Type()) {
      case 'M':
        sum_M_len_before_D += c->Length();
        break;
      case 'D':
        while (md.Next(unit)) {
          if (unit.first == MDType::MATCH) {
            match_len = stoi(unit.second);
            ref_pos += match_len;
            ref_consumer += match_len;
          }
          if (unit.first == MDType::SUB) {
            snvs_MDTags.emplace_back(bam.ChrID(), ref_pos, unit.second.length(), unit.second);
            ref_consumer += unit.second.length();
            ref_pos += unit.second.length();
          }
          if (ref_consumer == sum_M_len_before_D) break;
        }
        sum_M_len_before_D = 0;
        ref_pos += c->Length();
        break;
      default:
        break;
    }
  }

  ref_consumer = 0; 
  while (md.Next(unit)) {
    if (unit.first == MDType::MATCH) {
      match_len = stoi(unit.second);
      ref_pos += match_len;
      ref_consumer += match_len;
    }
    if (unit.first == MDType::SUB) {
      snvs_MDTags.emplace_back(bam.ChrID(), ref_pos, unit.second.length(), unit.second);
      ref_consumer += unit.second.length();
      ref_pos += unit.second.length();
    }
    if (ref_consumer == sum_M_len_before_D) break;
  }

  ref_pos = bam.Position();
  read_pos = bam.AlignmentPosition();
  auto b = cigar.begin();
  if (b->Type() == 'H') {
    if (b == cigar.begin()){
    read_pos -= b->Length();
    }
  }
  for (const auto& t : snvs_MDTags) {
    const int32_t ref_consume_goal = t.left() -  ref_pos;
    int32_t ref_consume = 0;
    int32_t read_consume = 0;
    while (ref_consume < ref_consume_goal) {
      if (b->Type() == 'M') {
        ref_consume += b->Length();
        read_consume += b->Length();
      } else if (b->Type() == 'D') {
        ref_consume += b->Length();
      } else if (b->Type() == 'I') {
        read_consume += b->Length();
      }
      ++b;
    }
    int32_t read_shift = read_consume - (ref_consume - ref_consume_goal);
    if (b != cigar.end() && b->Type() == 'I' && ref_consume == ref_consume_goal) { // SUB after an INS
      read_shift += b->Length();
    }
    std::string s = seq.substr(read_pos + read_shift, t.allele().length()); 
    result.emplace_back(bam.ChrID(), t.left(), t.length(), s);
    read_pos +=  read_consume;
    ref_pos += ref_consume;
  }

  return result;
}

template<typename BamRecord, typename Contig>
std::vector<neusomatic::bio::Variant<std::string, Contig>> GetVars(const BamRecord& bam) {
  auto results = GetIndels<BamRecord, Contig>(bam);
  auto snvs = GetSNVs<BamRecord, Contig>(bam);
  results.insert(results.end(), snvs.begin(), snvs.end());
  //std::sort(results.begin(), results.end());
  return results;
}
}/* ALLAMERICAN_VARIANT_HPP */
#endif 
