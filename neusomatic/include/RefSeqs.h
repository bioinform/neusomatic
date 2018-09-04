#ifndef CPPUTIL_BIO_REF_SEQS_H
#define CPPUTIL_BIO_REF_SEQS_H

#include <vector>
#include <map>
#include <string>
#include <seqan/seq_io.h>

namespace neusomatic { namespace bio {

class RefSeqs {
public:
  using Seq = seqan::Dna5String; // can templatize for other strings for the future
  explicit RefSeqs(std::string const& file_name, const bool load_seq = false) {
    load(file_name);
  }

  decltype(auto) operator[](size_t rid) const {
    return (rid_seq_[rid]);
  }

  decltype(auto) GetName(size_t rid) const {
    return (rid_str_[rid]);
  }

  decltype(auto) GetNames() const {
    return (rid_str_);
  }

  template<class Str>
  decltype(auto) GetId(Str const& str) const {
    const auto& itr = str_rid_.find(str);
    if (itr == str_rid_.end()) {
      throw std::runtime_error( str + " not found in reference database");
    }
    return itr->second;
  }

  auto size() const {
    return rid_seq_.size();
  }

  template<class Str>
  decltype(auto) operator[](Str const& str) const {
    const auto itr = str_rid_.find(str);
    if (itr == str_rid_.end()) {
      throw std::runtime_error( str + " not found in reference database");
    }
    return rid_seq_[ itr->second ];
  }

  void load(std::string const& file_name, const bool load_seq = false) {
    seqan::FaiIndex fai;
    if (not open(fai, file_name.c_str())) {
      throw std::runtime_error("cannot find index file for "+file_name);
    };

    rid_seq_.resize(numSeqs(fai));
    rid_str_.resize(numSeqs(fai));
    str_rid_.clear();
    for (size_t rid = 0; rid < numSeqs(fai); ++rid) {
      if (load_seq) {
        readSequence(rid_seq_[rid], fai, rid); // because this is not thread-safe
      }
      rid_str_[rid] = seqan::sequenceName(fai, rid);
      str_rid_[rid_str_[rid]] = rid;
    }
  }

private:
  // can use on-demand load, but seqan's fai is not thread safe
  std::vector<Seq> rid_seq_;
  std::vector<seqan::CharString> rid_str_;
  std::map<seqan::CharString, size_t> str_rid_;
};

} //namespace bio
} //namespace neusomatic


#endif
