#ifndef NEU_SOMATIC_VCF_H
#define NEU_SOMATIC_VCF_H

#include <seqan/vcf_io.h>

#include <string>

#include "RefSeqs.h"

namespace neusomatic {
namespace bio {
class VCFWriter {
public:
  VCFWriter(const std::string& vcf_file, const std::string& ref_file,
            const std::string& program_name) {
    if (!open(vcf_, vcf_file.c_str())) {
      throw std::runtime_error("Coult not open " + vcf_file + " for writing.");
    }
    WriteHeader(ref_file, program_name);
  }

  void Write(seqan::VcfRecord /*const&*/ record) {  // because seqan messed up constness
    writeRecord(vcf_, record);
  }

  decltype(auto) vcf_context() { return context(vcf_); }

private:
  seqan::VcfFileOut vcf_;
  void WriteHeader(const std::string& ref_file, const std::string& program_name) {
    neusomatic::bio::RefSeqs ref_seqs(ref_file);
    for (size_t i = 0; i < ref_seqs.size(); ++i) {
      appendValue(contigNames(context(vcf_)), ref_seqs.GetName(i));
    }
    seqan::VcfHeader header;
    appendValue(header, seqan::VcfHeaderRecord("fileformat", "VCFv4.2"));
    appendValue(header, seqan::VcfHeaderRecord("source", program_name));
    appendValue(header, seqan::VcfHeaderRecord("reference", ref_file));
    // appendValue(header, seqan::VcfHeaderRecord("INFO",
    // "<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples with Data\">"));
    appendValue(
        header,
        seqan::VcfHeaderRecord(
            "INFO", "<ID=AF,Number=A,Type=Float,Description=\"Allel Frequency\">"));
    // appendValue(header, seqan::VcfHeaderRecord("FORMAT",
    // "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    writeHeader(vcf_, header);
  }
};

}  // namespace bio
}  // namespace neusomatic

#endif
