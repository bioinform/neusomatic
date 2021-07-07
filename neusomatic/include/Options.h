
#ifndef LIB_INCLUDE_OPTIONS_H_
#define LIB_INCLUDE_OPTIONS_H_

#include <getopt.h>

#include <limits>
#include <string>
#include <thread>

namespace neusomatic {

static struct option long_options[] = {
    {"verbosity", required_argument, 0, 'v'},
    {"bam", required_argument, 0, 'b'},
    {"bed", required_argument, 0, 'L'},
    {"ref", required_argument, 0, 'r'},
    {"calculate_qual_stat", no_argument, 0, 'c'},
    {"min_mapq", required_argument, 0, 'q'},
    {"snp_min_bq", required_argument, 0, 'n'},
    {"snp_min_af", required_argument, 0, 'a'},
    {"ins_min_af", required_argument, 0, 'i'},
    {"del_min_af", required_argument, 0, 'l'},
    {"snp_min_ao", required_argument, 0, 'M'},
    {"out_vcf_file", required_argument, 0, 'f'},
    {"out_count_file", required_argument, 0, 'o'},
    {"fully_contained", no_argument, 0, 'y'},
    {"window_size", required_argument, 0, 'w'},
    {"num_threads", required_argument, 0, 't'},
    {"max_depth", required_argument, 0, 'd'},
    {"min_depth", required_argument, 0, 'h'},
    {"include_secondary", no_argument, 0, 's'},
    {"filter_duplicate", no_argument, 0, 'D'},
    {"filter_QCfailed", no_argument, 0, 'Q'},
    {"filter_improper_pair", no_argument, 0, 'E'},
    {"filter_mate_unmapped", no_argument, 0, 'F'},
    {"filter_improper_orientation", no_argument, 0, 'G'},
    {"report_all_alleles", no_argument, 0, 'A'},
    {"report_count_for_all_positions", no_argument, 0, 'C'},
    {0, 0, 0, 0}  // terminator
};

const char* short_options = "v:b:L:r:q:f:o:w:a:t:d:cysDQEFG";

void print_help() {
  std::cerr << "---------------------------------------------------\n";
  std::cerr << "Usage: scan_alignments [options] -b BAM -L BED -r REF\n";
  std::cerr << "General Options:\n";
  std::cerr << "-v/--verbosity,                       verbosity level from 0 to 9,       "
               "                                    default 0.\n";
  std::cerr << "-b/--bam,                             bam file path,                     "
               "                                    required.\n";
  std::cerr << "-L/--bed,                             targeted region in bam format,     "
               "                                    required.\n";
  std::cerr << "-r/--ref,                             reference file path,               "
               "                                    required.\n";
  std::cerr << "-c/--calculate_qual_stat,             calculating base quality and other "
               "stats,                              default False.\n";
  std::cerr << "-q/--min_mapq,                        minimum mapping quality,           "
               "                                    default 0.\n";
  std::cerr << "-n/--snp_min_bq,                        SNP minimum base quality,        "
               "                                       default 10.\n";
  std::cerr << "-a/--snp_min_af,                          SNP minimum allele freq,       "
               "                                            default 0.01.\n";
  std::cerr << "-i/--ins_min_af,                          INS minimum allele freq,       "
               "                                            default 0.01.\n";
  std::cerr << "-l/--del_min_af,                          DEL minimum allele freq,       "
               "                                            default 0.01.\n";
  std::cerr << "-M/--snp_min_ao,                          SNP minimum alternate count "
               "for low AF candidates, default 3.\n";
  std::cerr << "-f/--out_vcf_file,                    output vcf file path,              "
               "                                    required.\n";
  std::cerr << "-o/--out_count_file,                  output count file path,            "
               "                                    required.\n";
  std::cerr << "-w/--window_size,                     window size to scan the variants,  "
               "                                    default is 15.\n";
  std::cerr << "-y/--fully_contained,                 if this option is on. A read has "
               "to be fully contained in the region,  default is False.\n";
  // std::cerr<< "-t/--num_threads,                     number or thread used for building
  // the count matrix,                   default is 4.\n";
  std::cerr << "-d/--max_depth,                       maximum depth for building the "
               "count matrix,                           default is 40,000.\n";
  std::cerr << "-h/--min_depth,                       minimum depth for building the "
               "count matrix,                           default is 0.\n";
  std::cerr << "-s/--include_secondary,               consider secondary alignments,     "
               "                                    default is False.\n";
  std::cerr << "-D/--filter_duplicate,                filter duplicate reads if the flag "
               "is set,                             default is False.\n";
  std::cerr << "-Q/--filter_QCfailed,                 filter QC failed reads if the flag "
               "is set,                             default is False.\n";
  std::cerr << "-E/--filter_improper_pair,            filter improper pairs if the flag "
               "is set,                              default is False.\n";
  std::cerr << "-F/--filter_mate_unmapped,            filter reads with unmapeed mates "
               "if the flag is set,                   default is False.\n";
  std::cerr << "-G/--filter_improper_orientation,     filter reads with improper "
               "orientation (not FR) or different chrom,    default is False.\n";
  std::cerr << "-A/--report_all_alleles,     report all alleles per column,    default "
               "is False.\n";
  std::cerr << "-C/--report_count_for_all_positions,     report counts for all "
               "positions,    default is False.\n";
}

int parseInt(const char* optarg, int lower, const char* errmsg, void (*print_help)()) {
  long l;
  char* endPtr = NULL;
  l = strtol(optarg, &endPtr, 10);
  if (endPtr != NULL) {
    if (l < lower) {
      std::cerr << errmsg << std::endl;
      print_help();
      exit(1);
    }
    return (int32_t)l;
  }
  std::cerr << errmsg << std::endl;
  print_help();
  exit(1);
  return -1;
}

float parseFloat(const char* optarg, float lower, float upper, const char* errmsg,
                 void (*print_help)()) {
  float l;
  l = (float)atof(optarg);

  if (l < lower) {
    std::cerr << errmsg << std::endl;
    print_help();
    exit(1);
  }

  if (l > upper) {
    std::cerr << errmsg << std::endl;
    print_help();
    exit(1);
  }

  return l;

  std::cerr << errmsg << std::endl;
  print_help();
  exit(1);
  return -1;
}

template <typename Options>
int parse_options(int argc, char* argv[], Options& opt) {
  int option_index = 0;
  int next_option;
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
    switch (next_option) {
    case -1:
      break;
    case 'v':
      opt.verbosity() =
          parseInt(optarg, 0, "-v/--verbosity must be at least 0", print_help);
      break;
    case 'b':
      opt.bam_in() = optarg;
      break;
    case 'L':
      opt.target_region_in() = optarg;
      break;
    case 'r':
      opt.ref() = optarg;
      break;
    case 'c':
      opt.calculate_qual_stat() = true;
      break;
    case 'q':
      opt.min_mapq() =
          parseInt(optarg, 0, "-q/--min_mapq must be at least 0", print_help);
      break;
    case 'n':
      opt.snp_min_bq() =
          parseInt(optarg, 0, "-n/--snp_min_bq must be at least 0", print_help);
      break;
    case 'f':
      opt.vcf_out() = optarg;
      break;
    case 'd':
      opt.max_depth() =
          parseInt(optarg, 1, "-d/--max_depth must be at least 1", print_help);
      break;
    case 'h':
      opt.min_depth() =
          parseInt(optarg, 1, "-h/--min_depth must be at least 1", print_help);
      break;
    case 'o':
      opt.count_out() = optarg;
      break;
    case 'w':
      opt.window_size() =
          parseInt(optarg, 10, "-k/--window_size must be at least 10", print_help);
      break;
    case 'y':
      opt.fully_contained() = true;
      break;
    case 's':
      opt.include_secondary() = true;
      break;
    case 'D':
      opt.filter_duplicate() = true;
      break;
    case 'Q':
      opt.filter_QCfailed() = true;
      break;
    case 'E':
      opt.filter_improper_pair() = true;
      break;
    case 'F':
      opt.filter_mate_unmapped() = true;
      break;
    case 'G':
      opt.filter_improper_orientation() = true;
      break;
    case 'A':
      opt.report_all_alleles() = true;
      break;
    case 'C':
      opt.report_count_for_all_positions() = true;
      break;
    case 'a':
      opt.snp_min_allele_freq() = parseFloat(
          optarg, 0.0, 1.0, "-a/--snp_min_af must be between 0 and 1", print_help);
      break;
    case 'i':
      opt.ins_min_allele_freq() = parseFloat(
          optarg, 0.0, 1.0, "-i/--ins_min_af must be between 0 and 1", print_help);
      break;
    case 'l':
      opt.del_min_allele_freq() = parseFloat(
          optarg, 0.0, 1.0, "-l/--del_min_af must be between 0 and 1", print_help);
      break;
    case 'M':
      opt.snp_min_ao() =
          parseInt(optarg, 1, "-M/--snp_min_ao must be at least 1", print_help);
      break;
    case 't':
      // opt.num_threads() = parseInt(optarg, 1, "-t/--num_threads must be at least 1",
      // print_help);
      break;
    default:
      return 1;
    }
  } while (next_option != -1);

  return 0;
}

struct Options {
  Options(const int argc, char* argv[]) {
    std::string cmdline;
    for (int i = 0; i < argc; i++) {
      cmdline += argv[i];
      cmdline += " ";
    }
    std::cerr << "#" + cmdline << std::endl;

    int parse_ret = parse_options(argc, argv, *this);
    if (parse_ret) {
      print_help();
      exit(0);
    }

    if (this->bam_in_.empty()) {
      std::cerr << "no bam file is given. Exit..\n";
      print_help();
      exit(0);
    } else if (this->target_region_in_.empty()) {
      std::cerr << "no bed file is given. Exit..\n";
      print_help();
      exit(0);
    } else if (this->ref_.empty()) {
      std::cerr << "no reference file is given. Exit..\n";
      print_help();
      exit(0);
    } else if (this->vcf_out_.empty()) {
      std::cerr << "vcf output file required. Exit..\n";
      print_help();
      exit(0);
    } else if (this->count_out_.empty()) {
      std::cerr << "count output file required. Exit..\n";
      print_help();
      exit(0);
    }
  }

  decltype(auto) verbosity() { return (verbosity_); }

  decltype(auto) bam_in() { return (bam_in_); }

  decltype(auto) num_threads() { return (num_threads_); }

  // decltype(auto) bam_out() const {
  // return bam_out_;
  //}

  decltype(auto) vcf_out() { return (vcf_out_); }

  decltype(auto) count_out() { return (count_out_); }

  decltype(auto) target_region_in() { return (target_region_in_); }

  // const bool is_output_bam() const {
  // if (bam_out_.empty()) return false;
  // else return true;
  //}

  decltype(auto) ref() { return (ref_); }

  decltype(auto) min_mapq() const { return (min_mapq_); }

  decltype(auto) snp_min_bq() const { return (snp_min_bq_); }

  decltype(auto) min_mapq() { return (min_mapq_); }

  decltype(auto) snp_min_bq() { return (snp_min_bq_); }

  decltype(auto) calculate_qual_stat() { return (calculate_qual_stat_); }

  decltype(auto) window_size() { return (window_size_); }

  decltype(auto) snp_min_allele_freq() { return (snp_min_allele_freq_); }

  decltype(auto) ins_min_allele_freq() { return (ins_min_allele_freq_); }

  decltype(auto) del_min_allele_freq() { return (del_min_allele_freq_); }

  decltype(auto) snp_min_ao() { return (snp_min_ao_); }

  decltype(auto) fully_contained() { return (fully_contained_); }

  decltype(auto) max_depth() { return (max_depth_); }

  decltype(auto) min_depth() { return (min_depth_); }

  decltype(auto) include_secondary() { return (include_secondary_); }

  decltype(auto) include_secondary() const { return (include_secondary_); }

  decltype(auto) filter_duplicate() const { return (filter_duplicate_); }

  decltype(auto) filter_duplicate() { return (filter_duplicate_); }

  decltype(auto) filter_QCfailed() const { return (filter_QCfailed_); }

  decltype(auto) filter_QCfailed() { return (filter_QCfailed_); }

  decltype(auto) filter_improper_pair() const { return (filter_improper_pair_); }

  decltype(auto) filter_improper_pair() { return (filter_improper_pair_); }

  decltype(auto) filter_mate_unmapped() const { return (filter_mate_unmapped_); }

  decltype(auto) filter_mate_unmapped() { return (filter_mate_unmapped_); }

  decltype(auto) filter_improper_orientation() const {
    return (filter_improper_orientation_);
  }

  decltype(auto) filter_improper_orientation() { return (filter_improper_orientation_); }

  decltype(auto) report_all_alleles() const { return (report_all_alleles_); }

  decltype(auto) report_all_alleles() { return (report_all_alleles_); }

  decltype(auto) report_count_for_all_positions() const {
    return (report_count_for_all_positions_);
  }

  decltype(auto) report_count_for_all_positions() {
    return (report_count_for_all_positions_);
  }

private:
  unsigned verbosity_ = 0;
  std::string bam_in_;
  std::string bam_out_;
  std::string vcf_out_;
  std::string count_out_;
  std::string target_region_in_;
  std::string ref_;
  bool calculate_qual_stat_ = false;
  bool fully_contained_ = false;
  float snp_min_allele_freq_ = 0.01;
  float ins_min_allele_freq_ = 0.01;
  float del_min_allele_freq_ = 0.01;
  int snp_min_ao_ = 3;
  int min_mapq_ = 0;
  int snp_min_bq_ = 0;
  int window_size_ = 500;
  int num_threads_ = 1;
  int max_depth_ = 5000000;
  int min_depth_ = 0;
  bool include_secondary_ = false;
  bool filter_duplicate_ = false;
  bool filter_QCfailed_ = false;
  bool filter_improper_pair_ = false;
  bool filter_mate_unmapped_ = false;
  bool filter_improper_orientation_ = false;
  bool report_all_alleles_ = false;
  bool report_count_for_all_positions_ = false;
};
}  // namespace neusomatic

#endif /* LIB_INCLUDE_OPTIONS_H_ */
