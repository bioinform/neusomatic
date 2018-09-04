/*
 * Options.h
 *
 *  Created on: Feb 20, 2017
 *      Author: liur34
 */

#ifndef SEQAN_MSA_OPTIONS_H_
#define SEQAN_MSA_OPTIONS_H_

#include <limits>
#include <string>
#include <thread>
#include <getopt.h>

namespace neusomatic {

  static struct option  long_options[] = {
    {"input",                           required_argument,      0,        'i'},
    {"output",                          required_argument,      0,        'o'},
    {"match_score",                     required_argument,      0,        'A'},
    {"mismatch_penalty",                required_argument,      0,        'B'},
    {"gap_open_penalty",                required_argument,      0,        'O'},
    {"gap_ext_penalty",                 required_argument,      0,        'E'},
    {0, 0, 0, 0} // terminator
  };

  const char *short_options = "i:o:A:B:O:E:";

  void print_help()
  {
    std::cerr<< "---------------------------------------------------\n";
    std::cerr<< "Usage: msa [options] -i INPUT_FASTA -o OUTPUT_FASTA\n";
    std::cerr<< "General Options:\n";

    std::cerr<< "Required Options:\n";
    std::cerr<< "-i/--input,                           input fasta file path,                                                               required.\n";
    std::cerr<< "-o/--output,                          output fasta file path,                                                              required.\n";
    std::cerr<< "\n";
    std::cerr<< "Other Options:\n";
    std::cerr<< "-A/--match_score                      match score,                                                                         default is 10\n";
    std::cerr<< "-B/--mismatch_penalty                 penalty for having a mismatch,                                                       default is 8\n";
    std::cerr<< "-O/--gap_open_penalty                 penalty for opening a gap,                                                           default is 8\n";
    std::cerr<< "-E/--gap_ext_penalty                  penalty for extending a gap,                                                         default is 6\n";
  }

  int parseInt(const char* optarg, int lower, const char *errmsg, void (*print_help)()) {
    long l;
    char *endPtr= NULL;
    l = strtol(optarg, &endPtr, 10);
    if (endPtr != NULL) {
        if (l < lower ) {
            std::cerr << errmsg << std::endl;
            print_help();
            exit(1);
        }
        return (int32_t)l;
    }
    std::cerr<< errmsg <<std::endl;
    print_help();
    exit(1);
    return -1;
  }

  float parseFloat(const char* optarg, float lower, float upper, const char *errmsg, void (*print_help)()) {
    float l;
    l = (float)atof(optarg);

    if (l < lower) {
        std::cerr << errmsg <<std::endl;
        print_help();
        exit(1);
    }

    if (l > upper)
    {
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

  template<typename Options>
  int parse_options(int argc, char* argv[], Options& opt) {
    int option_index = 0;
    int next_option;
    do {
      next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
      switch(next_option) {
        case -1:
          break;
        case 'i':
          opt.input() = optarg;
          break;
        case 'o':
          opt.output() = optarg;
          break;
        case 'A':
          opt.match_score() = parseInt(optarg, 0, "match score must >= 0", print_help); 
          break;
        case 'B':
          opt.mismatch_score() = parseInt(optarg, 0, "mismatch penalty must >= 0", print_help); 
          break;
        case 'O':
          opt.gap_open_score() = parseInt(optarg, 0, "gap open penalty must >= 0", print_help); 
          break;
        case 'E':
          opt.gap_extension_score() = parseInt(optarg, 0, "gap extension penalty must >= 0", print_help); 
          break;
        default:
          print_help();
          return 1;
      }
    } while (next_option != -1);
    return 0;
  }

struct Options {
  Options(const int argc, char* argv[]) {

    std::string cmdline;
    for(int i=0; i<argc; i++){
      cmdline += argv[i];
      cmdline+=" ";
    }
    std::cerr << "#" + cmdline << std::endl;

    int parse_ret =  parse_options(argc, argv, *this);
    if(parse_ret){
       exit(0);
    }

    if (this->input_.empty()) {
      std::cerr << "no input file is given. Exit..\n";
      print_help();
      exit(0);
    } else if (this->output_.empty()) {
      std::cerr << "output file required. Exit..\n";
      print_help();
      exit(0);
    }
  }


  decltype(auto) input() {
    return (input_);
  }
  decltype(auto) output() {
    return (output_);
  }
  decltype(auto) match_score() {
    return (match_score_);
  }

  decltype(auto) mismatch_score() {
    return (mismatch_penalty_);
  }

  decltype(auto) gap_open_score() {
    return (gap_open_penalty_);
  }

  decltype(auto) gap_extension_score() {
    return (gap_ext_penalty_);
  }

  decltype(auto) input() const{
    return (input_);
  }
  decltype(auto) output() const{
    return (output_);
  }
  decltype(auto) match_score() const{
    return (match_score_);
  }

  decltype(auto) mismatch_score() const{
    return (mismatch_penalty_);
  }

  decltype(auto) gap_open_score() const{
    return (gap_open_penalty_);
  }

  decltype(auto) gap_extension_score() const{
    return (gap_ext_penalty_);
  }
private:
  std::string input_;
  std::string output_;
  int match_score_ = 10;
  int mismatch_penalty_ = 8;
  int gap_open_penalty_ = 8; 
  int gap_ext_penalty_ = 6;

};
}//namespace neusomatic



#endif /* LIB_INCLUDE_OPTIONS_H_ */
