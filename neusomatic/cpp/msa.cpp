#include <iostream>
#include <iostream>
#include <fstream>

#include <seqan/store.h>
#include <seqan/realign.h>

#include <msa_options.h>

void MSA(neusomatic::Options const& opt)
{    
    using namespace seqan;
    CharString seqFileName = opt.input().c_str();
    SeqFileIn seqFileIn(toCString(seqFileName));
    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;
    try
    {
        readRecords(ids, seqs, seqFileIn);
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return;
    }
    int n = length(ids);
    Align<Dna5String> align;
    resize(rows(align), n);
    for (int i = 0; i < n; ++i){
        assignSource(row(align, i), seqs[i]);
    }

    globalMsaAlignment(align, SimpleScore(opt.match_score(), -opt.mismatch_score(), -opt.gap_extension_score(), -opt.gap_open_score()));

    std::ofstream outfile;
    outfile.open(opt.output().c_str());
    for (int i = 0; i < n; ++i){
	    outfile << ">" << ids[i] << "\n";
	    outfile << row(align, i) << "\n";
	}
    outfile.close();
}

int main(int argc, char *argv []) {
  neusomatic::Options opt(argc, argv);
  MSA(opt);
  return 0;
}
