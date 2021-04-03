/** 
  * File Name :
  * Purpose :
  * Creation Date : 28-08-2017
  * Last Modified :
  * Created By : Mohammad Sahraeian  
  */





#include "SeqLib/BamReader.h"
#include "SeqLib/RefGenome.h"

#ifdef BGZF_MAX_BLOCK_SIZE
#pragma push_macro("BGZF_MAX_BLOCK_SIZE")
#undef BGZF_MAX_BLOCK_SIZE
#define BGZF_MAX_BLOCK_SIZE_BAK
#endif

#ifdef BGZF_BLOCK_SIZE
#pragma push_macro("BGZF_BLOCK_SIZE")
#undef BGZF_BLOCK_SIZE
#define BGZF_BLOCK_SIZE_BAK
#endif



#include "vcf.h"
#include "bedio.hpp"
#include "targeted_loading.hpp"
#include "msa.hpp"
#include "msa_utils.hpp"
#include "Options.h"

#ifdef BGZF_MAX_BLOCK_SIZE_BAK
#undef BGZF_MAX_BLOCK_SIZE_BAK
#pragma pop_macro("BGZF_MAX_BLOCK_SIZE")
#endif

#ifdef BGZF_BLOCK_SIZE_BAK
#undef BGZF_BLOCK_SIZE_BAK
#pragma pop_macro("BGZF_BLOCK_SIZE")
#endif


int main(int argc, char **argv) {
  neusomatic::Options opts(argc, argv);
  const std::string& bam_path = opts.bam_in();
  const std::string& bed_in = opts.target_region_in();
  const std::string& ref_in = opts.ref();
  const std::string& vcf_out = opts.vcf_out();
  const std::string& count_out = opts.count_out();
  int window_size = opts.window_size();
  float min_af = opts.min_allele_freq();
  float ins_min_af=min_af;
  float del_min_af=min_af;
  float snp_min_af=min_af;
  const bool calculate_qual_stat = opts.calculate_qual_stat(); 
  const bool report_all_alleles = opts.report_all_alleles(); 
  const bool report_count_for_all_positions = opts.report_count_for_all_positions(); 



  auto start_main = std::chrono::system_clock::now();

  //const std::map<char, int> empty_pileup_counts = {{'-', 0}, {'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};
  static const std::vector<char> nuc_code_char = {'A', 'C', 'G', 'T', '-', 'N'};
  const int matrix_base_pad = 7;

  using GInvStdString = neusomatic::bio::GenomicInterval<std::string>;
  using GInvInt = neusomatic::bio::GenomicInterval<int>;
  using MSA = typename neusomatic::MSABuilder<SeqLib::BamRecord, GInvInt>;
  using ContigGaps = typename MSA::ContigGaps;
  neusomatic::BedIO<GInvStdString> bed_reader(bed_in);
  std::vector<GInvStdString> bed_regions = bed_reader.ReadBed3_windowed(window_size);
  SeqLib::RefGenome ref_seqs;
  ref_seqs.LoadIndex(ref_in);

  neusomatic::bio::VCFWriter vcf_writer(vcf_out, ref_in, "VarCal");
  const unsigned contig_counts = seqan::length(seqan::contigNames(vcf_writer.vcf_context()));
  std::map<std::string, unsigned> chr_name_to_idx;
  for (unsigned i = 0; i < contig_counts; ++i) {
    chr_name_to_idx[seqan::toCString(seqan::contigNames(vcf_writer.vcf_context())[i])] = i; 
  }
  std::ofstream count_out_writer;
  count_out_writer.open(count_out);
  std::string count_out_str = "";

  neusomatic::CaptureLayout<GInvInt, SeqLib::BamRecord, SeqLib::BamReader> capture_layout(bam_path, bed_regions, opts);
  SeqLib::BamHeader bam_header = capture_layout.Reader().Header(); 
  
  try {
    for (auto it = chr_name_to_idx.cbegin(); it != chr_name_to_idx.cend(); ++it) {
      if (bam_header.Name2ID(it->first) != it->second) {
        std::cerr << " the bam file and the reference do not match.\n exit..\n please check the bam header and the reference file.\n";
        exit(1);
      }
    }
  } catch (const std::out_of_range& oor) {
        std::cerr << " the reference contains chromosome/contig name(s) that is/are not in the bam file.\n exit..\n please check the bam header and the reference file.\n";
    exit(1);
  } catch (const std::invalid_argument& ia) {
        std::cerr << " the reference contains chromosome/contig name(s) that is/are not in the bam file.\n exit..\n please check the bam header and the reference file.\n";
    exit(1);
  }

  auto end_main_0 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds_main_0 = end_main_0-start_main;
  std::time_t end_time_main_0 = std::chrono::system_clock::to_time_t(end_main_0);
  std::cout << "elapsed_main_0 time: " << elapsed_seconds_main_0.count() << "s\n";      


  int last_var_pos=-1;

  int cnt_region=0;
  auto timestamp = opts.verbosity() > 0;
  while (capture_layout.NextRegion(opts.fully_contained())) {

    auto start_0_while = std::chrono::system_clock::now();

    // a map from genomic interval -> a vector of alignReadIds
    for (auto targeted_region : capture_layout.ginv_records()) {


      auto start_0_for = std::chrono::system_clock::now();
      auto start_0 = std::chrono::system_clock::now();

      auto ginv = targeted_region.first;
      auto& records = targeted_region.second;
      GInvStdString ginvstr(bam_header.IDtoName(ginv.contig()), ginv.left(), ginv.right());

      if (opts.verbosity() > 2 || (cnt_region % 100 == 0 || timestamp)) {
        std::cerr<<"#On region "<<ginvstr<<"\n";
        std::cerr<<"#Aligned read number: "<<records.size()<<std::endl;

      }
      ++cnt_region;
      if (records.empty()) continue; 
      if (records.size() > opts.max_depth()) {
        records.resize(opts.max_depth());
      }
      const auto non_gapped_ref = ref_seqs.QueryRegion(ginvstr.contig(), ginvstr.left(), ginvstr.right() - 1);
      if (ginv.length() > non_gapped_ref.size())  {
        ginv.right() = ginv.left() + non_gapped_ref.size();
      }

      auto end_0 = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds_0 = end_0-start_0;
      std::time_t end_time_0 = std::chrono::system_clock::to_time_t(end_0);
      if (timestamp){
        std::cout << "elapsed_0 time: " << elapsed_seconds_0.count() << "s\n";      
      }
      auto start_1 = std::chrono::system_clock::now();

      neusomatic::CondensedArray<int> condensed_array;

      auto end_1 = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds_1 = end_1-start_1;
      std::time_t end_time_1 = std::chrono::system_clock::to_time_t(end_1);
      if (timestamp){
        std::cout << "elapsed_1 time: " << elapsed_seconds_1.count() << "s\n";      
      }
      auto start_2 = std::chrono::system_clock::now();

      MSA msa(ginv, records, non_gapped_ref);

      auto end_2 = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds_2 = end_2-start_2;
      std::time_t end_time_2 = std::chrono::system_clock::to_time_t(end_2);
      if (timestamp){
        std::cout << "elapsed_2 time: " << elapsed_seconds_2.count() << "s\n";      
      }
      auto start_3 = std::chrono::system_clock::now();
      
      condensed_array = neusomatic::CondensedArray<int>(msa, calculate_qual_stat);

      
      auto end_3 = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds_3 = end_3-start_3;
      std::time_t end_time_3 = std::chrono::system_clock::to_time_t(end_3);
      if (timestamp){
        std::cout << "elapsed_3 time: " << elapsed_seconds_3.count() << "s\n";      
      }
      // auto start_4 = std::chrono::system_clock::now();
      
      // auto cols = condensed_array.GetColSpace();
      // auto cols = condensed_array.GetColSpaceMQual();
      // auto cols = condensed_array.GetColSpaceStrand();
      // auto cols = condensed_array.GetColSpaceLSC();
      // auto cols = condensed_array.GetColSpaceRSC();
      // auto cols = condensed_array.GetColSpaceTag();
      // auto ncol = cols.size();
      // const auto ref = condensed_array.GetGappedRef();
      // const auto cc = neusomatic::ChangeCoordinates(ref);


      // auto end_4 = std::chrono::system_clock::now();
      // std::chrono::duration<double> elapsed_seconds_4 = end_4-start_4;
      // std::time_t end_time_4 = std::chrono::system_clock::to_time_t(end_4);
      // if (timestamp){
      //   std::cout << "elapsed_4 time: " << elapsed_seconds_4.count() << "s\n";      
      // }

      auto start_4 = std::chrono::system_clock::now();
      
      auto ncol = condensed_array.ncol();
      const auto ref = condensed_array.GetGappedRef();
      const auto cc = neusomatic::ChangeCoordinates(ref);

      auto end_4 = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds_4 = end_4-start_4;
      std::time_t end_time_4 = std::chrono::system_clock::to_time_t(end_4);
      if (timestamp){
        std::cout << "elapsed_4 time: " << elapsed_seconds_4.count() << "s\n";      
      }

      auto start_5 = std::chrono::system_clock::now();

      std::vector<int> var_cols;

      if (( cnt_region == 1 ) or (last_var_pos >= 0)) {
        var_cols.push_back(0);
      }

      last_var_pos = -1;

      for (size_t i = 0; i < ncol; i++) {

        // auto start_4p0 = std::chrono::system_clock::now();

        if (ginv.left() + cc.UngapPos(i) >=ginv.right()){ break;}
        auto ref_base = ref[i];
        ref_base = std::toupper(ref_base);
        auto ref_code = neusomatic::CondensedArray<int>::DnaCharToDnaCode(ref_base);

        if (ref_base == 'N') {
          ref_base = '-';
        }

        auto base_freq_ = condensed_array.GetBaseFreq(i);

        if (opts.verbosity()>1){
          std::cout<<"col "<<i<<": ";
          std::cout<<"(ref= "<< ref_base << ") ";
          for (int base = 0; base < (int) base_freq_.size(); ++base) {
            std::cout<<"("<<nuc_code_char[base] <<"): "<< base_freq_[base] <<"; ";
          }
          std::cout<< std::endl;
        }


        auto nrow = condensed_array.nrow()-base_freq_[5];
        base_freq_.erase(base_freq_.begin() + 5);

        std::vector<int> pileup_counts(base_freq_.size());
        int total_count=0;
        for (int base = 0; base < (int) base_freq_.size(); ++base) {
          pileup_counts[base] = base_freq_[base];
          total_count+=base_freq_[base];
        }


        if (total_count==0) {continue;}

        auto start_pos=ginv.left() + cc.UngapPos(i);
        if (ref_base!='-'){
          start_pos++;
        }

        // auto end_4p0 = std::chrono::system_clock::now();
        // std::chrono::duration<double> elapsed_seconds_4p0 = end_4p0-start_4p0;
        // std::time_t end_time_4p0 = std::chrono::system_clock::to_time_t(end_4p0);
        // std::cout << "elapsed_4p0 time: " << elapsed_seconds_4p0.count() << "s\n";      
        // auto start_4p1 = std::chrono::system_clock::now();


        std::map<int, int> alt_counts;
        auto ref_count = base_freq_[ref_code];
        auto var_code = ref_code; 
        int var_count = 0;
        int dp = ref_count;
        if (report_all_alleles and ref_base != '-'){
          for (int row = 0;  row < base_freq_.size(); ++row) {
            auto alt_cnt = base_freq_[row];
            if (( row != ref_code) and (alt_cnt > 0)){
              auto af = alt_cnt/float(alt_cnt+ref_count);
              if ((alt_cnt >= ref_count) or ((row == 4 and  af > del_min_af ) or
                                              (row != 4 and ref_base != '-' and af > snp_min_af ) or
                                              (ref_base =='-' and af > ins_min_af))){
                alt_counts.insert(std::pair<int, int>(row, alt_cnt));
                dp += alt_cnt;
              }
            }
          }
        }else{
          int major = -1;
          int major_count = 0;
          int minor = -1;
          int minor_count = 0;
          int minor2 = -1;
          int minor2_count = 0;

          for (int row = 0;  row < base_freq_.size(); ++row) {
            if (base_freq_[row] > major_count) {
              minor2 = minor;
              minor2_count = minor_count;
              minor_count = major_count;
              minor = major;
              major_count = base_freq_[row];
              major = row;
            } else if (base_freq_[row] > minor_count) {
              minor2 = minor;
              minor2_count = minor_count;
              minor_count = base_freq_[row];
              minor = row;
            } else if (base_freq_[row] > minor2_count) {
              minor2_count = base_freq_[row];
              minor2 = row;
            }
          }

          if (minor != -1 and major != -1){
            if (minor2 != -1 and ref_code == major and minor == 4 and ref_code != 4 ){
              if (minor2_count>0.5*minor_count){
                minor = minor2;
                minor_count = minor2_count;
              }
            }
          }
          auto af = minor_count/float(major_count+minor_count);
          if (major != ref_code){
            var_code = major;
            var_count = major_count;
          } else if (minor != ref_code and ( (minor == 4 and  af > del_min_af ) or
                                          (minor != 4 and ref_base != '-' and af > snp_min_af ) or
                                          (ref_base =='-' and af > ins_min_af))){
            var_code = minor;
            var_count = minor_count;
          }
          if (var_count > 0) { 
            alt_counts.insert(std::pair<int, int>(var_code,var_count));
            dp += var_count;
          }
        }
        // for(auto it = alt_counts.cbegin(); it != alt_counts.cend(); ++it)
        // {
        //     std::cout << it->first << " " << it->second << std::endl;
        // }          



        // auto end_4p1 = std::chrono::system_clock::now();
        // std::chrono::duration<double> elapsed_seconds_4p1 = end_4p1-start_4p1;
        // std::time_t end_time_4p1 = std::chrono::system_clock::to_time_t(end_4p1);
        // std::cout << "elapsed_4p1 time: " << elapsed_seconds_4p1.count() << "s\n";      
        // auto start_4p2 = std::chrono::system_clock::now();

        if ( alt_counts.size() > 0 ){
          var_cols.push_back(i);
          last_var_pos = ginv.left() + cc.UngapPos(i);
        }

        for(auto it = alt_counts.cbegin(); it != alt_counts.cend(); ++it)
        {
          auto var_code_ = it->first;
          auto var_count_ = it->second;
          auto record_info = "AF="+std::to_string((var_count_)/float(dp))+";DP="+std::to_string(nrow)+";RO="+std::to_string(ref_count)+";AO="+std::to_string(var_count_);
          auto gtinfo = "0/1:"+std::to_string(nrow)+":"+std::to_string(ref_count)+":"+std::to_string(var_count_);
          if (calculate_qual_stat){

            auto bqual_mean_ = condensed_array.GetBQMean(i);
            auto mqual_mean_ = condensed_array.GetMQMean(i);
            auto strand_mean_ = condensed_array.GetStrandMean(i);
            auto lsc_mean_ = condensed_array.GetLSCMean(i);
            auto rsc_mean_ = condensed_array.GetRSCMean(i);
            auto tag_mean_ = condensed_array.GetTagMean(i);
            
            int rsc_counts=0;
            int lsc_counts=0;
            for(auto it = lsc_mean_.cbegin(); it != lsc_mean_.cend(); ++it) {
              rsc_counts+=*it;
            }
            for(auto it = rsc_mean_.cbegin(); it != rsc_mean_.cend(); ++it) {
              lsc_counts+=*it;
            }


            record_info += ";ST="+std::to_string(int(round(ref_count*(strand_mean_[ref_code]/100))))+ \
                           ","+std::to_string(int(round(var_count_*(strand_mean_[var_code_]/100))))+ \
                           ";LS="+std::to_string(lsc_counts)+\
                           ";RS="+std::to_string(rsc_counts)+\
                           ";NM="+std::to_string(int(round(tag_mean_[var_code_][0])))+\
                           ";AS="+std::to_string(int(round(tag_mean_[var_code_][1])))+ \
                           ";XS="+std::to_string(int(round(tag_mean_[var_code_][2])))+ \
                           ";PR="+std::to_string(int(round(tag_mean_[var_code_][3])))+ \
                           ";CL="+std::to_string(int(round(tag_mean_[var_code_][4])))+ \
                           ";MQ="+std::to_string(int(round(mqual_mean_[var_code_])))+ \
                           ";BQ="+std::to_string(int(round(bqual_mean_[var_code_])));
            gtinfo += ":"+std::to_string(int(round(ref_count*(strand_mean_[ref_code]/100))))+","+ \
                      std::to_string(int(round(var_count_*(strand_mean_[var_code_]/100))))+":"+\
                      std::to_string(lsc_counts)+":"+\
                      std::to_string(rsc_counts)+":"+\
                      std::to_string(int(round(tag_mean_[var_code_][0])))+":"+\
                      std::to_string(int(round(tag_mean_[var_code_][1])))+":"+\
                      std::to_string(int(round(tag_mean_[var_code_][2])))+":"+\
                      std::to_string(int(round(tag_mean_[var_code_][3])))+":"+\
                      std::to_string(int(round(tag_mean_[var_code_][4])))+":"+\
                      std::to_string(int(round(mqual_mean_[var_code_])))+":"+\
                      std::to_string(int(round(bqual_mean_[var_code_])));
          }
          auto var_base = nuc_code_char[var_code_];  
          if (ref_base == '-') {ref_base = 'N';}
          if (var_base == '-') {var_base = 'N';}
          auto var_ref_pos=ginv.left() + cc.UngapPos(i);
          seqan::VcfRecord record;
          record.rID = ginv.contig();
          record.beginPos = var_ref_pos;
          record.id = ".";
          record.ref = ref_base;
          record.alt = var_base;
          record.qual = 100;
          record.filter = ".";
          record.info = record_info;
          char * info;
          if (calculate_qual_stat){
            record.format = "GT:DP:RO:AO:ST:LS:RS:NM:AS:XS:PR:CL:MQ:BQ";
          }else{
            record.format = "GT:DP:RO:AO";
          }
          appendValue(record.genotypeInfos, gtinfo);
          vcf_writer.Write(record);
          if (opts.verbosity()>1){
            std::cout<<"var: " << i << "," << var_ref_pos << ","<< ref_base << "," << var_base<<","<<nrow<<":"<<ref_count<<":"<<var_count_<<std::endl;
            std::cout<<"col "<<i<<": ";
            std::cout<<"(ref= "<< ref_base << ") ";
            for (size_t row = 0; row < base_freq_.size(); ++row) {
              std::cout<<"("<<row<<"): "<<base_freq_[row]<<",";
            }
            std::cout<< std::endl;
          }
        }

        // auto end_4p2 = std::chrono::system_clock::now();
        // std::chrono::duration<double> elapsed_seconds_4p2 = end_4p2-start_4p2;
        // std::time_t end_time_4p2 = std::chrono::system_clock::to_time_t(end_4p2);
        // std::cout << "elapsed_4p2 time: " << elapsed_seconds_4p2.count() << "s\n";      

      }

      if (last_var_pos>=0){
        if (ginv.right()>(last_var_pos+matrix_base_pad+2)){
          last_var_pos=-1;
        }
      }

      var_cols.push_back(std::max(0,int(ncol)-1));
                // for (auto i=std::max(processed_col+1,j-matrix_base_pad-1); i < std::min(int(ncol),j+matrix_base_pad+1); ++i){

      std::vector<int> count_cols;
      if (report_count_for_all_positions){
        for (size_t i = 0; i < ncol; i++) {
          count_cols.push_back(i);
        }
      }else{
        int processed_col = -1;
        for (auto& j : var_cols){
          // std::cout << j << "," << processed_col << std::endl;
          std::vector<int> count_cols_0;
          int cnt_0 = 0;
          for (auto i=j; i >= (processed_col+1); --i){
            if (ginv.left() + cc.UngapPos(i) >=ginv.right()){ break;}
            auto ref_base = ref[i];
            ref_base = std::toupper(ref_base);
            auto ref_code = neusomatic::CondensedArray<int>::DnaCharToDnaCode(ref_base);
            if (ref_base == 'N') {
              ref_base = '-';
            }
            count_cols_0.push_back(i);
            if (ref_base == '-') {
              continue;
            }
            cnt_0++;
            if (cnt_0 >= matrix_base_pad+2){
              break;
            }
          }
          std::reverse(count_cols_0.begin(),count_cols_0.end());
          count_cols.insert(count_cols.end(), count_cols_0.begin(), count_cols_0.end());

          std::vector<int> count_cols_1;
          int cnt_1 = 0;
          for (auto i=j+1; i < ncol; ++i){
            if (ginv.left() + cc.UngapPos(i) >=ginv.right()){ break;}
            auto ref_base = ref[i];
            ref_base = std::toupper(ref_base);
            auto ref_code = neusomatic::CondensedArray<int>::DnaCharToDnaCode(ref_base);
            if (ref_base == 'N') {
              ref_base = '-';
            }
            if (i >= (processed_col+1)){
              count_cols_1.push_back(i);
            }
            if (ref_base == '-') {
              continue;
            }
            cnt_1++;
            if (cnt_1 >= matrix_base_pad+1){
              break;
            }
          }
          count_cols.insert(count_cols.end(), count_cols_1.begin(), count_cols_1.end());          
          if (! count_cols.empty()) {
            processed_col = count_cols.back();
            // std::cout << ginv.left() + cc.UngapPos(processed_col) << std::endl;
          }
          // for (auto& iii : count_cols_0){
          //   std::cout << "count_cols_0, " << iii << std::endl;
          // }
          // for (auto& iii : count_cols_1){
          //   std::cout << "count_cols_1, " << iii << std::endl;
          // }

        }
      }

      auto end_5 = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds_5 = end_5-start_5;
      std::time_t end_time_5 = std::chrono::system_clock::to_time_t(end_5);
      if (timestamp){
        std::cout << "elapsed_5 time: " << elapsed_seconds_5.count() << "s\n";      
      }
      auto start_6 = std::chrono::system_clock::now();

      for (auto& i : count_cols){
        auto ref_base = ref[i];
        ref_base = std::toupper(ref_base);
        auto ref_code = neusomatic::CondensedArray<int>::DnaCharToDnaCode(ref_base);

        if (ref_base == 'N') {
          ref_base = '-';
        }
        auto start_pos=ginv.left() + cc.UngapPos(i);
        if (ref_base!='-'){
          start_pos++;
        }

        auto base_freq_ = condensed_array.GetBaseFreq(i);

        std::vector<int> pileup_counts(base_freq_.size());
        int total_count=0;
        for (int base = 0; base < (int) base_freq_.size(); ++base) {
          pileup_counts[base] = base_freq_[base];
          total_count+=base_freq_[base];
        }


        if (total_count==0) {continue;}

        if (calculate_qual_stat){ 
  
          auto bqual_mean_ = condensed_array.GetBQMean(i);
          auto mqual_mean_ = condensed_array.GetMQMean(i);
          auto strand_mean_ = condensed_array.GetStrandMean(i);
          auto lsc_mean_ = condensed_array.GetLSCMean(i);
          auto rsc_mean_ = condensed_array.GetRSCMean(i);
          auto tag_mean_ = condensed_array.GetTagMean(i);
          

          // count_out_writer<<bam_header.IDtoName(ginv.contig())<<"\t"<<start_pos<<"\t" \
          // << start_pos+1<<"\t"  << ref_base << "\t" \
          // << neusomatic::add_qual_col(pileup_counts, true)<<"\t" \
          // << neusomatic::add_qual_col(bqual_mean_)<<"\t" \
          // << neusomatic::add_qual_col(mqual_mean_)<<"\t" \
          // << neusomatic::add_qual_col(strand_mean_)<<"\t" \
          // << neusomatic::add_qual_col(lsc_mean_)<<"\t" \
          // << neusomatic::add_qual_col(rsc_mean_)<<"\t" \
          // << neusomatic::add_tag_col(tag_mean_, false, 0)<<"\t" \
          // << neusomatic::add_tag_col(tag_mean_, false, 1)<<"\t" \
          // << neusomatic::add_tag_col(tag_mean_, false, 2)<<"\t" \
          // << neusomatic::add_tag_col(tag_mean_, false, 3)<<"\t" \
          // << neusomatic::add_tag_col(tag_mean_, false, 4) \
          // << std::endl;


          count_out_writer<<bam_header.IDtoName(ginv.contig())<<"\t"<<start_pos<<"\t" \
          << start_pos+1<<"\t"  << ref_base << "\t";
          neusomatic::add_qual_col(count_out_writer, pileup_counts, true, "\t");
          neusomatic::add_qual_col(count_out_writer, bqual_mean_, false, "\t");
          neusomatic::add_qual_col(count_out_writer, mqual_mean_, false, "\t");
          neusomatic::add_qual_col(count_out_writer, strand_mean_, false, "\t");
          neusomatic::add_qual_col(count_out_writer, lsc_mean_, false, "\t");
          neusomatic::add_qual_col(count_out_writer, rsc_mean_, false, "\t");
          neusomatic::add_tag_col(count_out_writer, tag_mean_, false, 0, "\t");
          neusomatic::add_tag_col(count_out_writer, tag_mean_, false, 1, "\t");
          neusomatic::add_tag_col(count_out_writer, tag_mean_, false, 2, "\t");
          neusomatic::add_tag_col(count_out_writer, tag_mean_, false, 3, "\t");
          neusomatic::add_tag_col(count_out_writer, tag_mean_, false, 4, "\n");

          // count_out_str += bam_header.IDtoName(ginv.contig())+"\t"+std::to_string(start_pos)+"\t"+std::to_string(start_pos+1)+"\t"+ref_base+"\t";
          // count_out_str += neusomatic::add_qual_col(pileup_counts, true, "\t");
          // count_out_str += neusomatic::add_qual_col(bqual_mean_, false, "\t");
          // count_out_str += neusomatic::add_qual_col(mqual_mean_, false, "\t");
          // count_out_str += neusomatic::add_qual_col(strand_mean_, false, "\t");
          // count_out_str += neusomatic::add_qual_col(lsc_mean_, false, "\t");
          // count_out_str += neusomatic::add_qual_col(rsc_mean_, false, "\t");
          // count_out_str += neusomatic::add_tag_col(tag_mean_, false, 0, "\t");
          // count_out_str += neusomatic::add_tag_col(tag_mean_, false, 1, "\t");
          // count_out_str += neusomatic::add_tag_col(tag_mean_, false, 2, "\t");
          // count_out_str += neusomatic::add_tag_col(tag_mean_, false, 3, "\t");
          // count_out_str += neusomatic::add_tag_col(tag_mean_, false, 4, "\n");


          // std::cout<<bam_header.IDtoName(ginv.contig())<<"\t"<<start_pos<<std::endl;
          // std::cout<< start_pos+1<<"\t"  << ref_base << std::endl;
          // std::cout<< neusomatic::add_qual_col(pileup_counts, true)<<std::endl;
          // std::cout<< neusomatic::add_qual_col(bqual_mean_)<<std::endl;
          // std::cout<< neusomatic::add_qual_col(mqual_mean_)<<std::endl;
          // std::cout<< neusomatic::add_qual_col(strand_mean_)<<std::endl;
          // std::cout<< neusomatic::add_qual_col(lsc_mean_)<<std::endl;
          // std::cout<< neusomatic::add_qual_col(rsc_mean_)<<std::endl;
          // std::cout<< neusomatic::add_tag_col(tag_mean_, false, 0)<<std::endl;
          // std::cout<< neusomatic::add_tag_col(tag_mean_, false, 1)<<std::endl;
          // std::cout<< neusomatic::add_tag_col(tag_mean_, false, 2)<<std::endl;
          // std::cout<< neusomatic::add_tag_col(tag_mean_, false, 3)<<std::endl;
          // std::cout<< neusomatic::add_tag_col(tag_mean_, false, 4)<< std::endl;


        }else{
          count_out_writer<<bam_header.IDtoName(ginv.contig())<<"\t"<<start_pos<<"\t" \
          << start_pos+1<<"\t"  << ref_base << "\t" \
          <<pileup_counts[4]<<":"<<pileup_counts[0]<<":"<<pileup_counts[1]<<":"<<pileup_counts[2] \
          <<":"<<pileup_counts[3]<<std::endl;
            // count_out_str += bam_header.IDtoName(ginv.contig())+"\t"+std::to_string(start_pos)+"\t"+std::to_string(start_pos+1)+"\t"+ref_base+"\t";
          // count_out_str += neusomatic::add_qual_col(pileup_counts, true, "\n");
        }


      }


      auto end_6 = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds_6 = end_6-start_6;
      std::time_t end_time_6 = std::chrono::system_clock::to_time_t(end_6);
      if (timestamp){
        std::cout << "elapsed_6 time: " << elapsed_seconds_6.count() << "s\n";      
      }

      auto end_0_for = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds_0_for = end_0_for-start_0_for;
      std::time_t end_time_0_for = std::chrono::system_clock::to_time_t(end_0_for);
      if (timestamp){
        std::cout << "elapsed_0_for time: " << elapsed_seconds_0_for.count() << "s\n";      
      }
    } //end for
    auto end_0_while = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds_0_while = end_0_while-start_0_while;
    std::time_t end_time_0_while = std::chrono::system_clock::to_time_t(end_0_while);
    if (timestamp){
      std::cout << "elapsed_0_while time: " << elapsed_seconds_0_while.count() << "s\n";      
    }
  } //end while


  auto start_w = std::chrono::system_clock::now();

  // count_out_writer << count_out_str;

  auto end_w = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds_w = end_w-start_w;
  std::time_t end_time_w = std::chrono::system_clock::to_time_t(end_w);
  std::cout << "elapsed_w time: " << elapsed_seconds_w.count() << "s\n";      



  auto start_w2 = std::chrono::system_clock::now();

  count_out_writer.close();

  auto end_w2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds_w2 = end_w2-start_w2;
  std::time_t end_time_w2 = std::chrono::system_clock::to_time_t(end_w2);
  std::cout << "elapsed_w2 time: " << elapsed_seconds_w2.count() << "s\n";      




  auto end_main = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds_main = end_main-start_main;
  std::time_t end_time_main = std::chrono::system_clock::to_time_t(end_main);
  std::cout << "elapsed_main time: " << elapsed_seconds_main.count() << "s\n";      

  return 0;
}
