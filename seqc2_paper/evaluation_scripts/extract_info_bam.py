import multiprocessing
import argparse
import os
import shutil
import traceback
import logging

import pybedtools

import sys, os, re, pysam
import scipy.stats as stats
import numpy as np
import gzip

nan = float('nan')

def get_chromosomes_order(reference=None, bam=None):
    chroms_order = {}
    if reference:
        with pysam.FastaFile(reference) as fd:
            chroms_order = {chrom: chrom_index for chrom_index,
                            chrom in enumerate(fd.references)}

    if bam:
        with pysam.AlignmentFile(bam, "rb") as fd:
            chroms_order = {chrom: chrom_index for chrom_index,
                            chrom in enumerate(fd.references)}

    return chroms_order

## Dedup test for BAM file
def dedup_test(read_i, remove_dup_or_not=True):
    '''
    Return False (i.e., remove the read) if the read is a duplicate and if the user specify that duplicates should be removed.
    Else return True (i.e, keep the read)
    '''
    if read_i.is_duplicate and remove_dup_or_not:
        return False
    else:
        return True

### PYSAM ###
def position_of_aligned_read(read_i, target_position):
    '''
    Return the base call of the target position, or if it's a start of insertion/deletion.
    This target position follows pysam convension, i.e., 0-based.
    In VCF files, deletions/insertions occur AFTER the position.
    Return (Code, seq_i, base_at_target, indel_length, nearest insertion/deletion)
    The first number in result is a code:
    1) Match to reference, which is either a reference read or a SNV/SNP
    2) Deletion after the target position
    3) Insertion after the target position
    0) The target position does not match to reference, and may be discarded for "reference/alternate" read count purposes, but can be kept for "inconsistent read" metrics.
    '''


    for i, align_i in enumerate(read_i.get_aligned_pairs()):

        # If find a match:
        if align_i[1] == target_position:
            seq_i = align_i[0]
            break


    # If the target position is aligned:
    try:
        if seq_i is not None:
            base_at_target = read_i.seq[seq_i]

            # Whether if it's a Deletion/Insertion depends on what happens after this position:
            # If the match (i.e., i, seq_i) is the final alignment, then you cannot know if it's an indel
            # if "i" is NOT the final alignment:
            if i != len(read_i.get_aligned_pairs()) - 1:

                indel_length = 0
                # If the next alignment is the next sequenced base, then the target is either a reference read of a SNP/SNV:
                if read_i.get_aligned_pairs()[i+1][0] == seq_i+1 and read_i.get_aligned_pairs()[i+1][1] == target_position + 1:

                    code = 1 # Reference read for mismatch

                # If the next reference position has no read position to it, it is DELETED in this read:
                elif read_i.get_aligned_pairs()[i+1][0] == None and read_i.get_aligned_pairs()[i+1][1] == target_position + 1:

                    code = 2 # Deletion

                    for align_j in read_i.get_aligned_pairs()[ i+1:: ]:
                        if align_j[0] == None:
                            indel_length -= 1
                        else:
                            break

                # Opposite of deletion, if the read position cannot be aligned to the reference, it can be an INSERTION.
                # Insertions sometimes show up wit soft-clipping at the end, if the inserted sequence is "too long" to align on a single read. In this case, the inserted length derived here is but a lower limit of the real inserted length.
                elif read_i.get_aligned_pairs()[i+1][0] == seq_i+1 and read_i.get_aligned_pairs()[i+1][1] == None:

                    code = 3 # Insertion or soft-clipping

                    for align_j in read_i.get_aligned_pairs()[ i+1:: ]:
                        if align_j[1] == None:
                            indel_length += 1
                        else:
                            break

            # If "i" is the final alignment, cannt exam for indel:
            else:
                code = 1           # Assuming no indel
                indel_length = nan # Would be zero if certain no indel, but uncertain here

        # If the target position is deleted from the sequencing read (i.e., the deletion in this read occurs before the target position):
        else:
            code = 0
            base_at_target, indel_length = None, None
            

        # See if there is insertion/deletion within 5 bp of "i":
        if isinstance(indel_length, int):
            left_side_start = seq_i
            right_side_start = seq_i + abs(indel_length) + 1
            switch = 1
            for j in (3,2,1):
                for indel_seeker_i in left_side_start, right_side_start:

                    switch = switch * -1
                    displacement = j * switch
                    seq_j = indel_seeker_i + displacement

                    if 0 <= seq_j < len(read_i.get_aligned_pairs()):

                        # If the reference position has no base aligned to it, it's a deletion.
                        # On the other hand, if the base has no reference base aligned to it, it's an insertion.
                        if read_i.get_aligned_pairs()[ seq_j ][1] == None or read_i.get_aligned_pairs()[ seq_j ][0] == None:
                            break

        return code, seq_i, base_at_target, indel_length

    # The target position does not exist in the read
    except UnboundLocalError:
        return None, None, None, None


def from_bam(bam, my_coordinate, ref_base, first_alt, min_mq=1, min_bq=10):

    '''
    bam is the opened file handle of bam file
    my_coordiate is a list or tuple of 0-based (contig, position)
    '''
    
    indel_length = len(first_alt) - len(ref_base)
    reads = bam.fetch( my_coordinate[0], my_coordinate[1]-1, my_coordinate[1] )
    
    ref_for = ref_rev = alt_for = alt_rev = dp = 0
    noise_read_count = poor_read_count  = 0
        
    qname_collector = {}
    
    for read_i in reads:
        if not read_i.is_unmapped and dedup_test(read_i):
            
            dp += 1
            
            code_i, ith_base, base_call_i, indel_length_i = position_of_aligned_read(read_i, my_coordinate[1]-1 )
            
            if read_i.mapping_quality < min_mq and np.mean(read_i.query_qualities) < min_bq:
                poor_read_count += 1

            # Reference calls:
            if code_i == 1 and base_call_i == ref_base[0]:

                try:
                    qname_collector[read_i.qname].append(0)
                except KeyError:
                    qname_collector[read_i.qname] = [0]
            
                                                
                # Orientation
                if (not read_i.is_reverse) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                    ref_for += 1
                elif    read_i.is_reverse  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                    ref_rev += 1
                
            # Alternate calls:
            # SNV, or Deletion, or Insertion where I do not check for matching indel length
            elif (indel_length == 0 and code_i == 1 and base_call_i == first_alt) or \
                 (indel_length < 0  and code_i == 2 and indel_length == indel_length_i) or \
                 (indel_length > 0  and code_i == 3):

                try:
                    qname_collector[read_i.qname].append(1)
                except KeyError:
                    qname_collector[read_i.qname] = [1]

                # Orientation
                if (not read_i.is_reverse) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                    alt_for += 1
                elif    read_i.is_reverse  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                    alt_rev += 1
                            
            # Inconsistent read or 2nd alternate calls:
            else:
                
                try:
                    qname_collector[read_i.qname].append(2)
                except KeyError:
                    qname_collector[read_i.qname] = [2]
                
                noise_read_count += 1
    ro=ref_for+ref_rev
    ao=alt_for+alt_rev
    af=np.round(ao/float(dp+0.00001),4)

    return dp,ro,ao,af

def extract_split_info((bam,reference,tbc,min_mq,min_bq)):
    bam_ = pysam.AlignmentFile(bam, reference_filename=reference)
    my_info={}
    for chrom,pos,ref,alt,id_ in tbc:
        info=from_bam(bam_,[chrom,int(pos)],ref,alt,min_mq,min_bq)
        my_info[id_] = ":".join(map(str,info))
    return my_info


def extract_af(variants, output, bam, reference, min_mq, min_bq, num_threads):
    chroms_order=get_chromosomes_order(reference)
    chroms=pysam.FastaFile(reference).references
    info={}
    if os.path.exists(output):
        with open(output) as a_e:
            for line in a_e:
                if line[0]=="#":
                    continue
                chrom,pos,_,ref,alt,_,_,_,_,info_=line.strip().split()
                id_="-".join([str(chroms_order[chrom]),pos,ref,alt])
                info[id_]=info_

    to_be_checked=[]
    all_var_ids=[]
    if variants.endswith(".gz"):
        with gzip.open(variants,"rt") as a_e:
            for line in a_e:
                if line[0]=="#":
                    continue
                chrom,pos,_,ref,alt=line.strip().split()[0:5]
                id_="-".join([str(chroms_order[chrom]),pos,ref,alt])
                all_var_ids.append(id_)
                if id_ not in info:
                    to_be_checked.append([chrom,pos,ref,alt,id_])
    else:
        with open(variants) as a_e:
            for line in a_e:
                if line[0]=="#":
                    continue
                chrom,pos,_,ref,alt=line.strip().split()[0:5]
                id_="-".join([str(chroms_order[chrom]),pos,ref,alt])
                all_var_ids.append(id_)
                if id_ not in info:
                    to_be_checked.append([chrom,pos,ref,alt,id_])

    print "TB Checkecd:", len(to_be_checked)
    if to_be_checked:
        split=len(to_be_checked)/num_threads
        map_args=[]
        for i in range(num_threads):
            e=min((i+1)*split,len(to_be_checked))
            if i==num_threads-1:
                e=len(to_be_checked)
            map_args.append([bam, reference,to_be_checked[i*split:e],min_mq,min_bq])

        pool = multiprocessing.Pool(num_threads)
        try:
            extracted_info = pool.map_async(
                extract_split_info, map_args).get()
            pool.close()
        except Exception as inst:
            print inst
            pool.close()
            traceback.print_exc()
            raise Exception

        for i in extracted_info:
            info.update(i)

    print "Checked:", len(info)
    ids=sorted(all_var_ids,key=lambda x: map(int,x.split("-")[0:2]))
    with open(output,"w") as o_f:
        o_f.write("##fileformat=VCFv4.2\n")
        o_f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for id_ in ids:
            chrom_,pos,ref,alt=id_.split("-")
            chrom=chroms[int(chrom_)]
            dp,ro,ao,af=info[id_].split(":")
            o_f.write("\t".join([chrom,pos,".",ref,alt,af,ao,dp,"DP:RO:AO:AF",info[id_]])+"\n")




if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='extract bam info')
    parser.add_argument('--variants', type=str, help='variants vcf (or vcf.gz) file',
                        required=True)
    parser.add_argument('--bam', type=str,
                        help='bam file', required=True)
    parser.add_argument('--reference', type=str,
                        help='reference file', required=True)
    parser.add_argument('--min_mq', type=int,
                        help='minimum mapping quality', default=1)
    parser.add_argument('--min_bq', type=float,
                        help='minimum base quality', default=5)
    parser.add_argument('--num_threads', type=int,
                        help='number of threads', default=1)
    parser.add_argument('--output', type=str,
                        help='output vcf', required=True)
    args = parser.parse_args()
    print args

    try:
        extract_af(args.variants, args.output, args.bam, args.reference, args.min_mq, args.min_bq,args.num_threads)
    except Exception as e:
        print traceback.format_exc()
        print "Aborting!"
        print "extract_bam.py failure on arguments: {}".format(args)
        raise e
