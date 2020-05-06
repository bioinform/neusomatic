#!/usr/bin/env python3

import sys
import os
import re
import pysam
import scipy.stats as stats
import genomic_file_handlers as genome
from read_info_extractor import *
from collections import defaultdict

nan = float('nan')


class AlignmentFeatures:
    def __init__(self, bam, my_coordinate, ref_base, first_alt, min_mq=1, min_bq=10):
        '''
        bam is the opened file handle of bam file
        my_coordiate is a list or tuple of 0-based (contig, position)
        '''

        indel_length = len(first_alt) - len(ref_base)
        reads = bam.fetch(my_coordinate[0], my_coordinate[1] - 1, my_coordinate[1])

        # index 0 for ref, 1 for alt
        read_mq = [[], []]
        read_bq = [[], []]
        edit_distance = [[], []]
        flanking_indel = [[], []]
        pos_from_end = [[], []]
        concordance_counts = [[0, 0], [0, 0]]
        orientation_counts = [[0, 0], [0, 0]]
        soft_clip_counts = [[0, 0], [0, 0]]
        dp = 0
        MQ0 = 0

        noise_read_count = poor_read_count = 0

        qname_collector = defaultdict(list)

        for read_i in reads:
            if read_i.is_unmapped or not dedup_test(read_i):
                continue

            dp += 1

            code_i, ith_base, base_call_i, indel_length_i, flanking_indel_i = position_of_aligned_read(
                read_i, my_coordinate[1] - 1)

            if read_i.mapping_quality < min_mq and mean(read_i.query_qualities) < min_bq:
                poor_read_count += 1

            if read_i.mapping_quality == 0:
                MQ0 += 1

            is_ref_call = code_i == 1 and base_call_i == ref_base[0]
            is_alt_call = (indel_length == 0 and code_i == 1 and base_call_i == first_alt) or (
                        indel_length < 0 and code_i == 2 and indel_length == indel_length_i) or (
                                      indel_length > 0 and code_i == 3)

            # inconsistent read or second alternate calls
            if not (is_ref_call or is_alt_call):
                qname_collector[read_i.qname].append(2)
                noise_read_count += 1
                continue

            index = 1 if is_alt_call else 0

            qname_collector[read_i.qname].append(index)

            read_mq[index].append(read_i.mapping_quality)
            read_bq[index].append(read_i.query_qualities[ith_base])

            try:
                edit_distance[index].append(read_i.get_tag('NM'))
            except KeyError:
                pass

            if read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                concordance_counts[0 if read_i.is_proper_pair else 1][index] += 1
                orientation_counts[1 if read_i.is_reverse else 0][index] += 1

            is_soft_clipped = read_i.cigar[0][0] == cigar_soft_clip or read_i.cigar[-1][0] == cigar_soft_clip
            soft_clip_counts[1 if is_soft_clipped else 0][index] += 1

            # Distance from the end of the read:
            if ith_base is not None:
                pos_from_end[index].append(min(ith_base, read_i.query_length - ith_base))

            flanking_indel[index].append(flanking_indel_i)

        # unpack to get the ref and alt values
        ref_pos_from_end, alt_pos_from_end = pos_from_end
        self.ref_concordant_reads, self.alt_concordant_reads = concordance_counts[0]
        self.ref_discordant_reads, self.alt_discordant_reads = concordance_counts[1]
        self.ref_for, self.alt_for = orientation_counts[0]
        self.ref_rev, self.alt_rev = orientation_counts[1]
        self.ref_notSC_reads, self.alt_notSC_reads = soft_clip_counts[0]
        self.ref_SC_reads, self.alt_SC_reads = soft_clip_counts[1]

        # Done extracting info from BAM. Now tally them:
        ref_read_mq, alt_read_mq = read_mq
        self.ref_mq = mean(ref_read_mq)
        self.alt_mq = mean(alt_read_mq)
        self.z_ranksums_mq = stats.ranksums(alt_read_mq, ref_read_mq)[0]

        ref_read_bq, alt_read_bq = read_bq
        self.ref_bq = mean(ref_read_bq)
        self.alt_bq = mean(alt_read_bq)
        self.z_ranksums_bq = stats.ranksums(alt_read_bq, ref_read_bq)[0]

        ref_edit_distance, alt_edit_distance = edit_distance
        self.ref_NM = mean(ref_edit_distance)
        self.alt_NM = mean(alt_edit_distance)
        self.z_ranksums_NM = stats.ranksums(alt_edit_distance, ref_edit_distance)[0]
        self.NM_Diff = self.alt_NM - self.ref_NM - abs(indel_length)

        self.concordance_fet = stats.fisher_exact(concordance_counts)[1]
        self.strandbias_fet = stats.fisher_exact(orientation_counts)[1]
        self.clipping_fet = stats.fisher_exact(soft_clip_counts)[1]

        self.z_ranksums_endpos = stats.ranksums(alt_pos_from_end, ref_pos_from_end)[0]

        ref_flanking_indel, alt_flanking_indel = flanking_indel
        self.ref_indel_1bp = ref_flanking_indel.count(1)
        self.ref_indel_2bp = ref_flanking_indel.count(2) + self.ref_indel_1bp
        self.ref_indel_3bp = ref_flanking_indel.count(3) + self.ref_indel_2bp
        self.alt_indel_1bp = alt_flanking_indel.count(1)
        self.alt_indel_2bp = alt_flanking_indel.count(2) + self.alt_indel_1bp
        self.alt_indel_3bp = alt_flanking_indel.count(3) + self.alt_indel_2bp

        self.consistent_mates = self.inconsistent_mates = 0
        for one_count in map(lambda x: x.count(1), filter(lambda y: len(y) == 2, qname_collector.values())):
            # Both are alternative calls:
            if one_count == 2:
                self.consistent_mates += 1

            # One is alternate call but the other one is not:
            elif one_count == 1:
                self.inconsistent_mates += 1

        self.nref = self.ref_for + self.ref_rev
        self.nalt = self.alt_for + self.alt_rev
        self.dp = dp
        self.MQ0 = MQ0
        self.noise_read_count = noise_read_count
        self.poor_read_count = poor_read_count


def from_genome_reference(ref_fa, my_coordinate, ref_base, first_alt):
    '''
    ref_fa is the opened reference fasta file handle
    my_coordiate is a list or tuple of 0-based (contig, position)
    '''

    # Homopolymer eval (Make sure to modify for INDEL):
    # The min and max is to prevent the +/- 20 bases from exceeding the ends
    # of the reference sequence
    lseq = ref_fa.fetch(my_coordinate[0], max(
        0, my_coordinate[1] - 20), my_coordinate[1])
    rseq = ref_fa.fetch(my_coordinate[0], my_coordinate[
                        1] + 1, min(ref_fa.get_reference_length(my_coordinate[0]) + 1, my_coordinate[1] + 21))

    # This is to get around buy in old version of pysam that reads the
    # reference sequence in bytes instead of strings
    lseq = lseq.decode() if isinstance(lseq, bytes) else lseq
    rseq = rseq.decode() if isinstance(rseq, bytes) else rseq

    seq41_ref = lseq + ref_base + rseq
    seq41_alt = lseq + first_alt + rseq

    ref_counts = genome.count_repeating_bases(seq41_ref)
    alt_counts = genome.count_repeating_bases(seq41_alt)

    homopolymer_length = max(max(ref_counts), max(alt_counts))

    # Homopolymer spanning the variant site:
    ref_c = 0
    alt_c = 0
    for i in rseq:
        if i == ref_base:
            ref_c += 1
        else:
            break

    for i in lseq[::-1]:
        if i == ref_base:
            ref_c += 1
        else:
            break

    for i in rseq:
        if i == first_alt:
            alt_c += 1
        else:
            break

    for i in lseq[::-1]:
        if i == first_alt:
            alt_c += 1
        else:
            break

    site_homopolymer_length = max(alt_c + 1, ref_c + 1)

    return homopolymer_length, site_homopolymer_length


def somaticOddRatio(n_ref, n_alt, t_ref, t_alt, max_value=100):

    # Odds Ratio just like VarDict's output
    sor_numerator = n_alt * t_ref
    sor_denominator = n_ref * t_alt
    if sor_numerator == 0 and sor_denominator == 0:
        sor = nan
    elif sor_denominator == 0:
        sor = max_value
    else:
        sor = min(sor_numerator / sor_denominator, max_value)

    return sor

def max_sub_vocabularies(seq_length, max_subseq_length):
    # According to:
    # https://doi.org/10.1093/bioinformatics/18.5.679
    # capping the length of sub_string as an input parameter
    assert max_subseq_length <= seq_length
    
    counts = 0
    k = 1
    while k <= max_subseq_length:
        
        if 4**k < (seq_length - k + 1):
            counts = counts + 4**k
        else:
            counts = counts + (2*seq_length - k - max_subseq_length + 2) * (max_subseq_length - k + 1)/2
            break
        
        k += 1
                
    return counts


def subLC(sequence, max_substring_length=20):
    # Calculate linguistic sequence complexity according to
    # https://doi.org/10.1093/bioinformatics/18.5.679
    # Cut off substring at a fixed length
    sequence = sequence.upper()
    
    if not 'N' in sequence:
        
        number_of_subseqs     = 0
        seq_length            = len(sequence)
        max_number_of_subseqs = max_sub_vocabularies(seq_length, max_substring_length)
        
        set_of_seq_n = set()
        for i in range(1, min(max_substring_length+1, seq_length+1) ):
            set_of_seq_n.update((sequence[n: n+i] for n in range(len(sequence) - i + 1)))
        
        number_of_subseqs  = len(set_of_seq_n)
        lc = number_of_subseqs/max_number_of_subseqs
    
    else:
        lc = float('nan')

    return lc
