#-------------------------------------------------------------------------
# postprocess.py
# A wrapper that
# 1- Extract variants that need postprocessing (call to 'extract_postprocess_targets.py')
# 2- Resolve exact variants for target variants (call to
#    'resolve_variants.py' or 'long_read_indelrealign.py'/'resolve_scores.py')
# 3- Merge resolved variants and other predicted variants and
#    output the final NeuSomatic .vcf (call to 'merge_post_vcfs.py')
#-------------------------------------------------------------------------

import argparse
import os
import traceback
import logging
import shutil

import pybedtools
import pysam
import numpy as np

from extract_postprocess_targets import extract_postprocess_targets
from merge_post_vcfs import merge_post_vcfs
from resolve_variants import resolve_variants
from utils import concatenate_files, get_chromosomes_order
from long_read_indelrealign import long_read_indelrealign
from resolve_scores import resolve_scores
from _version import __version__


def add_vcf_info(work, reference, merged_vcf, candidates_vcf, ensemble_tsv,
                 output_vcf, pass_threshold, lowqual_threshold):
    merged_vcf = pybedtools.BedTool(merged_vcf)
    candidates_vcf = pybedtools.BedTool(candidates_vcf)
    ensemble_candids_vcf = []
    if ensemble_tsv:
        ensemble_candids_vcf = os.path.join(work, "ensemble_candids.vcf")
        with open(ensemble_tsv) as e_f:
            with open(ensemble_candids_vcf, "w") as c_f:
                c_f.write("##fileformat=VCFv4.2\n")
                c_f.write(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
                for line in e_f:
                    if "T_REF_FOR" in line:
                        header = line.strip().split()
                        chrom_id = header.index("CHROM")
                        pos_id = header.index("POS")
                        ref_id = header.index("REF")
                        alt_id = header.index("ALT")
                        dp_id = header.index("T_DP")
                        ref_fw_id = header.index("T_REF_FOR")
                        ref_rv_id = header.index("T_REF_REV")
                        alt_fw_id = header.index("T_ALT_FOR")
                        alt_rv_id = header.index("T_ALT_REV")
                        continue
                    fields = line.strip().split()
                    chrom = fields[chrom_id]
                    pos = fields[pos_id]
                    ref = fields[ref_id]
                    alt = fields[alt_id]
                    dp = int(fields[dp_id])
                    ro_fw = int(fields[ref_fw_id])
                    ro_rv = int(fields[ref_rv_id])
                    ao_fw = int(fields[alt_fw_id])
                    ao_rv = int(fields[alt_rv_id])
                    ro = ro_fw + ro_rv
                    ao = ao_fw + ao_rv
                    af = np.round(ao / float(ao + ro + 0.0001), 4)
                    c_f.write(
                        "\t".join(map(str, [chrom, pos, ".", ref, alt, ".", ".", ".", "GT:DP:RO:AO:AF", ":".join(map(str, ["0/1", dp, ro, ao, af]))])) + "\n")

    ensemble_candids_vcf = pybedtools.BedTool(ensemble_candids_vcf)
    in_candidates = merged_vcf.window(candidates_vcf, w=5)
    notin_candidates = merged_vcf.window(candidates_vcf, w=5, v=True)
    in_ensemble = merged_vcf.window(ensemble_candids_vcf, w=5)
    notin_any = notin_candidates.window(ensemble_candids_vcf, w=5, v=True)
    chroms_order = get_chromosomes_order(reference=reference)
    with pysam.FastaFile(reference) as rf:
        chroms = rf.references

    scores = {}
    tags_info = {}
    for s_e, dd in [0, in_candidates], [1, in_ensemble]:
        for x in dd:
            tag = "-".join([str(chroms_order[x[0]]), x[1], x[3], x[4]])
            scores[tag] = [x[5], x[6], x[7], x[9]]
            if tag not in tags_info:
                tags_info[tag] = []
            info = x[19].split(":")
            dp, ro, ao = map(int, info[1:4])
            af = float(info[4])
            is_same = x[1] == x[11] and x[3] == x[13] and x[4] == x[14]
            is_same_type = np.sign(
                len(x[3]) - len(x[13])) == np.sign(len(x[4]) - len(x[14]))
            dist = abs(int(x[1]) - int(x[11]))
            len_diff = abs((len(x[3]) - len(x[13])) - (len(x[4]) - len(x[14])))
            tags_info[tag].append(
                [~is_same, ~is_same_type, dist, len_diff, s_e, dp, ro, ao, af])
    fina_info_tag = {}
    for tag, hits in tags_info.iteritems():
        hits = sorted(hits, key=lambda x: x[0:5])
        fina_info_tag[tag] = hits[0][5:]

    for x in notin_any:
        tag = "-".join([str(chroms_order[x[0]]), x[1], x[3], x[4]])
        fina_info_tag[tag] = [0, 0, 0, 0]
        scores[tag] = [x[5], x[6], x[7], x[9]]

    tags = sorted(fina_info_tag.keys(), key=lambda x: map(int, x.split("-")[0:2]
                                                          ))
    with open(output_vcf, "w") as o_f:
        o_f.write("##fileformat=VCFv4.2\n")
        o_f.write("##NeuSomatic Version={}\n".format(__version__))
        o_f.write(
            "##FORMAT=<ID=SCORE,Number=1,Type=Float,Description=\"Prediction probability score\">\n")
        o_f.write("##FILTER=<ID=PASS,Description=\"Accept as a higher confidence somatic mutation calls with probability score value at least {}\">\n".format(
            pass_threshold))
        o_f.write("##FILTER=<ID=LowQual,Description=\"Less confident somatic mutation calls with probability score value at least {}\">\n".format(
            lowqual_threshold))
        o_f.write("##FILTER=<ID=REJECT,Description=\"Rejected as a confident somatic mutation with probability score value below {}\">\n".format(
            lowqual_threshold))
        o_f.write(
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        o_f.write(
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth in the tumor\">\n")
        o_f.write(
            "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count in the tumor\">\n")
        o_f.write(
            "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count in the tumor\">\n")
        o_f.write(
            "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele fractions of alternate alleles in the tumor\">\n")
        o_f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for tag in tags:
            chrom_id, pos, ref, alt = tag.split("-")
            qual, filter_, score, gt = scores[tag]
            dp, ro, ao, af = fina_info_tag[tag]
            info_field = "{};DP={};RO={};AO={};AF={};".format(
                score, dp, ro, ao, af)
            gt_field = "{}:{}:{}:{}:{}".format(gt, dp, ro, ao, af)
            o_f.write("\t".join(map(str, [chroms[int(chrom_id)], str(
                pos), ".", ref, alt, qual, filter_, info_field, "GT:DP:RO:AO:AF", gt_field])) + "\n")


def postprocess(work, reference, pred_vcf_file, output_vcf, candidates_vcf, ensemble_tsv,
                tumor_bam, min_len,
                postprocess_max_dist, long_read,
                lr_pad, lr_chunck_size, lr_chunck_scale,
                lr_snp_min_af, lr_ins_min_af, lr_del_min_af, lr_match_score, lr_mismatch_penalty,
                lr_gap_open_penalty, lr_gap_ext_penalty,
                pass_threshold, lowqual_threshold,
                msa_binary, num_threads):
    logger = logging.getLogger(postprocess.__name__)

    logger.info("----------------------Postprocessing-----------------------")
    if not os.path.exists(work):
        os.mkdir(work)

    candidates_preds = os.path.join(work, "candidates_preds.vcf")
    ensembled_preds = os.path.join(work, "ensembled_preds.vcf")
    pred_vcf = pybedtools.BedTool(pred_vcf_file)
    pred_vcf.window(candidates_vcf, w=5, v=True).saveas(ensembled_preds)
    pred_vcf.window(candidates_vcf, w=5, u=True).saveas(candidates_preds)

    logger.info("Extract targets")
    postprocess_pad = 1 if not long_read else 10
    extract_postprocess_targets(
        candidates_preds, min_len, postprocess_max_dist, postprocess_pad)

    no_resolve = os.path.join(work, "candidates_preds.no_resolve.vcf")
    target_vcf = os.path.join(work, "candidates_preds.resolve_target.vcf")
    target_bed = os.path.join(work, "candidates_preds.resolve_target.bed")
    resolved_vcf = os.path.join(work, "candidates_preds.resolved.vcf")

    logger.info("Resolve targets")
    if not long_read:
        resolve_variants(tumor_bam, resolved_vcf,
                         reference, target_vcf, target_bed, num_threads)
    else:
        work_lr_indel_realign = os.path.join(work, "work_lr_indel_realign")
        if os.path.exists(work_lr_indel_realign):
            shutil.rmtree(work_lr_indel_realign)
        os.mkdir(work_lr_indel_realign)
        ra_resolved_vcf = os.path.join(
            work, "candidates_preds.ra_resolved.vcf")
        long_read_indelrealign(work_lr_indel_realign, tumor_bam, None, ra_resolved_vcf, target_bed,
                               reference, num_threads, lr_pad,
                               lr_chunck_size, lr_chunck_scale, lr_snp_min_af,
                               lr_del_min_af, lr_ins_min_af,
                               lr_match_score, lr_mismatch_penalty, lr_gap_open_penalty,
                               lr_gap_ext_penalty, msa_binary)
        resolve_scores(tumor_bam, ra_resolved_vcf, target_vcf, resolved_vcf)

    all_no_resolve = concatenate_files(
        [no_resolve, ensembled_preds], os.path.join(work, "no_resolve.vcf"))

    logger.info("Merge vcfs")
    merged_vcf = os.path.join(work, "merged_preds.vcf")
    merge_post_vcfs(reference, resolved_vcf,
                    all_no_resolve, merged_vcf,
                    pass_threshold, lowqual_threshold)
    add_vcf_info(work, reference, merged_vcf,
                 candidates_vcf, ensemble_tsv, output_vcf,
                 pass_threshold, lowqual_threshold)

    logger.info("Output NeuSomatic prediction at {}".format(output_vcf))
    logger.info("Postprocessing is Done.")
    return output_vcf

if __name__ == '__main__':

    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description='Preprocess predictions for call mode')
    parser.add_argument('--reference', type=str,
                        help='reference fasta filename', required=True)
    parser.add_argument('--tumor_bam', type=str,
                        help='tumor bam', required=True)
    parser.add_argument('--pred_vcf', type=str,
                        help='predicted vcf', required=True)
    parser.add_argument('--output_vcf', type=str,
                        help='output final vcf', required=True)
    parser.add_argument('--candidates_vcf', type=str,
                        help='filtered candidate vcf', required=True)
    parser.add_argument('--ensemble_tsv', type=str,
                        help='Ensemble annotation tsv file (only for short read)', default=None)
    parser.add_argument('--min_len', type=int,
                        help='minimum INDEL len to resolve', default=4)
    parser.add_argument('--postprocess_max_dist', type=int,
                        help='max distance to neighboring variant', default=5)
    parser.add_argument(
        '--long_read', help='Enable long_read (high error-rate sequence) indel realignment', action="store_true")
    parser.add_argument(
        '--lr_pad', type=int, help='long_read indel realign: #base padding to the regions', default=1)
    parser.add_argument('--lr_chunck_size', type=int,
                        help='long_read indel realign: chuck split size for high depth', default=600)
    parser.add_argument('--lr_chunck_scale', type=float,
                        help='long_read indel realign: chuck scale size for high depth', default=1.5)
    parser.add_argument('--lr_snp_min_af', type=float,
                        help='long_read indel realign: SNP min allele freq', default=0.05)
    parser.add_argument('--lr_ins_min_af', type=float,
                        help='long_read indel realign: INS min allele freq', default=0.05)
    parser.add_argument('--lr_del_min_af', type=float,
                        help='long_read indel realign: DEL min allele freq', default=0.05)
    parser.add_argument('--lr_match_score', type=int,
                        help='long_read indel realign: match score', default=10)
    parser.add_argument('--lr_mismatch_penalty', type=int,
                        help='long_read indel realign: penalty for having a mismatch', default=8)
    parser.add_argument('--lr_gap_open_penalty', type=int,
                        help='long_read indel realign: penalty for opening a gap', default=8)
    parser.add_argument('--lr_gap_ext_penalty', type=int,
                        help='long_read indel realign: penalty for extending a gap', default=6)
    parser.add_argument('--pass_threshold', type=float,
                        help='SCORE for PASS (PASS for score => pass_threshold)', default=0.7)
    parser.add_argument('--lowqual_threshold', type=float,
                        help='SCORE for LowQual (PASS for lowqual_threshold <= score < pass_threshold)',
                        default=0.4)
    parser.add_argument('--msa_binary', type=str,
                        help='MSA binary', default="../bin/msa")
    parser.add_argument('--num_threads', type=int,
                        help='number of threads', default=1)
    parser.add_argument('--work', type=str,
                        help='work directory', required=True)
    args = parser.parse_args()
    logger.info(args)

    try:
        output_vcf = postprocess(args.work, args.reference, args.pred_vcf, args.output_vcf,
                                 args.candidates_vcf, args.ensemble_tsv,
                                 args.tumor_bam, args.min_len,
                                 args.postprocess_max_dist, args.long_read,
                                 args.lr_pad, args.lr_chunck_size, args.lr_chunck_scale,
                                 args.lr_snp_min_af, args.lr_ins_min_af, args.lr_del_min_af,
                                 args.lr_match_score, args.lr_mismatch_penalty,
                                 args.lr_gap_open_penalty,
                                 args.lr_gap_ext_penalty, args.pass_threshold, args.lowqual_threshold,
                                 args.msa_binary, args.num_threads)

    except Exception as e:
        logger.error(traceback.format_exc())
        logger.error("Aborting!")
        logger.error(
            "postprocess.py failure on arguments: {}".format(args))
        raise e
