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

from extract_postprocess_targets import extract_postprocess_targets
from merge_post_vcfs import merge_post_vcfs
from resolve_variants import resolve_variants
from utils import concatenate_files
from long_read_indelrealign import long_read_indelrealign
from resolve_scores import resolve_scores


FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logFormatter = logging.Formatter(FORMAT)
logger = logging.getLogger(__name__)
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)
logging.getLogger().setLevel(logging.INFO)


def postprocess(work, reference, pred_vcf_file, output_vcf, candidates_vcf, tumor_bam, min_len,
                postprocess_max_dist, long_read,
                lr_pad, lr_chunck_size, lr_chunck_scale,
                lr_snp_min_af, lr_ins_min_af, lr_del_min_af, lr_match_score, lr_mismatch_penalty,
                lr_gap_open_penalty, lr_gap_ext_penalty,
                pass_threshold, lowqual_threshold,
                msa_binary, samtools_binary, num_threads):

    logger.info("-----------------------------------------------------------")
    logger.info("Postprocessing")
    logger.info("-----------------------------------------------------------")

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
                               lr_gap_ext_penalty, msa_binary, samtools_binary)
        resolve_scores(tumor_bam, ra_resolved_vcf, target_vcf, resolved_vcf)

    all_no_resolve = concatenate_files(
        [no_resolve, ensembled_preds], os.path.join(work, "no_resolve.vcf"))

    logger.info("Merge vcfs")
    merge_post_vcfs(reference, resolved_vcf,
                    all_no_resolve, target_vcf, output_vcf,
                    pass_threshold, lowqual_threshold)

    logger.info("Output NeuSomatic prediction at {}".format(output_vcf))
    logger.info("Done.")
    return output_vcf

if __name__ == '__main__':
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
    parser.add_argument('--samtools_binary', type=str,
                        help='samtools binary', default="samtools")
    parser.add_argument('--num_threads', type=int,
                        help='number of threads', default=1)
    parser.add_argument('--work', type=str,
                        help='work directory', required=True)
    args = parser.parse_args()
    logger.info(args)

    try:
        output_vcf = postprocess(args.work, args.reference, args.pred_vcf, args.output_vcf,
                                 args.candidates_vcf, args.tumor_bam, args.min_len,
                                 args.postprocess_max_dist, args.long_read,
                                 args.lr_pad, args.lr_chunck_size, args.lr_chunck_scale,
                                 args.lr_snp_min_af, args.lr_ins_min_af, args.lr_del_min_af,
                                 args.lr_match_score, args.lr_mismatch_penalty,
                                 args.lr_gap_open_penalty,
                                 args.lr_gap_ext_penalty, args.pass_threshold, args.lowqual_threshold,
                                 args.msa_binary, args.samtools_binary, args.num_threads)

    except:
        traceback.print_exc()
        logger.error("Aborting!")
        raise Exception("postprocess.py failure on arguments: {}".format(args))
