#!/usr/bin/env python
#-------------------------------------------------------------------------
# read_callers_vcf.py
# read callers vcf files and generate ensemble tsv
#-------------------------------------------------------------------------
import argparse
import traceback
import logging
import re
import gzip

import genomic_file_handlers as genome
from read_info_extractor import rescale
from utils import skip_empty, get_chromosomes_order

import numpy as np

# Normal/Tumor index in the Merged VCF file, or any other VCF file that
# puts NORMAL first.
idxN, idxT = 0, 1
nan = float('nan')


def get_info_value(info_field, variable, ith_alt=None):
    logger = logging.getLogger(get_info_value.__name__)
    key_item = re.search(
        r'\b{}=([^;\s]+)([;\W]|$)'.format(variable), info_field)

    # The key has a value attached to it, e.g., VAR=1,2,3
    if key_item:
        if ith_alt is None:
            return key_item.groups()[0]
        else:
            return key_item.groups()[0].split(",")[ith_alt]

    # Perhaps it's simply a flag without "="
    else:
        key_item = info_field.split(';')
        return True if variable in key_item else False


def get_sample_value(fields, samples, variable, idx=0):

    var2value = dict(zip(fields.split(':'), samples[idx].split(':')))
    try:
        return var2value[variable]
    except KeyError:
        return None


def get_mutect2_info(filters, info, ith_alt):

    mutect_classification = 1 if (get_info_value(info,
                                                 'SOMATIC') or 'PASS' in filters) else 0

    # MuTect2 has some useful information:
    nlod = get_info_value(info, 'NLOD', ith_alt)
    nlod = float(nlod) if nlod else nan

    tlod = get_info_value(info, 'TLOD', ith_alt)
    tlod = float(tlod) if tlod else nan

    tandem = 1 if get_info_value(info, 'STR') else 0

    ecnt = get_info_value(info, 'ECNT')
    if ecnt:
        try:
            ecnt = int(ecnt)
        except ValueError:
            ecnt = nan
    else:
        ecnt = nan
    return mutect_classification, nlod, tlod, tandem, ecnt


def get_varscan2_info(info):
    varscan_classification = 1 if get_info_value(info,
                                                 'SOMATIC') else 0
    return varscan_classification


def get_somaticsniper_info(fields, samples, idxT):
    somaticsniper_classification = 1 if get_sample_value(fields, samples,
                                                         'SS', idxT) == '2' else 0
    if somaticsniper_classification == 1:
        score_somaticsniper = get_sample_value(fields, samples,
                                               'SSC', idxT)
        score_somaticsniper = int(
            score_somaticsniper) if score_somaticsniper else nan
    else:
        score_somaticsniper = nan

    return somaticsniper_classification, score_somaticsniper


def get_vardict_info(filters, info, fields, samples):

    if (filters == 'PASS') and ('Somatic' in info):
        vardict_classification = 1
    elif 'Somatic' in info:
        vardict_filters = filters.split(';')

        disqualifying_filters = \
            ('d7'      in vardict_filters  or 'd5' in vardict_filters) or \
            ('DIFF0.2' in vardict_filters) or \
            ('LongAT'  in vardict_filters) or \
            ('MAF0.05' in vardict_filters) or \
            ('MSI6'    in vardict_filters) or \
            ('NM4'     in vardict_filters  or 'NM4.25' in vardict_filters) or \
            ('pSTD'    in vardict_filters) or \
            ('SN1.5'   in vardict_filters) or \
            ( 'P0.05'  in vardict_filters  and float(get_info_value(info, 'SSF') ) >= 0.15 ) or \
            (('v3' in vardict_filters or 'v4' in vardict_filters)
             and int(get_sample_value(fields, samples, 'VD', 1)) < 3)

        no_bad_filter = not disqualifying_filters
        filter_fail_times = len(vardict_filters)

        if no_bad_filter and filter_fail_times <= 2:
            vardict_classification = 0.5
        else:
            vardict_classification = 0

    else:
        vardict_classification = 0

    # Somatic Score:
    score_vardict = get_info_value(info, 'SSF')
    if score_vardict:
        score_vardict = float(score_vardict)
        score_vardict = genome.p2phred(score_vardict, max_phred=100)
        score_vardict = rescale(score_vardict,       'phred', None, 1001)
    else:
        score_vardict = nan

    # MSI, MSILEN, and SHIFT3:
    msi = get_info_value(info, 'MSI')
    if msi:
        msi = float(msi)
    else:
        msi = nan
    msilen = get_info_value(info, 'MSILEN')
    if msilen:
        msilen = float(msilen)
    else:
        msilen = nan
    shift3 = get_info_value(info, 'SHIFT3')
    if shift3:
        shift3 = float(shift3)
    else:
        shift3 = nan

    return vardict_classification, msi, msilen, shift3, score_vardict


def get_muse_info(filters):
    if filters == 'PASS':
        muse_classification = 1
    elif filters == 'Tier1':
        muse_classification = 0.9
    elif filters == 'Tier2':
        muse_classification = 0.8
    elif filters == 'Tier3':
        muse_classification = 0.7
    elif filters == 'Tier4':
        muse_classification = 0.6
    elif filters == 'Tier5':
        muse_classification = 0.5
    else:
        muse_classification = 0
    return muse_classification


def get_strelka2_info(filters, info):
    strelka_classification = 1 if 'PASS' in filters else 0
    somatic_evs = get_info_value(info, 'SomaticEVS')
    qss = get_info_value(info, 'QSS')
    tqss = get_info_value(info, 'TQSS')
    return strelka_classification, somatic_evs, qss, tqss


def open_textfile(file_name):

    # See if the input file is a .gz file:
    if file_name.lower().endswith('.gz'):
        return gzip.open(file_name, 'rt')

    else:
        return open(file_name)


def read_callers_vcf(reference,
                     output_tsv,
                     mutect2_vcfs,
                     strelka2_vcfs,
                     varscan2_vcfs,
                     muse_vcfs,
                     vardict_vcfs,
                     somaticsniper_vcfs,
                     min_caller):

    logger = logging.getLogger(read_callers_vcf.__name__)

    logger.info(
        "----------------------Read Callers VCF------------------------")

    mutect2_info = {}
    if mutect2_vcfs:
        for mutect2_vcf in mutect2_vcfs:
            i_f = open_textfile(mutect2_vcf)
            for line in skip_empty(i_f):
                x = line.strip().split()
                chrom, pos, _, ref, alts, _, filters, info = x[0:8]
                for ith_alt, alt in enumerate(alts.split(",")):
                    if ref != alt:
                        mutect_classification, nlod, tlod, tandem, ecnt = get_mutect2_info(
                            filters, info, ith_alt)
                        var_id = "-".join([chrom, pos, ref, alt])
                        mutect2_info[var_id] = [
                            mutect_classification, nlod, tlod, tandem, ecnt]
            i_f.close()
    strelka2_info = {}
    if strelka2_vcfs:
        for strelka2_vcf in strelka2_vcfs:
            i_f = open_textfile(strelka2_vcf)
            for line in skip_empty(i_f):
                x = line.strip().split()
                chrom, pos, _, ref, alts, _, filters, info = x[0:8]
                strelka_classification, somatic_evs, qss, tqss = get_strelka2_info(
                    filters, info)
                for alt in alts.split(","):
                    if ref != alt:
                        var_id = "-".join([chrom, pos, ref, alt])
                        strelka2_info[var_id] = [
                            strelka_classification, somatic_evs, qss, tqss]
            i_f.close()
    vardict_info = {}
    if vardict_vcfs:
        for vardict_vcf in vardict_vcfs:
            i_f = open_textfile(vardict_vcf)
            for line in skip_empty(i_f):
                x = line.strip().split()
                chrom, pos, _, ref, alts, _, filters, info, fields = x[0:9]
                samples = x[9:]

                # In the REF/ALT field, non-GCTA characters should be
                # changed to N to fit the VCF standard:
                ref = re.sub(r'[^GCTA]', 'N', ref, flags=re.I)
                alts = re.sub(r'[^GCTA]', 'N', alts, flags=re.I)

                vardict_classification, msi, msilen, shift3, score_vardict = get_vardict_info(
                    filters, info, fields, samples)
                for alt in alts.split(","):
                    if ref != alt:
                        if 'TYPE=SNV' in info or 'TYPE=Deletion' in info or 'TYPE=Insertion' in info:
                            var_id = "-".join([chrom, pos, ref, alt])
                            vardict_info[var_id] = [
                                vardict_classification, msi, msilen, shift3, score_vardict]
                        elif 'TYPE=Complex' in info and (len(ref) == len(alt)):
                            for i, (ref_i, alt_i) in enumerate(zip(ref, alt)):
                                if ref_i != alt_i:
                                    var_id = "-".join([chrom,
                                                       str(int(pos) + i), ref_i, alt_i])
                                    vardict_info[var_id] = [
                                        vardict_classification, msi, msilen, shift3, score_vardict]
            i_f.close()
    varscan2_info = {}
    if varscan2_vcfs:
        for varscan2_vcf in varscan2_vcfs:
            i_f = open_textfile(varscan2_vcf)
            for line in skip_empty(i_f):
                x = line.strip().split()
                chrom, pos, _, ref, alts, _, filters, info = x[0:8]
                varscan_classification = get_varscan2_info(info)

                # Replace the wrong "G/A" with the correct "G,A" in ALT
                # column:
                alts = alts.replace('/', ',')

                # multiple sequences in the REF, as is the case in
                # VarScan2's indel output:
                ref = re.sub(r'[^\w].*$', '', ref)

                # Get rid of non-compliant characters in the ALT column:
                alts = re.sub(r'[^\w,.]', '', alts)

                # Eliminate dupliate entries in ALT:
                alts = re.sub(r'(\w+),\1', r'\1', alts)

                # VarScan2 output a line with REF allele as "M"
                if re.search(r'[^GCTAU]', ref, re.I):
                    continue

                for alt in alts.split(","):
                    if ref != alt:
                        var_id = "-".join([chrom, pos, ref, alt])
                        varscan2_info[var_id] = varscan_classification
            i_f.close()

    muse_info = {}
    if muse_vcfs:
        for muse_vcf in muse_vcfs:
            i_f = open_textfile(muse_vcf)
            for line in skip_empty(i_f):
                x = line.strip().split()
                chrom, pos, _, ref, alts, _, filters, info = x[0:8]
                muse_classification = get_muse_info(filters)
                for alt in alts.split(","):
                    if ref != alt:
                        var_id = "-".join([chrom, pos, ref, alt])
                        muse_info[var_id] = muse_classification
            i_f.close()

    somaticsniper_info = {}
    if somaticsniper_vcfs:
        for somaticsniper_vcf in somaticsniper_vcfs:
            i_f = open_textfile(somaticsniper_vcf)
            for line in skip_empty(i_f):
                x = line.strip().split()
                chrom, pos, _, ref, alts, _, filters, info, fields = x[0:9]
                samples = x[9:]
                ref = re.sub(r'[^GCTA]', 'N', ref, flags=re.I)
                somaticsniper_classification, score_somaticsniper = get_somaticsniper_info(
                    fields, samples, idxT)
                for alt in alts.split(","):
                    if ref != alt:
                        var_id = "-".join([chrom, pos, ref, alt])
                        somaticsniper_info[var_id] = [
                            somaticsniper_classification, score_somaticsniper]
            i_f.close()

    features = {}
    for var_id in (set(mutect2_info.keys()) | set(strelka2_info.keys()) | set(vardict_info.keys()) |
                   set(varscan2_info.keys()) | set(somaticsniper_info.keys()) | set(muse_info.keys())):
        num_callers = 0
        if var_id in mutect2_info:
            mutect_classification, nlod, tlod, tandem, ecnt = mutect2_info[
                var_id]
            num_callers += mutect_classification
        else:
            mutect_classification = 0
            nlod = tlod = tandem = ecnt = nan

        if var_id in strelka2_info:
            strelka_classification, somatic_evs, qss, tqss = strelka2_info[
                var_id]
            num_callers += strelka_classification
        else:
            strelka_classification = 0
            somatic_evs = qss = tqss = nan

        if var_id in vardict_info:
            vardict_classification, msi, msilen, shift3, score_vardict = vardict_info[
                var_id]
            num_callers += vardict_classification
        else:
            vardict_classification = 0
            msi = msilen = shift3 = score_vardict = nan

        if var_id in varscan2_info:
            varscan_classification = varscan2_info[var_id]
            num_callers += varscan_classification
        else:
            varscan_classification = 0

        if var_id in muse_info:
            muse_classification = muse_info[var_id]
            num_callers += muse_classification
        else:
            muse_classification = 0

        if var_id in somaticsniper_info:
            somaticsniper_classification, score_somaticsniper = somaticsniper_info[
                var_id]
            num_callers += somaticsniper_classification
        else:
            somaticsniper_classification = 0
            score_somaticsniper = nan

        if num_callers >= min_caller:
            features[var_id] = [mutect_classification, nlod, tlod, tandem, ecnt,
                                strelka_classification, somatic_evs, qss, tqss,
                                vardict_classification, msi, msilen, shift3, score_vardict,
                                varscan_classification,
                                muse_classification,
                                somaticsniper_classification, score_somaticsniper]

    chrom_order = get_chromosomes_order(reference)
    ordered_vars = sorted(features.keys(), key=lambda x: [
                          chrom_order["-".join(x.split("-")[:-3])], int(x.split("-")[1])])
    n_variants = len(ordered_vars)
    logger.info("Number of variants: {}".format(n_variants))
    header = ["CHROM", "POS", "ID", "REF", "ALT", "if_MuTect", "if_VarScan2", "if_SomaticSniper", "if_VarDict", "MuSE_Tier",
              "if_Strelka", "Strelka_Score", "Strelka_QSS",
              "Strelka_TQSS", "Sniper_Score", "VarDict_Score",
              "M2_NLOD", "M2_TLOD", "M2_STR", "M2_ECNT", "MSI", "MSILEN", "SHIFT3"]

    with open(output_tsv, "w") as o_f:
        o_f.write("\t".join(header) + "\n")
        for var_id in ordered_vars:
            mutect_classification, nlod, tlod, tandem, ecnt, \
                strelka_classification, somatic_evs, qss, tqss, \
                vardict_classification, msi, msilen, shift3, score_vardict, \
                varscan_classification, \
                muse_classification, \
                somaticsniper_classification, score_somaticsniper = features[
                    var_id]

            f = [mutect_classification, varscan_classification, somaticsniper_classification,
                 vardict_classification, muse_classification, strelka_classification,
                 somatic_evs, qss, tqss,
                 score_somaticsniper, score_vardict,
                 nlod, tlod, tandem, ecnt,
                 msi, msilen, shift3]
            chrom = "-".join(var_id.split("-")[:-3])
            pos, ref, alt = var_id.split("-")[-3:]
            o_f.write(
                "\t".join([chrom, pos, ".", ref, alt] + list(map(lambda x: str(x).replace("nan", "0"), f))) + "\n")

    logger.info("Done Reading Callers' Features.")
    return output_tsv


if __name__ == '__main__':
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description='extract extra features for standalone mode')
    parser.add_argument('--reference', type=str, help='reference fasta filename',
                        required=True)
    parser.add_argument('--output_tsv', type=str, help='output features tsv',
                        required=True)
    parser.add_argument('--mutect2_vcfs', type=str, nargs="*",
                        help='MuTect2 VCFs',
                        default=None)
    parser.add_argument('--strelka2_vcfs', type=str, nargs="*",
                        help='Strelka2 VCFs',
                        default=None)
    parser.add_argument('--varscan2_vcfs', type=str, nargs="*",
                        help='VarScan2 VCFs',
                        default=None)
    parser.add_argument('--muse_vcfs', type=str, nargs="*",
                        help='MuSE VCFs',
                        default=None)
    parser.add_argument('--vardict_vcfs', type=str, nargs="*",
                        help='VarDict VCFs',
                        default=None)
    parser.add_argument('--somaticsniper_vcfs', type=str, nargs="*",
                        help='SomaticSniper VCFs',
                        default=None)
    parser.add_argument('--min_caller', type=float,
                        help='Number of minimum callers support a call',
                        default=0.5)
    args = parser.parse_args()
    logger.info(args)

    try:
        output = read_callers_vcf(args.reference,
                                  args.output_tsv,
                                  args.mutect2_vcfs,
                                  args.strelka2_vcfs,
                                  args.varscan2_vcfs,
                                  args.muse_vcfs,
                                  args.vardict_vcfs,
                                  args.somaticsniper_vcfs,
                                  args.min_caller,
                                  )
        if output is None:
            raise Exception("read_callers_vcf failed!")
    except Exception as e:
        logger.error(traceback.format_exc())
        logger.error("Aborting!")
        logger.error(
            "read_callers_vcf.py failure on arguments: {}".format(args))
        raise e
