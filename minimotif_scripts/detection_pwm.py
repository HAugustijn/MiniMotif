"""  """

import pandas as pd
import seqlogo
from collections import Counter
import os
import subprocess


def create_pwm(pfm, pscount, outdir, reg_name):
    """ creates a position wight matrix of a position frequency matrix"""
    pwm_file = f"{outdir}/{reg_name}_PWM.tsv"
    pwm = (seqlogo.pfm2pwm(pd.DataFrame.from_dict(pfm), background=None, pseudocount=pscount)).T
    pwm.to_csv(pwm_file, sep="\t", index=False, index_label=False, header=False)
    return pwm_file, pwm


def get_sequence_gc_content(sequence):
    """ calculate the GC content of the bgc region's sequence """
    counter = Counter(sequence)
    return sum(counter[char] for char in "GC") / len(sequence)


def get_bg_distribution(sequence):
    """ calculate background distribution based on GC percentage of the
    query sequence """
    gc_percentage = get_sequence_gc_content(sequence) / 2
    at_percentage = 0.5 - gc_percentage
    percentage = ' '.join(map(str, (at_percentage, gc_percentage,
                                    gc_percentage, at_percentage)))
    return percentage


def run_moods(fasta_file, pwm, reg_type, reg_name, outdir, bg_dist,
              gb_name, threshold, batch):
    """ Run moods to scan sequences with a PWM """

    outfile = f"{outdir}/{reg_name}_{gb_name}_{reg_type}.moods"

    if batch:
        cmd_moods = f"moods-dna.py -S {pwm} -s {fasta_file} -o {outfile} " \
                f"-p {threshold} --bg {bg_dist} --batch"
    else:
        cmd_moods = f"moods-dna.py -S {pwm} -s {fasta_file} -o {outfile} " \
                f"-p {threshold} --bg {bg_dist}"

    try:
        if not os.path.exists(outfile):
            res_map = subprocess.check_output(cmd_moods, shell=True,
                                              stderr=subprocess.STDOUT)
    except(subprocess.CalledProcessError):
        print(f'Unable to run moods with command: {cmd_moods}')
    return outfile


def parse_moods(moods_results, reg_name, gb_name, reg_type, thres, outdir):
    """ parse results from MOODS output """
    out_dict = {}
    with open(moods_results, "r") as moods_output:
        for line in moods_output:
            ID, PWM_name, loc, strand, score, seq, no = line.split(",")
            region = ID.split("~")[0]
            start_loc = ID.split("~")[1].split(":")[0]
            full_start_loc = int(start_loc) + int(loc)
            full_end_loc = int(start_loc) + int(loc) + int(len(seq))
            full_loc = f"{full_start_loc}:{full_end_loc}"
            conf = set_confidence(float(score), thres)
            if region not in out_dict.keys():
                out_dict[region] = []
                out_dict[region].append([full_loc, strand, score, conf, seq])
            else:
                out_dict[region].append([full_loc, strand, score, conf, seq])
    write_output(out_dict, reg_name, gb_name, reg_type, outdir)
    return


def set_confidence(score, thresholds):
    """ From the moods results score and pwm threshold, determine the confidence """
    med_str_thres = (thresholds[1] + thresholds[0]) / 2
    if score <= thresholds[1]:
        return "weak"
    if score >= med_str_thres:
        return "strong"
    return "medium"


def set_threshold(pwm):
    """ set a additional threshold to separate top hits from medium """
    max_pwm = sum(pwm.max(numeric_only=True))
    strict_threshold = 0
    for val in (pwm.max().to_frame()).values:
        if val >= 1.9:
            strict_threshold = strict_threshold + val
    return [max_pwm, strict_threshold[0]]


def write_output(results_dict, reg_name, gb_name, reg_type, outdir):
    """ Writes a tsv output file """
    with open(f"{outdir}/{reg_name}_{gb_name}_{reg_type}_pwm_results.tsv", "w") as outfile:
        outfile.write("Region\tLocation\tStrand\tScore\tConfidence\tSequence\n")
        for key in results_dict.keys():
            for i in results_dict[key]:
                outfile.write(f"{key}\t{i[0]}\t{i[1]}\t{i[2]}\t{i[3]}\t{i[4]}\n")
    return


def run_pwm_detection(genbank_file, pfm, pseudocount, reg_name, gbk_regions, coding, pvalue, batch, outdir):
    """ Run MOODS to detect TFBS occurrences """
    gb_name = genbank_file.split("/")[-1].split(".")[0]
    bg_dis = get_bg_distribution(gbk_regions)  # bg distribution of the full genome

    pwm_file, pwm = create_pwm(pfm, pseudocount, outdir, reg_name)
    thresholds_pwm = set_threshold(pwm)
    if coding:
        co_fasta = f"{outdir}/{gb_name}_co_region.fasta"
        moods_results = run_moods(co_fasta, pwm_file, "co", reg_name, outdir, bg_dis, gb_name, pvalue, batch)
        parse_moods(moods_results, reg_name, gb_name, "co", thresholds_pwm, outdir)

    reg_fasta = f"{outdir}/{gb_name}_reg_region.fasta"
    moods_results = run_moods(reg_fasta, pwm_file, "reg", reg_name, outdir, bg_dis, gb_name, pvalue, batch)
    parse_moods(moods_results, reg_name, gb_name, "reg", thresholds_pwm, outdir)
