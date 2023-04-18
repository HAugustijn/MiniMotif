#!/usr/bin/env python3

from MiniMotif.minimotif_scripts import *
import argparse
import os
import json as jsonlib

parser = argparse.ArgumentParser(description="", usage='''
-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
Generic command: 
     minimotif.py -I [binding_profiles] -G [genome] -O [outdir]
_________________________________________________________________________
Mandatory arguments:
    -I  Provide the binding profiles in fasta format
    -G  Provide a genbank file of the genome of interest for TFBS scanning
    -O  Put the path to the output folder where results should be deposited

Optional arguments:

-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
''')

parser.add_argument("-i", "--input", nargs='+', help=argparse.SUPPRESS,
                    required=False)
parser.add_argument("-O", "--outdir", help=argparse.SUPPRESS,
                    required=True)
parser.add_argument("-G", "--genbank", nargs='+', help=argparse.SUPPRESS,
                    required=True)
parser.add_argument("-w", "--min_width", help=argparse.SUPPRESS,
                    required=False, default=10)
parser.add_argument("-ps", "--pseudocount", help=argparse.SUPPRESS,
                   required=False, default=0.1)
parser.add_argument("-l", "--logo", help=argparse.SUPPRESS,
                    required=False, action='store_true')
parser.add_argument("-co", "--coding", help=argparse.SUPPRESS,
                    required=False, action='store_true')
parser.add_argument("-r", "--regregion", nargs='+', help=argparse.SUPPRESS,
                    required=False, default=[-350, 50])
parser.add_argument("-c", "--coregion", nargs='+', help=argparse.SUPPRESS,
                    required=False, default=[-50, 40])
parser.add_argument("-p", "--pvalue", help=argparse.SUPPRESS,
                    required=False, default=0.00001)
parser.add_argument("-pc", "--precal", help=argparse.SUPPRESS,
                    required=False, action='store_true')
parser.add_argument("-b", "--batch", help=argparse.SUPPRESS,
                    required=False, default=True)
parser.add_argument("-m", "--mode", help=argparse.SUPPRESS,
                    required=False, default="spacer_masking")
parser.add_argument("-ic", "--ic_threshold", type=float,
                    help=argparse.SUPPRESS, required=False, default=1.0)
parser.add_argument("-ladj", "--adjust_length", type=float,
                    help=argparse.SUPPRESS, required=False, default=1)

args = parser.parse_args()

GB_REGIONS = {}

if not os.path.exists(args.outdir):
    os.mkdir(args.outdir)


for genbank_file in args.genbank:
    # functions are part of parse_genbank.py
    gb_name = genbank_file.split("/")[-1].split(".")[0]
    if is_gbk(genbank_file):  # check gbk format
        print(f"Extracting sequences from genbank file: {gb_name}")
        reg_region, co_region, complete_seq = extract_regions(genbank_file, args.coregion, args.regregion)
        GB_REGIONS[genbank_file] = complete_seq
        if args.coding:
            write_fastas(co_region, "co", gb_name, args.outdir)
        write_fastas(reg_region, "reg", gb_name, args.outdir)

if args.precal:  # Run detection on precalculated PWMs
    base = os.path.dirname(os.path.abspath(__file__))
    path_pwm = os.path.join(base, 'data/pwms.json')
    with open(path_pwm, "r", encoding="utf-8") as file:
        regulators = jsonlib.load(file)
    for reg in regulators.keys():
        pwm_file = f"{args.outdir}/{reg}_PWM.tsv"
        pwm = pd.DataFrame(regulators[reg]['pwm'])
        pwm.to_csv(pwm_file, sep="\t", index=False, index_label=False, header=False)

        thresholds_pwm = [regulators[reg]['max_score'], regulators[reg]['min_score']]

        for genbank_file in args.genbank:
            gb_name = genbank_file.split("/")[-1].split(".")[0]
            bg_dis = get_bg_distribution(GB_REGIONS[genbank_file])
            reg_fasta = f"{args.outdir}/{gb_name}_reg_region.fasta"
            if args.coding:
                co_fasta = f"{args.outdir}/{gb_name}_co_region.fasta"
                moods_results = run_moods(reg_fasta, pwm_file, "co", reg, args.outdir, bg_dis,
                                          gb_name, args.pvalue, args.batch)
                parse_moods(moods_results, reg, gb_name, "co", thresholds_pwm, args.outdir)
            moods_results = run_moods(reg_fasta, pwm_file, "reg", reg, args.outdir, bg_dis,
                                      gb_name, args.pvalue, args.batch)
            parse_moods(moods_results, reg, gb_name, "reg", thresholds_pwm, args.outdir)

else:
    if args.input:
        for input_file in args.input:
            # functions are part of preprocessing_bindingsites.py
            reg_name = input_file.split("/")[-1].split(".")[0]
            if is_fasta(input_file):  # check fasta format
                print(f"Analysing regulator: {reg_name}")
                # run MEME to extract the binding motif from the input sequences
                meme_resutls = run_meme(input_file, args.outdir, args.min_width)
                con_motif, motifs = parse_meme(meme_resutls)

                # functions are part of generate_sequence_logo.py
                pfm = create_pfm(motifs)  # construct position frequency matrix (PFM)
                shan_ic = get_ic(pfm)  # calculate information content (IC)
                if args.logo:
                    create_logo(shan_ic, reg_name, args.outdir)

                # functions are part of method_selection.py
                mean_ic = calc_mean_ic(shan_ic)
                for genbank_file in args.genbank:
                    if mean_ic < 0.8:  # functions are part of detection_hmm.py
                        print(f"Running HMM detection for a gapped sequence motif")
                        hmm_models = prep_hmm_detection(input_file, reg_name, args.mode,
                                                        args.ic_threshold, args.outdir)
                        run_hmm_detection(genbank_file, reg_name, hmm_models, args.coding,
                                          args.adjust_length, args.outdir)

                    elif 0.8 <= mean_ic <= 1.2:  # functions are part of both detection_pwm/hmm.py
                        print(f"Running both HMM and PWM detection")
                        ic_threshold = 1.4
                        mode = "positional_masking"
                        run_pwm_detection(genbank_file, pfm, args.pseudocount, reg_name,
                                      GB_REGIONS[genbank_file], args.coding, args.pvalue, args.batch, args.outdir)
                        hmm_models = prep_hmm_detection(input_file, reg_name, mode, ic_threshold, args.outdir)
                        run_hmm_detection(genbank_file, reg_name, hmm_models, args.coding,
                                          args.adjust_length, args.outdir)

                    elif mean_ic >= 1.2:  # functions are part of detection_pwm.py
                        print(f"Running PWM detection for a ungapped sequence motif")
                        run_pwm_detection(genbank_file, pfm, args.pseudocount, reg_name,
                                      GB_REGIONS[genbank_file], args.coding, args.pvalue, args.batch, args.outdir)
        extensions_to_move = [".sto", ".hmm", ".fasta", ".meme", ".moods", "PWM.tsv"]
        movetodir(args.outdir + os.sep, "bin", extensions_to_move)
    else:
        print('Please provide the input binding profiles in fasta format or use the -pc flag to run precalculated PWMs')
