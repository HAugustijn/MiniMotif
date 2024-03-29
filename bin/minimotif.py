#!/usr/bin/env python3

from minimotif_scripts import *
import argparse
import os
import json as jsonlib
from rich.console import Console
from datetime import datetime


parser = argparse.ArgumentParser(description="", usage='''
-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
Generic command: 
     minimotif.py -i [binding_profiles] -G [genome] -O [outdir]
_____________________________________________________________________________________________
Mandatory arguments:
    -G  Provide a genbank file of the genome of interest for TFBS scanning
    -O  Put the path to the output folder where results should be deposited

Optional arguments:
    -i  Provide the binding profiles in fasta format
    -w  Minimal width of the meme detection module. Default: 10
    -ps Pseudocount used to generate the PWM matrices. Default: 0.1
    -l  Use this flag to output .png sequence logo files
    -co Include this flag to detect TFBSs occurrences in coding regions
    -r  Range of the regulatory region. Default: -350 50
    -c  Range between genes that are considered to be co-regulated. Default: -50 40
    -p  P-value threshold used for the PWM detection module. Default: 0.00001
    -pc Add this flag to run on pre-calculated PWM matrices
    -b  Run MOODS in batch mode. In this mode, the p-value is not separately 
        calculated which increases the run speed. Default: True
    -m  Mode for the HMM detection module. Options: spacer_masking or positional_masking.
        Positional_masking masks nucleotides individually, if their information content 
        is over the given threshold. Spacer_masking assumes that nucleotides belonging 
        to -10 and -35 regions are significantly more conserved than the spacer nucleotides.
        Default: spacer_masking
    -ic Information content threshold. Default: 1.0
    -la Adjust the length of the alignments that are outputted from the 
        script, in comparison with full alignments. The default is 1 nucleotide
        less than the global alignment between pHMM models and the query
        sequence. Default: 1
    -am Analysis mode. Default: auto (gapped, ungapped, both)
-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
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
parser.add_argument("-sm", "--spacermode", help=argparse.SUPPRESS, choices=['spacer_masking', 'positional_masking',
                                                                            "no_masking"],
                    required=False, default="spacer_masking")
parser.add_argument("-ic", "--ic_threshold", type=float,
                    help=argparse.SUPPRESS, required=False, default=1.0)
parser.add_argument("-la", "--adjust_length", type=float,
                    help=argparse.SUPPRESS, required=False, default=1)
parser.add_argument('-am', '--analysis_mode', help= argparse.SUPPRESS, choices=['auto', 'gapped', 'ungapped', 'both'],
                    default='auto', required=False)


args = parser.parse_args() # parse arguments
GB_REGIONS = {}  # dictionary to store the genbank regions
GB_STRAND = []  # dict to store the strand orientation
PROD = []  # store the product dict
console = Console() # Create console instance

if not os.path.exists(args.outdir):
    os.mkdir(args.outdir)

# TODO: This is optimised to receive multiple genbank files, shall the documentation be updated?

for genbank_file in args.genbank:
    # functions are part of parse_genbank.py
    gb_name = genbank_file.split("/")[-1].split(".")[0]
    if is_gbk(genbank_file):  # check gbk format
        console.print(f"[bold cyan]{ datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Extracting sequences from genbank"
                      f" file: {gb_name}[/bold cyan]")
        reg_region, co_region, up_region, complete_seq, gene_strand_dict, product_dict = parse_gb(genbank_file, args.coregion,
                                                                                       args.regregion)
        GB_REGIONS[genbank_file] = complete_seq
        GB_STRAND.append(gene_strand_dict)
        PROD.append(product_dict)
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
                moods_results = run_moods(co_fasta, pwm_file, "co", reg, args.outdir, bg_dis,
                                          gb_name, args.pvalue, args.batch)
                parse_moods(moods_results, reg, gb_name, "co", thresholds_pwm, args.outdir)
            moods_results = run_moods(reg_fasta, pwm_file, "reg", reg, args.outdir, bg_dis,
                                      gb_name, args.pvalue, args.batch)
            parse_moods(moods_results, reg, gb_name, "reg", thresholds_pwm, args.outdir)
    extensions_to_move = [".sto", ".hmm", ".fasta", ".meme", ".moods", "PWM.tsv"]
    movetodir(args.outdir + os.sep, "bin", extensions_to_move)


else:
    if args.input:
        for input_file in args.input:
            # functions are part of preprocessing_bindingsites.py
            reg_name = input_file.split("/")[-1].split(".")[0]
            if is_fasta(input_file):  # check fasta format
                console.print(
                    f"[bold cyan]{ datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Analysing regulator:"
                    f" {reg_name}[/bold cyan]")

                # run MEME to extract the binding motif from the input sequences
                meme_results = run_meme(input_file, args.outdir, args.min_width)
                con_motif, motifs = parse_meme(meme_results)
                if not con_motif:
                    pass
                else:
                    # functions are part of generate_sequence_logo.py
                    pfm = create_pfm(motifs)  # construct position frequency matrix (PFM)
                    shan_ic = get_ic(pfm)  # calculate information content (IC)
                    pfm_file = f"{args.outdir}/{reg_name}_PFM.tsv"
                    pfm2 = (pd.DataFrame.from_dict(pfm)).T
                    pfm2.to_csv(pfm_file, sep="\t", index=False, index_label=False, header=False)
                    if args.logo:
                        create_logo(shan_ic, reg_name, args.outdir)

                    # functions are part of method_selection.py
                    mean_ic = calc_mean_ic(shan_ic)
                    for i, genbank_file in enumerate(args.genbank):
                        if args.analysis_mode == 'auto':
                            # TODO: Maybe also set an upper limit for the number of positions per motif?
                            if len(shan_ic) >= 10:
                                if mean_ic < 0.8:  # functions are part of detection_hmm.py
                                    console.print(
                                        f"[bold cyan]{ datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Running HMM"
                                        f" detection for a gapped sequence motif[/bold cyan]")
                                    hmm_models = prep_hmm_detection(input_file, reg_name, args.spacermode,
                                                                    args.ic_threshold, args.outdir)
                                    run_hmm_detection(genbank_file, reg_name, hmm_models, args.coding,
                                                      args.adjust_length, args.outdir, GB_STRAND[i], PROD[i])

                                elif 0.8 <= mean_ic <= 1.2:  # functions are part of both detection_pwm/hmm.py
                                    console.print(
                                        f"[bold cyan]{ datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Running both HMM"
                                        f" and PWM detection[/bold cyan]")
                                    ic_threshold = 1.4
                                    mode = "positional_masking"
                                    run_pwm_detection(genbank_file, pfm, args.pseudocount, reg_name,
                                                      GB_REGIONS[genbank_file], args.coding, args.pvalue, args.batch,
                                                      args.outdir)
                                    hmm_models = prep_hmm_detection(input_file, reg_name, mode, ic_threshold,
                                                                    args.outdir)
                                    run_hmm_detection(genbank_file, reg_name, hmm_models, args.coding,
                                                      args.adjust_length, args.outdir, GB_STRAND[i], PROD[i])

                                elif mean_ic >= 1.2:  # functions are part of detection_pwm.py
                                    console.print(
                                        f"[bold cyan]{ datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Running PWM "
                                        f"detection for a ungapped sequence motif[/bold cyan]")
                                    run_pwm_detection(genbank_file, pfm, args.pseudocount, reg_name,
                                                      GB_REGIONS[genbank_file], args.coding, args.pvalue, args.batch,
                                                      args.outdir)

                            else:  # run PWM detection for short sequence motifs
                                console.print(
                                    f"[bold cyan]{ datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Running PWM detection"
                                    f" for a short sequence motif[/bold cyan]")
                                run_pwm_detection(genbank_file, pfm, args.pseudocount, reg_name,
                                                  GB_REGIONS[genbank_file], args.coding, args.pvalue, args.batch,
                                                  args.outdir)
                        elif args.analysis_mode == 'gapped':
                            console.print(
                                f"[bold cyan]{ datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Running HMM detection "
                                f"for a gapped sequence motif[/bold cyan]")
                            hmm_models = prep_hmm_detection(input_file, reg_name, args.spacermode, args.ic_threshold,
                                                            args.outdir)
                            run_hmm_detection(genbank_file, reg_name, hmm_models, args.coding, args.adjust_length,
                                              args.outdir, gene_strand_dict, product_dict)
                        elif args.analysis_mode == 'ungapped':
                            console.print(
                                f"[bold cyan]{ datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Running PWM detection "
                                f"for a ungapped sequence motif[/bold cyan]")
                            run_pwm_detection(genbank_file, pfm, args.pseudocount, reg_name, GB_REGIONS[genbank_file],
                                              args.coding, args.pvalue, args.batch, args.outdir)
                        elif args.analysis_mode == 'both':
                            console.print(
                                f"[bold cyan]{ datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Running both HMM and PWM"
                                f" detection[/bold cyan]")

                            ic_threshold = 1.4 # TODO: How is this defined?
                            mode = "positional_masking"
                            run_pwm_detection(genbank_file, pfm, args.pseudocount, reg_name,
                                              GB_REGIONS[genbank_file], args.coding, args.pvalue, args.batch,
                                              args.outdir)
                            hmm_models = prep_hmm_detection(input_file, reg_name, mode, ic_threshold, args.outdir)
                            run_hmm_detection(genbank_file, reg_name, hmm_models, args.coding,
                                              args.adjust_length, args.outdir, gene_strand_dict, product_dict)

        extensions_to_move = [".sto", ".hmm", ".fasta", ".meme", ".moods", "PWM.tsv"]
        movetodir(args.outdir + os.sep, "bin", extensions_to_move)

    else:
        console.print(
            f"[bold red]{ datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Please provide the input binding profiles in"
            f" fasta format or use the -pc flag to run precalculated PWMs[/bold red]")
