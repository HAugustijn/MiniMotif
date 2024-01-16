import os
import sys
import time
import logging
from rich.console import Console
from datetime import datetime
import argparse

class Logger:
    def __init__(self, log_file):
       self.console = Console()
       self.log_file = log_file
       self.logger = logging.getLogger('Logger')
       self.logger.setLevel(logging.INFO)
       handler = logging.FileHandler(self.log_file)
       handler.setLevel(logging.INFO)
       formatter = logging.Formatter('%(asctime)s - %(message)s')
       handler.setFormatter(formatter)
       self.logger.addHandler(handler)

    def log(self, message):
       self.console.log(message)
       stripped_message = message.replace("[bold cyan]", "").\
           replace("[/bold cyan]", "").replace("[bold red]", "").replace("[/bold red]", "")
       self.logger.info(stripped_message)

    def log_only(self, message):
         stripped_message = message.replace("[bold cyan]", "").\
              replace("[/bold cyan]", "").replace("[bold red]", "").replace("[/bold red]", "")
         self.logger.info(stripped_message)

    def print_only(self, message):
        self.console.log(message)

    def log_error(self, message):
       self.logger.info(message)

    def start_timer(self):
       self.start_time = time.time()

    def stop_timer(self):
       end_time = time.time()
       execution_time = end_time - self.start_time
       self.log(f"Execution time: {execution_time} seconds")


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


args = parser.parse_args()
if not os.path.exists(f"{args.outdir}"):
    os.mkdir(f"{args.outdir}")
logger = Logger(f"{args.outdir}/{args.outdir}_minimotif_log.txt")