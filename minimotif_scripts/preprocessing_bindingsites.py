"""  """

from Bio import SeqIO
import os
import subprocess
import sys
from rich.console import Console
from datetime import datetime
from minimotif_scripts.logger import logger

def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        if not any(fasta):
            logger.log(
                f"[bold red]Please provide the input in fasta"
                f" format of file: {filename}[/bold red]")
            return False
        else:
            return True

def check_fasta_length(fasta_file, min_width):
    """ Check the length of the input sequences
    """
    with open(fasta_file, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        for record in fasta:
            if len(record.seq) < min_width:
                return False
    return True


def run_meme(fasta_file, outdir, min_width):
    """ Run MEME to identify motifs in binding sequences """
    reg_name = fasta_file.split("/")[-1].split(".")[0]
    out_file = f"{outdir}/{reg_name}.meme"
    cmd_meme = f"meme {fasta_file} -o {outdir} -text -dna -minw " \
               f"{min_width} -nmotifs 3 -evt 0.05 -revcomp > {out_file}"
    try:
        if not os.path.exists(out_file):
            subprocess.check_output(cmd_meme, shell=True,
                                    stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        logger.log(
            f"[bold red]Unable to run meme with"
            f" command: {cmd_meme}[/bold red]")
        sys.exit()
    return out_file


def parse_motif_coordinates(meme_out):
    """ Extract the coordinates and consensus of the motif from MEME results """
    consensus_motif = ""
    start_motif = 0
    end_motif = 0
    for num, line in enumerate(meme_out, 1):
        if "E-value = " in line:
            consensus_motif = line.split()[1]
        elif line.startswith("BL   MOTIF"):
            start_motif = num
        elif line.startswith("//"):
            end_motif = num
    return consensus_motif, start_motif, end_motif


def parse_meme(meme_results):
    """ parses the meme output to lists of motifs and consensus sequences """
    motif_list = []
    with open(meme_results, "r") as meme_out:
        consensus_motif, start_motif, end_motif = parse_motif_coordinates(meme_out)
        if not consensus_motif:
            logger.log(
                f"[bold red]Could not parse results from the "
                f"MEME output. Please check the MEME output[/bold red]")
            return False, False

        else:
            # seek the index in the file and extract the motifs
            meme_out.seek(0)
            txt = meme_out.readlines()
            motifs = txt[start_motif: end_motif - 1]
            for motif in motifs:
                motif_list.append(motif.split(" ")[-4])

            return consensus_motif, motif_list
