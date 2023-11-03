"""  """
from Bio import AlignIO
from operator import itemgetter
import os
import subprocess
import math
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import re
import pandas as pd
from rich.console import Console
from datetime import datetime


console = Console()

def align_sequences(input_fasta, reg_name, outdir):
    """Align the sequences in a FASTA file using the MAFFT algorithm.

    input_fasta: str, the path to the input FASTA file
    reg_name: str, the name of the regulator
    outdir: str, the directory path to save the output files
    """

    outfile = f"{outdir}/{reg_name}.fasta"
    cmd_mafft = f"mafft --auto {input_fasta} > {outfile}"
    try:
        if not os.path.exists(outfile):
            subprocess.check_output(cmd_mafft, shell=True, stderr=subprocess.STDOUT)

    except Exception as e:
        raise Exception(console.print(f"[bold red]{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Unable to align "
                                      f"sequences with mafft: {str(e)}[/bold red]"))


def convert_to_stockholm(reg_name, outdir, al_format='stockholm'):
    """Converts a given FASTA file to a Stockholm format alignment file and returns the parsed MSA object.

    reg_name: str, the name of the regulator
    outdir: str, the directory path to save the output files
    al_format: str, the format of the alignment file (default: 'stockholm')

    Returns:
    - An AlignIO.MSA object representing the alignment.
    """

    AlignIO.convert(f'{outdir}/{reg_name}.fasta', 'fasta', f'{outdir}/{reg_name}.sto', 'stockholm')

    return AlignIO.read(f'{outdir}/{reg_name}.sto', al_format)


def calculate_gap_freq(msa_obj):
    """ Calculate gap frequency per position of an MSA object

    Arguments:
        msa_obj: An alignIO.MSA object

    Returns:
        A list of floats representing the alignment and gap frequencies
    """
    gap_freq_list = []

    for col_no in range(len(list(msa_obj[0]))):
        list_input = list(msa_obj[:, col_no])
        no_of_gaps = list_input.count('-')
        length = len(list_input)
        proportion_of_gaps = no_of_gaps / length
        gap_freq_list.append(proportion_of_gaps)

    return gap_freq_list


def gap_checker(gap_freq_list, gap_threshold=0):
    """ Checks if an MSA position has a gap frequency above a given threshold

    Arguments:
        gap_freq_list: List of floats, each float per pos gap frequency
        gap_threshold: float, default = 0

    Returns:
        A list of tuples, were each tuple has the gap index as first element
        and the gap freq as a second element. i.e. (22, 0.455)
    """
    gap_summ_list = []

    for index in range(0, len(gap_freq_list)):
        if gap_freq_list[index] > gap_threshold:
            gap_pos_freq = index, gap_freq_list[index]
            gap_summ_list.append(gap_pos_freq)
    # Sorted, at a latter stage we want to trim gap positions from last to
    # first, so the indexes of the trimmed position do not alter the nex
    # iteration
    gap_summ_list = sorted(gap_summ_list, key=itemgetter(0), reverse=True)

    return gap_summ_list


def shannon_entropy(column_nuc_list):
    """Calculate Shannon's Entropy of a column of the alignment

    :param column_nuc_list: A list of nucleotides(str) for one MSA position
    :return: sh_entropy: A position's Shannon's entropy
    The position's Shannon's Entropy calculated according to:
    Entropy =  - sum(px*log(px)), where px = frequency of each nucleotide
    """
    sh_entropy = 0
    all_chars = set(column_nuc_list)
    # total per column residues
    col_length = len(column_nuc_list)
    # Number of residues in column
    for char in all_chars:
        # p_i = No of nucleotides of type i/ No of nucleotides per column
        p_i = column_nuc_list.count(char) / float(col_length)
        # log2 tranformation of Pi & sum
        sh_entropy += - (p_i * (math.log(p_i, 2)))

    return sh_entropy


def sh_max_finder(msa_obj):
    """Finds ungapped position with max Shannon's entropy in an MSA

    :param msa_obj: AlignIO MSA object
    :return: max_ent_pos_sum: tuple of two elements: max_ent_pos[0] is the
    index of the first instance of a position with the higher information
    content in the MSA, max_ent_pos[1] is the corresponding information content
    value. In case that spacer characters are Ns, then the first column that
    contains N characters is outputted. and max_ent_pos[1] is "N case".
    """
    shannon_entropy_list_ungapped = []
    shannon_entropy_list_gapped = []

    # To avoid considering a position that contains gaps, for a max information
    # content position, we calculate the information content of ungapped
    # positions
    for col_no in range(len(list(msa_obj[0]))):
        column_nuc_list = list(msa_obj[:, col_no])
        shannon_entropy_list_gapped.append(
            shannon_entropy(column_nuc_list))
        if "N" not in column_nuc_list:
            if "-" not in column_nuc_list:
                shannon_entropy_list_ungapped.append(
                    shannon_entropy(column_nuc_list))
            # Extract the ungapped position with the max Shannon's entropy
            max_ent_pos_value = max(shannon_entropy_list_ungapped)
            max_ent_pos_index = shannon_entropy_list_gapped.index(
                max_ent_pos_value)
            max_ent_pos_sum = max_ent_pos_index, max_ent_pos_value
        else:
            if "-" not in column_nuc_list:
                # For this case, we only need to extract a column that is full
                # of N characters.
                max_ent_pos_value = 'N case'
                max_ent_pos_index = col_no
                max_ent_pos_sum = max_ent_pos_index, max_ent_pos_value
                break

    return max_ent_pos_sum



def create_msa_obj(np_array, seq_ids):
    """Creates an AlignIO MSA object from a numpy array & list of ids

    :param np_array: Numpy array that represents a DNA alignment
    :param seq_ids: List of IDs of the aligned sequences(str)
    :return: An AlignIO MSA object
    """
    seq_record_list = []
    counter = 0

    for x in np_array:
        dna_str = "".join(x)
        record_temp = SeqRecord(Seq(dna_str), id=seq_ids[counter],
                                description="")
        seq_record_list.append(record_temp)
        counter += 1

    return AlignIO.MultipleSeqAlignment(seq_record_list)


def msa_mod_length_editor(msa_obj, max_ent_pos, gap_summ_list):
    """Artificially create MSA alignments of variable spacer length

    :param msa_obj: AlignIO MSA object
    :param max_ent_pos: Tuple, (position with max information content,
    information content value)
    :param gap_summ_list: List of tuples,Each tuple-> (gap_index,gap_freq)
    :return: sequences_list, a list of MSA objects of variable spacer
    length.
    """
    sequences_list = []
    msa_obj.sort()  # Sort MSA object, more robust approach
    align_array = np.array([list(rec) for rec in msa_obj], dtype=str)
    # Store the IDs of the sequences to a list
    align_seq_ids = [rec.id for rec in msa_obj]
    # Parse the column with the highest information content in the alignment
    extra_column = align_array[:, max_ent_pos[0]]
    # Create a column with the highest information content, in a format that
    # can be concatenated with the rest of the sequence
    extra_column_array = np.array([list(nuc) for nuc in extra_column],
                                  dtype=str)
    # First, turn all gaps to high information content positions
    for gap in gap_summ_list:
        position_to_slice = gap[0]
        first_part = align_array[:, :position_to_slice]
        last_part = align_array[:, position_to_slice + 1:]
        align_array = np.concatenate((first_part, extra_column_array,
                                      last_part), axis=1)
    # Turn array to MSA object
    full_alignment = create_msa_obj(align_array, align_seq_ids)
    # Append the ungapped sequence that has the gaps replaced to a list
    sequences_list.append(full_alignment)
    # Now iteratively trim the gap positions
    for gap in gap_summ_list:
        position_to_slice = gap[0]
        first_part = align_array[:, :position_to_slice]
        last_part = align_array[:, position_to_slice + 1:]
        align_array = np.concatenate((first_part, last_part), axis=1)
        trimmed_alignment = create_msa_obj(align_array, align_seq_ids)
        sequences_list.append(trimmed_alignment)

    return sequences_list


def processed_msa_file_writer(sequences_list, reg_name, outdir):
    """Writes alignment files based on given alignment objects.

    :param sequences_list: List of MSA objects
    :param reg_name: name of regulator as string
    :param outdir: output directory
    :return: alignment_filenames, list of strs, each str an alignment
    filename that is created
    Naming convention: Processed_alignment_length_{alignment_length}.sto
    """
    alignment_filenames = []
    for msa_obj in sequences_list:
        al_length = msa_obj.get_alignment_length()
        filename = f"{reg_name}_Processed_alignment_length_{al_length}.sto"
        out_file = os.path.join(outdir, filename)
        with open(out_file, "w") as handle:
            AlignIO.write(msa_obj, handle, "stockholm")
            alignment_filenames.append(filename)
    return alignment_filenames


def run_alimask(msa_file, range_to_mask, outdir):
    """Run  HMMER flavor 'alimask' to mask a part of an MSA

    :param msa_file: AlignIO MSA object
    :param outdir: path to output directory
    :param range_to_mask: A str in the format "Number-Number"
    :return: alimask_out, Name of the alimask output file
    """
    filename = msa_file.split('.')
    alimask_out = f"{outdir}/{filename[0]}_masked.sto"
    cmd_alimask = f"alimask --dna --alirange {range_to_mask}" \
                  f" {outdir}/{msa_file} {alimask_out}"

    try:
        if not os.path.exists(alimask_out):
            subprocess.check_output(cmd_alimask, shell=True, stderr=subprocess.STDOUT)

    except Exception as e:
        raise Exception(console.print(
            f"[bold red]{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Failed to run Alimask: {str(e)}[/bold red]"))

    return alimask_out


def calculate_n_freq(msa_obj):
    """Calculates the frequency of N characters per MSA position

    :param msa_obj: An alignIO.MSA object
    :return: n_freq_list: A list of floats
    Each of the list elements corresponds to an alignment position,
    and each float to the position's N character frequency. i.e. (22, 1.0)
    """
    n_freq_list = []
    for col_no in range(len(list(msa_obj[0]))):
        list_input = list(msa_obj[:, col_no])
        no_of_ns = list_input.count('N')
        length = len(list_input)
        proportion_of_ns = no_of_ns / length
        n_freq_list.append(proportion_of_ns)

    return n_freq_list


def n_checker(n_freq_list):
    """Checks which positions of MSA contain 'N' characters and their frequency

    :param n_freq_list: A list of floats. Per position N frequency
    :return: n_index_list: A list of ints, Each int -> index of position
    with N character.
    """
    n_index_list = []
    for index in range(len(n_freq_list)):
        if n_freq_list[index] > 0:
            n_index_list.append(index)

    return n_index_list


def spacer_tracker(msa_obj, high_information_content_index_list):
    """Finds the positions in an MSA that belong to the spacer region.

    :param msa_obj: AlignIO MSA object
    :param high_information_content_index_list: Indices of high information
    content pos
    :return: spacer_list, List of ints, each int corresponds to an
    alignment position that will be masked at a latter stage.
    """
    gap_freq_list = calculate_gap_freq(msa_obj)
    gap_summ_list = gap_checker(gap_freq_list)
    n_freq_list = calculate_n_freq(msa_obj)
    n_index_list = n_checker(n_freq_list)

    # if there is a gapped position, add it to the list.
    if len(gap_summ_list) > 0:
        for gap in gap_summ_list:
            if gap[0] not in high_information_content_index_list:
                high_information_content_index_list.append(gap[0])
        high_information_content_index_list.sort()
    else:
        high_information_content_index_list.sort()
    # Do the same if there are positions with Ns
    if len(n_index_list) > 0:
        for n_index in n_index_list:
            if n_index not in high_information_content_index_list:
                high_information_content_index_list.append(n_index)
        high_information_content_index_list.sort()
    else:
        high_information_content_index_list.sort()
    # The resulting list contains the range of positions to be masked

    return high_information_content_index_list


def shannon_entropy_msa(msa_obj):
    """Calculates Shannnon's entropy across an MSA.

    :param msa_obj: AlignIO MSA object
    :return: shannon_entropy_list: List of floats <- Each element
    corresponds to a position of MSA, each float to per position Shannon's
    entropy
    """
    shannon_entropy_list = []

    for col_no in range(len(list(msa_obj[0]))):
        column_nuc_list = list(msa_obj[:, col_no])
        shannon_entropy_list.append(shannon_entropy(column_nuc_list))

    return shannon_entropy_list


def positions_with_high_information_content(shannon_entropy_list, spacer_en_threshold=1):
    """Checks which positions have information content over a given threshold.

    :param shannon_entropy_list:List of floats,per pos Shannon's entropy
    :param spacer_en_threshold: float, information content threshold, positions
    with higher information content are outputted.
    :return:high_information_content_index_list: List of ints, each int
    represents the index of the high information content position.
    """
    high_information_content_index_list = []

    for index in range(len(shannon_entropy_list)):
        if shannon_entropy_list[index] > spacer_en_threshold:
            high_information_content_index_list.append(index)

    return high_information_content_index_list


def spacer_masker(sequences_list, alignment_filenames, mode, outdir, information_content=1):
    """Masks the spacer region of an MSA.

    :param sequences_list: List of strings, each string is a sequence
    :param alignment_filenames: List of strings, each string is a filename
    :param mode: String, either "spacer_masking", "positional_masking" or "no_masking"
    :param outdir: String, path to output directory
    :param information_content: Float, information content threshold
    :return: masked_seq_filenames: Set of strings, each string is a filename
    """
    masked_seq_filenames = set()
    full_masked_string = ""
    range_to_mask = ""

    # Find the positions with high information content in the alignment
    for index, seq in enumerate(sequences_list):
        msa_s_information_content = shannon_entropy_msa(seq)
        spacer_high_ent_pos = positions_with_high_information_content(msa_s_information_content, information_content)
        spacer_pos = spacer_tracker(seq, spacer_high_ent_pos)
        if mode == "spacer_masking":
            # find smaller and largest index with high information content, mask everything in between them
            range_to_mask = f"{min(spacer_pos) + 1}-{max(spacer_pos) + 1}"
            masked_seq = run_alimask(alignment_filenames[index], range_to_mask, outdir)
            masked_seq_filenames.add(masked_seq)
        elif mode == "positional_masking":
            # Mode: positional_masking, mask every position that is over the threshold individually.
            full_masked_string = ",".join(f"{pos + 1}-{pos + 1}" for pos in spacer_pos)
            masked_seq = run_alimask(alignment_filenames[index], full_masked_string, outdir)
            masked_seq_filenames.add(masked_seq)
        elif not spacer_pos or mode == "no_masking":
            # Nothing to be masked: Alignment passes to the next step as it is
            masked_seq_filenames.add(f"{outdir}/{alignment_filenames[index]}")

    return list(masked_seq_filenames)


def masked_seq_name_editor(masked_seq_filenames):
    """Add an ID line to a  group of masked MSA files, parsed from filename

    :param masked_seq_filenames: List of masked  MSA files
    :return: None
    If there is not a GF ID line in the MSA description, hmmbuild fails to
    recognise that the input MSAs are different from each other.
    """
    for file in masked_seq_filenames:
        with open(file, 'r', encoding='utf-8') as mask_seq_file:
            lines = mask_seq_file.readlines()
            # Name the alignment according to the file name
            lines[1] = f"#=GF ID {file.split('/')[-1]}\n"

        with open(file, 'w', encoding='utf-8') as mask_seq_file:
            mask_seq_file.writelines(lines)

    return None


def msa_file_merger(masked_seq_filenames, reg_name, outdir):
    """Merges the masked MSA files to one file

    :param masked_seq_filenames: List of masked  MSA files
    :param reg_name: Regulator name
    :param outdir: output directory
    :return: full_msa_filename, the name of the merged MSA file
    MSAs are merged to one file, so hmmbuild can build multiple models from
    one input file.
    """
    full_msa_filename = f"{outdir}/{reg_name}_msa_db.sto"
    with open(full_msa_filename, 'w') as msa_db_file:
        for file in masked_seq_filenames:
            with open(file, "r") as masked_msa_file:
                msa_db_file.write(masked_msa_file.read())
            msa_db_file.write('\n\n')

    return full_msa_filename


def run_hmmbuild(ali_out, hmmbuild_out):
    """Creates hmm models by running  HMMER flavor 'hmmbuild'

    :param ali_out: A masked MSA file or a multiple MSA file
    :param hmmbuild_out: The file that the output will be saved(.hmm)
    :return: hmmbuild_out:  The file that the output will be saved(.hmm)
    """

    cmd_hmmbuild = f"hmmbuild --dna {hmmbuild_out} {ali_out}"

    try:
        if not os.path.exists(hmmbuild_out):
            subprocess.check_output(cmd_hmmbuild, shell=True, stderr=subprocess.STDOUT)
    except Exception as e:
        raise Exception(console.print(
            f"[bold red]{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} -Unable to run hmmbuild: {str(e)}[/bold red]"))


    return hmmbuild_out


def model_creator(full_msa_filename, reg_name, outdir):
    """Wrapper function that creates and names HMM profile models

    :param full_msa_filename: File with multiple masked MSAs
    :param reg_name: name of regulator as string
    :param outdir: output directory
    :return: hmm_full_file, filename of hmmbuild output
    """
    hmm_models = run_hmmbuild(full_msa_filename, f"{outdir}/{reg_name}.hmm")
    return hmm_models

#TODO: Note: The length of each hmm model is also parsed from detection_hmm.py> motif_length_parser() function
def parse_hmm(hmm_models):
    """parses the sigma_factor.hmm to extract the length of the sequence and
    the sequence itself

    :param hmm_models: file in .hmm format
    :return:
    sequence_with_n: list of strings with the sigma factor sequence with N's
    on the masked area
    short_sequence: list of strings with the sigma factor sequence with
    -(total amount of N)- to show the size and location of the masked area
    length: list of integers that are the total length of the sequence
    """

    data_boolean = False
    length = []
    hmm = []
    sequence_with_n = []
    short_sequence = []
    # extract the length of the sequence and the masked or unmasked sequence
    with open(hmm_models, "r") as hmm_file:
        for lines in hmm_file:
            if lines.startswith('LENG'):
                leng = (lines.split())[1]
                length.append(leng)
            line = lines.split()
            for item in line:
                if item == "COMPO":
                    data_boolean = True
                if data_boolean:
                    hmm.append(item)
                if item == "//":
                    data_boolean = False
    # remove the non letters from the sequence. remove COMPO from the
    # sequence and make all the letters capitals
    hmm_sequence = list([letter for letter in hmm if letter.isalpha()])
    hmm_sequence = [index.strip() for index in
                    ' '.join(hmm_sequence).split("COMPO")]
    del hmm_sequence[0]
    hmm_sequence = [character.upper() for character in hmm_sequence]
    # remove all the M characters and remove the spaces between characters.
    # count the amount of N in the sequence.
    for character in hmm_sequence:
        remove_m = character.replace('M', '')
        sequence_with_n.append(remove_m.replace(" ", ""))
    # remove all the N except for the middle 3 N
    for seq in sequence_with_n:
        cleaner_sequence = []
        # split the sequence before and after a group or a single N
        sequences = re.split("(N{1,})", seq)
        for characters in sequences:
            # if there are more then 2 N's in the sequence
            if characters.count("N") > 2:
                # count the amount of N in the sequence
                number_of_n = characters.count("N")
                # replace the first N with -
                characters = characters.replace("N", "-", 1)
                # replace the "-" with "- number of N's -"
                characters = characters.replace("-", "-" + str(number_of_n) + "-")
                # remove all the N's that are still in the sequence
                characters = characters.replace("N", "")
                # if there are more then 0 N's
            if characters.count("N") > 0:
                # replace N with -
                characters = characters.replace("N", "-")
            cleaner_sequence.append(characters)
        # join the split sequence to become one whole sequence again
        short_sequence.append("".join(cleaner_sequence))

    return sequence_with_n, short_sequence, length


def create_hmm_database_output(sequence_with_n, short_sequence, length, reg_name):
    """ create table from hmm output that can be used in the database. convert
    it to a list so it can be combined with multiple sigma factors into one
    output file when the tool is run with multiple input files.

    :param sequence_with_n: list of strings with the regulator sequence
    with N's on the masked area
    :param short_sequence: ist of strings with the sigma factor sequence with
    -(total amount of N)- to show the size and location of the masked area
    :param length: list of integers that are the total length of the sequence
    :param reg_name: string of the name of the regulator
    :return:
    table_list: list of the hmm table
    """
    table = pd.DataFrame(sequence_with_n, index=None, columns=['sequence'])
    column_short_sequence = pd.Series(short_sequence, dtype='object')
    table.insert(loc=0, column='short_sequence', value=column_short_sequence)
    table['masked'] = length
    table['sig_factor_name'] = reg_name
    table = table[['sig_factor_name', 'masked', 'sequence', 'short_sequence']]
    table = table.rename(columns={col: "" for col in table})
    table_list = table.values.tolist()

    return table_list


def run_hmmpress(hmm_full_file):
    """Presses a multiple pHMM file to a pHMM database

    :param hmm_full_file: file with multiple pHMMs
    :return: hmm_full_file, Pressed pHMM database file
    """
    cmd_hmmpres = f"hmmpress -f {hmm_full_file}"
    try:
        subprocess.check_output(cmd_hmmpres, shell=True, stderr=subprocess.STDOUT)
    except Exception as e:
        raise Exception(console.print(
            f"[bold red]{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} -Unable to run hmmpress: {str(e)}[/bold red]"))

    return hmm_full_file


def prep_hmm_detection(reg_fasta, reg_name, mask_mode, ic_thres, outdir):
    """Prepares the input files for the hmm detection step

    :param reg_fasta: fasta file with the regulator sequence
    :param reg_name: name of the regulator
    :param mask_mode: mode of masking
    :param ic_thres: information content threshold
    :param outdir: output directory

    :return: Pressed pHMM database file
      """
    align_sequences(reg_fasta, reg_name, outdir)  # to handle unequal input sequence lengths
    alignment = convert_to_stockholm(reg_name, outdir)
    gap_freq_list = calculate_gap_freq(alignment)
    gap_summary = gap_checker(gap_freq_list)
    # Find position with max information content and no gap characters
    max_ent_pos = sh_max_finder(alignment)
    # Create MSAs with variable spacer lengths
    sequences_list = msa_mod_length_editor(alignment, max_ent_pos, gap_summary)
    alignment_filenames = processed_msa_file_writer(sequences_list, reg_name, outdir)
    # Mask sequence spacers using alimask
    if mask_mode in ["spacer_masking", "positional_masking"]:
        masked_seq_filenames = spacer_masker(sequences_list, alignment_filenames,
                                             mask_mode, outdir, ic_thres)
    elif mask_mode == "no_masking":
        masked_seq_filenames = spacer_masker(sequences_list, alignment_filenames, mask_mode, reg_name, outdir)
    else:
        raise ValueError(
            "Invalid mode. Try -m spacer_masking(default), -m positional_masking, or -m no_masking")
    # Add an annotation line to the MSA files
    masked_seq_name_editor(masked_seq_filenames)
    msa_db_filename = msa_file_merger(masked_seq_filenames, reg_name, outdir)
    # Create pHMM models from the alignment file
    hmm_models = model_creator(msa_db_filename, reg_name, outdir)
    sequence_with_n, short_sequence, length = parse_hmm(hmm_models)
    # combine the parsed data into output that can be used for the database
    db_output = create_hmm_database_output(sequence_with_n, short_sequence, length, reg_name)
    run_hmmpress(hmm_models)

    return hmm_models
