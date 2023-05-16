""" Script containing all functions for TFBS detection using HMMs """
import os
import subprocess
from Bio import SeqIO


def nhmmscan_wrapper(fasta_file, reg_name, gb_name, reg_type, hmm_db_file, outdir):
    """Runs nhmmscan against multiple DNA query sequences

    :param outdir: path to the output directory
    :param reg_type: string, regulatory type
    :param gb_name: string, name of the gb file
    :param reg_name: string, name of the regulator
    :param fasta_file: DNA query fasta file
    :param hmm_db_file:Pressed pHMM model database file
    :return:tab_out_list: nhmmscan output in tabular format
    """
    tabular_out = f"{outdir}/nhmmscan_{reg_name}_{gb_name}_{reg_type}_tab_out.txt"
    full_out = f"{outdir}/nhmmscan_{reg_name}_{gb_name}_{reg_type}_full_out.txt"
    run_nhmmscan(hmm_db_file, fasta_file, tabular_out, full_out)
    return tabular_out


def run_nhmmscan(hmm_full_file, query_seq_db, tabular_outfile,
                 regular_outfile):
    """Runs  HMMER flavor 'nhmmscan', search pHMMs db against sequence queries

    :param hmm_full_file: Pressed pHMM database file
    :param query_seq_db: Sequences to search for matches with the pHMMs
    :param tabular_outfile: Output in tabular format file
    :param regular_outfile: Output in regular HMMEr format file
    :return: tabular_outfile: Output in tabular format file
    """
    cmd_nhmmscan = f"nhmmscan --tblout {tabular_outfile} --max -T 0.1"  \
                   f" {hmm_full_file} {query_seq_db}> {regular_outfile}"
    try:
        if not os.path.exists(tabular_outfile):
            subprocess.check_output(cmd_nhmmscan, shell=True, stderr=subprocess.STDOUT)

    except subprocess.CalledProcessError:
        print('Unable to run nhmmscan')

    return tabular_outfile


def motif_length_parser(hmm_full_file):
    """Parses hmm_db file and extracts motif length and name per record.

    :param hmm_full_file: hmm_db file
    :return:motif_length_list: A list of tuples,(model_name,model_length)
    """
    motif_length_list = []
    with open(hmm_full_file) as hmm_file:
        for line in hmm_file:
            if line.startswith("NAME"):
                model_name = line.strip().split()[1]
            if line.startswith("LENG"):
                motif_length = line.strip().split()[1]
                mod_name_length_tup = model_name, int(motif_length)
                motif_length_list.append(mod_name_length_tup)

    return motif_length_list


def nhmmscan_out_file_writer(outfile, metadata, real_hits):
    """Formats processed nhmmscan output and writes it to a file

    :param outfile: Name of the output file
    :param metadata: List of lines, describing nhmmscan metadata
    :param real_hits: List of lists, each list a line <--str
    :return: None
    """
    initial_metadata = metadata[:2]
    final_metadata = metadata[2:]

    with open(outfile, 'w') as out_file:

        out_file.write('Processed nhmmscan output file '
                       'to include only hits properly aligned '
                       'to the HMM profiles'
                       ' \n')
        for line in initial_metadata:
            out_file.write(line + '\n')
        for line in real_hits:
            line = '\t  '.join([col for col in line])
            out_file.write(line + '\n')
        for line in final_metadata:
            out_file.write(line)
    return None


def nhmmscan_out_parser(nhmmscan_tab, motif_length_list, adjust_length,
                        outdir, reg_type, gb_name, reg_name):
    """Parses hits that are adequately aligned with the model from nhmmscan out

    :param nhmmscan_tab: nhmmscan output tabular format
    :param motif_length_list:List of tuples, (model_name,model_length)
    :param adjust_length: int, defines the allowed difference between
    the length of the model and the alignment between a hit and the model.
    If 0 only full alignments are allowed. If 1, alignments that are
    model_length - 1 are also included to the output.
    :param outdir: path to the output directory
    :param reg_type: string, regulatory type
    :param gb_name: string, name of the gb file
    :param reg_name: string, name of the regulator
    :return:proc_output_filenames, parsed nhmmscan output filenames
    full_length_real_hits: List of hits that are fully aligned with the
    model
    partial_len_real_hits: List of hits that are partially aligned with
    the model
    """
    metadata = []
    real_hits = []
    track_dup = {}
    full_len_real_hits = []
    partial_len_real_hits = []

    outfile = f"{outdir}/nhmmscan_{reg_name}_{gb_name}_{reg_type}_region_tab_out.txt"

    with open(nhmmscan_tab) as inp_file:
        for line in inp_file:
            if line.startswith('#'):
                metadata.append(line)
            else:
                line = line.strip().split()
                for mot_name, mot_length in motif_length_list:
                    # If pHMM alignment starts from 1 and ends to
                    # motif length - adjust_length, put the entry to
                    # an intermediate dict
                    if int(line[4]) <= 3 and line[0] == mot_name \
                            and int(line[5]) >= mot_length - adjust_length:

                        if line[2] not in track_dup.keys():
                            track_dup[line[2]] = [line]
                        else:
                            track_dup[line[2]].append(line)

    for value in track_dup.values():
        value = sorted(value, key=lambda x: (float(x[6]), -float(x[13])))
        # If there are multiple hits for the same regulatory region, sort
        # them based on the position that they start highest nhmmscan score
        # (- sign reverses the sorting, so higher nhmmscan values are
        # outputted first.)
        if len(value) > 1:
            for i, j in enumerate(value[:-1]):
                # Pairwise comparison between the hits of the same region
                # if multiple models are detecting the same binding site,
                # modify the annotation of all the duplicates.
                if j[6] == value[i + 1][6] or j[7] == value[i + 1][7]:
                    value[i + 1][2] = f"Duplicate_{value[i + 1][2]}"

    for value in track_dup.values():
        for hit in value:
            # Filter out duplicates and put the final hits to a list
            if hit[2].startswith("Duplicate"):
                continue
            else:
                real_hits.append(hit)

    for real_hit in real_hits:
        for mot_name, mot_length in motif_length_list:
            # Distinguish between full hits and partial hits
            if int(real_hit[5]) == mot_length and real_hit[0] == mot_name and int(real_hit[4]) == 1:
                full_len_real_hits.append(
                    f"{real_hit[2]}_{real_hit[6]}")
            else:
                partial_len_real_hits.append(
                    f"{real_hit[2]}_{real_hit[6]}")

    nhmmscan_out_file_writer(outfile, metadata, real_hits)
    real_hits = []
    metadata = []
    track_dup = {}

    return outfile, full_len_real_hits, partial_len_real_hits


def gene_strand_prod_parser(gb_file):
    """ Parses gene strand and product from a .gb genome annotation file

    :param gb_file: a .gb Genbank genome annotation file
    :return:gene_strand_dict: dict with gene names or locus tags as keys and
    '+' or '-' as values, depending on the gene's orientation.
    product_dict: A dict with gene names or locus tags as keys and predicted
    products as values.
    *** If a CDS does not have a product annotation, an N/A value is put in as
    a qualifier in the annotation. Those genes are also stored in a list so
    in a latter version they can be shown to the user for interpretation
    reasons.
    """
    product_dict = {}
    gene_strand_dict = {}
    no_product_genes_list = []
    for rec in SeqIO.parse(gb_file, "genbank"):
        for f in rec.features:
            if f.type == "CDS":
                if "locus_tag" in f.qualifiers:
                    gene_name = f.qualifiers["locus_tag"]
                    gene_name = ''.join(gene_name)
                    gene_name = gene_name.replace('-', '.')
                    gene_name = gene_name.split()
                else:
                    gene_name = f.qualifiers["gene"]
                    gene_name = ''.join(gene_name)
                    gene_name = gene_name.replace('-', '.')
                    gene_name = gene_name.split()
                if "product" not in f.qualifiers.keys():
                    f.qualifiers["product"] = "N/A"
                    no_product_genes_list.append(gene_name[0])
                gene_strand = f.strand
                product_dict[gene_name[0]] = f.qualifiers["product"]
                if gene_strand == 1:
                    gene_strand_dict[gene_name[0]] = "+"
                else:
                    gene_strand_dict[gene_name[0]] = "-"

    return gene_strand_dict, product_dict


def genes_with_bs_file_writer(output_filename, sfbs_list, gene_strand_dict, product_dict):
    """Writes a summary  file of gene targets for the Sigma factor

    :param output_filename: Name of the output file
    :param sfbs_list: a list of binding site summary. In detail, the format:
    bs_gene,bs genomic location ,bs strand, bs sequence, nhmmscan score,
    nhmmscan e-value, length of bs sequence, alignment type, alignment_length,
    model length.
    In cases that there are two genes in the form of Gene1-Gene2, Then the
    first element bs_gene is replaced by 2 flanking genes, Gene1, Gene2 and
    the length of the list is 11 in contrast with the first case that it is
    10
    :param gene_strand_dict: A dictionary with gene IDs are keys and '+' or '-'
    as values.
    :param product_dict: A dict with gene IDs as keys and the annotated gene
    product as value
    :return: None
    *** This function also outputs if a gene has the right orientation to be
    regulated by the binding site. If it is not, then an N/A is appended to the
    regulated gene list. otherwise, the gene name is appended to the list and
    the product is also written to the output file.
    """
    regulated_gene_list = []

    with open(output_filename, 'w') as out_fn:
        out_fn.write(
            f'#Target Genes\tBinding_site_genomic_coordinates\tBinding_site_sequence\tReported_bs_length'
            f'\tAlignment_length\tModel_length\tAl_type\tnhmmscan_eval\tnhmmscan_score\tStrand'
            f'\tGene_strand\tRegulated_genes\tProduct\n')

        for bs_summary in sfbs_list:
            if len(bs_summary) == 11:
                left_gene_strand = gene_strand_dict[bs_summary[0]]
                right_gene_strand = gene_strand_dict[bs_summary[1]]
                left_gene_product = product_dict[bs_summary[0]]
                right_gene_product = product_dict[bs_summary[1]]
                flanking_genes_strand = \
                    f"{left_gene_strand},{right_gene_strand}"
                if bs_summary[3] == "-" and left_gene_strand == "-":
                    regulated_gene_list.append(bs_summary[0])
                    product_string = f"{bs_summary[0]}:{left_gene_product[0]}"
                elif bs_summary[3] == "+" and right_gene_strand == "+":
                    product_string = f"{bs_summary[1]}:{right_gene_product[0]}"
                    regulated_gene_list.append(bs_summary[1])
                else:
                    regulated_gene_list.append("N/A")
                    product_string = 'N/A'
                genes_string = f"{bs_summary[0]}-{bs_summary[1]}"
                gene_strand_string = f"{genes_string}\t{bs_summary[2]}\t{bs_summary[4]}\t{bs_summary[7]}" \
                                     f"\t{bs_summary[9]}\t{bs_summary[10]}\t{bs_summary[8]}\t{bs_summary[6]}" \
                                     f"\t{bs_summary[5]}\t{bs_summary[3]}\t{flanking_genes_strand}" \
                                     f"\t{','.join(regulated_gene_list)}\t{product_string}"

            else:
                gene1 = bs_summary[0].split("-")[0]
                gene_orient = gene_strand_dict[gene1]
                gene_product = product_dict[gene1]

                if bs_summary[2] == "-" and gene_orient == "-":
                    regulated_gene_list.append(bs_summary[0])
                    gene_product_string = f"{gene1}:{gene_product[0]}"
                elif bs_summary[2] == "+" and gene_orient == "+":
                    regulated_gene_list.append(bs_summary[0])
                    gene_product_string = f"{gene1}:{gene_product[0]}"
                else:
                    regulated_gene_list.append(["N/A"])
                    gene_product_string = "N/A"
                gene_strand_string = f"{gene1}\t{bs_summary[1]}\t{bs_summary[3]}\t{bs_summary[6]}" \
                                     f"\t{bs_summary[8]}\t{bs_summary[9]}\t{bs_summary[7]}\t{bs_summary[5]}" \
                                     f"\t{bs_summary[4]}\t{bs_summary[2]}\t{gene_orient}\t{regulated_gene_list[0][0]}" \
                                     f"\t{gene_product_string}"
            regulated_gene_list = []
            out_fn.write(f"{gene_strand_string}\n")
    return


def binding_site_seq_parser(fasta_file, gene_name, start_pos, end_pos, strand, al_type, model_length):
    """ Parses the sequence, genomic location & length from given fasta file

    :param fasta_file:str, File that contains genomic regions in fasta format.
    :param gene_name:str, Name of the gene that is under investigation
    :param start_pos:str, position that the genomic region starts
    :param end_pos:str, position that the genomic region ends
    :param strand:str, '+' or '-"
    :param al_type: str, full or partial depending on the input length
    :param model_length: int, length of the model that is aligned with the
    sequence of interest.
    If 0 only full alignments are allowed. If 1, alignments with length
    model_length - 1 are also included to the output.
    :return: binding_site_seq: Seq Object, sequence of the binding site
    gen_loc_bs: Genomic location of the binding site
    len(binding_site_seq): Length of the reported sequence
    """
    record_iterator = SeqIO.parse(fasta_file, "fasta")
    binding_site_seq = ""
    gen_loc_bs = ""

    for seq_record in record_iterator:
        seq_rec_list = seq_record.id.strip().split('~')
        seq_rec_name = seq_rec_list[0]

        if seq_rec_name == gene_name:
            seq_rec_start_loc, seq_rec_end_loc = seq_rec_list[1].split(":")
            binding_site_loc = sorted([int(start_pos), int(end_pos)])

            if al_type == "Full":
                gen_loc_bs = f"{int(seq_rec_start_loc) + binding_site_loc[0]}"\
                             f"-{int(seq_rec_start_loc) + binding_site_loc[1]}"
                if strand == "+":
                    binding_site_seq = seq_record.seq[
                             binding_site_loc[0] - 1: binding_site_loc[1]]
                elif strand == "-":
                    binding_site_seq = seq_record.seq[
                                   binding_site_loc[0] - 1:
                                   binding_site_loc[1]].reverse_complement()
                break
            elif al_type == "Partial":
                gen_loc_bs = f"{int(seq_rec_start_loc) + binding_site_loc[0]}"\
               f"-{int(seq_rec_start_loc)+ binding_site_loc[0] + model_length}"
                if strand == "+":
                    binding_site_seq = seq_record.seq[
                                   binding_site_loc[0] - 1: binding_site_loc[
                                                    0] - 1 + model_length]
                elif strand == "-":
                    binding_site_seq = seq_record.seq[
                                   binding_site_loc[1] - model_length:
                                   binding_site_loc[1]].reverse_complement()
                break

    return binding_site_seq, gen_loc_bs, len(binding_site_seq)


def genes_with_bs_finder(processed_output_filenames, gene_orient_dict, product_dict, fasta_file,
                         full_len, partial_len, reg_name, gb_name, reg_type, outdir):
    """Wrapper:extracts Sigma factor target genes from nhmmscan processed
    output

    :param processed_output_filenames: nhmmscan processed output files
    :param gene_orient_dict: A dict with genes as keys and '+' or '-' as values
    :param product_dict:  A dict with genes as keys and the annotated gene
    gene product as value
    :param fasta_file: DNA query fasta files
    :param full_len: List of hits that are fully aligned with the model
    :param partial_len: List of hits that are partially aligned with the model
    :param outdir: path to the output directory
    :param reg_type: string, regulatory type
    :param gb_name: string, name of the gb file
    :param reg_name: string, name of the regulator
    """
    sfbs_list = []
    output_filename = f"{outdir}/{reg_name}_{gb_name}_{reg_type}_hmm_results.tsv"

    with open(processed_output_filenames) as p_outfile:
        for line in p_outfile:
            if line.startswith('#') or line.startswith('Processed nhmmscan'):
                continue
            elif len(line) == 1:
                continue  # Avoid irrelevant lines
            else:
                line = line.strip().split()
                al_length = line[5]
                binding_start = line[6]
                binding_end = line[7]
                model_length = line[10]
                strand = line[11]
                nhmmscan_eval = line[12]
                nhmmscan_score = line[13]

                if len(line) > 1:
                    # Parse gene name: Either Gene or Gene1-Gene2
                    bs_region = f"{line[2].split('~')[0]}"
                    if f"{line[2]}_{line[6]}" in full_len:
                        al_type = "Full"
                        bs_seq, gen_bs_loc, len_bs_seq = \
                            binding_site_seq_parser(fasta_file, bs_region, binding_start,
                                                    binding_end, strand, al_type, int(model_length))
                    elif f"{line[2]}_{line[6]}" in partial_len:
                        al_type = 'Partial'
                        bs_seq, gen_bs_loc, len_bs_seq = \
                            binding_site_seq_parser(fasta_file, bs_region, binding_start,
                                                    binding_end, strand, al_type, int(model_length))

                    sfbs_list.append((bs_region, gen_bs_loc, strand, bs_seq, nhmmscan_score,
                                      nhmmscan_eval, len_bs_seq, al_type, al_length, model_length))

        genes_with_bs_file_writer(output_filename, sfbs_list, gene_orient_dict, product_dict)
        sfbs_list = []
    return


def run_hmm_detection(gbk_file, reg_name, hmm_models, coding, adj_len, outdir):
    """ Wrapper function for running the HMM detection

    :param hmm_models: hmm models created by the prep_hmm_detection.py script
    :param coding: True/False, detect for coding regions
    :param adj_len: float, adjust the length of the alignments
    :param outdir: path to the output directory
    :param gbk_file: string, name of the gb file
    :param reg_name: string, name of the regulator

     """
    gb_name = (gbk_file).split("/")[-1].split(".")[0]
    reg_fasta = f"{outdir}/{gb_name}_reg_region.fasta"

    tab_out = nhmmscan_wrapper(reg_fasta, reg_name, gb_name, "reg", hmm_models, outdir)
    motif_sum_list = motif_length_parser(hmm_models)
    nhmmscan_file, full_len, partial_len = nhmmscan_out_parser(tab_out, motif_sum_list,
                                                               adj_len, outdir, "reg", gb_name,reg_name)

    gene_orient_dict, product_dict = gene_strand_prod_parser(gbk_file)
    genes_with_bs_finder(nhmmscan_file, gene_orient_dict, product_dict, reg_fasta,
         full_len, partial_len, reg_name, gb_name, "reg", outdir)

    if coding:
        co_fasta = f"{outdir}/{gb_name}_co_region.fasta"

        tab_out = nhmmscan_wrapper(co_fasta, reg_name, gb_name, "co", hmm_models, outdir)
        motif_sum_list = motif_length_parser(hmm_models)
        nhmmscan_file, full_len, partial_len = nhmmscan_out_parser(tab_out, motif_sum_list,
                                                                   adj_len, outdir, "co", gb_name,reg_name)

        gene_orient_dict, product_dict = gene_strand_prod_parser(gbk_file)
        genes_with_bs_finder(nhmmscan_file, gene_orient_dict, product_dict, co_fasta,
             full_len, partial_len, reg_name, gb_name, "co", outdir)
