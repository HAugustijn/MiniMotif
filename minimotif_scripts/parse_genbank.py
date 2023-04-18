"""  """

from Bio import SeqIO


def is_gbk(filename):
    """ Determines if the input file is in genbank format """
    with open(filename, "r") as handle:
        gbk = SeqIO.parse(handle, "genbank")
        if not any(gbk):
            print(f'Please provide the genome input (-G) in GenBank format for file: {filename}')
            return False
        else:
            return True


def extract_regions(genbank_file, cotrans_region, reg_region):
    """ Extracts coding and regulatory regions for genbank files """
    regulatory_region = []  # region in which TFBSs mostly occur
    coding_region = []  # CDS sequences of genes
    CDS_region = []  # used to store intermediate information
    counter = 0

    with open(genbank_file, "r") as gb_file:
        # Parse genbank file and extract CDS regions of genome of interest
        for rec in SeqIO.parse(gb_file, "genbank"):
            genome_length = len(rec)
            complete_seq = rec.seq
            for f in rec.features:
                if f.type == "CDS":
                    if "locus_tag" in f.qualifiers:
                        gene_name = f.qualifiers["locus_tag"]
                        if type(gene_name) == list:
                            gene_name = gene_name[0]
                    elif "gene" in f.qualifiers:
                        gene_name = f.qualifiers["gene"]
                        if type(gene_name) == list:
                            gene_name = gene_name[0]
                    else:
                        print(f"Couldn't identify the CDS locus tag in the following qualifiers: {f.qualifiers}")
                    if str(f.location).startswith("join"):
                        start = min(f.location)
                        end = max(f.location) + 1
                    else:
                        start = f.location._start.position
                        end = f.location._end.position
                    CDS_region.append([gene_name, start, end, f.strand])
                    intergene_seq = rec.seq[start:end]
                    seq_region = f"{start}:{end}"
                    counter += 1
                    coding_region.append([gene_name, seq_region, f.strand, intergene_seq])
                    if counter == 1 and start >= cotrans_region[1]:
                        intergene_seq = rec.seq[start + reg_region[0]:start + reg_region[1]]
                        seq_region = f"{start + reg_region[0]}:{start + reg_region[1]}"
                        regulatory_region.append([gene_name, seq_region, f.strand, intergene_seq])
            # Based on the location of the CDS genes, identify regulatory regions
            for i, pos in enumerate(CDS_region[1:]):
                # Compare current start position to previous end position
                last_end = CDS_region[i][2]
                this_start = pos[1]
                gene_region_name = '-'.join([CDS_region[i][0], pos[0]])
                strand_pos = [CDS_region[i][3], pos[3]]
                if this_start - last_end >= cotrans_region[1] and strand_pos != [1, -1]:
                    # Remove regions which are likely co-transcribed and exclude terminator regions
                    intergene_seq = rec.seq[(this_start) + reg_region[0]:(this_start) + reg_region[1]]
                    seq_region = f"{this_start + reg_region[0]}:{this_start + reg_region[1]}"
                    regulatory_region.append([gene_region_name, seq_region, strand_pos, intergene_seq])
                elif 1 < this_start - last_end < cotrans_region[1] and strand_pos == [-1, 1]:
                    # Include shorter regions between two start codons
                    intergene_seq = rec.seq[(last_end):(this_start)]
                    seq_region = f"{last_end}:{this_start}"
                    regulatory_region.append([gene_region_name, seq_region, strand_pos, intergene_seq])
            if CDS_region[-1][3] == 1:  # Add the region after the last CDS
                intergene_seq = rec.seq[(CDS_region[-1][2]) + reg_region[0]:genome_length]
                seq_region = f"{CDS_region[-1][2] + reg_region[0]}:{genome_length}"
                regulatory_region.append([gene_region_name, seq_region, strand_pos, intergene_seq])
    return regulatory_region, coding_region, complete_seq


def write_fastas(region, region_type, genome_name, outdir):
    """ Writes fasta files of the four previously identified regions """
    file_to_write = f"{outdir}/{genome_name}_{region_type}_region.fasta"
    with open(file_to_write, "w") as outfile:
        for i, seq in enumerate(region):
            fasta_header = f">{region[i][0]}~{region[i][1]}~{genome_name}"
            sequence = region[i][3]
            outfile.write(f"{fasta_header}\n{sequence}\n")
    return
