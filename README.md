# MiniMotif (beta version)

## Introduction

Minimotif detects transcription factor binding sites in a given genome, by combining the power of Position Weight Matrices (PWMs) and profile Hidden Markov Models (pHMMs). If the binding site of interest is gapless, then a PWM is created and the tool MOODS is used to find any occurences of the motif within the genome. Alternatively, if the binding sites contain gaps (i.e sigma factor binding sites wih variable spacer length), then MiniMotif constructs pHMMs and interrogates the genome with the nhmmscan flavor of HMMer.

## Instalation

MiniMotif downloaded with the following command:

```
git clone https://github.com/HAugustijn/MiniMotif.git
```

Then install all the dependencies from the minimotif.yml file with:

```
conda env create -f minimotif.yml minimotif
conda activate minimotif
```

## Quick usage 

In general, MiniMotif can be used with the following command:

```
minimotif [optional arguments] -G [genome_file] -O [output_directory]
```
Please read the following paragraphs for further information.

### 1) Query a genome with precalculated PWMs

Minimotif requires a genome file in .gbk format and allows the automated search of a enome by a set of precalculated PWMs from transcription factors of Streptomyces coelicolor, using the following command:

```
minimotf -G [genome_file] -O [output_directory] -pc

```
Notes: 
1. The genome filename has to be formatted as: [organism]_genome.gbk. i.e. scoe_genome.gbk.
2. Specifying an output directory is mandatory.

### 2) Query a genome using custom binding site sequences

The user can specify a binding site file in a .fasta format, for a given transcription factor. Each sequence in the multi-fasta file corresponds to one binding site. The sequences are used to construct a binding site profile/

Example: test.fasta
```
>1
ACTGGTCTAGACAACT
>2
ACTGGTCTAGACAAGA
>3
ACTGGTCTACACCAGT
>4
ACAGGTCTACACCACT
>5
AGTGGTGTAGACCACC
>6
ATTGGTCTAAACCACA
```
Then, using the following command the tool decides if the profile is gapped or ungapped, based on an Information content  heuristic:

```
minimotif -i test.fasta -G [genome_file] -O [output_directory]
```
If the user knows that the motif is gapped, ungapped or wants both the PWM and pHMM branches to be used, then the flag -am (--analysis-mode) allows it:

```
minimotif -i test.fasta -am gapped -G [genome_file] -O [output_directory]
```
Notes: -am can be set to "ungapped" (PWMs), "gapped" (pHMMs), "both" (PWMs and pHMMs) and "auto"( Default)

Here's a full description of all the optional arguments:

```
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
```
# License


