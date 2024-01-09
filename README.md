# MiniMotif (beta version)

## Project Description

Minimotif is a tool that detects transcription factor binding sites in a given genome. 

Minimotif detects transcription factor binding sites (TFBS) in a given genome, by combining the power of Position Weight Matrices (PWMs) and profile Hidden Markov Models (pHMMs). If the binding site of interest is gapless, then a Position Weight Matrix (PWM) is created and the tool MOODS is used to find any occurrences of the motif within the genome. Alternatively, if the binding site contain gaps (i.e sigma factor binding sites with variable spacer length), then MiniMotif constructs profile Hidden Markov Models (pHMMs) and interrogates the genome with the nhmmscan flavor of HMMER. In addition, it allows scanning a genome with a premade set of TFBSs. 

## Requirements

The following instructions require the installation of the following on your machine (links to installation guidelines) : 
1. Git ( https://github.com/git-guides/install-git )
2. conda ( https://conda.io/projects/conda/en/latest/user-guide/install/index.html )



## Instalation

Download MiniMotif with the following command:

```
git clone https://github.com/HAugustijn/MiniMotif.git
```
Note: Requires installation of git.

Then install all the dependencies from the minimotif.yml file with:

```
cd MiniMotif
conda env create -f minimotif.yml minimotif
conda activate MiniMotif 
```
Note: Requires installation of conda.

Note 2: Remember to activate the MiniMotif environment everytime you use MiniMotif!

## Quick usage 

Generally, MiniMotif can be used with the following command:

```
python3 minimotif.py [optional arguments] -G [genome_file] -O [output_directory]
```
Please read the following paragraphs for further information.

### 1) Query a genome with precalculated PWMs

Minimotif requires a genome file in .gbk format and allows the automated search of a genome by a set of precalculated PWMs from transcription factors of Streptomyces coelicolor, using the following command:

```
python3 minimotif.py -pc -G [genome_file] -O [output_directory] 

```
Notes: 
1. The genome filename has to be formatted as: [organism]_genome.gbk. i.e. scoe_genome.gbk.
2. Specifying an output directory is mandatory.

### 2) Query a genome using custom binding site sequences

The user can specify a binding site file in a .fasta format, for a given transcription factor. Each sequence in the multi-fasta file corresponds to one binding site. The sequences are used to construct a binding site profile.

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
Then, using the following command the tool decides if the profile is gapped or ungapped, based on a Shannon Information content heuristic:

```
python3 minimotif.py -i test.fasta -G [genome_file] -O [output_directory]
```
If the user knows that the motif is gapped, ungapped or wants both the PWM and pHMM branches to be used, then the flag -am (--analysis-mode) allows it:

```
python3 minimotif.py -i test.fasta -am gapped -G [genome_file] -O [output_directory]
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

## References



# License


