# ancseq
#### Version 1.0.0

## Table of contents
- [What is ancseq?](#What-is-ancseq)
- [Installation](#Installation)
  + [Dependencies](#Dependencies)
  + [Installation using bioconda](#Installation-using-bioconda)
- [Usage](#Usage)
  + [Example 1 : Running ancseq for nucleotide sequence alignment](#Example-1--Running-ancseq-for-nucleotide-sequence-alignment)
  + [Example 2 : Running ancseq for amino acid sequence alignment](#Example-2--Running-ancseq-for-amino-acid-sequence-alignment)
  + [Example 3 : Running ancseq for codon sequence alignment](#Example-3--Running-ancseq-for-codon-sequence-alignment)
  + [Example 4 : Running ancseq with ```--fast``` option](#Example-4--Running-ancseq-with---fast-option)
- [Outputs](#Outputs)


## What is ancseq?
Ancestral sequence reconstruction is a powerfull technique to reconstruct ancestral states from a multiple sequence alignment. ancseq is a wrapper tool to do ancestral sequence reconstruction using [IQ-TREE](http://www.iqtree.org).

#### Citation
- Under preparation...


## Installation
### Dependencies
#### Softwares
- [IQ-TREE](http://www.iqtree.org)

#### Python (>=3.5) libraries
- [biopython](https://biopython.org)
### Installation using conda
You can install ancseq using [anaconda](https://www.anaconda.com).
```
git clone https://github.com/YuSugihara/ancseq.git
cd ancseq
conda env create -f ancseq.yml
conda activate ancseq
```

## Usage

```
$ ancseq -h
usage: ancseq -s <ALIGNED_FASTA> -m <MODE> -o <OUT_DIR> [-t <INT>] [-f]

ancseq version 1.0.0

options:
  -h, --help           show this help message and exit
  -s , --seq           Sequence alignment in FASTA format.
  -m , --mode          Sequence type. [NT/AA/CODON]
  -o , --out           Output directory. The given name must not exist.
  -t , --threads       Number of threads. [8]
  -b , --bootstrap     Replicate for bootstrap. [1000]
  --max-report         Maximum number of ambiguous sites to report at the same position. [5]
  --min-prob           Minimum probability of being reported as an ambiguous site. [0.05]
  --min-gap-prob       Minimum probability of replacing the ancestral state with a gap. [0.5]
  --fast               Use -fast option in IQ-TREE [FLASE]
  --model              Input substitution model for IQ-TREE. [None]
  --stop-pseudo-codon  Stop calculation of pvalues of pseudo-codon [FLASE]
  --asr-only           Skip building tree and reconstruct ancestral states only [FLASE]
  -v, --version        show program's version number and exit
```

+ [Example 1 : Running ancseq for nucleotide sequence alignment](#Example-1--Running-ancseq-for-nucleotide-sequence-alignment)
+ [Example 2 : Running ancseq for amino acid sequence alignment](#Example-2--Running-ancseq-for-amino-acid-sequence-alignment)
+ [Example 3 : Running ancseq for codon sequence alignment](#Example-3--Running-ancseq-for-codon-sequence-alignment)
+ [Example 4 : Running ancseq with ```--fast``` option](#Example-4--Running-ancseq-with---fast-option)


### Example 1 : Running ancseq for nucleotide sequence alignment
```
ancseq -s test_nuc.fasta \
       -m DNA \
       -o out_dir
```

`-s` :  Nucleotide sequence alignment in fasta format.

`-m` :  Sequence type

`-o` : Name of the output directory. The given name should not exist.

### Example 2 : Running ancseq for amino acid sequence alignment
```
ancseq -s test_nuc.fasta \
       -m AA \
       -o out_dir
```

`-s` :  Amino acid sequence alignment in fasta format.

`-m` :  Sequence type

`-o` : Name of the output directory. The given name should not exist.

### Example 3 : Running ancseq for codon sequence alignment
**The codon mode can take a very long time!**
```
ancseq -s test_nuc.fasta \
       -m CODON \
       -o out_dir
```

`-s` :  Codon sequence alignment in fasta format.

`-m` :  Sequence type

`-o` : Name of the output directory. The given name should not exist.

### Example 4 : Running ancseq with ```--fast``` option
```
ancseq -s test_nuc.fasta \
       -m DNA \
       -o out_dir \
       --fast
```

`-s` :  Nucleotide sequence alignment in fasta format.

`-m` :  Sequence type

`-o` : Name of the output directory. The given name should not exist.

`--fast` : use ```-fast``` option in IQ-TREE.


## Outputs
Inside of `OUT_DIR` is like below.
```
├── 00_tree
│  ├── 00_iqtree.err
│  ├── 00_iqtree.out
│  ├── test_nuc.fasta
│  ├── test_nuc.fasta.bionj
│  ├── test_nuc.fasta.ckp.gz
│  ├── test_nuc.fasta.contree
│  ├── test_nuc.fasta.iqtree
│  ├── test_nuc.fasta.log
│  ├── test_nuc.fasta.mldist
│  ├── test_nuc.fasta.model.gz
│  ├── test_nuc.fasta.splits.nex
│  └── test_nuc.fasta.treefile
├── 10_asr
│  ├── 10_iqtree.err
│  ├── 10_iqtree.out
│  ├── test_nuc.fasta
│  ├── test_nuc.fasta.ckp.gz
│  ├── test_nuc.fasta.iqtree
│  ├── test_nuc.fasta.log
│  ├── test_nuc.fasta.state.gz
│  └── test_nuc.fasta.treefile
├── 20_indels
│  ├── 20_iqtree.err
│  ├── 20_iqtree.out
│  ├── test_nuc.fasta.binary
│  ├── test_nuc.fasta.binary.ckp.gz
│  ├── test_nuc.fasta.binary.iqtree
│  ├── test_nuc.fasta.binary.log
│  ├── test_nuc.fasta.binary.state.gz
│  └── test_nuc.fasta.binary.treefile
└── 30_result
   ├── ancestral_state_result.treefile
   ├── ancestral_state_result.fasta
   ├── ancestral_state_result_with_gap.fasta
   ├── ancestral_state_result.sort.tsv
   ├── ancestral_state_result.tsv.gz
   └── ancestral_state_result.pseudo_codon.tsv.gz
```
- The phylogenetic tree reconstructed by IQ-TREE can be found in `00_tree`.
- The results of the ancestral sequence reconstruction can be found in `30_result`.
  + `ancestral_state_result.treefile`: Phylogenetic tree with the node labels.
  + `ancestral_state_result.fasta`: FASTA file of the ancestral sequences without gaps.
  + `ancestral_state_result_with_gap.fasta`: FASTA file of the ancestral sequences with gaps.
  + `ancestral_state_result.sort.tsv` : Probabilities of the ancestral states. 
