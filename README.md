# ancseq
#### Version 1.1.0

## Table of contents
- [What is ancseq?](#what-is-ancseq)
- [Installation](#installation)
  + [Dependencies](#dependencies)
  + [Installation using conda](#installation-using-conda)
- [Usage](#usage)
  + [Example 1 : Running ancseq for nucleotide sequence alignment](#example-1--running-ancseq-for-nucleotide-sequence-alignment)
  + [Example 2 : Running ancseq for amino acid sequence alignment](#example-2--running-ancseq-for-amino-acid-sequence-alignment)
  + [Example 3 : Running ancseq for codon sequence alignment](#example-3--running-ancseq-for-codon-sequence-alignment)
  + [Example 4 : Running ancseq with ```--fast``` option](#example-4--Running-ancseq-with---fast-option)
- [Outputs](#outputs)
- [Workflow of ancseq](#workflow-of-ancseq)
  + [IQ-TREE command 1 : Building a phylogenetic tree](#iq-tree-command-1)
  + [IQ-TREE command 2 : Reconstructing ancestral sequences](#iq-tree-command-2)
  + [IQ-TREE command 3 : Reconstructing inserstions and deletions (INDELs)](#iq-tree-command-3)
- [What is a pseudo-codon probability in ancseq?](#what-is-a-pseudo-codon-probability)


## What is ancseq?
<img src="https://github.com/YuSugihara/ancseq/blob/main/images/ancseq_workflow.png" width=400>

Ancestral sequence reconstruction is a technique to reconstruct ancestral states from a multiple sequence alignment. ancseq is a wrapper tool to reconstruct ancestral sequences using [IQ-TREE](http://www.iqtree.org). See more detail workflow of ancseq [here](#Workflow-of-ancseq).

#### Citation
- Under preparation...


## Installation
### Dependencies
#### Softwares
- [IQ-TREE](http://www.iqtree.org)

#### Python (>=3.5) libraries
- [biopython](https://biopython.org)

### Installation using conda
You can install ancseq with the dependencies using [anaconda](https://www.anaconda.com).
```bash
git clone https://github.com/YuSugihara/ancseq.git
cd ancseq
conda env create -f ancseq.yml
conda activate ancseq
pip install .
```

## Usage

```bash
$ ancseq -h
usage: ancseq -s <ALIGNED_FASTA> -m <MODE> -o <OUT_DIR> [-t <INT>]

ancseq version 1.1.0

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
  --model              Specify substitution model for IQ-TREE. IQ-TREE searches the best substitution
                       model using ModelFinder in default [MFP]
  --outgroup           Specify outgroup for IQ-TREE. [None]
  --stop-pseudo-codon  Stop calculation of pvalues of pseudo-codon [FLASE]
  --asr-only           Skip building tree and reconstruct ancestral states only [FLASE]
  -v, --version        show program's version number and exit
```

**We recoomend to specify the outgroup to avoid misinterpretation of the ancestral states of nodes. See more detail [here](#Example-4--Running-ancseq-specifing-outgroup).**

+ [Example 1 : Running ancseq for nucleotide sequence alignment](#example-1--running-ancseq-for-nucleotide-sequence-alignment)
+ [Example 2 : Running ancseq for amino acid sequence alignment](#example-2--running-ancseq-for-amino-acid-sequence-alignment)
+ [Example 3 : Running ancseq for codon sequence alignment](#example-3--running-ancseq-for-codon-sequence-alignment)
+ [Example 4 : Running ancseq specifing outgroup](#example-4--running-ancseq-specifing-outgroup)
+ [Example 5 : Running ancseq with ```--fast``` option](#example-5--running-ancseq-with---fast-option)


### Example 1 : Running ancseq for nucleotide sequence alignment
```
ancseq -s test_nuc.fasta \
       -m DNA \
       -o out_dir
```

`-s` : Nucleotide sequence alignment in fasta format.

`-m` : Sequence type.

`-o` : Name of the output directory. The given name should not exist.

### Example 2 : Running ancseq for amino acid sequence alignment

```bash
ancseq -s test_nuc.fasta \
       -m AA \
       -o out_dir
```

`-s` : Amino acid sequence alignment in fasta format.

`-m` : Sequence type.

`-o` : Name of the output directory. The given name should not exist.

### Example 3 : Running ancseq for codon sequence alignment

**!!!WARNING!!!** The codon mode can take a very long time to build a phylogenetic tree. Therefore, we would recommend to run ancseq in DNA mode even if your alignmment is codon-based.

```bash
ancseq -s test_codon.fasta \
       -m CODON \
       -o out_dir
```

`-s` : Codon sequence alignment in fasta format.

`-m` : Sequence type.

`-o` : Name of the output directory. The given name should not exist.

### Example 4 : Running ancseq specifing outgroup

You can reconstruct the ancestral states without specifying the outgroup. However, the ancestral states of the node may be misinterpreted when you visualize the tree. Therefore, we recommend to specify the outgroup to avoid misinterpretation of ancestral states of the node. IQ-TREE converts the rooted tree to the unrooted tree in defalt.

```
ancseq -s test_nuc.fasta \
       -m DNA \
       --outgroup seq_id \
       -o out_dir
```

`-s` : Nucleotide sequence alignment in fasta format.

`-m` : Sequence type.

`--outgroup` : Sequence ID of outgroup.

`-o` : Name of the output directory. The given name should not exist.

### Example 5 : Running ancseq with ```--fast``` option

```bash
ancseq -s test_nuc.fasta \
       -m DNA \
       -o out_dir \
       --fast
```

`-s` : Nucleotide sequence alignment in fasta format.

`-m` : Sequence type.

`-o` : Name of the output directory. The given name should not exist.

`--fast` : Use ```-fast``` option in IQ-TREE.


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

## Workflow of ancseq

<img src="https://github.com/YuSugihara/ancseq/blob/main/images/ancseq_workflow.png" width=400>

- [IQ-TREE command 1](#iq-tree-command-1) : Building a phylogenetic tree.
- [IQ-TREE command 2](#iq-tree-command-2) : Reconstructing ancestral sequences.
- [IQ-TREE command 3](#iq-tree-command-3) : Reconstructing inserstions and deletions (INDELs).

### IQ-TREE command 1

```bash
iqtree -s ${INPUT_FASTA} \
       -st ${SEQ_TYPE} \
       -T ${NUM_THREADS} \
       -B ${NUM_BOOTSTRAP} \
       -m MFP \
       1> /OUT_DIR/00_tree/00_iqtree.out \
       2> /OUT_DIR/00_tree/00_iqtree.err
```

`-s` : Sequence alignment in fasta format.

`-st` : Sequence type.

`-T` : Number of threads.

`-B` : Replicates for ultrafast bootstrap.

`-m MFP` : Extended model selection followed by tree inference.

If you specify the ```--fast``` option in ancseq, the IQ-TREE command 1 will change as follows.

```bash
iqtree -s ${INPUT_FASTA} \
       -st ${SEQ_TYPE} \
       -T ${NUM_THREADS} \
       --alrt ${NUM_BOOTSTRAP} \
       -m MFP \
       --fast \
       1> /OUT_DIR/00_tree/00_iqtree.out \
       2> /OUT_DIR/00_tree/00_iqtree.err
```

`-s` : Sequence alignment in fasta format.

`-st` : Sequence type.

`-T` : Number of threads.

`--alrt` : Replicates for SH approximate likelihood ratio test.

`-m MFP` : Extended model selection followed by tree inference.

`--fast` : Fast search to resemble FastTree.

### IQ-TREE command 2

```bash
iqtree -asr \
       -s ${INPUT_FASTA} \
       -te /OUT_DIR/00_tree/${INPUT_FASTA}.treefile \
       -st ${SEQ_TYPE} \
       -T ${NUM_THREADS} \
       -m ${MODEL} \
       -o ${OUTGROUP} \
       -keep_empty_seq \
       1> /OUT_DIR/10_asr/10_iqtree.out \
       2> /OUT_DIR/10_asr/10_iqtree.err
```

`-asr` : Ancestral state reconstruction by empirical Bayes.

`-s` : Sequence alignment in fasta format.

`-te` : Tree file.

`-st` : Sequence type.

`-T` : Number of threads.

`-m` : Model name.

`-o` : Sequence ID of the outgroup if you specify the outgroup with `--outgroup`.

`-keep_empty_seq` : Keep empty sequences in the alignment.


### IQ-TREE command 3

```bash
iqtree -asr \
       -s ${INPUT_FASTA}.binary \
       -te /OUT_DIR/00_tree/${INPUT_FASTA}.treefile \
       -st BIN \
       -T ${NUM_THREADS} \
       -blfix \
       -m JC2 \
       -o ${OUTGROUP} \
       -keep_empty_seq \
       1> /OUT_DIR/20_indels/20_iqtree.out \
       2> /OUT_DIR/20_indels/20_iqtree.err
```

`-asr` : Ancestral state reconstruction by empirical Bayes.

`-s` : Sequence alignment in fasta format.

`-te` : Tree file.

`-st BIN` : Binary sequence type.

`-T` : Number of threads.

`-blfix` : Fix branch lengths of tree passed via `-t` or `-te`.

`-m JC2` : Jukes-Cantor type model for binary data.

`-o` : Sequence ID of the outgroup if you specify the outgroup with `--outgroup`.

`-keep_empty_seq` : Keep empty sequences in the alignment.

#### References

1. Aadland K, Pugh C, Kolaczkowski B. 2019. High-Throughput Reconstruction of Ancestral Protein Sequence, Structure, and Molecular Function. In: Sikosek T ed. Computational Methods in Protein Evolution. Methods in Molecular Biology. New York, NY: Springer, 135–170. DOI: [10.1007/978-1-4939-8736-8_8](https://doi.org/10.1007/978-1-4939-8736-8_8).

2. VanAntwerp J, Finneran P, Dolgikh B, Woldring D. 2022. Ancestral Sequence Reconstruction and Alternate Amino Acid States Guide Protein Library Design for Directed Evolution. In: Traxlmayr MW ed. Yeast Surface Display. Methods in Molecular Biology. New York, NY: Springer US, 75–86. DOI: [10.1007/978-1-0716-2285-8_4](https://doi.org/10.1007/978-1-0716-2285-8_4).

## What is a pseudo-codon probability in ancseq?
