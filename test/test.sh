# Example 1: running ancseq in DNA mode

ancseq -s data/test_nt.fasta \
       -m DNA \
       -o test_nt

# Example 2: running ancseq in CODON mode

ancseq -s data/test_nt.fasta \
       -m CODON \
       -o test_codon

# Example 3: running ancseq in AA mode

ancseq -s data/test_aa.fasta \
       -m AA \
       -o test_aa

# Example 4: running ancseq specifing outgroup

ancseq -s data/test_nt.fasta \
       -m DNA \
       -o test_outgroup \
       --outgroup AVR-Mgk5_GE16_2

# Example 5: running ancseq with --fast option

ancseq -s data/test_nt.fasta \
       -m DNA \
       -o test_fast \
       --fast \
       --outgroup AVR-Mgk5_GE16_2

# Example 6: halting calculation of codon probabilities

ancseq -s data/test_nt.fasta \
       -m DNA \
       -o test_no_codon_prob \
       --stop-codon-prob

# Example 7: running the full MLP alignment with its original tree and model

ancseq -s data/test_mlp.fasta \
       -m DNA \
       -o test_mlp \
       --asr-only \
       --model GTR+F+I+G4 \
       --outgroup GCA_036507275.1_Cucurbitales_Sicyos_edulis_Helixer_CM071850.1_001929.1
