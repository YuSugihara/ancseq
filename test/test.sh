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

# Example 6: running a larger DNA alignment with a fixed tree and model

ancseq -s data/test_nt_large.fasta \
       -m DNA \
       -o test_nt_large \
       --asr-only \
       --model GTR+F+I+G4 \
       --outgroup seq157
