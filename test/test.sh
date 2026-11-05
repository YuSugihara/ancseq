# Example 1: running ancseq in DNA mode

ancseq -s test_nt.fasta \
       -m DNA \
       -o test_nt

# Example 2: running ancseq in CODON mode

ancseq -s test_nt.fasta \
       -m CODON \
       -o test_codon

# Example 3: running ancseq in AA mode

ancseq -s test_aa.fasta \
       -m AA \
       -o test_aa

# Example 4: running ancseq specifing outgroup

ancseq -s test_nt.fasta \
       -m DNA \
       -o test_outgroup \
       --outgroup AVR-Mgk5_GE16_2

# Example 5: running ancseq with --fast option

ancseq -s test_nt.fasta \
       -m DNA \
       -o test_fast \
       --fast \
       --outgroup AVR-Mgk5_GE16_2

# Example 6: halting calculation of codon probabilities

ancseq -s test_nt.fasta \
       -m DNA \
       -o test_no_codon_prob \
       --stop-codon-prob
       