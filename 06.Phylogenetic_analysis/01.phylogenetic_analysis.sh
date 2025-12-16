mafft --auto --thread 4 --anysymbol --maxiterate 1000 my_sequence.fasta > my_sequence.mafft
trimal -in my_sequence.mafft -out my_sequence.trimal -automated1
fasttree -nt -gtr my_sequence.trimal > my_sequence.tree