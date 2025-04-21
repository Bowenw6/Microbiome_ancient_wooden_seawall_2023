
# perform alignment with muscle
nohup muscle -align filtered_CK_fasta_16S.fasta -output aligned_filtered_CK_fasta_16S.fasta &
# construct tree with iqtree
nohup iqtree -s aligned_filtered_CK_fasta_16S.fasta -m MFP -B 1000 --bnni -cptime 6000 -T AUTO &

# perform alignment with muscle
nohup muscle -align filtered_Bla_fasta_16S.fasta -output aligned_filtered_Bla_fasta_16S.fasta &
#  construct tree with iqtree
nohup iqtree -s aligned_filtered_Bla_fasta_16S.fasta -m MFP -B 1000 --bnni -cptime 6000 -T AUTO &


# perform alignment with muscle
nohup muscle -align filtered_Gre_fasta_16S.fasta -output aligned_filtered_Gre_fasta_16S.fasta &
# construct tree with iqtree
nohup iqtree -s aligned_filtered_Gre_fasta_16S.fasta -m MFP -B 1000 --bnni -cptime 6000 -T AUTO &


# perform alignment with muscle
nohup muscle -align filtered_Whi_fasta_16S.fasta -output aligned_filtered_Whi_fasta_16S.fasta &
# construct tree with iqtree
nohup iqtree -s aligned_filtered_Whi_fasta_16S.fasta -m MFP -B 1000 --bnni -cptime 6000 -T AUTO &


# perform alignment with muscle
nohup muscle -align filtered_All_fasta_16S.fasta -output aligned_filtered_All_fasta_16S.fasta &
# construct tree with iqtree
nohup iqtree -s aligned_filtered_All_fasta_16S.fasta -m MFP -B 1000 --bnni -cptime 6000 -T AUTO &
