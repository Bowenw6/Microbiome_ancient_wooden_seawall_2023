####################################
#############    16S    ############
####################################

conda activate /Users/bowenwang/miniconda3/envs/qiime2-2020.11

asv_table_16S="ASV_filtered_16S.tsv"
asv_table_ITS="ASV_filtered_ITS.tsv"


biom convert -i $asv_table_16S -o asv_table_16S.biom --table-type="OTU table" --to-json
biom convert -i $asv_table_ITS -o asv_table_ITS.biom --table-type="OTU table" --to-json







conda activate picrust2-2.6.0-github


picrust2_pipeline.py -s filtered_ASVs_16S.fa -i asv_table_16S.biom -o picrust2_out_pipeline_16S -p 4
picrust2_pipeline.py -s filtered_ASVs_ITS.fa -i asv_table_ITS.biom -o picrust2_out_pipeline_ITS -p 4












# ####################################
# #############    ITS    ############
# ####################################


conda create -n funguild python=3.8 -y
conda activate funguild

pip install requests

git clone https://github.com/UMNFuN/FUNGuild.git
cd FUNGuild








conda activate /Users/bowenwang/miniconda3/envs/qiime2-2020.11

asv_table_ITS="ITS_for_FUNGuild.tsv"

python ./FUNGuild/Guilds_v1.1.py -otu $asv_table_ITS -db fungi -m -u

