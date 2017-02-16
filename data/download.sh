


#=======================================================================
# Capture Hi-C data from Mifsud2015
#=======================================================================
mkdir -p Mifsud2015
wget -P Mifsud2015 http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2323/E-MTAB-2323.additional.1.zip
unzip Mifsud2015/E-MTAB-2323.additional.1.zip -d Mifsud2015


# reformat promoter-promoter intaction to have each gene pair as separate line.
python ../python/reformat_CaptureHiC.py \
  -i Mifsud2015/TS5_GM12878_promoter-promoter_significant_interactions.txt \
  -o Mifsud2015/TS5_GM12878_promoter-promoter_significant_interactions.txt.genePairs
