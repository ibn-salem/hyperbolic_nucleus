


#=======================================================================
# Capture Hi-C data from Mifsud2015
#=======================================================================
mkdir -p Mifsud2015
wget -P Mifsud2015 http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2323/E-MTAB-2323.additional.1.zip
unzip Mifsud2015/E-MTAB-2323.additional.1.zip -d Mifsud2015


# reformat promoter-promoter intaction to have each gene pair as separate line.
python ./python/reformat_CaptureHiC.py \
  -i Mifsud2015/TS5_GM12878_promoter-promoter_significant_interactions.txt \
  -o Mifsud2015/TS5_GM12878_promoter-promoter_significant_interactions.txt.genePairs

#=======================================================================
# UCSC lift over chains
#=======================================================================

mkdir -p UCSC
wget -P UCSC http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
gunzip UCSC/*.gz

# download liftOver tool from UCSC:
mkdir -p bin
wget -P bin http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod u+x bin/liftOver


#=======================================================================
# TADs in hESC from Dixon et al. 2012
#=======================================================================
mkdir -p Dixon2012

C="hESC"

wget -P Dixon2012 http://chromosome.sdsc.edu/mouse/hi-c/${C}.domain.tar.gz 
# extract and rename
tar xvfz Dixon2012/${C}.domain.tar.gz -C Dixon2012
cp Dixon2012/${C}/combined/total.combined.domain Dixon2012/${C}.hg18.bed
    
# liftover to hg19
bin/liftOver \
      Dixon2012/${C}.hg18.bed \
      UCSC/hg18ToHg19.over.chain \
      Dixon2012/${C}.hg18.bed.hg19.bed \
      Dixon2012/${C}.hg18.bed.hg19.bed_unmapped.bed
