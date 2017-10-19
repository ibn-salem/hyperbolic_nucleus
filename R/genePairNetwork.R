#!/usr/bin/Rscript
#=======================================================================
#
#   build a gene-pair network based on Caputer Hi-C based chromatin 
#	interaction frequencies.
#
#=======================================================================

require(biomaRt)        # to retrieve human paralogs from Ensembl
require(TxDb.Hsapiens.UCSC.hg19.knownGene)

#-------------------------------------------------------------------
# 1. get ENSG genes with positions
#-------------------------------------------------------------------
seqInfo <- seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)

ensemblGRCh37 <- useMart(host = "grch37.ensembl.org", 
                         biomart = "ENSEMBL_MART_ENSEMBL", 
                         dataset = "hsapiens_gene_ensembl", 
                         verbose = FALSE)

geneAttributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", 
                   "start_position", "end_position", "strand", 
                   "transcription_start_site","gene_biotype")
geneFilters = "chromosome_name"

# read "normal" human chromosome names (without fixes and patches)
geneValues = c(1:22, "X", "Y")
allGenes = getBM(attributes = geneAttributes, mart = ensemblGRCh37, 
                 filters = geneFilters, values = geneValues)

# filter for known genes only and only for protein coding genes
#knownCodingGenes = allGenes[allGenes$status=="KNOWN" & allGenes$gene_biotype== "protein_coding",]

# unique gene entry by ENSG ID symbol:
#genes = knownCodingGenes[!duplicated(knownCodingGenes$ensembl_gene_id),]
genes = allGenes[!duplicated(allGenes$ensembl_gene_id),]
#~ rownames(genes) <- genes$ensembl_gene_id

# make GRanges object for all known prot coding genes

tssGR = GRanges(
        paste0("chr", genes$chromosome_name),
        IRanges(genes$transcription_start_site, genes$transcription_start_site),
        strand = ifelse(genes$strand == 1, '+', '-'), 
        names = genes$ensembl_gene_id, 
        genes[,c("hgnc_symbol", "gene_biotype")],
        seqinfo = seqInfo
        )
names(tssGR) = genes$ensembl_gene_id
tssGR <- sort(tssGR)

#-------------------------------------------------------------------
# 2. parse Caputure Hi-C data 
#-------------------------------------------------------------------
# promoter-promoter interaction from Capture Hi-C (Mifsud2015a)
# CAPTURC_FILE="data/Mifsud2015/TS5_GM12878_promoter-promoter_significant_interactions.txt"
CAPTURC_FILE = "data/Mifsud2015/TS5_GM12878_promoter-promoter_significant_interactions.txt.genePairs"

# Check if file exists, if not, download and format via shell script
if (!file.exists(CAPTURC_FILE)) {
  system("sh data/download.sh")
  system("mv Mifsud2015 data")
}

# due to sparse matrix data structure non available pairs will get 0 counts. This needs to be addressed in downstream analysis
# see function parseCaptureHiC(inFile=CAPTURC_FILE, tssGR) in paralog_regulation project

# parse input file as data frame
classes <- sapply(read.delim(CAPTURC_FILE, nrows = 5, header=TRUE, stringsAsFactors=FALSE ), class)
inData <- read.delim(CAPTURC_FILE, header=TRUE, colClasses = classes, stringsAsFactors=FALSE)

#-------------------------------------------------------------------
# 3. Transform in gene pair file
#-------------------------------------------------------------------


#-------------------------------------------------------------------
# 4. Annotate with genomic distance
#-------------------------------------------------------------------

s1 <- start(tssGR)[match(inData[,1], names(tssGR))]
s2 <- start(tssGR)[match(inData[,2], names(tssGR))]

chr1 <- as.vector(seqnames(tssGR))[match(inData[,1], names(tssGR))]
chr2 <- as.vector(seqnames(tssGR))[match(inData[,2], names(tssGR))]

inData$dist <- abs(s2 - s1) / 1000
inData$dist[chr1 != chr2] <- NA

#-------------------------------------------------------------------
# 5. Output tab-separated file
#-------------------------------------------------------------------
dir.create("results")

message("INFO: Network has ", nrow(inData), " interactions between ", 
        length(unique(c(inData[,1], inData[,2]))), " genes")

write.table(inData, file = "results/Mifsud2015_GM12787_with_dist.tsv",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


# filter for pairs with log(exp/obs) >= 10 
fltDF <- inData[inData$log.obs.exp. >= 10,]

message("INFO: Filtered network has ", nrow(fltDF), " interactions between ", 
        length(unique(c(fltDF[,1], fltDF[,2]))), " genes")

write.table(inData, 
            file = "results/Mifsud2015_GM12787_with_dist.lnObsExp_10.tsv",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#Save the data frames
save(inData, fltDF, file = "results/edge_lists.RData")

system("rm -r data/Mifsud2015")

#-------------------------------------------------------------------
# 6. Write node annotation table (chromsome for each gene):
#-------------------------------------------------------------------
require(dplyr) # for joining nodes with annotation table

# take the unique list of ENSG ids that are present in the network
genesInNetwork <- unique(c(inData[,1], inData[,2]))

# build data.frame from node list
netGenesDF <- data.frame(
  id = genesInNetwork,
  stringsAsFactors = FALSE
)

# build a data.frame for the annotations
geneDF <- as.data.frame(mcols(tssGR))
geneDF$chr <- as.character(seqnames(tssGR))
geneDF$names <- as.character(geneDF$names)

# join node list and annotation by using NAs for not available gene IDs.
nodeDF <- netGenesDF %>%
  left_join(geneDF, by = c("id" = "names"))

# out put note table
write.table(nodeDF, file = "results/Mifsud2015_GM12787_with_dist.nodeDF.tsv",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#Save the data frames as .RData file
save(nodeDF, file = "results/nodeDF.RData")

