# ==== Mutation data cleanup ====
#rm(list = ls(pattern = "ccle"))
require(data.table)
path = "C:/Users/danbeltran/Documents/ECEN766/PQD_Data/"

ccle_mut <- fread(paste0(path, "CCLE_mutations.csv"))
table(ccle_mut$isCOSMIChotspot)
table(ccle_mut$isTCGAhotspot)
table(ccle_mut$Variant_Type)
length(unique(ccle_mut$DepMap_ID))

dim(ccle_mut)
ccle_mut[1,]
colnames(ccle_mut)

# Calculate number of mutations per cell line
temp <- ccle_mut[, c("Variant_Type", "DepMap_ID")]
temp[, nMut := .N, by = "DepMap_ID"]
temp
unique(temp$Variant_Type)

#fwrite(unique(temp$Variant_Type),paste0(path, "OUTPUT_DepMap_21Q2_Mutations_by_Cell.csv"), sep = ',')

# For simplicity, extract only SNP data for now: this discards ~90,000 mutations
# ccle_mut <- ccle_mut[Variant_Type == "SNP"]
dim(ccle_mut)
t <- ccle_mut[, c("DepMap_ID", "Chromosome", "Strand", "Start_position", "End_position")]
dim(unique(t))
length(unique(ccle_mut$DepMap_ID))
# Keep relevant columns/features
# Aside: Should the sequence change be provided, or just whether the SNP is deleterious or not?
ccle_mut <- ccle_mut[, c("DepMap_ID", "Hugo_Symbol", "Chromosome", "Start_position", "End_position", "Strand",
                         "Variant_Classification", "Variant_Type", "isDeleterious",
                         "isTCGAhotspot", "isCOSMIChotspot", "Genome_Change", "cDNA_Change")]
dim(ccle_mut)
length(unique(ccle_mut$DepMap_ID))
table(ccle_mut$isDeleterious)
table(ccle_mut$isTCGAhotspot)
table(ccle_mut$isCOSMIChotspot)

# Daniel B ^ works

# ==== CCLE Mut Overlap with COSMIC CGC ====
# Perhaps it's best to use the mutations in genes that COSMIC considers important, like another paper in
# the field (~500 genes)
# Or, we can use a binary vector for genes and whether they have a deleterious mutation: this will result in 
# ~20,000 parameters
length(unique(ccle_mut$Hugo_Symbol))

length(unique(ccle_mut[isCOSMIChotspot == T]$Hugo_Symbol))
length(unique(ccle_mut[isTCGAhotspot == T]$Hugo_Symbol))
length(unique(ccle_mut[isDeleterious == T]$Hugo_Symbol))

tcga_hotspot_genes <- unique(ccle_mut[isTCGAhotspot == T]$Hugo_Symbol)

# Read COSMIC Cancer Gene Census data
cgc <- fread(paste0(path,"cancer_gene_census.csv"))

dim(cgc)
cgc[1:5, 1:20]
length(unique(cgc$`Gene Symbol`))
length(unique(cgc$HGVSG))
# Get Genes in this census
cgc_genes <- unique(cgc$`Gene Symbol`)
cgc[Tier == 1]
length(unique(cgc$`Genome Location`))  # 922,732
# rm(cgc)

# Subset DepMap mutations based on the CGC genes
sum(unique(ccle_mut$Hugo_Symbol) %in% unique(cgc_genes))
ccle_mut <- ccle_mut[Hugo_Symbol %in% cgc_genes]
length(unique(ccle_mut$DepMap_ID))

sum(ccle_mut$isDeleterious)
ccle_mut[Variant_Classification == "Missense_Mutation"]
length(unique(ccle_mut[isDeleterious == T]$Hugo_Symbol))
ccle_mut[isDeleterious == T]


# TODO: Use CGC to check for overlap with CCLE cell lines, then collapse to whether each of the 700 genes for
# that cell line has a mutation listed in the CGC
length(unique(cgc$`Mutation genome position`))  # ~922,000 unique mutations
unique(ccle_mut$NCBI_Build)  # CCLE is with GRCh 37
unique(cgc$GRCh)  # CGC has GRCh 38
# We must "lift over" the mutations from 37 to 38 before checking for overlap
if (!require(liftOver)) {
  BiocManager::install("liftOver")
  require(liftOver)
  require(rtracklayer)
}

# liftOver requires a chain file to convert 37 to 38: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/

chain_path <- paste0(path, "hg19ToHg38.over.chain")
grch_37_38_chain <- import.chain(chain_path)

# Must add "chr" to start of chromosome names
ccle_mut$Chromosome <- paste0("chr", ccle_mut$Chromosome)
# Must convert positions to GRanges
ccle_mut_gr <- makeGRangesFromDataFrame(df = ccle_mut, keep.extra.columns = T,
                                        seqnames.field = "Chromosome", start.field = "Start_position",
                                        end.field = "End_position", strand.field = "Strand")
length(unique(ccle_mut_gr$DepMap_ID))

# Lift over
lifted_ccle_mut <- liftOver(x = ccle_mut_gr, chain = grch_37_38_chain)
# Convert GRangesList to GRanges
lifted_ccle_mut <- unlist(lifted_ccle_mut)
# Convert back to data.table
lifted_ccle_mut <- as.data.table(lifted_ccle_mut)
# Note: Genome_Change is now out of date!
# Remove chr from seqnames
lifted_ccle_mut$seqnames <- gsub("chr", "", lifted_ccle_mut$seqnames)
# Can find the overlap of Mutation genome position in CGC with a newly created column based on CCLE positions
lifted_ccle_mut[, Mutation_Position := paste0(seqnames, ':', start, '-', end)]

ccle_mut$seqnames <- gsub("chr", "", ccle_mut$Chromosome)
ccle_mut[, Mutation_Position := paste0(seqnames, ':', as.character(Start_position), '-', as.character(End_position))]


length(unique(lifted_ccle_mut$DepMap_ID))

sum(ccle_mut$Mutation_Position %in% unique(cgc$`Genome Location`))

# Now find the overlap with CGC (which already has GRCh38)
subset <- lifted_ccle_mut[Mutation_Position %in% unique(cgc$`Genome Location`)]
table(subset$Variant_Type)
length(unique(subset$DepMap_ID))
# IMPORTANT! There is a loss of 8 cell lines (which do not have a mutation that is in
# CGC) using the Tier 1 data only

# Alternative (March 2021) ====
# Take those mutations that are COSMIC or TCGA hotspots, ignoring CGC
subset <- ccle_mut[isTCGAhotspot | isCOSMIChotspot]

uniqueN(subset$Hugo_Symbol)
table(ccle_mut[isTCGAhotspot == T]$Variant_Classification)
uniqueN(ccle_mut[isTCGAhotspot == T]$Hugo_Symbol)
uniqueN(ccle_mut[isCOSMIChotspot == T]$Hugo_Symbol)

uniqueN(subset$isTCGAhotspot)

### Create a vector of mutations for each cell line with the CGC genes
length(unique(subset$Hugo_Symbol))
sub_dcast <- dcast.data.table(data = subset[, c("DepMap_ID", "Hugo_Symbol")],
                              formula = DepMap_ID ~ Hugo_Symbol, fun.aggregate = length, value.var = "DepMap_ID")
dim(sub_dcast)
sub_dcast[1:5, 1:50]
sum(sub_dcast$A1BG)
sum(sub_dcast$A1CF)


fwrite(sub_dcast, paste0(path, "FILTERD_DepMap_21Q2_Mutations_by_Cell.csv"), sep = ',')

depmap_samples <- fread(paste0(path, "Data/DRP_Training_Data/DepMap_21Q2_Line_Info.csv"))
sub_dcast <- merge(sub_dcast, depmap_samples[, c("DepMap_ID", "stripped_cell_line_name", "primary_disease")],
                   by = "DepMap_ID")
setcolorder(sub_dcast, c("DepMap_ID", "stripped_cell_line_name", "primary_disease"))
sub_dcast[1:5, 1:50]

# Save
fwrite(sub_dcast, paste0(path, "Data/DRP_Training_Data/DepMap_21Q2_Mutations_by_Cell.csv"), sep = ',')
dim(cgc_muts)
cgc_muts[1:5, 1:5]
typeof(cgc_muts[1,2])


