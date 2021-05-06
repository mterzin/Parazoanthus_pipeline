###############################
# Stats and graph compilation #
###############################
# Tutorials: 
# 1. https://cran.r-project.org/web/packages/poppr/vignettes/mlg.html
# 2. This one seems to have everything I need: https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html
# Combine 2 with this one: http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
# 4. Suggestions are here: https://juliehopper.wordpress.com/tag/poppr/

# Install the packages I need
# R will first look if the packages are installed, 
# and then later install them if not  

list.of.packages <- c("adegenet", "ape", "ggplot2","gtools", "devtools", "dplyr", "hierfstat",
  "igraph", "pegas", "phangorn", "poppr", "viridis", "pheatmap", "mmod",
  "RColorBrewer", "reshape2", "vcfR", "vegan", "fsthet", "lattice", "treemap", "magrittr", "diveRsity")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Installing adegenet was a bit tricky, done from the terminal using this tutorial in the end: 
# https://zoomadmin.com/HowToInstall/UbuntuPackage/r-cran-adegenet
library(adegenet)
library(ape)
library(ggplot2)
library(gtools)
library(igraph)
library(pegas)
library(phangorn)
library(poppr)
library(RColorBrewer)
library(reshape2)
library(vcfR)
library(vegan)
library(fsthet)
library(lattice)
library(dplyr)
library(treemap)
library(magrittr)
library(pheatmap)
library(viridis)
library(mmod)
library(hierfstat)
# Followed this to install diveRsity: 
# https://github.com/kkeenan02/diveRsity
library(diveRsity)

# Setting work directory
setwd("/home/markoterzin/Documents/Parazoanthus_2bRAD/Marko/Script_validation_neutral_loci/O_and_Y")

# Let's load the VCF file
C.VCF <- read.vcfR("populations.snps.vcf")
C.VCF 

# Since VCF data does not typically include any sort of population information, I have to load this data 
# separately from the text-delimited file (my popmap file)
pop.data <- read.table("popmap.tsv", sep = "\t", header = TRUE)

# We can now check that all the samples in the VCF and the population data frame are included:
all(colnames(C.VCF@gt)[-1] == pop.data$Sample_ID)

# Converting vcf to genlight format
C.genlight <- vcfR2genlight(C.VCF)

# Set ploidy to 2, because Antipathella subpinnata is a diploid
ploidy(C.genlight) <- 2

# Our biological question requires predetermined populations, so we add
# the Population column from our pop.data data frame to the pop slot of our genlight object:
pop(C.genlight) <- pop.data$Population

# Let's check out this genlight object:
C.genlight@gen # Tells me which samples have a lot of missing data
C.genlight@n.loc # Tells me how many SNPs I have
C.genlight@ind.names # Returns labels for each individual
C.genlight@loc.names # Returns names for each locus/SNP
C.genlight@position # Position of every SNP
C.genlight@ploidy # returns/sets ploidy of the individuals
C.genlight@pop # returns/sets a factor grouping individuals
C.genlight@other # returns/sets misc information stored as a list.
C.genlight@hierarchy # not set to anything
C.genlight@strata # still not set

# Filtering out missing data (at locus level)
# and replacing with mean allele frequency
C.genlight_no_NAs <- tab(C.genlight, NA.method = "mean")
# Converting it back to genlight
C.genlight_no_NAs <- as.genlight(C.genlight_no_NAs)
C.genlight_no_NAs

# Ploidy and population information need to be reset now
ploidy(C.genlight_no_NAs) <- 2
pop(C.genlight_no_NAs) <- pop.data$Population
# Checking
C.genlight_no_NAs

# Plot to see the coverage
# Open pdf
pdf("01_Coverage.pdf", 
    width = 7, height = 5)
# Plot
plot(C.genlight, main="With missing data")
plot(C.genlight_no_NAs, main="NAs replaced with mean allele frequency")
# Close the pdf file
dev.off()

#########################
###########################################
###**** Genotype accumulation curve ****###
###########################################
# Now that I have my object, I can try this tutorial: https://grunwaldlab.github.io/Population_Genetics_in_R/First_Steps.html
# Good to have a first look at my data

# The genotype accumulation curve: a tool that allows you to assess how much power you have to 
# discriminate between unique individuals given a random sample of n loci. This analysis is 
# particularly important for clonal organisms to confirm that a plateau has been reached in the 
# number of loci necessary to discriminate individuals.

#############
# IMPORTANT # 
# I cannot do this on a genlight object, so import data_genepop
data_genepop <- read.genepop("populations.snps.gen")
# Include Population information
pop(data_genepop) <- pop.data$Population
data_genepop

# Open pdf
pdf("02A_data_Genotype_accummulation_curve.pdf", 
width = 200, height = 8)
# Plot
gac <- genotype_curve(data_genepop, sample = 1000, quiet = TRUE)
# Close the pdf file
dev.off()

# We specified sample = 1000 in our function call. This means that for each boxplot, n loci were 
# randomly sampled 1000 times in order to create the distribution.

# It is expected that a big number of loci will be needed to
# discriminate between clonal organisms.

#########################
##########################################
###**** Hardy-Weinberg equilibrium ****###
##########################################
# Next, let’s determine if our populations are in Hardy-Weinberg equilibrium. We will test for HWE 
# using the function hw.test() from the pegas package. This will compute the χ2 statistic over the 
# entire data set and compute two P-values, one analytical and one derived from permutations.
data_genepop.full <- hw.test(data_genepop, B = 1000) 
# performs 1000 permuatations (Monte Carlo tests basically)
data_genepop.full

# If we wanted to check what the HWE statistic for each population is, we should first separate the 
# populations with the function seppop().
data_genepop.pop <- seppop(data_genepop) %>% lapply(hw.test, B = 1000)
data_genepop.pop

# Visualization of population-wise p-values

# Now we have one matrix per sample, but all we care about are the p-values, which are in the 
# third column. We can use the functions sapply and [ to loop to create a matrix that only contains 
# populations in columns and loci in rows.
data_genepop.full.mat <- sapply(data_genepop.pop, "[", i = TRUE, j = 3) # Take the third column with all rows
data_genepop.full.mat

# This output is still hard to sift through. An easy way to analyze this is by visualizing this as a 
# heatmap. Since we only care whether or not a given locus is in or out of HWE, we will thus define 
# an α value and set everything above that value to 1 so that we can visually detect candidate loci 
# where we might reject the Ho
alpha  <- 0.05
data_genepop.mat <- data_genepop.full.mat
data_genepop.mat[data_genepop.mat > alpha] <- 1
data_genepop.mat

# Now we can create a simple heatmap with levelplot.
# Open pdf
pdf("03A_HWE_heatmap_no_FDR_correction.pdf", 
width = 400, height = 5)
# Plot
levelplot(data_genepop.mat)
# Close the pdf file
dev.off()
# All loci shown in pink are loci suspected of not being in HWE with p ≤ 0.05

# But this is based on non-adjusted p-values, and we want to remove false positives!!!
# This will be done using the 'stats' R package: 
# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html

# Let's first take a look at which p adjustment methods are available
p.adjust.M <- p.adjust.methods[p.adjust.methods != "fdr"]
p.adjust.M # So these are all the methods

# My p vals can now be adjusted using all 6 available methods
p.adj <- sapply(p.adjust.M, simplify="matrix", function(meth) p.adjust(data_genepop.full.mat, 
                                                                       meth))
p.adj

# Export table
write.csv(p.adj, file = "padj.csv", quote = F)

# We can try to correct using any of the 6 correction methods available in R stats package:
# "holm"       "hochberg"   "hommel"     "bonferroni" "BH"         "BY"         "none" 

# Let's check if any of the correction methods gave significant p-adjusted values (<0.05)

# No correction
None_pval <- p.adj[,7,drop=FALSE] # extracting the 7th column (no correction)
None <- None_pval <= 0.05 # We get a logical index by comparing with 0.05
sum(None) # Get the sum of the TRUE values to find the number of rows that have p.value less than 0.05

# Holm
Holm_pval <- p.adj[,1,drop=FALSE] # extracting the 1st column (Holm correction)
Holm <- Holm_pval <= 0.05 # We get a logical index by comparing with 0.05
sum(Holm) # Get the sum of the TRUE values to find the number of rows that have p.value less than 0.05

# Hochberg
Hochberg_pval <- p.adj[,2,drop=FALSE] # extracting the 2nd column (Hochberg correction)
Hochberg <- Hochberg_pval <= 0.05 # We get a logical index by comparing with 0.05
sum(Hochberg) # Get the sum of the TRUE values to find the number of rows that have p.value less than 0.05

# Hammel
Hammel_pval <- p.adj[,3,drop=FALSE] # extracting the 3rd column (Hammel correction)
Hammel <- Hammel_pval <= 0.05 # We get a logical index by comparing with 0.05
sum(Hammel) # Get the sum of the TRUE values to find the number of rows that have p.value less than 0.05

# Bonferroni
Bonferroni_pval <- p.adj[,4,drop=FALSE] # extracting the 4th column (Bonferroni correction)
Bonferroni <- Bonferroni_pval <= 0.05 # We get a logical index by comparing with 0.05
sum(Bonferroni) # Get the sum of the TRUE values to find the number of rows that have p.value less than 0.05

# Benjamini Hochberg
BH_pval <- p.adj[,5,drop=FALSE] # extracting the 5th column (Benjamini Hochberg correction)
BH <- BH_pval <= 0.05 # We get a logical index by comparing with 0.05
sum(BH) # Get the sum of the TRUE values to find the number of rows that have p.value less than 0.05

# Benjamini Yekutieli
BY_pval <- p.adj[,6,drop=FALSE] # extracting the 6th column (Benjamini Yekutieli)
BY <- BY_pval <= 0.05 # We get a logical index by comparing with 0.05
sum(BY) # Get the sum of the TRUE values to find the number of rows that have p.value less than 0.05

# Now we can remove loci out of HWE
# For instance: trying with Benjamini Yekutieli correction - BY
data_BY_pvals <- p.adj[,6,drop=FALSE] # Take the sixth (BY) column with all rows
data_BY_pvals <- as.data.frame(data_BY_pvals)
# Adding loci names to the data frame
data_genepop_loci.names <- data_genepop@loc.fac
data_genepop_loci.names <- as.data.frame(data_genepop_loci.names)
# Removing duplicated rows
data_genepop_loci.names <- distinct(data_genepop_loci.names, .keep_all = TRUE)
# Check how many populations you have in your dataset
data_genepop@pop
# Multiplying by the number of populations in times=
data_genepop_loci.names <- data_genepop_loci.names[rep(
  row.names(data_genepop_loci.names), times=14), 1]
data_genepop_loci.names <- as.data.frame(data_genepop_loci.names)
# Combining the 2 dataframes to get loci names and BY corrected p vals 
P_loci.names_and_BY <- c(data_genepop_loci.names, data_BY_pvals)
P_loci.names_and_BY <- as.data.frame(P_loci.names_and_BY)

# Finding all loci where BY corrected p val is significant
P_loci.names_and_BY_below_0.05 <- P_loci.names_and_BY[!rowSums
                                       (P_loci.names_and_BY[-1] > 0.05),]
# Removing empty fields
P_loci.names_and_BY_below_0.05 <- P_loci.names_and_BY_below_0.05[complete.cases
                                                                 (P_loci.names_and_BY_below_0.05), ]

# Keeping only the loci names
P_loci_out_of_HWE_BY <- P_loci.names_and_BY_below_0.05
P_loci_out_of_HWE_BY$BY <- NULL
# Removing duplicated names
P_loci_out_of_HWE_BY <- distinct(P_loci_out_of_HWE_BY, .keep_all = TRUE)

# Removal of loci out of HWE
##############################################
# create vector of all loci
all_loci <- locNames(data_genepop)
# create vector containing loci to remove
removeloc <- as.matrix(P_loci_out_of_HWE_BY)
# create list of loci to keep
keeploc <- setdiff(all_loci, removeloc)
# filter loci in genind object
data_genepop_no_HWE <- data_genepop[loc = keeploc]

# check number of loci in genind obj
length(locNames(data_genepop))
length(locNames(data_genepop_no_HWE))

#########################
############################
###**** Missing data ****###
############################
pdf("03B_Average_missing_loci_per_population.pdf", 
width = 350, height = 6)
# Plot
info_table(data_genepop_no_HWE, type = "missing", plot = TRUE)
# Close the pdf file
dev.off()

######################
######################
###**** Ploidy ****###
############################
pdf("03C_Ploidy_plot.pdf", 
width = 200, height = 10)
# Plot
data_genepop.ploidy <- info_table(data_genepop_no_HWE, 
                                          type = "ploidy", 
                                          plot = TRUE, 
                                          low = "black", 
                                          high = "orange")
# Close the pdf file
dev.off()

######################
###############################################
###**** Calculating genotypic diversity ****###
###############################################
# Let us get a first impression of the diversity found in this data using the summary function, poppr:
data_summary_stats <- poppr(data_genepop_no_HWE)
write.csv(data_summary_stats, 
          file = "Spreadsheet_02A_summary_stats_data_genepop.csv", 
          row.names = T)
# Abbreviation 	Statistic
# Pop 	Population name.
# N 	Number of individuals observed.
# MLG 	Number of multilocus genotypes (MLG) observed.
# eMLG 	The number of expected MLG at the smallest sample size ≥ 10 based on rarefaction
# SE 	Standard error based on eMLG.
# H 	Shannon-Wiener Index of MLG diversity (Shannon, 2001).
# G 	Stoddart and Taylor’s Index of MLG diversity (Stoddart & Taylor, 1988).
# lambda 	Simpson’s Index (Simpson, 1949).
# E.5 	Evenness, E5 (Pielou, 1975; Ludwig & Reynolds, 1988; Grünwald et al., 2003).
# Hexp 	Nei’s unbiased gene diversity (Nei, 1978).
# Ia 	The index of association, IA (Brown, Feldman & Nevo, 1980; Smith et al., 1993).
# rbarD 	The standardized index of association, r¯d[@].

#########################
########################################
###**** Deal with missing values ****###
########################################
# Tutorial here: https://grunwaldlab.github.io/Population_Genetics_in_R/Population_Strata.html
#######################################################

# Removing low coverage individuals first: more than 30% missing data
data_no_low_coverage <- missingno(data_genepop_no_HWE, type = "geno", cutoff = 0.3)
# Replacing NAs in missing loci with mean
data_no_NAs <- tab(data_no_low_coverage, NA.method="mean")
data_final <- as.genind(data_no_NAs)

# Removing the EST individual, as it's only one replicate
data_final <- data_final[indNames(data_final) != "T_Pax_EST_4_concatenated"]

# Check for missing data now
missingno(data_final, type = "loci")

# Modification of the popmap file
pop.data_final2 <- indNames(data_final)
pop.data_final2 <- as.data.frame(pop.data_final2)
colnames(pop.data_final2) <- "Sample_ID"
pop.data_final <- left_join(pop.data_final2, pop.data)
pop(data_final) <- pop.data_final$Population
# Checking the object
data_final@pop

# Taking a look at the object
data_final
# This is the genind object for the 'ALL LOCI' dataset, containing neutral loci and outliers

#########################
#############################################################
###**** Importing strata to use a hierarchical design ****###
#############################################################
# I want test the effect of 'Morphotype' and 'Population' on genetic structuring

# Import strata
strata <- read.delim("strata.tsv")
# Now removing all the individuals that are not in the genind object
strata <- left_join(pop.data_final2, strata)
# And removing the SampleID column
strata$Sample_ID <- NULL

# Now add to my object
strata(data_final) <- strata
nameStrata(data_final) <- ~MorphotypePopulation/Morphotype/Population
# let's look at the dataset now
data_final
# We see that we have three vectors of different names in the 'strata' slot for microbov (MorphotypePopulation, 
# Morphotype, Population)

#########################
########################################
###**** Looking for OUTLIER loci ****###
########################################
###############
# in BayeScan #
###############

# Getting the BayeScan format
# I need to say which populations I want
pop.select <- c("ALA", "ALAOR", "CAMP", "CUB", "GIA", "PTF", "PUG", "PVEN", "PVENOR", "SAR", "GIAOR", "SAROR")

# Converting genind to bayescan
data_final_bayescan <- radiator::genomic_converter(
  data = data_final,
  strata = "strataBayeScan_P.tsv",
#  whitelist.markers = "whitelist.filtered.markers.tsv",
#  blacklist.id = "blacklist.id.tsv",
  output = "bayescan",
  filename = "bayescan.haplotypes"
)

# to run BayeScan:
# use the bayescan.haplotypes.txt file and run from the terminal (in the right directory)

# This is the command I used
# export PATH=/home/markoterzin/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin${PATH} && /home/markoterzin/Programs_Bioinformatics/BayeScan/BayeScan2.1/binaries/BayeScan2.1_linux64bits bayescan.haplotypes.txt -o bayescanLongRunOD10 -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 10

plot_bayescan<-function(res,FDR=0.05,size=1,pos=0.35,highlight=NULL,name_highlighted=F,add_text=T)
{
  if (is.character(res))
    res=read.table(res)
  
  colfstat=5
  colq=colfstat-2
  
  highlight_rows=which(is.element(as.numeric(row.names(res)),highlight))
  non_highlight_rows=setdiff(1:nrow(res),highlight_rows)
  
  outliers=as.integer(row.names(res[res[,colq]<=FDR,]))
  
  ok_outliers=TRUE
  if (sum(res[,colq]<=FDR)==0)
    ok_outliers=FALSE;
  
  res[res[,colq]<=0.0001,colq]=0.0001
  
  # plot
  plot(log10(res[,colq]),res[,colfstat],xlim=rev(range(log10(res[,colq]))),xlab="log10(q value)",ylab=names(res[colfstat]),type="n")
  points(log10(res[non_highlight_rows,colq]),res[non_highlight_rows,colfstat],pch=19,cex=size)
  
  if (name_highlighted) {
    if (length(highlight_rows)>0) {
      text(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],row.names(res[highlight_rows,]),col="red",cex=size*1.2,font=2)
    }
  }
  else {
    points(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],col="red",pch=19,cex=size)
    # add names of loci over p and vertical line
    if (ok_outliers & add_text) {
      text(log10(res[res[,colq]<=FDR,][,colq])+pos*(round(runif(nrow(res[res[,colq]<=FDR,]),1,2))*2-3),res[res[,colq]<=FDR,][,colfstat],row.names(res[res[,colq]<=FDR,]),cex=size)
    }
  }
  lines(c(log10(FDR),log10(FDR)),c(-1,1),lwd=2)
  
  return(list("outliers"=outliers,"nb_outliers"=length(outliers)))
}

  # Visualize the outliers!
  
  pdf("06_Bayescan_outliers.pdf",
      width = 10, height = 4)
  # Plotting
  bayescan_results <- plot_bayescan("03_radiator_genomic_converter_20210112@2137/bayescanLongRunOD10_fst.txt", # bayescan file
                                    0,
                                    FDR=0.05)
  # Close pdf
  dev.off()
  
  bayescan_results$outliers # the list of outliers
  bayescan_results$nb_outliers # number of outliers
  
  # Now removing those from the data_final object
  bayescan_outliers <- as.data.frame(bayescan_results$outliers)
  loci_no_HWE <- locNames(data_final)
  loci_no_HWE <- as.data.frame(loci_no_HWE)
  # Adding number_ID column
  loci_no_HWE <- loci_no_HWE %>% mutate(id = row_number())
  names(bayescan_outliers)[names(bayescan_outliers) == "bayescan_results$outliers"] <- "id"
  remove_loci <- left_join(bayescan_outliers,loci_no_HWE)
  remove_loci$id <- NULL
  remove_loci # These are the ones that need to be removed! 
  remove_loci <- data.matrix(remove_loci, rownames.force = NA)
  
  data_neutral_bayescan=data_final[loc=-remove_loci] # It worked!
  
  # Checking the number of loci
  data_final
  data_neutral_bayescan
  
  ###################################################################
  # Now trying the analysis with all loci and outliers/neutral only #
  ###################################################################
  # The second R script will be used for the 'ALL LOCI' dataset
  
  ###################################################
  ###**** Multi-locus genotype (MLG) analysis ****###
  ###################################################
  # Using the naive algorithm
  pdf("05A_MLGs_plot.pdf",
      width = 8, height = 6)
  # Plot
  data_genepop.tab <- mlg.table(data_neutral_bayescan)
  # Close the pdf file
  dev.off()
  
  ###################################################
  ###**** Multi-locus genotype (MLG) analysis ****###
  ######******    with a 3% threshold    ******######
  ###################################################
  # This tutorial: https://grunwaldlab.github.io/poppr/reference/mlg.filter.html
  pc <- as.genclone(data_neutral_bayescan, threads = 1L) # convert to genclone object
  data_MLGs_thr_0.03 <- mlg.filter(pc, 
                                   threshold = 0.03, # 3 %
                                   distance = "nei.dist", # Nei distances
                                   threads = 1L,
                                   stats = 'ALL')
  
  # Visualizing in ggplot2
  library(ggplot2)
  
  # But let's set the colors first before making a plot
  cols <- c("khaki", # ALA 
            "orange", # ALAOR
            "yellow1", # CAMP
            "gold2", # CUB
            "lightgoldenrod", # GIA
            "moccasin", # PTF
            "yellow3", # PUG
            "snow3", # PVEN
            "orangered", # PVENOR
            "peachpuff3", # SAR
            "orange2", # GIAOR
            "orange3") # SAROR
  
  # Just to see how their data set looks like
  # data(mpg, package="ggplot2")
  # mpg <- read.csv("http://goo.gl/uEeRGu")
  
  # Scatterplot
  theme_set(theme_bw())  # pre-set the bw theme.
  ind.namesMLG <- indNames(data_neutral_bayescan)
  ind.namesMLG <- as.character(ind.namesMLG) # These are my individuals
  MLGs_0.03 <- data_MLGs_thr_0.03$MLGS # And the MLGs constructed based on 0.03 % threshold
  popinfoMLG <- data_neutral_bayescan@pop # This is the information on populations
  
  # Combining into one data frame
  MLG_plot <- cbind(ind.namesMLG, popinfoMLG, MLGs_0.03)
  MLG_plot <- as.data.frame(MLG_plot)
  MLG_plot
  
  # Let's visualize!
  pdf("05_MLG_0.03_threshold.pdf",
      width = 16, height = 9)
  # Plot
  g <- ggplot(MLG_plot, 
              aes(x=ind.namesMLG, 
                  y=MLGs_0.03, 
                  colour=data_neutral_bayescan@pop))
  # g + geom_count(col="lightblue4", show.legend=F) + # choosing the color of the points
  g + geom_point(size=2, stat='identity') + 
    scale_color_manual(values = cols) +
    facet_grid(~data_neutral_bayescan@pop, # placing the ind within their pops
               scales = "free_x", 
               space = "free") + 
    # facet_grid(cols = vars(Savaglia_final_neutral2@pop)) +
    labs(subtitle="based on a 3% similarity cut-off (Nei distances)", # putting the title / axis names
         y="MLGs",
         x="Individual names", 
         title="MLGs Plot") +
    theme(axis.text.x = element_text(angle = 90, # rotating the sample ID for 90 degrees 
                                     vjust = 0.5, 
                                     hjust=1)) 
  # Close pdf
  dev.off()
  
  #########################
  ################################################
  ###**** Locus stats, heterozygosity, HWE ****###
  ################################################
  
  # A rigorous population genetic analysis looks closely at the data to assess quality and identify 
  # outliers or problems in the data such as erroneous allele calls. This chapter focuses on analysis 
  # on a per-locus level. While there are statistics that analyze populations across loci, it is 
  # important to analyze each locus independently to make sure that one locus is not introducing bias 
  # or spurious errors into the analysis. Locus summary statistics: A quick way to assess quality of 
  # the data is to determine the number, diversity, expected heterozygosity, and evenness of the 
  # alleles at each locus.
  locus_table(data_neutral_bayescan)
  ## 
  ## allele = Number of observed alleles
  # I have 2 alleles within each locus because I used the --write_single_snp
  ## 
  ## 1-D = Simpson index
  ## 
  ## Hexp = Nei's 1978 gene diversity
  
  # We can also do this for each population
  locus_table(data_neutral_bayescan, pop = "ALA")
  locus_table(data_neutral_bayescan, pop = "CAMP")
  locus_table(data_neutral_bayescan, pop = "CUB")
  locus_table(data_neutral_bayescan, pop = "ALAOR")
  locus_table(data_neutral_bayescan, pop = "GIA")
  locus_table(data_neutral_bayescan, pop = "GIAOR")
  locus_table(data_neutral_bayescan, pop = "PTF")
  locus_table(data_neutral_bayescan, pop = "PUG")
  locus_table(data_neutral_bayescan, pop = "PVEN")
  locus_table(data_neutral_bayescan, pop = "PVENOR")
  locus_table(data_neutral_bayescan, pop = "SAR")
  locus_table(data_neutral_bayescan, pop = "SAROR")
  # But we can remove these monomorphic loci, and this is important because they are phylogenetically 
  # uninformative!
  nLoc(data_neutral_bayescan)  # Let's look at our data set, note how many loci we have.
  informative_data_neutral_bayescan <- informloci(data_neutral_bayescan)
  nLoc(informative_data_neutral_bayescan)
  
  #########################
  ###########################################################
  ###**** Genotypic richness, diversity, and evenness ****###
  ###########################################################
  # Tutorial here: https://grunwaldlab.github.io/Population_Genetics_in_R/Genotypic_EvenRichDiv.html
  
  # Let's have a look at populations only
  setPop(data_neutral_bayescan) <- ~Population
  data_neutral_bayescan
  
  # To calculate genotypic richness, diversity, and evenness, we can use the poppr function:
  data_final_diversity <- poppr(data_neutral_bayescan)
  data_final_diversity
  
  # Genotypic richness
  # The number of observed MLGs is equivalent to genotypic richness. A type of a rarefaction method... 
  # (See tutorial)
  
  # Open pdf
  pdf("08A_Genotypic_richness_rarefaction.pdf", 
  width = 10, height = 10)
  # Plot
  data_final.tab <- mlg.table(data_neutral_bayescan, plot = FALSE)
  min_sample <- min(rowSums(data_final.tab))
  rarecurve(data_final.tab, sample = min_sample, xlab = "Sample Size", ylab = "Expected MLGs")
  title("Rarefaction of all my populations: I get a linear correlation because I don't have any clones")
  # Close the pdf file
  dev.off()
  # I get a linear correlation because I don't have any clones!
  
  # Genotypic diversity
  
  # Diversity measures incorporate both genotypic richness and abundance. There are three measures of 
  # genotypic diversity employed by poppr, the Shannon-Wiener index (H), 
  # Stoddart and Taylor’s index (G), and Simpson’s index (lambda)
  N  <- data_final_diversity$N  # number of samples
  N
  lambda <- data_final_diversity$lambda # Simpson's index
  lambda
  # This is repetitive, and was already calculated as lambda before
  (N/(N - 1)) * lambda  # Corrected Simpson's index
  # I get 1 for everything after correction...
  
  # Genotypic evenness
  
  # Evenness is a measure of the distribution of genotype abundances, where in a population with 
  # equally abundant genotypes yields a value equal to 1 and a population dominated by a single 
  # genotype is closer to zero. That is why I am getting 1!
  
  # Open pdf
  pdf("08B_Genotypic_eveness.pdf", 
  width = 10, height = 10)
  # Plot
  data_final_tab <- mlg.table(data_neutral_bayescan)
  # Close pdf
  dev.off()

#########################
######################################
###**** Linkage disequilibrium ****###
######################################
# Tutorial here: https://grunwaldlab.github.io/Population_Genetics_in_R/Linkage_disequilibrium.html

# Linkage disequilibrium test is useful to determine if populations are clonal (where significant disequilibrium
# is expected due to linkage among loci) or sexual (where linkage among loci is not expected). The null 
# hypothesis tested is that alleles observed at different loci are not linked if populations are sexual while 
# alleles recombine freely into new genotypes during the process of sexual reproduction. In molecular ecology 
# we typically use the index of association or related indices to test this phenomenon.
# We will analyze all the populations with the index of association and use 999 permutations of the data in 
# order to get a p-value. Note that the p-value is calculated with the original observation included.

# Open pdf
pdf("09_Linkage_disequilibrium_PUG.pdf", 
    width = 5, height = 5)
# Plot
PUG <- popsub(data_neutral_bayescan, "PUG")
ia(PUG, sample = 999)
PUG
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_PTF.pdf", 
    width = 5, height = 5)
# Plot
PTF <- popsub(data_neutral_bayescan, "PTF")
ia(PTF, sample = 999)
PTF
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_CUB.pdf", 
    width = 5, height = 5)
# Plot
CUB <- popsub(data_neutral_bayescan, "CUB")
ia(CUB, sample = 999)
CUB
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_CAMP.pdf", 
    width = 5, height = 5)
# Plot
CAMP <- popsub(data_neutral_bayescan, "CAMP")
ia(CAMP, sample = 999)
CAMP
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_ALA.pdf", 
    width = 5, height = 5)
# Plot
ALA <- popsub(data_neutral_bayescan, "ALA")
ia(ALA, sample = 999)
ALA
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_ALAOR.pdf", 
width = 5, height = 5)
# Plot
ALAOR <- popsub(data_neutral_bayescan, "ALAOR")
ia(ALAOR, sample = 999)
ALAOR
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_GIA.pdf", 
    width = 5, height = 5)
# Plot
GIA <- popsub(data_neutral_bayescan, "GIA")
ia(GIA, sample = 999)
GIA
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_GIAOR.pdf", 
width = 5, height = 5)
# Plot
GIAOR <- popsub(data_neutral_bayescan, "GIAOR")
ia(GIAOR, sample = 999)
GIAOR
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_PVEN.pdf", 
    width = 5, height = 5)
# Plot
PVEN <- popsub(data_neutral_bayescan, "PVEN")
ia(PVEN, sample = 999)
PVEN
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_PVENOR.pdf", 
width = 5, height = 5)
# Plot
PVENOR <- popsub(data_neutral_bayescan, "PVENOR")
ia(PVENOR, sample = 999)
PVENOR
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_SAR.pdf", 
    width = 5, height = 5)
# Plot
SAR <- popsub(data_neutral_bayescan, "SAR")
ia(SAR, sample = 999)
SAR
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_SAROR.pdf", 
width = 5, height = 5)
# Plot
SAROR <- popsub(data_neutral_bayescan, "SAROR")
ia(SAROR, sample = 999)
SAROR
# Close pdf
dev.off()

#########################
#########################################
###**** Population structure: GST ****###
#########################################
# Tutorial here: https://grunwaldlab.github.io/Population_Genetics_in_R/Pop_Structure.html

# Gst

# Let's use Hendrick’s standardized GST to assess population structure among these populations 
# (Hedrick, 2005).

# This is where I used the mmod (Modern Methods of Differentiation) R package
Gst_Hedrick(data_neutral_bayescan)	

#########################
#############################################################
###**** Discriminant Analysis of Principal Components ****###
#############################################################

# The DAPC is a multivariate statistical approach that uses populations defined a priori to maximize the variance 
# among populations in the sample by partitioning it into between-population and within-population components

# How can we determine the optimal number of PCs to retain?
#############################################################################
# Using the alpha score and cross validation to find an ideal number of PCs #
#############################################################################
# Use the number of PCs as suggested by the calculated scores

# LET'S START WITH ALPHA SCORE
# The trade-off between power of discrimination and over-fitting can be measured by the ascore, which is simply 
# the difference between the proportion of successful reassignment of the analysis (observed discrimination) and 
# values obtained using random groups (random discrimination). It can be seen as the proportion of successful 
# reassignment corrected for the number of retained PCs. It is implemented by a.score, which relies on repeating 
# the DAPC analysis using randomized groups, and computing a-scores for each group, as well as the average 
# a-score:

C.dapc_alpha <- dapc(data_neutral_bayescan, n.da=100, n.pca=100)
temp_alpha <- a.score(C.dapc_alpha)
names(temp_alpha)
temp_alpha$tab[1:6,1:6]
temp_alpha$pop.score
temp_alpha$mean

# The number of retained PCs can be chosen so as to optimize the a-score; this is achived by optim.a.score:
# Open pdf
pdf("15_Optimizing_alpha_score.pdf", width = 10, height = 10)
# Plot
C.dapc_alpha <- dapc(data_neutral_bayescan, n.da=100, n.pca=100)
temp_alpha <- optim.a.score(C.dapc_alpha)
# Close the pdf file
dev.off()

# We perform the analysis again with the number of PCs as suggested by alpha score, and then map the membership
# probabilities as before:

# Cross-validation

# Cross-validation (carried out with the function xvalDapc) provides an objective optimisation procedure for 
# identifying the ’golidlocks point’ in the trade-off between retaining too few and too many PCs in the model. 
# In cross-validation, the data is divided into two sets: a training set (typically comprising 90% of the data) 
# and a validation set (which contains the remainder (by default, 10%) of the data). With xvalDapc, the 
# validation set is selected by stratified random sampling: this ensures that at least one member of each group 
# or population in the original data is represented in both training and validation sets. DAPC is carried out on 
# the training set with variable numbers of PCs retained, and the degree to which the analysis is able to 
# accurately predict the group membership of excluded individuals (those in the validation set) is used to 
# identify the optimal number of PCs to retain. At each level of PC retention, the sampling and DAPC procedures 
# are repeated n.rep times.

# Open pdf
pdf("15_Cross_validation.pdf", 
    width = 10, height = 10)
# Plot
mat <- as.matrix(data_neutral_bayescan)
grp <- pop(data_neutral_bayescan)
xval <- xvalDapc(mat, grp, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 100, xval.plot = TRUE)
# Close the pdf file
dev.off()

xval[2:6]

##############################################################################
# Now I can make the DAPC using the number of PCs as suggested by alpha
# optimization and cross validation scores

# Open pdf
pdf("15_DAPC.pdf", 
width = 10, height = 5)
# Plot
C.dapc <- dapc(data_neutral_bayescan, n.pca = 13, n.da = 3)
# To confirm that the DAPC is similar to the PCA we can plot the data in a scatter plot.
scatter(C.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
posi.pca = "topleft", cleg = 0.75)
# Close the pdf file
dev.off()

# And finally STRUCTURE-like PLOTS
# Open pdf
pdf("15_STRUCTURE-like_PLOT.pdf", width = 10, height = 4)
# Plot
#compoplot(pnw.dapc,col = function(x) cols, posi = 'top')
compoplot(C.dapc, col = cols, posi = 'top')
# Close the pdf file
dev.off()

# But we want to organise them within populations
# Open pdf
pdf("15_STRUCTURE-like_PLOT_PER_POPULATION.pdf", 
width = 15, height = 4)
# Plot
dapc.results <- as.data.frame(C.dapc$posterior)
dapc.results$pop <- pop(data_neutral_bayescan)
dapc.results$indNames <- rownames(dapc.results)

# ggplot2 has specific requirements for the structure of the data frame format, as it requires each 
# observation in rows, and all different values of these observations in columns 
# (i.e., a long format data frame). To transform the data frame we use the function melt from the package 
# reshape2. melt reorganizes the data frame into the required data frame format, where each membership 
# probability observation for a given population is a row with the sample name, original population, 
# and assigned population as columns.

dapc.results <- melt(dapc.results)

# I get this: Using pop, indNames as id variables
# BUT
# Ignore the prompt for now. Then, we rename the columns into more familiar terms:

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

# ggplot2 will plot the dapc.results data frame we reorganized using melt, using the samples on the X-axis 
# and membership probabilities on the Y-axis. The fill color will indicate the original population assignments. 
# Each facet represents the original population assignment for each sample:

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~Original_Pop, scales = "free_x", space = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p
# Close the pdf file
dev.off()

  # This bar plot shows us a more organized perspective of our data set by contrasting the population membership 
# probability assignments against their original populations.

# Open pdf
pdf("15_DAPC_no_centroids.pdf", width = 10, height = 5)
# Plot
scatter(C.dapc, col=cols, scree.da=T, # show eigenvalues
        scree.pca = T, # show how many PCs were retained
        legend = T, # show populations
        posi.leg = "topleft", # legend position
cell=1.5, cex=2, bg="white",cstar=0)
# Close the pdf file
dev.off()

#####################
#####################
###**** AMOVA ****###
#####################

# Doing AMOVA from this Tutorial: https://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html
# This manual is good in case of doubts: https://cran.r-project.org/web/packages/poppr/poppr.pdf
# Check this tutorial as well: https://github.com/marinegenomicslab/workflow/wiki/AMOVA-in-R

# AMOVA stands for Analysis of MOlecular VAriance and is a method to detect population differentiation 
# utilizing molecular markers (Excoffier, Smouse & Quattro, 1992). This procedure was initially implemented for 
# DNA haplotypes, but applies to any marker system. The implementation of AMOVA in poppr requires two very basic 
# components: (1) A distance matrix derived from the data and (2) a separate table used to partition the data 
# into different stratifications.

# The distance matrix can be calculated using any distance as long as it is euclidean.
###########################
# IMPORTANT, from the poppr manual: https://cran.r-project.org/web/packages/poppr/poppr.pdf
# On  Euclidean  Distances:: With  the ade4 implementation  of  AMOVA  (utilized  by poppr), distances must be 
# Euclidean (due to the nature of the calculations).  Unfortunately, many genetic distance  measures  are  not  
# always  euclidean  and  must  be  corrected  for  before  being  analyzed. Poppr automates this with three 
# methods implemented in ade4: quasieuclid(), lingoes(), andcailliez(). The correction of these distances 
# should not adversely affect the outcome of the analysis.

# Let's try this
data.euclidean_final <- bitwise.dist(data_neutral_bayescan, euclidean = TRUE)
###################################################################################
# Validated! The plots look the same if I refer myself to dist=C.euclidean_no_NAs #
#   and if I just plot without it, meaning that corrections mentioned above have  # 
#   been deployed. BUT STILL NEED TO USE within = FALSE, otherwise df are wrong!  #
###################################################################################
# We will do AMOVA in poppr package with the previously created C.genlight_no_NAs object. Strata has been created already

# In panmictic populations, we would expect to see most of the variance arise from within samples. If we see 
# that the most of the variance occurs among samples within populations or among populations, then there is 
# evidence that we have some sort of population structure. In the case of clonal organisms, this would help 
# support a hypothesis of clonal reproduction.

# Let’s invoke the AMOVA functions with and without clone correction:
# AMOVA
data_final_amova <- poppr.amova(data_neutral_bayescan, # genind object
                                        ~Morphotype/Population, # hierarchy 
                                        within = F) # this is to use populations as lowest level

# without clone correction
# within: logical. When this is set to TRUE (Default), variance within individuals are calculated as 
# well. If this is set to FALSE, the lowest level of the hierarchy will be the sample level (in my 
# case populations).
# IMPORTANT: This example from the tutorial was performed with a data set of dominant (AFLP) markers, 
# but it can also be performed on codominant markers such as SNPs. These provide more information 
# because within sample (individual) variance is also assessed
data_final_amova

data_final_amovacc <- poppr.amova(data_neutral_bayescan, # genind object
#                                  ~Population, # no hierarchy for Orange or Slender
                                ~Morphotype/Population, # Use this hierarchy for the Orange and Yellow dataset                                
                                  within = F,
                                  clonecorrect = TRUE) # with clone correction
# with clone correction
data_final_amovacc

# You can export data in the form of a table with write.table

# Significance testing
set.seed(1999)
  data_final_amova_signif   <- randtest(data_final_amova, nrepet = 999)
data_final_amova_cc_signif <- randtest(data_final_amovacc, nrepet = 999)
# This was done with no correction method though!

data_final_amova_signif 
data_final_amova_cc_signif

# Open pdf
pdf("16_AMOVA_significance.pdf", width = 10, height = 6)
# Plot
plot(data_final_amova_signif)
# Close the pdf file
dev.off()

pdf("16_AMOVA_significance_clone_correction.pdf", width = 10, height = 6)
# Plot
plot(data_final_amova_cc_signif)
# Close the pdf file
dev.off()
# By samples they mean populations!

#####################
#########################
###**** PERMANOVA ****###
#########################

# Looking at populations only
data.adonis_P <- adonis(data.euclidean_final ~ Population,
                                strata,
                                perm=200)
data.adonis_P

# But how can I see which populations differ significantly??
# Trying pairwise.adonis function now from here: https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
# Code by Pedro Martinez Arbizu
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 

# Preparing the data_final file for adonis
data.outliers.1 <- as.matrix(data_neutral_bayescan)
data.outliers.1

# Now add populations and Sample IDs
popinfo <- strata$Population
popinfo

morphotypeinfo <- strata$Morphotype
morphotypeinfo

ind.names <- indNames(data_neutral_bayescan)
ind.names <- as.factor(ind.names)
ind.names # But the order is not the same here as in data.outliers

data.outliers.2 <- cbind.data.frame(popinfo, data.outliers.1)
# Now it almost looks like their gpop file, just need to add ind.names somehow
data.adonis.outliers.1 <- cbind.data.frame(popinfo, morphotypeinfo, data.outliers.1)

# Convert row names into first column
data.outliers <- tibble::rownames_to_column(data.outliers.2, "ind.names")
data.adonis.outliers <- tibble::rownames_to_column(data.adonis.outliers.1, "ind.names")

# similarity euclidean from vegdist and Benjamini Hochberg correction
data_P_pairwise_adonis_eucl_holm <- pairwise.adonis(x=data.adonis.outliers[,4:2641], 
                                                          # Selecting loci
                                                factors=data.adonis.outliers$popinfo, 
                                                # Saying I want pairwise comparisons of populations
                                                sim.function='vegdist',
                                                sim.method='euclidean',
                                                p.adjust.m='holm') # the p.value correction method, 
# one of the methods supported by p.adjust(); default is 'bonferroni'
data_P_pairwise_adonis_eucl_holm

# They used Holm correction for Euclidean distances in the original script
data_M_pairwise_adonis_eucl_holm <- pairwise.adonis(x=data.adonis.outliers[,4:2641], 
                                                            # Selecting loci
                                                            factors=data.adonis.outliers$morphotypeinfo, 
                                                            # Saying I want pairwise comparisons of populations
                                                            sim.function='vegdist',
                                                            sim.method='euclidean',
                                                            p.adjust.m='holm') # the p.value correction method, 
# one of the methods supported by p.adjust(); default is 'bonferroni'
data_M_pairwise_adonis_eucl_holm

# Export now!
write.csv(data_P_pairwise_adonis_eucl_holm, 
          file = "data_P_pairwise_adonis.csv", 
          quote = F)
write.csv(data_M_pairwise_adonis_eucl_holm, 
          file = "data_M_pairwise_adonis.csv", 
          quote = F)

# Betadisper now!
# Different from what I do in adonis!
# This tests for pairwise differences in variance among populations

# Implements Marti Anderson's PERMDISP2 procedure for the analysis of multivariate homogeneity of 
# group dispersions (variances). betadisper is a multivariate analogue of Levene's test for 
# homogeneity of variances
data.mod_P <- with(strata, betadisper(data.euclidean_final, Population))
data.mod_P

anova(data.mod_P)
permutest(data.mod_P)
Tukey_P <- TukeyHSD(data.mod_P)
Tukey_P

pdf("16_Adonis_plots_Populations.pdf", 
width = 10, height = 4)
# Plot
plot(data.mod_P, col=c("khaki", # ALA 
                       "orange", # ALAOR
                       "yellow1", # CAMP
                       "gold2", # CUB
                       "lightgoldenrod", # GIA
                       "orange2", # GIAOR
                       "moccasin", # PTF
                       "yellow3", # PUG
                       "snow3", # PVEN
                       "orangered", # PVENOR
                       "peachpuff3", # SAR
                       "orange3")) # SAROR
boxplot(data.mod_P, col=c("khaki", # ALA 
                          "orange", # ALAOR
                          "yellow1", # CAMP
                          "gold2", # CUB
                          "lightgoldenrod", # GIA
                          "orange2", # GIAOR
                          "moccasin", # PTF
                          "yellow3", # PUG
                          "snow3", # PVEN
                          "orangered", # PVENOR
                          "peachpuff3", # SAR
                          "orange3")) # SAROR
# Close the pdf file
dev.off()

# Now Morphotype
data.mod_M <- with(strata, betadisper(data.euclidean_final, Morphotype))
data.mod_M

anova(data.mod_M)
permutest(data.mod_M)
Tukey_M <- TukeyHSD(data.mod_M)
Tukey_M

pdf("16_Adonis_plot_Morphotype.pdf", 
    width = 10, height = 4)
# Plot
plot(data.mod_M, col=c("yellow", "orange"))
boxplot(data.mod_M, col=c("yellow", "orange"))
# Close the pdf file
dev.off()

#####################
##############################
###**** Fst statistics ****###
##############################

# Try this tutorial: http://adegenet.r-forge.r-project.org/files/Barcelona2015/practical-MVAintro.1.0.pdf
# This one also seems great: https://cran.r-project.org/web/packages/hierfstat/vignettes/hierfstat.html
# Population structure is traditionally measured and tested using F statistics, in particular the Fst, which 
# measures population differentiation (as the proportion of allelic variance occuring between  groups). The  
# package hierfstat implements a wealth of F statistics and related tests, now designed to work natively 
# with genind objects.

# Let's look at the overall structuring first

# Just to have a quick look at the data, let's see the basic stats of our object first: 
# https://rdrr.io/cran/hierfstat/man/basic.stats.html
data_basic_stats <- basic.stats(data_neutral_bayescan[,-1])
data_basic_stats
write.csv(as.matrix(data_basic_stats$perloc), file = "data_basic_stats.csv", quote = F)
# Overall stats on locus level!

# But Ho, Hs and Fis are also given per population! (and are the only vals reported in
# Carreras et al 2019)
# So export them and calculate the means in excel
write.csv(data_basic_stats$Ho, file = "Spreadsheets_04_data_Ho_perPOP.csv", quote = F)
write.csv(data_basic_stats$Fis, file = "Spreadsheets_04_data_Fis_perPOP.csv", quote = F)
write.csv(data_basic_stats$Hs, file = "Spreadsheets_04_data_Hs_perPOP.csv", quote = F)

# We first compute overall F statistics, and then use Goudet’s G statistics to test the existence 
# of population structure
data.fst <- fstat(data_neutral_bayescan)
data.fst

# And if you want to look at Fst only
data.fst_only <- fstat(data_neutral_bayescan, fstonly=TRUE)
data.fst_only

# Goudet's test
data.gtest.M <- gstat.randtest(data_neutral_bayescan,
                                       nsim=999,
                                       method="global", # "global": tests for genetic 
                                       # structuring given 'pop'.
                                       pop=data_neutral_bayescan$strata$Morphotype)

data.gtest.P <- gstat.randtest(data_neutral_bayescan,
                                       nsim=999,
                                       method="global", # "global": tests for genetic 
                                       # structuring given 'pop'.
                                       pop=data_neutral_bayescan$strata$Population)

# Let's make a plot for the G test
pdf("17_Fst_g_test.pdf", 
width = 5, height = 5)
# Plot
plot(data.gtest.M, main = "Effect of Morphotype")
plot(data.gtest.P, main = "Effect of Population")
# Close the pdf file
dev.off()

# Looks like there is some structuring!

############################################
# PAIRWISE COMPARISONS here
# Script from Carreras et al. (2019) 
############################################

# Here is the script: (da ta_h is the data in hierfstat format
data_h <- genind2hierfstat(data_neutral_bayescan)
# WC_fst is the Fst matrix 
WC_fst <- genet.dist(data_h,method = "Nei87")
                     
# p-values for Fst
library(parallel)
mat.obs <- as.matrix(WC_fst) # To get a matrix
NBPERM <- 999 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(data_h[,1]),(data_h[,-1])),method = "Nei87"),mc.cores = 4)
                     
                     library(ade4)
                     allTests <- list()
                     for(i in 1:(nrow(mat.obs)-1)){
                       for(j in 2:nrow(mat.obs)){
                         allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
                       }
                     }
                     
                     pvals <- matrix(data=NA,nrow=12,ncol=12) # 12 because we have 12 localities
                     x=0
                     for (i in 1:(nrow(mat.obs)-1)){
                       for(j in 2:nrow(mat.obs)){
                         x=x+1
                         #     pvals[i,j] <- allTests[[x]][[6]]
                         pvals[i,j] <- allTests[[x]]$pvalue
                       }
                     }
                     
                     pvals[lower.tri(pvals)] <- NA # Putting NAs in the lower triangle
                     fst_pval <- mat.obs
                     fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)] # Putting p vals in the
                     # upper triangle
                     diag(fst_pval) <- NA # Putting NAs in the diagonal
                     write.csv(fst_pval,file="Fst-pvalues_NEI_999boot.csv",quote = F)

# The p vals need to be FDR corrected now! BY correction was used in Carreras et al

# This was done before                     
# p.adjust.M <- p.adjust.methods[p.adjust.methods != "fdr"]
# p.adjust.M # So these are all the methods
                     
# Here I am getting the non corrected p vals from 999 permutations
pvals_for_BY <- fst_pval[upper.tri(fst_pval)]
pvals_for_BY # checking the object to see how they were arranged

# Doing the p adjustment
p.adj <- sapply(p.adjust.M, simplify="matrix", function(meth) p.adjust(pvals_for_BY,meth))
p.adj # in the column none you can see the uncorrected ones, the order is still the same as 
# in pvals_for_BY!
                     
# Let's keep only the BY correction as in Carreras et al 2019
# Take the sixth column with all rows
A_fst_BY_pvals <- p.adj[,6,drop=FALSE] 
# Round to 3 digits
A_fst_BY_pvals <- round(A_fst_BY_pvals, digits = 3)
A_fst_BY_pvals

# Add these into the upper triangle of the matrix
fst_BYpval <- fst_pval
fst_BYpval[upper.tri(fst_BYpval)] <- A_fst_BY_pvals
fst_BYpval
# Putting NAs in the diagonal
diag(fst_BYpval) <- NA
# And export
write.csv(fst_BYpval,file="Fst-pvalues_NEI_999boot_BY_correction.csv",quote = F)

# Plotting a heatmap now from tutorial: https://slowkow.com/notes/pheatmap-tutorial/
# The viridis color palettes: https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html

# Split the columns into 12 groups/populations
col_groups <- substr(colnames(mat.obs), 1, 12) # Because we have 12 populations
table(col_groups)

# Data frame with column annotations
mat_col <- data.frame(Populations = col_groups)
rownames(mat_col) <- colnames(mat.obs)
mat_col

# List with colors for each annotation
mat_colors <- list(Population = brewer.pal(n = nPop(data_neutral_bayescan), name = "Dark2"))
names(mat_colors$Population) <- unique(col_groups)
mat_colors

pdf('19_Fst_values_heatmap.pdf', width=10,height=6)
pheatmap(
  mat   = mat.obs,
  color = viridis(1000),
  # color = brewer.pal(n = 9, name = "Reds"), # If using RcolorBrewer
  border_color  = NA,
  show_colnames = TRUE,
  show_rownames = TRUE,
  # annotation_row = mat_col,
  annotation_col = mat_col,
  annotation_colors = mat_colors,
  drop_levels   = TRUE,
  fontsize  = 14,
  main  = "Fst values",
  cutree_cols   = 3,
  cutree_rows   = 3,
  cluster_rows  = T,
  display_numbers   = T
)
dev.off()

####################################
# NOW MAKING THE PHYLOGENETIC TREE
###########################################
# For my data maybe best to use Provesti! #
###########################################

# Make a dendrogram now
data_neutral_bayescan_dist <- provesti.dist(data_neutral_bayescan)
data_neutral_bayescan_dist

# Open pdf
pdf("10_Provesti_tree.pdf", 
    width = 15, height = 15)
# Plot
# We will use this distance matrix to create a neighbor-joining tree.
theTree <- data_neutral_bayescan_dist %>%
  nj() %>%# calculate neighbor-joining tree
  ladderize() # organize branches by clade
plot(theTree)
add.scale.bar(length = 0.05) # add a scale bar showing 5% difference.
# Close pdf
dev.off()

# A tree is a hypothesis and one way of generating support is to bootstrap loci. 
# This can be achieved with the poppr function aboot.

# Open pdf
pdf("10_Provesti_tree_with_bootstrap.pdf", 
    width = 15, height = 15)
# Plot
set.seed(999)
aboot(data_neutral_bayescan, 
      dist = provesti.dist, 
      sample = 200, 
      tree = "upgma", 
      cutoff = 50, 
      quiet = TRUE)
# Close pdf
dev.off()

# If we wanted to analyze all of the populations against one another, it would be better to create 
# a bootstrapped dendrogram based on a genetic distance. To do this, we will add 3 stratifications 
# to our data: # MorphotypePopulation, Morphotype, and Population. 
# We will then set the population to Morphotype by Population, convert the data to a genpop object and 
# then create a tree using aboot with Nei’s genetic distance.

# Open pdf
pdf("10_Populations_tree_with_bootstrap_Nei_distances.pdf", 
    width = 6, height = 5)
# Plot
# Analysis
set.seed(999)
data_neutral_bayescan %>%
  genind2genpop(pop = ~Morphotype/Population) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = nei.dist)
# Close pdf
dev.off()
# Looks like those from the deep are separating!

# Now making the tree in colors!
# We will build a distance tree to obtain an initial assessment of the population structure of the 
# A. subpinnata samples in the Mediterranean Sea.
tree <- aboot(data_neutral_bayescan, 
              tree = "upgma", # UPGMA
              distance = prevosti.dist, 
              sample = 100, 
              showtree = F,
              cutoff = 50, 
              quiet = T)

# Next, we will color the tips of the tree based on the population of origin of the samples, and draw 
# conclusions from what we observe in the tree. We will save it as a pdf file

# Open pdf
pdf("10_NJ_tree_Provesti_distances.pdf", 
    width = 15, height = 15)
# Plot
# cols <- brewer.pal(n = nPop(data_neutral_bayescan), name = "Dark2")
plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(data_neutral_bayescan)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
#legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
legend('topleft', legend = c("ALA", "ALAOR", "CAMP", "CUB", "GIA", "PTF", "PUG", "PVEN", "PVENOR", "SAR", "GIAOR", "SAROR"), 
       fill = cols, 
       border = T, 
       bty = "n", 
       cex = 2)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")
# Close the pdf file
dev.off()

#########################
    ### THE END ###
#########################

# Now that the whole pipeline has been run, let's save the R environment file
save.image("O_and_Y_neutral.RData")

# Figure 2 was created using the following tutorial: https://github.com/Tom-Jenkins/admixture_pie_chart_map_tutorial
