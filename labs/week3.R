# Task 1 ----
# load packages
library(rhdf5) #provides functions for handling hdf5 file formats (kallisto outputs bootstraps in this format)
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(EnsDb.Hsapiens.v86) #replace with your organism-specific database package
library(beepr) #just for fun
library(datapasta) # great for copy/paste data into the R environment
library(biomaRt) # an alternative for annotation

listMarts() #default host is ensembl.org, and most current release of mammalian genomes
#listMarts(host="parasite.wormbase.org") #access to parasite worm genomes
#listMarts(host="protists.ensembl.org") #access to protozoan genomes

#choose the 'mart' you want to work with
myMart <- useMart(biomart="ENSEMBL_MART_ENSEMBL")
#take a look at all available datasets within the selected mart
available.datasets <- listDatasets(myMart)
#now grab the ensembl annotations for fer
fer.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "mpfuro_gene_ensembl")
fer.filters <- listFilters(fer.anno)
fer.attributes <- listAttributes(fer.anno)

Tx.fer <- getBM(attributes=c('ensembl_transcript_id_version',
                             'transcript_start',
                             'transcript_end',
                             'description',
                              'entrezgene_id',
                             'external_gene_name',
                             'pfam'),
                mart = fer.anno)

Tx.fer <- as_tibble(Tx.fer)

# Task 2 ----
antiviral_genes <- c('IFIT2', 'OAS2', 'IRF1', 'IFNAR1', 'MX1')

ferret_promoters <- getSequence(id = antiviral_genes, 
            type = 'external_gene_name',
            seqType = 'gene_flank',
            upstream = 1000,
            mart = fer.anno)

ferret_promoters <- as.tibble(ferret_promoters)


