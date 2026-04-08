# LAB 4 ------
# step1: the essentials ----
# this chunk contains the minimal essential code from this script. Simply uncomment the lines below and run the code.
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(EnsDb.Hsapiens.v86) #replace with your organism-specific database package
targets <- read_tsv("./schistosoma/studyDesign.txt")# read in your study design
path <- file.path("./schistosoma",targets$sample, "abundance.tsv") # set file paths to your mapped data
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = TRUE, #determines whether your data represented at transcript or gene level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

Txi_gene_real <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, #determines whether your data represented at transcript or gene level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)


# step 2: the essentials ----
library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)

sampleLabels <- targets$sample
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "transcriptID")
colnames(log2.cpm.df) <- c("transcriptID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = -1, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=5 #user defined
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "transcriptID")
colnames(log2.cpm.filtered.df) <- c("transcriptID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = -1, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "transcriptID")
colnames(log2.cpm.filtered.norm.df) <- c("transcriptID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = -1, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12) # part one answer

# step 3: the essentials ----
library(tidyverse)
library(DT)
library(gt)
library(plotly)

# factoring all the categorical variables
# sex

sex <- targets$sex
sex <- factor(sex)

strain <- targets$strain
strain <- factor(strain)

drugSensitivity <- targets$drugSensitivity
drugSensitivity <- factor(drugSensitivity)

timpoint <- targets$timpoint
timpoint <- factor(timpoint)

dataType <- targets$dataType
dataType <- factor(dataType)

runAccession <- targets$runAccession
runAccession <- factor(runAccession)

pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = sex) + # task 2 answer: sex
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot)

# female
female <- targets %>% filter(sex == "female" & timpoint == "12hr") %>% pull(sample)
female_controls <- targets %>% filter(sex == "female" & timpoint == "control") %>% pull(sample)

mydata.df <- mutate(log2.cpm.filtered.norm.df,
                    female.AVG = (F12h_LEPZQ_1 + F12h_NMRI_1 + F12h_LE_1 + F12h_NMRI_2 + F12h_LEPZQ_2 +
                                    F12h_LE_2 + F12h_NMRI_3 + F12h_LE_3 + F12h_LEPZQ_3)/9, 
                    female_controls.AVG = (FCtl_LEPZQ_1 + FCtl_NMRI_1 + FCtl_LE_1 + FCtl_NMRI_2 + FCtl_LEPZQ_2 +
                                             FCtl_LE_2 + FCtl_NMRI_3 + FCtl_LE_3 + FCtl_LEPZQ_3)/9,
                    #now make columns comparing each of the averages above that you're interested in
                    LogFC = (female.AVG - female_controls.AVG)) %>% #note that this is the first time you've seen the 'pipe' operator
  mutate_if(is.numeric, round, 2) %>% 
  arrange(desc(LogFC)) %>% 
  select(transcriptID, female.AVG, female_controls.AVG, LogFC) %>% 
  head(10)
  

datatable(mydata.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         #dom = "Blfrtip", 
                         #buttons = c("copy", "csv", "excel"),
                         lengthMenu = c("10", "25", "50", "100")))
