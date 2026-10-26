## Code for merging DADA2 output from embryos, larvae and tank water. Followed by code for Phyloseq. 

# merge sequence tables from embryos, larvae and tank DADA2

library(dada2)
setwd("/Volumes/Sarah_SSD/embryo_microbiome")

#read in sequence tables
seqtab_emb <- readRDS(file = "/Volumes/Sarah_SSD/embryo_microbiome/embryo_fastq/seqtab_emb.nochim.rds") #from dada2, see additional file 1
seqtab_lar_tank <- readRDS(file = "/Volumes/Sarah_SSD/embryo_microbiome/larvae_tank_fastq/seqtab_LT.nochim.rds") #from dada2, see additional file 2

#merge sequence tables
seqtab.nochim <- mergeSequenceTables(seqtab_emb, seqtab_lar_tank)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))

#save asv file
write.csv(seqtab.nochim, file="asv.csv") #save for later
saveRDS(seqtab.nochim, file = "asv.rds") #save for ease of use later

# assigning taxonomy to ASV table
taxa <- assignTaxonomy(seqtab.nochim, "/Users/jskoshak/Desktop/lobster_embryo2/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(taxa, file = "taxa.rds") #save rds for future use

taxa.print <- taxa # Removing sequence row names for display only
rownames(taxa.print) <- NULL
print(taxa.print) #sanity check
write.csv(taxa, file="taxa.csv") # write csv taxa file



#############Phyloseq###############
#continue phyloseq to make figures for manuscript. These are  embryo fastq reads results merged with tank and larvae dada2 taxa and asv files #


taxa <- readRDS("taxa.rds") # load in merged taxa object (if reopening)
seqtab.nochim <- readRDS("asv.rds") # load in merged asv object (if reopening)

# packages
library(phyloseq) 
library(ggplot2)
library(dplyr)
library(vegan)
library(tidyr)
library(tidyverse)

# build phyloseq object with asv and taxa table 
physeq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   tax_table(taxa))
physeq #check that it worked

map <- read.csv("map.csv") #csv to metadata
map <- sample_data(map) #add map to sample_data of ps object
rownames(map) <- map$Sample_ID #Sample_ID in map file should match Sample_ID in asv file
print(map) #sanity check

phyloseq_merge <- merge_phyloseq(physeq, map) #add metadata to ps object
print(phyloseq_merge) #sanity check, ps object should now have sample data 


####Subsetting####
phyob <- phyloseq_merge %>%  
  subset_taxa(
    Order!="Mitochondria" &  # remove taxa with mitochondria or chloroplasts
      Order !='Chloroplast' &
      Class != 'Cyanobacteria') #cyanobacteria sneaked through so removed manually

phyob <- filter_taxa(phyob, function(x) sum(x > 10) > (0*length(x)), prune=T) # I prune taxa that don't appear more than 10 times across all samples. 

#store sequences in refseq slot of phyloseq object but change to ASV1, ASV2 for readability. 
dna <- Biostrings::DNAStringSet(taxa_names(phyob))
names(dna) <- taxa_names(phyob)
phyob <- merge_phyloseq(phyob, dna)
taxa_names(phyob) <- paste0("ASV", seq(ntaxa(phyob)))
phyob # sanity check - reseq slot should have same number of seqs as otu table 


### More Subsetting ### 
# At this point, I make separate ps objects for various figures, statistics, and to check controls #

## Check controls ##
phyob_controls <- subset_samples(phyob, Location =="na") # 2 sets of pos/neg controls
phyob_controls <- prune_taxa(taxa_sums(phyob_controls)>0, phyob_controls) #prune

# separate further by positive and negative 
phyob_positive <- subset_samples(phyob_controls, Embryogenesis =="Positive")
phyob_positive <- prune_taxa(taxa_sums(phyob_positive)>0, phyob_positive) #prune, contains 8 taxa
tax_table(phyob_positive) # Contains all 8 expected in Zymo mock community
taxa_sums(phyob_positive) 

#check negatives
phyob_negative <- subset_samples(phyob_controls, Embryogenesis =="Negative")
phyob_negative <- prune_taxa(taxa_sums(phyob_negative)>0, phyob_negative) #prune, contains 1 taxon
tax_table(phyob_negative) # R. thalassococcus
taxa_sums(phyob_negative) # 5 reads

#which negative was potentially contaminated?
phyob_neg_emb <- subset_samples(phyob_negative, Sample_ID=="Negative")
phyob_neg_emb <- prune_taxa(taxa_sums(phyob_neg_emb)>0, phyob_neg_emb) #no reads

phyob_neg_TL <- subset_samples(phyob_negative, Sample_ID=="NC")
phyob_neg_TL <- prune_taxa(taxa_sums(phyob_neg_TL)>0, phyob_neg_TL) #1 taxon
tax_table(phyob_neg_TL) # R. thalassococcus
taxa_sums(phyob_neg_TL) # 5 reads. We move on. 


#remove controls from phyob 
phyob <- subset_samples(phyob, Sample_Type!= "na") 
phyob <- prune_taxa(taxa_sums(phyob)>0, phyob) #prune unwanted taxa. Pruning is not necessary but I like the cleanliness of it. 

#ps for individual variation
phyob_indiv <- subset_samples(phyob, Embryogenesis == "Very Late") # individual variation samples 
phyob_indiv <- prune_taxa(taxa_sums(phyob_indiv)>0, phyob_indiv) #prune

#Remove individual variation samples
phyob_novlate <- subset_samples(phyob, Embryogenesis!= "Very Late") #remove individual variation samples. Should contain 35 samples (early, late embryos. larvae and tank)
phyob_novlate <- prune_taxa(taxa_sums(phyob_novlate) > 0, phyob_novlate) #prune taxa
phyob_novlate <- subset_samples(phyob_novlate, Sample_ID!= "T4L1E2") #remove outlier
saveRDS(phyob_novlate, file = "phyob_novlate.rds")

#embryo only ps object (no indiv variation)
phyob_emb <- subset_samples(phyob_novlate, Sample_Type!= "Tank" & Sample_Type!= "Larvae") #remove tank and larvae 
phyob_emb <- prune_taxa(taxa_sums(phyob_emb) > 0, phyob_emb) #should contain 24 samples

#early embryos only
phyob_early <- subset_samples(phyob_emb, Embryogenesis!= "Late") #early embryos only
phyob_early <- prune_taxa(taxa_sums(phyob_early) > 0, phyob_early) #prune

#late embryos only
phyob_late <- subset_samples(phyob_emb, Embryogenesis!= "Early") #late embryos only
phyob_late <- prune_taxa(taxa_sums(phyob_late) > 0, phyob_late) #prune

#larvae and tank only 
phyob_lt <- subset_samples(phyob, Sample_Type!= "Late_Embryo" & Sample_Type!= "Early_Embryo" & Sample_Type!= "Very_Late_Embryo") #11 samples
phyob_lt <- prune_taxa(taxa_sums(phyob_lt) >0, phyob_lt) #prune

#larvae only 
phyob_larvae <- subset_samples(phyob_lt, Sample_Type == "Larvae")
phyob_larvae <- prune_taxa(taxa_sums(phyob_larvae) > 0, phyob_larvae)

#tank only
phyob_tank <- subset_samples(phyob_lt, Sample_Type== "Tank")
phyob_tank <- prune_taxa(taxa_sums(phyob_tank) > 0, phyob_tank)

########### Rarefaction Curves ###########
#Rarefaction checks sequencing depth to ensure ASVs were adequately sequenced to capture the diversity of the community ##
devtools::install_github("gauravsk/ranacapa")
BiocManager::install("multtest")
library("ranacapa")
library("multtest")

#I split it up for speed, and for visualizing sequencing depth by sample type. These all look good. 
rf_tank <- ggrare(phyob_tank, step = 50, se = TRUE, label = "Sample_ID")
rf_larvae <- ggrare(phyob_larvae, step = 50, se = TRUE, label = "Sample_ID")
rf_embryo <- ggrare(phyob_emb, step = 50, se = TRUE, label = "Sample_ID")
rf_indiv <- ggrare(phyob_indiv, step = 50, se = TRUE, label = "Sample_ID")


######### Bar Plots ###########

#factor and relevel samples for figures
sample_data(phyob_novlate)$Sample_Type <- factor(
  sample_data(phyob_novlate)$Sample_Type, levels = c('Early_Embryo', 'Late_Embryo', 'Larvae', 'Tank') )

levels(sample_data(phyob_novlate)$Sample_Type) <- c(Early_Embryo = "Early Embryo", Late_Embryo = "Late Embryo", Larvae = "Larvae", Tank = "Tank") #renaming to make figures look nice

#Class barplot
phyob_class <- phyob_novlate %>%
  tax_glom(taxrank='Order') %>% #select/combine at specified taxa level
  transform_sample_counts(function(x) {(x/sum(x))*100} ) %>% #convert to relative abundance
  psmelt() %>% #switch to ggplot usable format
  filter(Abundance > 5) %>% #remove low abundance taxa, this can be changed to fit your samples
  mutate(Order = fct_reorder(Order, Abundance, .fun = sum, .desc = TRUE))


class_bp <- ggplot(phyob_class, aes(x=Sample, y=Abundance, fill=Order)) + 
  facet_grid(~factor(Sample_Type, levels=c("Early Embryo", "Late Embryo", "Larvae", "Tank")), scales = "free", space = "free")+ #makes panels on your plot
  guides(x = guide_axis(angle = 90)) + #sample names turn 90
  geom_bar(stat='identity', color = "black") + # black around bars
  scale_fill_manual(values = set_colors40_2) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 2, keyheight = 1, ncol = 1)) +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(size = 11))+
  ylab('Relative Order Abundance (>5%) \n') +
  xlab('Sample ID')

class_bp

#### Order barplot #####
phyob_order <- phyob_novlate %>%
  tax_glom(taxrank='Order') %>% #select/combine at specified taxa level
  transform_sample_counts(function(x) {(x/sum(x))*100} ) %>% #convert to relative abundance
  psmelt() %>% #switch to ggplot usable format
  filter(Abundance >10) %>% #filter by taxa with relative abundance >1%
  mutate(Order = fct_reorder(Order, Abundance, .fun = sum, .desc = TRUE)) 

order_bp <- ggplot(phyob_order, aes(x=Sample, y=Abundance, fill=Order)) + 
  facet_grid(~factor(Sample_Type, levels=c("Early Embryo", "Late Embryo", "Larvae", "Tank")), scales = "free", space = "free")+ #makes panels on your plot
  guides(x = guide_axis(angle = 90)) + #sample names turn 90
  geom_bar(stat='identity', color = "black") + # black around bars
  scale_fill_manual(values = set_colors40_2) + #colors set below
  guides(fill = guide_legend(reverse = FALSE, keywidth = 2, keyheight = 1, ncol = 1)) + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(size = 11))+
  ylab('Relative Order Abundance (>10%) \n') +
  xlab('Sample ID')

order_bp

### To make it interactive for readability - just for fun ###
library(plotly)
ggplotly(class_bp)


############ Alpha Diversity ############
#### Non - Rarefied Table with Standard Error #####
## There is an additional table showing rarefied alpha diversity, but the trends are the same and there are no statistically significant differences.

#Estimate richness and combine into table
diversity.early <- estimate_richness(phyob_early, split = TRUE, measures = c("Observed", "Shannon", "Simpson", "InvSimpson")) #estimate richness
early.means <- t(colMeans(diversity.early)) #take the mean of each column
early.means <- data.frame(Mean = early.means) #coerce to data frame

early.means$Sample_Type <- c("Early Embryo") # add column with early embryo sample_type

calc.se <- function(x) {
  se <- sd(x) / sqrt(length(x))
  return(se)
} #function to calculate the standard error

standard_error_early <- t((apply(diversity.early, 2, calc.se))) #calculate se 
standard_error_early <- data.frame(SE = standard_error_early) #coerce to dataframe


early_div_calc <- cbind(early.means, standard_error_early) #combine mean and SE by column
order_col <- c("Sample_Type", "Mean.Observed", "SE.Observed", "Mean.Shannon", "SE.Shannon", "Mean.Simpson", "SE.Simpson", "Mean.InvSimpson", "SE.InvSimpson") #reorder columns
early_div_table <- early_div_calc[, order_col] # make table


#estimate late embryo richness and combine into table
diversity.late <- estimate_richness(phyob_late, split = TRUE, measures = c("Observed", "Shannon", "Simpson", "InvSimpson")) #estimate richness
late.means <- t(colMeans(diversity.late)) #take the mean
late.means <- data.frame(Mean = late.means) #make dataframe

late.means$Sample_Type <- c("Late Embryo") #add Late Embryo column

standard_error_late <- t((apply(diversity.late, 2, calc.se))) #calculate standard error
standard_error_late <- data.frame(SE = standard_error_late) #make dataframe

late_div_calc <- cbind(late.means, standard_error_late) #bind mean and se dataframes
div_table_late <- late_div_calc[, order_col] #reorder columns as above

div_tab_all <- rbind(early_div_table, div_table_late) #combine early and late tables

# estimate larvae richness, calculate standard error and combine into table
diversity.larvae <- estimate_richness(phyob_larvae, split = TRUE, measures = c("Observed", "Shannon", "Simpson", "InvSimpson"))
diversity.larvae.means <- t(colMeans(diversity.larvae))
diversity.larvae.means <- data.frame(Mean = diversity.larvae.means)
diversity.larvae.means$Sample_Type <- c("Larvae")

standard_error_larvae <- t((apply(diversity.larvae, 2, calc.se)))
standard_error_larvae <- data.frame(SE = standard_error_larvae)
larvae_div_calc <- cbind(diversity.larvae.means, standard_error_larvae)
div_table_larvae <- larvae_div_calc[, order_col]
div_tab_all <- rbind(div_tab_all, div_table_larvae) #combine embryo and larvae tables

# estimate tank richness, calculate standard error and combine into table
diversity.tank <- estimate_richness(phyob_tank, split = TRUE, measures = c("Observed", "Shannon", "Simpson", "InvSimpson"))
diversity.tank.means <- t(colMeans(diversity.tank))
diversity.tank.means <- data.frame(Mean = diversity.tank.means)
diversity.tank.means$Sample_Type <- c("Tank")

standard_error_tank <- t((apply(diversity.tank, 2, calc.se)))
standard_error_tank <- data.frame(SE = standard_error_tank)

tank_div_calc <- cbind(diversity.tank.means, standard_error_tank)
div_table_tank <- tank_div_calc[, order_col]

div_tab_all <- rbind(div_tab_all, div_table_tank) # add tank
write_csv(div_tab_all, file ="diversity_table_all.csv") #write it to a csv



############# Kruskal-Walis Test ##############
install.packages("FSA")
library(FSA)

set.seed(14)
rich.stats <- estimate_richness(phyob_novlate, split = TRUE, measures = c("Observed", "Shannon", "Simpson", "InvSimpson")) %>%
  rownames_to_column(var= "Sample_ID") %>%
  left_join(as.data.frame(sample_data(phyob_novlate)), by = "Sample_ID") #manipulate data for KW test

##KW by sample type - all samples
observed_kw <- kruskal.test(Observed ~ Sample_Type, data = rich.stats)
shannon_kw <- kruskal.test(Shannon ~ Sample_Type, data = rich.stats)
simpson_kw <- kruskal.test(Simpson ~ Sample_Type, data = rich.stats)
Invsimpson_kw <- kruskal.test(InvSimpson ~ Sample_Type, data = rich.stats)

## KW embryo by Embryogenesis ##
rich.stats.emb <- estimate_richness(phyob_emb, split = TRUE, measures = c("Observed", "Shannon", "Simpson", "InvSimpson")) %>%
  rownames_to_column(var= "Sample_ID") %>%
  left_join(as.data.frame(sample_data(phyob_emb)), by = "Sample_ID")

observed_kw_emb <- kruskal.test(Observed ~ Embryogenesis, data = rich.stats.emb)
shannon_kw_emb <- kruskal.test(Shannon ~ Embryogenesis, data = rich.stats.emb)
simpson_kw_emb <- kruskal.test(Simpson ~ Embryogenesis, data = rich.stats.emb)
Invsimpson_kw_emb <- kruskal.test(InvSimpson ~ Embryogenesis, data = rich.stats.emb)

## KW larvae vs. tank by Sample_type ##
rich.stats.lt <- estimate_richness(phyob_lt, split = TRUE, measures = c("Observed", "Shannon", "Simpson", "InvSimpson")) %>%
  rownames_to_column(var= "Sample_ID") %>%
  left_join(as.data.frame(sample_data(phyob_lt)), by = "Sample_ID")

observed_kw_lt <- kruskal.test(Observed ~ Sample_Type, data = rich.stats.lt)
shannon_kw_lt <- kruskal.test(Shannon ~ Sample_Type, data = rich.stats.lt)
simpson_kw_lt <- kruskal.test(Simpson ~ Sample_Type, data = rich.stats.lt)
Invsimpson_kw_lt <- kruskal.test(InvSimpson ~ Sample_Type, data = rich.stats.lt)



######### Dunn Test ##########
dunnO_shannon <- dunnTest(Shannon ~ Sample_Type,
                          data=rich.stats, #all sample types
                          method="bh") #p values adjusted with Benjamini-Hochberg method





#### Perkins eye index data x Shannon diversity ####
diversity.emb <- estimate_richness(phyob_emb, split = TRUE, measures = c("Shannon"))%>%
  rownames_to_column(var= "Sample_ID") #estimate shannon for all embryos and name Sample_ID column

diversity.emb.early <- estimate_richness(phyob_early, split = TRUE, measures = c("Shannon"))%>%
  rownames_to_column(var= "Sample_ID") 

diversity.emb.late <- estimate_richness(phyob_late, split = TRUE, measures = c("Shannon"))%>%
  rownames_to_column(var= "Sample_ID") 

# merge diversity table with PEI data by sampleID
map.emb.only <- read_csv("map.emb.only.csv") # read in map file with only sample data of embryos, including PEI data

merged <- merge(diversity.emb, map.emb.only, by = "Sample_ID") #merge by Sample_ID to add Shannon Diversity to embryo map file

### plot shannon diversity x PEI for all embryos
library(ggpubr)
install.packages("ggpmisc")
library(ggpmisc)

pei.shan <- ggplot(merged, aes(x=PEI, y=Shannon, color = Embryogenesis)) + 
  geom_point(size = 2) +
  scale_color_manual(values = c("#33A02C", "#E41A1C")) +
  geom_smooth(aes(color = Embryogenesis, fill = Embryogenesis), method = "lm", se = FALSE) +
  xlab('Perkins Eye Index') +
  ylab('Shannon Diversity')+
  stat_fit_glance(
    method = 'lm',
    method.args = list(formula = y ~ x),
    geom = 'text',
    aes(label = paste("p =", signif(..p.value.., digits = 3))),
    label.x = c(169, 472), # Adjust position as needed
    label.y = c(3.28, 3.28),   # Adjust position as needed
    size = 4.5
  ) +
  stat_cor(aes(label = paste(..rr.label..)), # adds R^2 value
           r.accuracy = 0.01,
           label.x = c(150, 450), label.y = c(3.4, 3.4) , size = 4.5) +
  stat_regline_equation(aes(label = ..eq.label..), # adds equation to linear regression
                        label.x = c(150, 450), label.y = c(3.5, 3.5), size = 4.5) +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10))
#title(main = "Early Embryo Shannon Diversity by Perkins Eye index")
pei.shan



### not sure how to separate text for formula and R2 so make 2 plots, early and late, then combine

sample_data(phyob_early)





########## Principal Components Analysis ############
min_lib <- min(sample_sums(phyob_novlate)) #find minimum number of reads in sample with lowest number of reads
min_lib #6612
#subsample
set.seed(3)
subphyob_all <- rarefy_even_depth(phyob_novlate, sample.size = min_lib, verbose=FALSE, replace=TRUE) #randomly subsample your ASVs so that all samples have the same number
subphyob_all #new ps object with rarefied and subsampled ASVs 

#factor and relevel
sample_data(subphyob_all)$Sample_Type <-factor(sample_data(subphyob_all)$Sample_Type, levels=c('Early Embryo', 'Late Embryo', 'Larvae', 'Tank'))

sample_data(subphyob_all)$Embryogenesis <-factor(sample_data(subphyob_all)$Embryogenesis, levels=c('Early Embryo', 'Late Embryo', 'Larvae', 'Tank')) 

sample_data(subphyob_all)$Location <-factor(sample_data(subphyob_all)$Location, levels=c('Maine', 'Massachusetts')) 

sample_data(subphyob_all)$Treatment <-factor(sample_data(subphyob_all)$Treatment, levels=c('Control', 'Warming', 'Acidification', "Multistressor")) 

subphyob_pcoa_all <- ordinate(physeq = subphyob_all, method ='PCoA', distance = 'bray') #ordinate using PCoA and Bray curtis distances

all_PCOA <- plot_ordination(physeq = subphyob_all, ordination = subphyob_pcoa_all, color='Sample_Type', shape = 'Sample_Type') +
  geom_point(aes(color = Sample_Type), alpha=1, size=4) +
  stat_ellipse(aes(group = Sample_Type), type = "norm", level = 0.95) +
  scale_color_manual(values = vals1, name = "") +
  scale_shape_manual(values = vals2, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.85, 0.2),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 

vals1 <- c('Early Embryo' = '#33A02C', 'Late Embryo' = '#E41A1C', 'Larvae' = '#BC80BD', 'Tank' = '#2166AC') #colors
vals2 <- c('Early Embryo' = 16, 'Late Embryo' = 16, 'Larvae' = 17, 'Tank' = 15) #shapes
vals3 <- c('Maine' = 16, 'Massachusetts' = 17)


##### PCoA Tank & Larvae #####
min_lib <- min(sample_sums(phyob_lt)) #find minimum number of reads in sample with lowest number of reads
min_lib #6612
#subsample
set.seed(4)
subhab_TL <- rarefy_even_depth(phyob_lt, sample.size = min_lib, verbose=FALSE, replace=TRUE) #randomly subsample your ASVs so that all samples have the same number
subhab_TL #new ps object with rarefied and subsampled ASVs 

sample_data(subhab_TL)$Sample_Type <-factor(sample_data(subhab_TL)$Sample_Type, levels=c('Larvae', 'Tank'))

hab_pcoa_TL <- ordinate(physeq = subhab_TL, method ='PCoA', distance = 'bray')

TL_PCOA <- plot_ordination(physeq = subhab_TL, ordination = hab_pcoa_TL, color='Sample_Type', shape = 'Sample_Type') +
  geom_point(aes(color = Sample_Type), alpha=1, size=4) +
  stat_ellipse(aes(group = Sample_Type), type = "norm", level = 0.95) +
  #ggtitle("Principal Coordinate Analysis", subtitle = "Larvae, and tank water") +
  scale_color_manual(values = vals3, name = "") +
  scale_shape_manual(values = vals4, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.85, 0.2),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 

vals3 <- c('Larvae' = '#BC80BD', 'Tank' = '#2166AC')
vals4 <- c('Larvae' = 17, 'Tank' = 15)

### Combine two PCoAs into one plot for manuscript ###
library(ggpubr)
ggarrange(all_PCOA, TL_PCOA, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2,
          font.label = list(size = 20))



### PCOA for Location - Embryos only ###

min_lib <- min(sample_sums(phyob_emb)) #find minimum number of reads in sample with lowest number of reads
min_lib #6612
#subsample
set.seed(4)
subphyob_all <- rarefy_even_depth(phyob_emb, sample.size = min_lib, verbose=FALSE, replace=TRUE) #randomly subsample your ASVs so that all samples have the same number
subphyob_all #new ps object with rarefied and subsampled ASVs 

#factor and relevel

sample_data(subphyob_all)$Location <-factor(sample_data(subphyob_all)$Location, levels=c('Maine', 'Massachusetts')) 

sample_data(subphyob_all)$Treatment <-factor(sample_data(subphyob_all)$Treatment, levels=c('Control', 'Warming', 'Acidification', "Multistressor")) 

subphyob_pcoa_all <- ordinate(physeq = subphyob_all, method ='PCoA', distance = 'bray') #ordinate using PCoA and Bray curtis distances

all_PCOA <- plot_ordination(physeq = subphyob_all, ordination = subphyob_pcoa_all, color='Location', shape = 'Sample_Type') +
  geom_point(aes(color = Location), alpha=1, size=4) +
  #stat_ellipse(aes(group = Sample_Type), type = "norm", level = 0.95) +
  scale_color_manual(values = vals4, name = "") +
  scale_shape_manual(values = vals3, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.45, 0.1),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 
all_PCOA
vals1 <- c('Early Embryo' = '#33A02C', 'Late Embryo' = '#E41A1C', 'Larvae' = '#BC80BD', 'Tank' = '#2166AC') #colors
vals2 <- c('Early Embryo' = 16, 'Late Embryo' = 16, 'Larvae' = 17, 'Tank' = 15) #shapes

vals3 <- c('Early Embryo' = 16, 'Late Embryo' = 17)
vals4 <- c('Maine' = '#33A02C', "Massachusetts" = "#E41A1C")





### PCOA for Treatment - Embryos only ###

min_lib <- min(sample_sums(phyob_emb)) #find minimum number of reads in sample with lowest number of reads
min_lib #6612
#subsample
set.seed(4)
subphyob_all <- rarefy_even_depth(phyob_emb, sample.size = min_lib, verbose=FALSE, replace=TRUE) #randomly subsample your ASVs so that all samples have the same number
subphyob_all #new ps object with rarefied and subsampled ASVs 

#factor and relevel

sample_data(subphyob_all)$Location <-factor(sample_data(subphyob_all)$Location, levels=c('Maine', 'Massachusetts')) 

sample_data(subphyob_all)$Treatment <-factor(sample_data(subphyob_all)$Treatment, levels=c('Control', 'Warming', 'Acidification', "Multistressor")) 

subphyob_pcoa_all <- ordinate(physeq = subphyob_all, method ='PCoA', distance = 'bray') #ordinate using PCoA and Bray curtis distances

tx_pcoa <- plot_ordination(physeq = subphyob_all, ordination = subphyob_pcoa_all, color='Treatment', shape = 'Location') +
  geom_point(aes(color = Treatment), alpha=1, size=4) +
  annotate("text", x= -0.35, y=0.25, label ="Early") +
  annotate("text", x= 0.5, y=0.25, label ="Late") +
  stat_ellipse(aes(group = Sample_Type), type = "norm", level = 0.95) +
  scale_color_manual(values = vals5, name = "") +
  scale_shape_manual(values = vals6, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.3, 0.17),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 
tx_pcoa


vals5 <- c('Control' = '#33A02C', 'Warming' = '#E41A1C', 'Acidification' = '#BC80BD', 'Multistressor' = '#2166AC') #colors
vals6 <- c('Maine' = 16, 'Massachusetts' = 17)



########## Principal Components Analysis - Individual Variation ############
min_lib <- min(sample_sums(phyob_indiv)) #find minimum number of reads in sample with lowest number of reads
min_lib #140434
#subsample
set.seed(4)
subphyob_all <- rarefy_even_depth(phyob_indiv, sample.size = min_lib, verbose=FALSE, replace=TRUE) #randomly subsample your ASVs so that all samples have the same number
subphyob_all #new ps object with rarefied and subsampled ASVs 

#factor and relevel

sample_data(subphyob_all)$Treatment <-factor(sample_data(subphyob_all)$Lobster_ID, levels=c('T6_L1', "T7_L1")) 

subphyob_pcoa_all <- ordinate(physeq = subphyob_all, method ='PCoA', distance = 'bray') #ordinate using PCoA and Bray curtis distances

all_PCOA <- plot_ordination(
  physeq = subphyob_all,
  ordination = subphyob_pcoa_all,
  color='Lobster_ID',
  shape='Lobster_ID'
) +
  geom_point(aes(color = Lobster_ID), alpha=1, size=4) +
  stat_ellipse(aes(group = Lobster_ID), type = "norm", level = 0.95) +
  scale_color_manual(values = vals1, labels = c("T6_L1" = "Individual 1","T7_L1" = "Individual 2"),name = "") +
  scale_shape_manual(values = vals2, labels = c("T6_L1" = "Individual 1","T7_L1" = "Individual 2"),name = "") +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.85, 0.2),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")
  )
vals1 <- c('T6_L1' = '#33A02C', 'T7_L1' = '#E41A1C') #colors
vals2 <- c('T6_L1' = 16, 'T7_L1' = 17) #shapes







######### PERMANOVA ############
set.seed(12)
library(vegan)

## Remove unpaired lobsters (remove samples with no early - late pairing)
phyob_emb_trimmed <- subset_samples(phyob_emb, Lobster_ID!= "T11_L1" & Lobster_ID!= "T11_L2" & Lobster_ID!= "T3_L1" & Lobster_ID!= "T4_L1" & Lobster_ID!= "T6_L2") 

bray_emb_trimmed <- phyloseq::distance(phyob_emb_trimmed, method = "bray") #calculate distance on paired embryos
sampledf7 <- data.frame(sample_data(phyob_emb_trimmed))

#are communities different on embryos due to embryo stage (early v late)?
adon_emb_trimmed <- adonis2(bray_emb_trimmed ~ Embryogenesis, data = sampledf7, strata = sampledf7$Lobster_ID) 
#pval = 0.007812 - difference between early and late is significant, but test dispersion

dispersion_emb <- betadisper(bray_emb_trimmed, sampledf7$Embryogenesis) #calculates the average distance of samples to their group centroid in multivariate space.
anova_dispersion_emb <- anova(dispersion_emb) #test for differences in dispersion. If this ANOVA is significant (p < 0.05), it means at least one group has a different level of dispersion which potentially violates PERMANOVA’s assumption. Pvalue is 0.006403 - dispersion between the groups are different, but this does not necessarily negate permanova results. 
write.csv(anova_dispersion_emb, file = "emb_disp.csv")
plot(dispersion_emb) #dispersions are different but centroids are separate






#### Differences in embryos from Treatment? ######
adon_emb_tx <- adonis2(bray_emb_trimmed ~ Treatment, data = sampledf7) #0.89
dispersion_emb_tx <- betadisper(bray_emb_trimmed, sampledf7$Treatment)
disp_emb_tx <- anova(dispersion_emb_tx) #p=0.6775. dispersions are not different
write.csv(disp_emb_tx, file = "emb_disp_tx.csv")
emb_tx_tukey <- TukeyHSD(dispersion_emb_tx)# all pvals >0.6
plot(dispersion_emb_tx)

## Early embryos treatment ##
bray_early <- phyloseq::distance(phyob_early, method = "bray")
meta_early <- data.frame(sample_data(phyob_early))
adon_emb_tx_early <- adonis2(bray_early ~ Treatment, data = meta_early) #0.875
dispersion_emb_early_tx <- betadisper(bray_early, meta_early$Treatment)
disp_emb_tx_early <- anova(dispersion_emb_early_tx) #p=0.1451. dispersions are not different
emb_tx_tukey_early <- TukeyHSD(dispersion_emb_early_tx)# all pvals >0.6
plot(dispersion_emb_early_tx)

## Late Embryos Treatment ##
bray_late <- phyloseq::distance(phyob_late, method = "bray")
meta_late <- data.frame(sample_data(phyob_late))
adon_emb_tx_late <- adonis2(bray_late ~ Treatment, data = meta_late) #0.433
dispersion_emb_late_tx <- betadisper(bray_late, meta_late$Treatment)
disp_emb_tx_late <- anova(dispersion_emb_late_tx) #p=0.4063. dispersions are not different
emb_tx_tukey_late <- TukeyHSD(dispersion_emb_late_tx)# all pvals >0.45
plot(dispersion_emb_late_tx)




### Differences in embryo community from Location? ####
adon_emb_loc <- adonis2(bray_emb_trimmed ~ Location, data = sampledf7) #p=0.776
dispersion_emb_loc <- betadisper(bray_emb_trimmed, sampledf7$Location)
disp_emb_loc <- anova(dispersion_emb_loc) #p=0.9667
write.csv(disp_emb_loc, file = "emb_disp_loc.csv")
plot(dispersion_emb_loc)

## Early embryo location ##
adon_emb_loc_early <- adonis2(bray_early ~ Location, data = meta_early) #0.093
dispersion_emb_early_loc <- betadisper(bray_early, meta_early$Location)
disp_emb_loc_early <- anova(dispersion_emb_early_loc) #p=0.8373. dispersions are not different
plot(dispersion_emb_early_loc)

## Late embryo location ##
adon_emb_loc_late <- adonis2(bray_late ~ Location, data = meta_late) #0.815
dispersion_emb_late_loc <- betadisper(bray_late, meta_late$Location)
disp_emb_loc_late <- anova(dispersion_emb_late_loc) #p=0.1253. dispersions are not different
plot(dispersion_emb_late_loc)



#### PERMANOVA for Larvae ####
#remove control larvae bc you cannot run permanova on a group with only one sample
phyob_larvae_trimmed <- subset_samples(phyob_larvae, Sample_ID!= "T11L1L")
bray_larvae <- phyloseq::distance(phyob_larvae_trimmed, method = "bray") #calculate distance on paired embryos
sampled_lar <- data.frame(sample_data(phyob_larvae_trimmed))

#are communities different on larvae due treatment?
adon_lar <- adonis2(bray_larvae ~ Treatment, data = sampled_lar) 
#pval = 0.7333 - no difference in community by treatment. run betadisper to check dispersions

dispersion_lar <- betadisper(bray_larvae, sampled_lar$Treatment) #calculates the average distance of samples to their group centroid in multivariate space — the multivariate analogue of variance.
anova_dispersion_lar <- anova(dispersion_lar) # dispersion cannot be tested, "“PERMANOVA indicated a significant effect of Treatment on larval microbiome composition. However, betadisper could not be reliably computed due to zero within-group variance (residual sum of squares = 0), so homogeneity of dispersions could not be assessed.”

plot(dispersion_lar) 









### Variable interactions?? ###
loc_tx_emb <- adonis2(bray_emb_trimmed ~ Location * Treatment, data = sampledf7, strata = "Lobster_ID")


loc_tx_early <- adonis2(bray_early ~ Location * Treatment, data = meta_early)
loc_tx_late <- adonis2(bray_late ~ Location * Treatment, data = meta_late)



###### Use all samples for permanova - paired early and late embryos - difference by sample type #####
phyob_pairs <- merge_phyloseq(phyob_emb_trimmed, phyob_lt)
sample_data(phyob_pairs)$PairID <- c("T3_L2", "T3_L2", "T4_L2", "T4_L2", "T5_L1", "T5_L2", "T5_L1", "T5_L2", "T6_L1", "T6_L1", "T7_L1", "T7_L1", "T7_L2", "T7_L2", "T8_L1", "T8_L1", "T8_L2", "T8_L2", "na", "na", "na", "na", "na", "na", "na", "na", "na", "na", "na") # add pairs to sample_data for stratification
phyob_bray_pair <- phyloseq::distance(phyob_pairs, method = "bray")
sampledf8 <- data.frame(sample_data(phyob_pairs))

sampledf8$Strata <- ifelse(
  sampledf8$Embryogenesis %in% c("Early", "Late"),
  sampledf8$PairID,         # keep the pair ID for paired samples
  sampledf8$Sample_ID        # unpaired samples get unique strata
)

adonis2(phyob_bray_pair ~ Sample_Type * Sample_Type, data = sampledf8, permutations = 999, strata = sampledf8$Strata)

dispersion_paired <- betadisper(phyob_bray_pair, sampledf8$Sample_Type)
anova(dispersion_paired)
TukeyHSD(dispersion_paired)
plot(dispersion_paired)


#### Colors Used #####

set_colors40_2 <- c("#1B9E77",
                    "#D95F02",
                    "#7570B3",
                    "#E7298A",
                    "#66A61E",
                    "#E6AB02",
                    "#A6761D",
                    "#666666",
                    "#8DD3C7",
                    "#FFFFB3",
                    "#BEBADA",
                    "#FB8072",
                    "#80B1D3",
                    "#FDB462",
                    "#B3DE69",
                    "#FCCDE5",
                    "#BC80BD",
                    "#CCEBC5",
                    "#FFED6F",
                    "#8C510A",
                    "#01665E",
                    "#5AB4AC",
                    "#F6E8C3",
                    "#D8B365",
                    "#F5F5F5",
                    "#B2182B",
                    "#2166AC",
                    "#4DAF4A",
                    "#984EA3",
                    "#FF7F00",
                    "#E41A1C",
                    "#FFFF33",
                    "#A6CEE3",
                    "#1F78B4",
                    "#B2DF8A",
                    "#33A02C",
                    "#FB9A99",
                    "#E31A1C",
                    "#CAB2D6",
                    "#C7EAE5",
                    "#FFFF10")


