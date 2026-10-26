##### Lobster sample core analysis - no tank water, no very late embryos ######

setwd("/Volumes/Sarah_SSD/embryo_microbiome")
phyob_novlate <- readRDS("phyob_novlate.rds")

library("phyloseq")
library("ggplot2")
library("microbiome")

phyob_ell <- subset_samples(phyob_novlate, Sample_Type != "Tank")
phyob_ell <- prune_taxa(taxa_sums(phyob_ell)>0, phyob_ell)

pseq.ell.genus <- tax_glom(phyob_ell, taxrank = "Genus") #collapse to genera


# line plot of num_genera X prevalence for absolute counts Early, Late & Larvae
list_prevalences <- seq(0, 1.0, 0.1)
prevalence <- list_prevalences*100

core(phyob_ell_genus, detection= 0, prevalence = 0.0) #check .0 to .9, record in list below (prevalence is strictly greater)
number_asvs <- c(501, 221, 155, 117, 91, 75, 53, 33, 19, 8, 3 )

p <- data.frame(prevalence, number_asvs) #make data frame

plot1 <- ggplot(data = p, mapping = aes(x = prevalence, y = number_asvs)) +
  geom_point()+
  scale_x_discrete("Prevalence %", limits = seq(0,100,10)) +
  scale_y_discrete("Number of ASVs", limits = seq(0,2800, 100)) +
  labs(title = "Lobster Core ASVs by Prevalence") +
  geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95)


###
phyob_early.genus <- tax_glom(phyob_early, taxrank = "Genus") #collapse early embryo ps object to genera
early.core.100 <- core_members(phyob_early.genus, prevalence = .99) #which taxa are on 100% of early embryo samples

phyob_late.genus <- tax_glom(phyob_late, taxrank = "Genus")
late.core.100 <- core_members(phyob_late.genus, prevalence = .99)

phyob_larvae.genus <- tax_glom(phyob_larvae, taxrank = "Genus")
larvae.core.100 <- core_members(phyob_larvae.genus, prevalence = .99)


#you can make core phyloseq objects - helpful for making bar plots, etc. 
pseq.core.early <- core(phyob_early.genus, detection = 0, prevalence = .99)
pseq.core.late <- core(phyob_late.genus, detection = 0, prevalence = .99)
pseq.core.larvae <- core(phyob_larvae.genus, detection = 0, prevalence = .99)




#merge core pseq objects for barplot

pseq.core.el <- merge_phyloseq(pseq.core.early, pseq.core.late)
pseq.core.ell <- merge_phyloseq(pseq.core.el, pseq.core.larvae)

### Core barplot ###
sample_data(pseq.core.ell)$Embryogenesis <- factor(
  sample_data(pseq.core.ell)$Embryogenesis, levels = c('Early', 'Late', 'Larvae') )

sample_data(pseq.core.ell)$Sample_Type <- factor(
  sample_data(pseq.core.ell)$Sample_Type, levels = c('Early Embryo', 'Late Embryo', 'Larvae') )

levels(sample_data(phyob_ell)$Embryogenesis) <- c(Early = "Early Embryo", Late = "Late Embryo", Larvae = "Larvae")


phyob_class <- pseq.core.ell %>%
  #tax_glom(taxrank='Genus') %>% #select/combine at specified taxa level
  transform_sample_counts(function(x) {(x/sum(x))*100} ) %>% #coverts everything into relative abundance
  psmelt() %>% #switch to ggplot usable format
  filter(Abundance > 0) %>% #remove low abundance taxa, this can be changed to fit your samples
  arrange(Class) #sort alphabetically by class


class_bp <- ggplot(phyob_class, aes(x=Sample, y=Abundance, fill= Family)) + 
  facet_grid(~factor(Sample_Type, levels=c("Early Embryo", "Late Embryo", "Larvae")), scales = "free", space = "free")+ #makes panels on your plot
  guides(x = guide_axis(angle = 90)) + #sample names turn 90
  geom_bar(stat='identity', color = "black") + # black around bars
  scale_fill_manual(values = set_colors52) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 2, keyheight = 1, ncol = 2)) +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(size = 11))+
  ylab('Relative Family Abundance (>0%) \n') +
  xlab('Sample ID')

class_bp

library(plotly)
ggplotly(class_bp)

#### Using interactive bar plot, I determined which genera are shared and made a Venn Diagram. I'm sure there is a way to do this in R but I did it by hand ####



set_colors52 <- c("#95ae74",
                  "#5668e6",
                  "#494b93",
                  "#6cd84e",
                  "#5762be",
                  "#c0dc47",
                  "#7d89e9",
                  "#73a72f",
                  "#366bb3",
                  "#59e186",
                  "#6196e2",
                  "#b7dd79",
                  "#646ea3",
                  "#49a65a",
                  "#98a6d7",
                  "#557e2c",
                  "#3f78a6",
                  "#8fdf9f",
                  "#646c89",
                  "#57e5c4",
                  "#3a718b",
                  "#d1e2ab",
                  "#55b0db",
                  "#69773d",
                  "#b6cded",
                  "#4a7f51",
                  "#63d6e0",
                  "#737b55",
                  "#c6e9db",
                  "#3a7883",
                  "#8fd3bb",
                  "#667560",
                  "#a3bfc6",
                  "#3f7a61",
                  "#7794ab",
                  "#40a488",
                  "#52756e",
                  "#8ba691",
                  "#59a1a0",
                  "#4071f3",
                  "#42a359",
                  "#f995aa",
                  "#7a4523",
                  "#c479a3",
                  "#965d39",
                  "#e991ac",
                  "#9c4e34",
                  "#a15276",
                  "#df8c93",
                  "#9e4b4f",
                  "#ce6174",
                  "#c36b65")


