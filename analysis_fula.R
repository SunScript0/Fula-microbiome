library(vegan)
library(RColorBrewer)
library(Maaslin2)
library(MASS)
library(gridExtra)
library(ape)
library(picante)
library(ggpubr)
library(GGally)
library(plotly)
library(rbiom)
library(sjmisc)
library(reshape2)
library(patchwork)
library(cowplot)
library(scales)
library(tidyverse)

setwd("/scratch/fmorandi/fula_test/Fula-microbiome")
source("./utils.R")

datapath = "./03_processed_data"
outpath = "./04_analysis_output/"

set.seed(1337)

##### SETUP OUTPUT #####

dat = Sys.time()
dat = gsub(" ", "_", dat)
dat = gsub(":", "-", dat)

# Create output folders
outpath = paste0(outpath, "/outputs_", dat)
dir.create(outpath)

# Redirect stdout to file
sink(file = paste0(outpath, "/main_script_res_", dat, ".txt"), type = "output")

# List to store all panels to be assembled at the end
plots = list()

##### GLOBAL PLOTTING SETTINGS #####

theme_set(theme_light(base_size = 8))
update_geom_defaults("text", list(size = 8*0.35))
pval_size = 8*0.35

cols = list()
cols$adults = c("#1F78B4", "#FF7F00", "#33A02C", "#E31A1C")
cols$infants = c("#5b9fb5", "#e6b950", "#8eb36f", "#d96e6c")
cols$all = c(cols$adults, cols$infants)[c(1,3,5,7,2,4,6,8)]
cols$times = hue_pal()(4)

##### READ DATA #####

load(file=paste0(datapath, "/fula_data.Rdata"))

##### SAMPLE COUNTS #####

cat("Pairs:", length(unique(meta$Index)), "\n")
cat("Total samples:\n")
meta %>%
  group_by(Group3) %>%
  summarize(length(Group3))
cat("Passing QC:\n")
meta %>%
  filter(PassesQC) %>%
  group_by(Group3) %>%
  summarize(length(Group3))

cat("T1 mother ages:", "\n")
meta %>%
  filter(Group4 == "Adult.1") %>%
  summarize(min = min(AgeM), mean = mean(AgeM), max = max(AgeM))
meta %>%
  filter(Group4 == "Adult.1") %>%
  t.test(AgeM ~ Location, .)

cat("T1 child ages:", "\n")
meta %>%
  filter(Group4 == "Child.1") %>%
  summarize(min = min(AgeC), mean = mean(AgeC), max = max(AgeC))
meta %>%
  filter(Group4 == "Child.1") %>%
  t.test(AgeC ~ Location, .)

##### SUMMARY TABLE (T1) #####

summary_mothers = meta %>%
  filter(AdultChild == "Adult") %>%
  group_by(Group3) %>%
  summarize(count = n(), age_mean = mean(AgeM), age_sd = sd(AgeM),
            height_mean = mean(Height, na.rm=T), height_sd = sd(Height, na.rm=T),
            weight_mean = mean(Weight, na.rm=T), weight_sd = sd(Weight, na.rm=T),
            BMI_mean = mean(BMI, na.rm=T), BMI_sd = sd(BMI, na.rm=T),
            WHR_mean = mean(WHR, na.rm=T), WHR_sd = sd(WHR, na.rm=T),
            entColi_pos = sum(EntColi, na.rm=T), entColi_perc = 100*sum(EntColi, na.rm=T) / n()) %>%
  mutate_at(c("age_mean", "age_sd", "height_mean", "height_sd", 
              "weight_mean", "weight_sd", "BMI_mean", "BMI_sd", "entColi_perc"), round, 1) %>%
  mutate_at(c("WHR_mean", "WHR_sd"), round, 2)

summary_infants = meta %>%
  filter(AdultChild == "Child") %>%
  group_by(Group3) %>%
  summarize(count = n(), age_mean = mean(AgeC), age_sd = sd(AgeC),
            male = sum(Gender == "M", na.rm=T), female = sum(Gender == "F", na.rm=T),
            height_mean = mean(Height, na.rm=T), height_sd = sd(Height, na.rm=T),
            weight_mean = mean(Weight, na.rm=T), weight_sd = sd(Weight, na.rm=T),
            BMI_mean = mean(BMI, na.rm=T), BMI_sd = sd(BMI, na.rm=T),
            headCirc_mean = mean(HeadCirc, na.rm=T), headCirc_sd = sd(HeadCirc, na.rm=T),
            cSec_n = sum(cSection, na.rm=T), cSec_perc = 100*sum(cSection, na.rm=T) / n(),
            sinceDiv_mean = mean(SinceDivers., na.rm=T), sinceDiv_sd = sd(SinceDivers., na.rm=T),
            weaned_n = sum(Weaned, na.rm=T), weaned_perc = 100*sum(Weaned, na.rm=T) / n()) %>%
  mutate_if(is.numeric, round, 1)

write.table(summary_mothers, paste0(outpath, "/summary_mothers.tsv"), sep = "\t", quote=F)
write.table(summary_infants, paste0(outpath, "/summary_infants.tsv"), sep = "\t", quote=F)

##### CHOSE RAREFACTION LEVEL #####

plots$qc = list()

plots$qc$coverage = meta %>%
  mutate(label = str_replace_all(Group3, "\\.", "\n")) %>%
  ggplot(., aes(x = label, y = FilteredPairs, color=label)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2), size=1)+
  geom_hline(yintercept = 14000, color = "red") +
  scale_color_manual(values = cols$all)+
  guides(color="none")+
  theme(axis.title.x = element_blank())+
  labs(y="Depth")

# Genus
rarec = rarecurve(ab$gen, step=100, label=F, tidy=T)
plots$qc$raref_gen = 
  ggplot(rarec, aes(x=Sample, y=Species, group=Site))+
  geom_line(alpha=0.5)+
  geom_vline(xintercept=14000, color = "red")+
  labs(x="Depth", y="Genera")

# Species
rarec = rarecurve(ab$spe, step=100, label=F, tidy=T)
plots$qc$raref_spe = 
  ggplot(rarec, aes(x=Sample, y=Species, group=Site))+
  geom_line(alpha=0.5)+
  geom_vline(xintercept=14000, color = "red")+
  labs(x="Depth", y="Species")

# ASVs
rarec = rarecurve(ab$asv, step=100, label=F, tidy=T)
plots$qc$raref_asv = 
  ggplot(rarec, aes(x=Sample, y=Species, group=Site))+
  geom_line(alpha=0.5)+
  geom_vline(xintercept=14000, color = "red")+
  labs(x="Depth", y="ASVs")

##### SAMPLE COVERAGE #####

cat("Raw pairs:\n")
meta %>%
  summarize(min=min(RawPairs), max=max(RawPairs), mean=mean(RawPairs), sd=sd(RawPairs),)
cat("Filtered pairs:\n")
meta %>%
  summarize(min=min(FilteredPairs), max=max(FilteredPairs), mean=mean(FilteredPairs), sd=sd(FilteredPairs),)

##### REMOVE QC FAILED SAMPLES #####

meta = subset(meta, PassesQC)
for (level in names(ab)) {
  ab[[level]] = subset(ab[[level]], rownames(ab[[level]]) %in% rownames(meta))
}

##### NORMALIZE ABUNDANCE #####

rarefaction_level = 14000
ab$gen_rare = rrarefy(ab$gen, rarefaction_level)
ab$spe_rare = rrarefy(ab$spe, rarefaction_level)
ab$asv_rare = rrarefy(ab$asv, rarefaction_level)

ab$phy_norm = (1e6 * ab$phy) / rowSums(ab$phy)

##### PHYLUM BARPLOT #####

plots$phy_barplots = list()

# Combine phyla with mean abundance < 1% as "Other"
rare_phyla = which((100 * sapply(ab$phy_norm, mean) / 1e6) < 1)
rare_phyla_ab = rowSums(ab$phy_norm[, rare_phyla])
ab_phy_sorted = ab$phy_norm[, -rare_phyla] # Will be sorted later
colnames(ab_phy_sorted) = str_replace_all(taxa$phy[colnames(ab_phy_sorted), "Phylum"], "p__", "")
ab_phy_sorted$Other = rare_phyla_ab
phyla = colnames(ab_phy_sorted)
ab_phy_sorted = merge(ab_phy_sorted, meta, by="row.names", suffixes = c("", ".y"))

cols_phyla = c("#8DD3C7", "#FFFFB3", "#FB8072", "#CCCCCC", "#72B2FB")

plots$phy_barplots$p1 = ab_phy_sorted %>%
  filter(AdultChild == "Adult" & TimePoint == 1) %>%
  group_by(Location) %>%
  mutate(Order = rank(-Firmicutes)) %>%
  mutate_at("Order", as.factor) %>%
  pivot_longer(., phyla) %>%
  ggplot(., aes(x = Order, y = value, fill = name)) +
    geom_bar(position = 'fill', stat = 'identity', width=1) +
    scale_fill_manual(values = cols_phyla) +
    theme_classic() +
    theme(axis.text.x=element_blank(), 
          text = element_text(size = 8)) +
    facet_grid(~ Location) +
    labs(x = "T1 mothers", y="Relative abundance", fill = "Phylum")
  
plots$phy_barplots$p2 = ab_phy_sorted %>%
  filter(AdultChild == "Child" & TimePoint == 1) %>%
  group_by(Location) %>%
  mutate(Order = rank(-Firmicutes)) %>%
  mutate_at("Order", as.factor) %>%
  pivot_longer(., phyla) %>%
  ggplot(., aes(x = Order, y = value, fill = name)) +
    geom_bar(position = 'fill', stat = 'identity', width=1) +
    scale_fill_manual(values = cols_phyla) +
    theme_classic() + 
    theme(axis.text.x=element_blank(), 
          text = element_text(size = 8)) +
    facet_grid(~ Location) +
    labs(x = "T1 infants", y="Relative abundance", fill = "Phylum")

plots$phy_barplots$p3 = ab_phy_sorted %>%
  filter(AdultChild == "Adult" & TimePoint == 2) %>%
  group_by(Location) %>%
  mutate(Order = rank(-Firmicutes)) %>%
  mutate_at("Order", as.factor) %>%
  pivot_longer(., phyla) %>%
  ggplot(., aes(x = Order, y = value, fill = name)) +
    geom_bar(position = 'fill', stat = 'identity', width=1) +
    scale_fill_manual(values = cols_phyla) +
    theme_classic() + 
    theme(axis.text.x=element_blank(), 
          text = element_text(size = 8)) +
    facet_grid(~ Location) +
    labs(x = "T2 mothers", y="Relative abundance", fill = "Phylum")

plots$phy_barplots$p4 = ab_phy_sorted %>%
  filter(AdultChild == "Child" & TimePoint == 2) %>%
  group_by(Location) %>%
  mutate(Order = rank(-Firmicutes)) %>%
  mutate_at("Order", as.factor) %>%
  pivot_longer(., phyla) %>%
  ggplot(., aes(x = Order, y = value, fill = name)) +
    geom_bar(position = 'fill', stat = 'identity', width=1) +
    scale_fill_manual(values = cols_phyla) +
    theme_classic() + 
    theme(axis.text.x=element_blank(), 
          text = element_text(size = 8)) +
    facet_grid(~ Location) +
    labs(x = "T2 infants", y="Relative abundance", fill = "Phylum")

plots$phy_barplots$combined = ggarrange(plotlist = plots$phy_barplots[c("p1", "p2", "p3", "p4")], 
                                        ncol=2, nrow = 2,  common.legend = TRUE, legend="right")

phylum_summary = ab_phy_sorted %>%
  dplyr::select(2:6, "Group3") %>% 
  group_by(Group3) %>%
  summarize(across(everything(), list(mean=mean, std=sd)))
phylum_summary = as.data.frame(phylum_summary)
rownames(phylum_summary) = phylum_summary$Group3
phylum_summary = dplyr::select(phylum_summary, -Group3)
phylum_summary = as.data.frame(t(phylum_summary))
phylum_summary = phylum_summary / 1e4
phylum_summary

phylum_summary_pooled = ab_phy_sorted %>%
  dplyr::select(2:6, "Group4") %>% 
  group_by(Group4) %>%
  summarize(across(everything(), list(mean=mean, std=sd)))
phylum_summary_pooled = as.data.frame(phylum_summary_pooled)
rownames(phylum_summary_pooled) = phylum_summary_pooled$Group4
phylum_summary_pooled = dplyr::select(phylum_summary_pooled, -Group4)
phylum_summary_pooled = as.data.frame(t(phylum_summary_pooled))
phylum_summary_pooled = phylum_summary_pooled / 1e4
phylum_summary_pooled

phy_comp_rVu = my_compare_means(
  data = ab_phy_sorted, 
  dvars = c("Firmicutes", "Bacteroidota", "Proteobacteria", "Actinobacteriota"), 
  indvar = "Location", 
  grouping = c("AdultChild", "TimePoint"))
print(phy_comp_rVu)

phy_comp_1V2 = my_compare_means(
  data = ab_phy_sorted, 
  dvars = c("Firmicutes", "Bacteroidota", "Proteobacteria", "Actinobacteriota"), 
  indvar = "TimePoint", 
  grouping = c("AdultChild", "Location"),
  pairing = "Index")
print(phy_comp_1V2)

rm(ab_phy_sorted, phy_comp_1V2, phy_comp_rVu, phylum_summary, phylum_summary_pooled, 
   cols_phyla, phyla, rare_phyla, rare_phyla_ab)

##### REMOVE HIGH PROTEOBACTERIA SAMPLES #####

cat("How many infants are high on proteo?\n")
sum(meta$Proteobacteria > 50 & meta$AdultChild == "Child")

cat("Is high proteo more common in rural?\n")
meta %>%
  filter(AdultChild == "Child") %>%
  mutate(HighProteo = Proteobacteria > 50) %>%
  select(Location, HighProteo) %>%
  table() %>%
  fisher.test()

cat("Are high proteo younger?\n")
meta %>%
  filter(AdultChild == "Child") %>%
  filter(TimePoint == 1) %>%
  mutate(HighProteo = Proteobacteria > 50) %>%
  wilcox.test(AgeC ~ HighProteo, .)

mds = list()
mds$all = list()

mds$all$bray_with_proteo = ab$asv_rare %>%
  vegdist() %>%
  as.matrix() %>%
  cmdscale(., eig = TRUE, k = 2)
mds$all$unifracU_with_proteo = ab$asv_rare %>%
  t() %>%
  rbiom::unifrac(., weighted=F, tree=tree) %>%
  as.matrix() %>%
  cmdscale(., eig = TRUE, k = 2)
mds$all$unifracW_with_proteo = ab$asv_rare %>%
  t() %>%
  rbiom::unifrac(., weighted=T, tree=tree) %>%
  as.matrix() %>%
  cmdscale(., eig = TRUE, k = 2)

mds_points = list()
mds_points$all = list()
for (type in names(mds$all)) {
  mds_points$all[[type]] = as.data.frame(mds$all[[type]]$points)
  colnames(mds_points$all[[type]]) = paste0("PC", 1:ncol(mds_points$all[[type]]))
  mds_points$all[[type]] = merge(mds_points$all[[type]], meta, by = 0)
  rownames(mds_points$all[[type]]) = mds_points$all[[type]]$Row.names
  mds_points$all[[type]] = mds_points$all[[type]][,-c(1)]
}

plots$beta = list()
point_size = 1.5

plots$beta$all_pcoa = ggplot(mds_points$all$unifracW, aes(x=PC1, y=PC2, color=Group4, shape=Proteobacteria >50 & AdultChild == "Child"))+
  geom_point(size=point_size, alpha=0.5)+
  scale_shape_manual("Proteobacteria > 50% & Infant", values=c(19, 8))+
  scale_color_discrete("Group", labels = c("Mothers T1", "Infants T1", "Mothers T2", "Infants T2"))

meta = subset(meta, Proteobacteria < 50 | AdultChild == "Adult")
for (level in names(ab)) {
  ab[[level]] = subset(ab[[level]], rownames(ab[[level]]) %in% rownames(meta))
}

##### ALPHA DIVERSITY: compare metrics #####

plots$alpha = list()

alpha_shan = diversity(ab$asv)
alpha_shan_rare = diversity(ab$asv_rare)
alpha_chao = estimateR(ab$asv)[2,]
alpha_chao_rare = estimateR(ab$asv_rare)[2,]
alpha_ace = estimateR(ab$asv)[4,]
alpha_ace_rare = estimateR(ab$asv_rare)[4,]
alpha_pd = pd(ab$asv, tree)[,1]
alpha_pd_rare = pd(ab$asv_rare, tree)[,1]

alpha_div_comparison = data.frame("shannon" = alpha_shan, "shannon_rare" = alpha_shan_rare,
                                  "chao" = alpha_chao, "chao_rare" = alpha_chao_rare,
                                  "ace" = alpha_ace, "ace_rare" = alpha_ace_rare,
                                  "pd" = alpha_pd, "pd_rare" = alpha_pd_rare)

plots$alpha$metric_comp = ggpairs(alpha_div_comparison)

# Chao1 and ACE are basically the same, will drop ACE
# Rarefaction or not does not make much difference, will go for rarefied
alpha = data.frame("Shannon" = alpha_shan_rare,
                   "Chao1" = alpha_chao_rare,
                   "PD" = alpha_pd_rare) # Use rarefied
alpha = merge(alpha, meta, by="row.names")

##### ALPHA DIVERSITY: compare location #####

point_size = 1

# Reorder to compare locations easily
alpha$Group3 = factor(alpha$Group3, levels=c(
  "Rural.Adult.1", "Urban.Adult.1",
  "Rural.Adult.2", "Urban.Adult.2",
  "Rural.Child.1", "Urban.Child.1",
  "Rural.Child.2", "Urban.Child.2")
)

alpha_comp_location = my_compare_means(
  data = alpha, 
  dvars = c("Shannon", "Chao1", "PD"), 
  indvar = "Location", 
  grouping = c("AdultChild", "TimePoint"))

alpha_comp_location$group1 = rep(c("Rural.Adult.1", "Rural.Adult.2",
                                   "Rural.Child.1", "Rural.Child.2"), each = 3)
alpha_comp_location$group2 = rep(c("Urban.Adult.1", "Urban.Adult.2",
                                   "Urban.Child.1", "Urban.Child.2"), each = 3)
alpha_comp_location$p.format = formatC(alpha_comp_location$p, 2)

# Chao1 adults
plots$alpha$adults_chao = alpha %>%
  filter(AdultChild == "Adult") %>%
  ggplot(., aes(x=Group3, y=Chao1, color=Group3))+
    scale_color_manual(values = cols$adults) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size=point_size)+
    stat_pvalue_manual(subset(alpha_comp_location, Dependent == "Chao1" & AdultChild == "Adult"),
                              y.position = 480, label = "p.format", size = pval_size) +
    ylim(0, 520) +
    theme(legend.position="none")+
    labs(x="") +
    scale_x_discrete(labels=c(
      "Rural\nmothers\nT1", "Urban\nmothers\nT1",
      "Rural\nmothers\nT2", "Urban\nmothers\nT2"))

# Chao1 infants
plots$alpha$infants_chao = alpha %>%
  filter(AdultChild == "Child") %>%
  ggplot(., aes(x=Group3, y=Chao1, color=Group3))+
    scale_color_manual(values = cols$infants) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size=point_size)+
    stat_pvalue_manual(subset(alpha_comp_location, Dependent == "Chao1" & AdultChild == "Child"),
                       y.position = 350, label = "p.format", size = pval_size) +
    ylim(0, 380) +
    theme(legend.position="none")+
    labs(x="") +
    scale_x_discrete(labels=c(
      "Rural\ninfants\nT1", "Urban\ninfants\nT1",
      "Rural\ninfants\nT2", "Urban\ninfants\nT2"))
  
# Shannon adults
plots$alpha$adults_shan = alpha %>%
  filter(AdultChild == "Adult") %>%
  ggplot(., aes(x=Group3, y=Shannon, color=Group3))+
    scale_color_manual(values = cols$adults) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size=point_size)+
    stat_pvalue_manual(subset(alpha_comp_location, Dependent == "Shannon" & AdultChild == "Adult"),
                       y.position = 6, label = "p.format", size = pval_size) +
    ylim(0, 6.5) +
    theme(legend.position="none")+
    labs(x="") +
    scale_x_discrete(labels=c(
      "Rural\nmothers\nT1", "Urban\nmothers\nT1",
      "Rural\nmothers\nT2", "Urban\nmothers\nT2"))

# Shannon infants
plots$alpha$infants_shan = alpha %>%
  filter(AdultChild == "Child") %>%
  ggplot(., aes(x=Group3, y=Shannon, color=Group3))+
    scale_color_manual(values = cols$infants) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size=point_size)+
    stat_pvalue_manual(subset(alpha_comp_location, Dependent == "Shannon" & AdultChild == "Child"),
                       y.position = 4.6, label = "p.format", size = pval_size) +
    ylim(0, 5.2) +
    theme(legend.position="none")+
    labs(x="") +
    scale_x_discrete(labels=c(
      "Rural\ninfants\nT1", "Urban\ninfants\nT1",
      "Rural\ninfants\nT2", "Urban\ninfants\nT2"))

# PD adults
plots$alpha$adults_PD = alpha %>%
  filter(AdultChild == "Adult") %>%
  ggplot(., aes(x=Group3, y=PD, color=Group3))+
    scale_color_manual(values = cols$adults) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size=point_size)+
    stat_pvalue_manual(subset(alpha_comp_location, Dependent == "PD" & AdultChild == "Adult"),
                       y.position = 42, label = "p.format", size = pval_size) +
    ylim(0, 45) +
    theme(legend.position="none")+
    labs(x="") +
    scale_x_discrete(labels=c(
      "Rural\nmothers\nT1", "Urban\nmothers\nT1",
      "Rural\nmothers\nT2", "Urban\nmothers\nT2"))

# PD infants
plots$alpha$infants_PD = alpha %>%
  filter(AdultChild == "Child") %>%
  ggplot(., aes(x=Group3, y=PD, color=Group3))+
    scale_color_manual(values = cols$infants) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size=point_size)+
    stat_pvalue_manual(subset(alpha_comp_location, Dependent == "PD" & AdultChild == "Child"),
                       y.position = 27, label = "p.format", size = pval_size) +
    ylim(0, 31) +
    theme(legend.position="none")+
    labs(x="") +
    scale_x_discrete(labels=c(
      "Rural\ninfants\nT1", "Urban\ninfants\nT1",
      "Rural\ninfants\nT2", "Urban\ninfants\nT2"))

##### ALPHA DIVERSITY: compare timepoint #####

# Reorder to compare timepoints easily
alpha$Group3 = factor(alpha$Group3, levels=c(
  "Rural.Adult.1", "Rural.Adult.2",
  "Urban.Adult.1", "Urban.Adult.2",
  "Rural.Child.1", "Rural.Child.2",
  "Urban.Child.1", "Urban.Child.2")
)

alpha_comp_time = my_compare_means(
  data = alpha, 
  dvars = c("Shannon", "Chao1", "PD"), 
  indvar = "TimePoint", 
  grouping = c("AdultChild", "Location"),
  pairing = "Index")

alpha_comp_time$group1 = rep(c("Rural.Adult.1", "Rural.Child.1",
                               "Urban.Adult.1", "Urban.Child.1"), each = 3)
alpha_comp_time$group2 = rep(c("Rural.Adult.2", "Rural.Child.2",
                               "Urban.Adult.2", "Urban.Child.2"), each = 3)
alpha_comp_time$p.format = formatC(alpha_comp_time$p, 2)

# Chao1 adults
plots$alpha$adults_chao_time = alpha %>%
  filter(AdultChild == "Adult") %>%
  ggplot(., aes(x=Group3, y=Chao1, color=Group3))+
    scale_color_manual(values = cols$adults[c(1,3,2,4)]) +
    geom_boxplot(outlier.shape = NA) +
    geom_line(aes(group = TimePairs), color = "grey") +
    geom_jitter(width = 0, size=point_size) +
    stat_pvalue_manual(subset(alpha_comp_time, Dependent == "Chao1" & AdultChild == "Adult"),
                       y.position = 480, label = "p.format", size = pval_size) +
    ylim(0, 520) +
    theme(legend.position="none")+
    labs(x="") +
    scale_x_discrete(labels=c(
      "Rural\nmothers\nT1", "Rural\nmothers\nT2",
      "Urban\nmothers\nT1", "Urban\nmothers\nT2"))

# Chao1 infants
plots$alpha$infants_chao_time = alpha %>%
  filter(AdultChild == "Child") %>%
    ggplot(., aes(x=Group3, y=Chao1, color=Group3))+
    scale_color_manual(values = cols$infants[c(1,3,2,4)]) +
    geom_boxplot(outlier.shape = NA) +
    geom_line(aes(group = TimePairs), color = "grey") +
    geom_jitter(width = 0, size=point_size) +
    stat_pvalue_manual(subset(alpha_comp_time, Dependent == "Chao1" & AdultChild == "Child"),
                       y.position = 350, label = "p.format", size = pval_size) +
    ylim(0, 380) +
    theme(legend.position="none")+
    labs(x="") +
    scale_x_discrete(labels=c(
      "Rural\ninfants\nT1", "Rural\ninfants\nT2",
      "Urban\ninfants\nT1", "Urban\ninfants\nT2"))

# Shannon adults
plots$alpha$adults_shan_time = alpha %>%
  filter(AdultChild == "Adult") %>%
  ggplot(., aes(x=Group3, y=Shannon, color=Group3))+
    scale_color_manual(values = cols$adults[c(1,3,2,4)]) +
    geom_boxplot(outlier.shape = NA) +
    geom_line(aes(group = TimePairs), color = "grey") +
    geom_jitter(width = 0, size=point_size) +
    stat_pvalue_manual(subset(alpha_comp_time, Dependent == "Shannon" & AdultChild == "Adult"),
                       y.position = 6, label = "p.format", size = pval_size) +
    ylim(0, 6.5) +
    theme(legend.position="none")+
    labs(x="") +
    scale_x_discrete(labels=c(
      "Rural\nmothers\nT1", "Rural\nmothers\nT2",
      "Urban\nmothers\nT1", "Urban\nmothers\nT2"))

# Shannon infants
plots$alpha$infants_shan_time = alpha %>%
  filter(AdultChild == "Child") %>%
  ggplot(., aes(x=Group3, y=Shannon, color=Group3))+
    scale_color_manual(values = cols$infants[c(1,3,2,4)]) +
    geom_boxplot(outlier.shape = NA) +
    geom_line(aes(group = TimePairs), color = "grey") +
    geom_jitter(width = 0, size=point_size) +
    stat_pvalue_manual(subset(alpha_comp_time, Dependent == "Shannon" & AdultChild == "Child"),
                       y.position = 4.6, label = "p.format", size = pval_size) +
    ylim(0, 5.2) +
    theme(legend.position="none")+
    labs(x="") +
    scale_x_discrete(labels=c(
      "Rural\ninfants\nT1", "Rural\ninfants\nT2",
      "Urban\ninfants\nT1", "Urban\ninfants\nT2"))

# PD adults
plots$alpha$adults_PD_time = alpha %>%
  filter(AdultChild == "Adult") %>%
  ggplot(., aes(x=Group3, y=PD, color=Group3))+
    scale_color_manual(values = cols$adults[c(1,3,2,4)]) +
    geom_boxplot(outlier.shape = NA) +
    geom_line(aes(group = TimePairs), color = "grey") +
    geom_jitter(width = 0, size=point_size) +
    stat_pvalue_manual(subset(alpha_comp_time, Dependent == "PD" & AdultChild == "Adult"),
                       y.position = 42, label = "p.format", size = pval_size) +
    ylim(0, 45) +
    theme(legend.position="none")+
    labs(x="") +
    scale_x_discrete(labels=c(
      "Rural\nmothers\nT1", "Rural\nmothers\nT2",
      "Urban\nmothers\nT1", "Urban\nmothers\nT2"))

# PD infants
plots$alpha$infants_PD_time = alpha %>%
  filter(AdultChild == "Child") %>%
  ggplot(., aes(x=Group3, y=PD, color=Group3))+
    scale_color_manual(values = cols$infants[c(1,3,2,4)]) +
    geom_boxplot(outlier.shape = NA) +
    geom_line(aes(group = TimePairs), color = "grey") +
    geom_jitter(width = 0, size=point_size) +
    stat_pvalue_manual(subset(alpha_comp_time, Dependent == "PD" & AdultChild == "Child"),
                       y.position = 27, label = "p.format", size = pval_size) +
    ylim(0, 31) +
    theme(legend.position="none")+
    labs(x="") +
    scale_x_discrete(labels=c(
      "Rural\ninfants\nT1", "Rural\ninfants\nT2",
      "Urban\ninfants\nT1", "Urban\ninfants\nT2"))

rm(alpha_comp_location, alpha_comp_time, alpha_div_comparison,
   alpha_ace, alpha_ace_rare, alpha_chao, alpha_chao_rare,
   alpha_pd, alpha_pd_rare, alpha_shan, alpha_shan_rare)

##### ALPHA DIVERSITY: time effect #####
# Some time-related variables are highly collinear in the infants
# In particular infant age and time since diet diversification
# Here I check which is most important for alpha diversity and select one
# so they are not both entered in the same model

point_size = 1

plots$alpha_reg = list()

plots$alpha_reg$time_since_preg = alpha %>%
  filter(AdultChild == "Adult") %>%
  ggplot(., aes(x=AgeC, y=Chao1, color=TimePoint))+
    geom_point(size = point_size, alpha=0.5)+
    geom_smooth(method="lm")+
    stat_cor() +
    theme(legend.position = "left")+
    labs(x="Time since pregnancy [months]")+
    ylim(c(0,500))

plots$alpha_reg$time_child_age = alpha %>%
  filter(AdultChild == "Child") %>%
  ggplot(., aes(x=AgeC, y=Chao1, color=TimePoint))+
  geom_point(size = point_size, alpha=0.5)+
  geom_smooth(method="lm")+
  stat_cor() +
  labs(x="Infant age [months]")

plots$alpha_reg$time_since_divers = alpha %>%
  filter(AdultChild == "Child") %>%
  ggplot(., aes(x=SinceDivers., y=Chao1, color=TimePoint))+
  geom_point(size = point_size, alpha=0.5)+
  geom_smooth(method="lm")+
  stat_cor() +
  labs(x="Time since diversification [months]")

# Infant age is more relevant to alpha div than time since diet diversification

##### ALPHA DIVERSITY: other effects' exploration #####

# Adults
p1 = alpha %>%
  filter(AdultChild == "Adult") %>%
  ggplot(., aes(x = BMI, y = Chao1)) +
  geom_point()+
  stat_cor()

p2 = alpha %>%
  filter(AdultChild == "Adult") %>%
  ggplot(., aes(x = WHR, y = Chao1)) +
  geom_point()+
  stat_cor()

p3 = alpha %>%
  filter(AdultChild == "Adult") %>%
  ggplot(., aes(x = IngredN, y = Chao1)) +
  geom_point()+
  stat_cor()

p4 = alpha %>%
  filter(AdultChild == "Adult") %>%
  ggplot(., aes(x = SharingMealAvg, y = Chao1)) +
  geom_point()+
  stat_cor()

p5 = alpha %>%
  filter(AdultChild == "Adult") %>%
  ggplot(., aes(x = PSystAvg, y = Chao1)) +
  geom_point()+
  stat_cor()

p6 = alpha %>%
  filter(AdultChild == "Adult") %>%
  ggplot(., aes(x = PDiasAvg, y = Chao1)) +
  geom_point()+
  stat_cor()

p7 = alpha %>%
  filter(AdultChild == "Adult") %>%
  ggplot(., aes(x = BMI, y = WHR)) +
  geom_point()+
  stat_cor()

plots$alpha_reg$explo_adults =  ggarrange(p1, p2, p3, p4, p5, p6, p7, ncol = 2, nrow = 4)

# Infants
p1 = alpha %>%
  filter(AdultChild == "Child") %>%
  ggplot(., aes(x = BMI, y = Chao1, color=TimePoint)) +
  geom_point()+
  stat_cor()

p2 = alpha %>%
  filter(AdultChild == "Child") %>%
  ggplot(., aes(x = SharingMealAvg, y = Chao1)) +
  geom_point()+
  stat_cor()

p3 = alpha %>%
  filter(AdultChild == "Child") %>%
  ggplot(., aes(x = IngredN, y = Chao1)) +
  geom_point()+
  stat_cor()

p4 = alpha %>%
  filter(AdultChild == "Child") %>%
  ggplot(., aes(x = HeadCirc, y = Chao1, color=TimePoint)) +
  geom_point()+
  stat_cor()

plots$alpha_reg$explo_infants = ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

rm(p1, p2, p3, p4, p5, p6, p7)

##### ALPHA DIVERSITY: regress adults #####

point_size = 1

alpha_nona = alpha %>%
  filter(AdultChild == "Adult") %>%
  drop_na(BMI) %>% 
  mutate(Index = as.factor(Index))

# Formula with all fixed effects to explore
full_formula = Chao1 ~ Location*AgeC+AgeM+BMI+EntColi+IngredN+SharingMealAvg
  
# Step AIC selection
null_model = lm(Chao1 ~ 1, data = alpha_nona)
full_model = lm(full_formula, data = alpha_nona)
step_model = stepAIC(null_model, direction = "both", scope=list(lower=null_model, upper=full_model))
summary(step_model)

cat("Is ent coli colonization more common in rural?\n")
meta %>%
  filter(AdultChild == "Adult") %>%
  dplyr::select(Location, EntColi) %>%
  table() %>%
  fisher.test()

plots$alpha_reg$ent_coli = ggplot(alpha_nona, aes(x = Group3, y = Chao1, color=EntColi)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.2), size=point_size)+
  scale_color_manual(values = c("#888888", "purple"))+
  theme(legend.position = "left")+
  labs(x="") +
  scale_x_discrete(labels=c(
    "Rural\nmothers\nT1","Rural\nmothers\nT2","Urban\nmothers\nT1","Urban\nmothers\nT2"))

##### ALPHA DIVERSITY: regress infants #####

alpha_nona = alpha %>%
  filter(AdultChild == "Child") %>%
  drop_na(BMI, Gender, FirstMeal, AgeC, SinceDivers.) %>% 
  mutate(Index = as.factor(Index)) %>%
  group_by(TimePoint) %>%
  mutate(outlier = is_outlier(Chao1))

print(paste(sum(alpha_nona$outlier), "infant samples are outliers for chao1"))
plots$alpha_reg$outliers_infants=alpha_nona %>%
  ggplot(., aes(x=Group3, y=Chao1))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(color=outlier), width = 0.2, size=point_size) +
  ylim(0, 300) +
  labs(x="") +
  scale_x_discrete(labels=c(
    "Rural\ninfants\nt1", "Rural\ninfants\nt2",
    "Urban\ninfants\nt1", "Urban\ninfants\nt2"))
# These are pretty big outliers, best to remove them

# Before and after outlier removal
p1 = ggplot(alpha_nona, aes(x=AgeC, y=Chao1, color=outlier))+
  geom_point() +
  labs(x="Infant age [months]")
alpha_nona = alpha_nona %>%
  filter(!outlier)
p2 = ggplot(alpha_nona, aes(x=AgeC, y=Chao1, color=outlier))+
  geom_point() +
  labs(x="Infant age [months]")
plots$alpha_reg$outliers_before_after=ggarrange(p1, p2, nrow=2)

# Formula with all fixed effects to explore
full_formula = Chao1 ~ Location*AgeC+AgeM+BMI+Gender+FirstMeal

# Step AIC selection
null_model = lm(Chao1 ~ 1, data = alpha_nona)
full_model = lm(full_formula, data = alpha_nona)
step_model = stepAIC(null_model, direction = "both", scope=list(lower=null_model, upper=full_model))
summary(step_model)

cat("Regress only for infants at T2 to check effects of food related variables\n")
full_formula2 = Chao1 ~ Location*AgeC+AgeM+BMI+Gender+SharingMealAvg+IngredN
alpha_nona2 = filter(alpha_nona, TimePoint == 2)

null_model2 = lm(Chao1 ~ 1, data = alpha_nona2)
full_model2 = lm(full_formula2, data = alpha_nona2)
step_model2 = stepAIC(null_model2, direction = "both", scope=list(lower=null_model2, upper=full_model2))
summary(step_model2)
# Nothing new compared to the previous model: effect of environment and maybe BMI

##### DISEASE INCIDENCE: infants #####

cat("Are health conditions more common in urban or rural infants at T1?\n")

# These notes are in HealthStatusC for T1 but are not diseases
not_diseases = c("Amoxicillin treatment in the last 72h",
                 "Delayed cicatrization of belly button",
                 "Formula fed in the first 48h of life, breastfed after",
                 "Occasionally fed with cow milk",
                 "")
meta %>%
  filter(Group4 == "Child.1") %>%
  mutate(AnyCondition = !(HealthStatusC %in% not_diseases)) %>%
  dplyr::select(Location, AnyCondition) %>%
  table()

cat("Test at T1\n")
meta %>%
  filter(Group4 == "Child.1") %>%
  mutate(AnyCondition = !(HealthStatusC %in% not_diseases)) %>%
  dplyr::select(Location, AnyCondition) %>%
  table() %>%
  fisher.test()

cat("Are health conditions more common in urban or rural infants at T2?\n")
meta %>%
  filter(Group4 == "Child.2") %>%
  mutate(AnyCondition = HealthStatusC != "") %>%
  dplyr::select(Location, AnyCondition) %>%
  table()

cat("Test at T2\n")
meta %>% 
  filter(Group4 == "Child.2") %>%
  mutate(AnyCondition = HealthStatusC != "") %>%
  dplyr::select(Location, AnyCondition) %>%
  table() %>%
  fisher.test()

cat("Are allergy symtoms more common in urban or rural?\n")
meta %>%
  filter(AdultChild == "Child") %>%
  group_by(Group3) %>%
  summarise(Allergy = sum(grepl("allerg", HealthStatusC, ignore.case = T) | grepl("allerg", HealthStatusPriorC, ignore.case = T)))

cat("Test at T2\n")
meta %>% 
  filter(Group4 == "Child.2") %>%
  mutate(Allergy = grepl("allerg", HealthStatusC, ignore.case = T) | grepl("allerg", HealthStatusPriorC, ignore.case = T)) %>%
  dplyr::select(Location, Allergy) %>%
  table() %>%
  fisher.test()

##### BETA DIVERSITY: prep #####

beta = list()
beta$bray = vegdist(ab$asv_rare)
beta$unifracU = rbiom::unifrac(t(ab$asv_rare), weighted=F, tree=tree)
beta$unifracW = rbiom::unifrac(t(ab$asv_rare), weighted=T, tree=tree)

##### BETA DIVERSITY: compare metrics #####

plots$beta$metric_comp = ggpairs(as.data.frame(do.call(cbind, beta)), 
                                 lower = list(continuous = "density"))

beta$bray = as.matrix(beta$bray)
beta$unifracU = as.matrix(beta$unifracU)
beta$unifracW = as.matrix(beta$unifracW)

# Here I look at how well these beta div metrics disctiminate between groups
# i.e. adult vs child or individuals. This serves as sanity check

p1 = in_vs_out_group_dists(beta$bray, meta, "AdultChild") %>%
  ggplot(., aes(x=value, fill=type))+
  geom_density(alpha=0.5)+
  xlim(c(0,1))
p2 = in_vs_out_group_dists(beta$unifracU, meta, "AdultChild") %>%
  ggplot(., aes(x=value, fill=type))+
  geom_density(alpha=0.5)+
  xlim(c(0,1))
p3 = in_vs_out_group_dists(beta$unifracW, meta, "AdultChild") %>%
  ggplot(., aes(x=value, fill=type))+
  geom_density(alpha=0.5)+
  xlim(c(0,1))
plots$beta$hists_AdultChild = ggarrange(p1, p2, p3, nrow=3, common.legend = T)
# Obviously unifracW discriminates adults vs children much better

p1 = in_vs_out_group_dists(beta$bray, filter(meta, AdultChild == "Adult"), "Index") %>%
  ggplot(., aes(x=value, fill=type))+
  geom_density(alpha=0.5)+
  xlim(c(0,1))
p2 = in_vs_out_group_dists(beta$unifracU, filter(meta, AdultChild == "Adult"), "Index") %>%
  ggplot(., aes(x=value, fill=type))+
  geom_density(alpha=0.5)+
  xlim(c(0,1))
p3 = in_vs_out_group_dists(beta$unifracW, filter(meta, AdultChild == "Adult"), "Index") %>%
  ggplot(., aes(x=value, fill=type))+
  geom_density(alpha=0.5)+
  xlim(c(0,1))
plots$beta$hists_individuals_adults = ggarrange(p1, p2, p3, nrow = 3, common.legend = T)
# All 3 metrics seem equally good at recognizing that the same individuals are closest to themselves

p1 = in_vs_out_group_dists(beta$bray, filter(meta, AdultChild == "Child"), "Index", list("TimePoint.x != TimePoint.y")) %>%
  ggplot(., aes(x=value, fill=type))+
  geom_density(alpha=0.5)+
  xlim(c(0,1))
p2 = in_vs_out_group_dists(beta$unifracU, filter(meta, AdultChild == "Child"), "Index", list("TimePoint.x != TimePoint.y")) %>%
  ggplot(., aes(x=value, fill=type))+
  geom_density(alpha=0.5)+
  xlim(c(0,1))
p3 = in_vs_out_group_dists(beta$unifracW, filter(meta, AdultChild == "Child"), "Index", list("TimePoint.x != TimePoint.y")) %>%
  ggplot(., aes(x=value, fill=type))+
  geom_density(alpha=0.5)+
  xlim(c(0,1))
plots$beta$hists_individuals_children = ggarrange(p1, p2, p3, nrow = 3, common.legend = T)
# As expected, t1 carries no info on t2 for children

##### BETA DIVERSITY: adonis adults #####

plots$beta_ado = list()

samples_to_keep = rownames(filter(meta, AdultChild == "Adult"))

# Adonis data prep
adonis_data = list()
adonis_data$bray = prep_for_adonis(beta$bray, meta, samples_to_keep, c("AgeC", "AgeM", "EntColi", "BMI"))
adonis_data$unifracU = prep_for_adonis(beta$unifracU, meta, samples_to_keep, c("AgeC", "AgeM", "EntColi", "BMI"))
adonis_data$unifracW = prep_for_adonis(beta$unifracW, meta, samples_to_keep, c("AgeC", "AgeM", "EntColi", "BMI"))

nperms = 2000
set.seed(1337)

# Just location
adonis2(as.dist(adonis_data$bray[[1]]) ~ Location, 
        data = adonis_data$bray[[2]],
        by = "terms",
        permutations = nperms)
adonis2(as.dist(adonis_data$unifracU[[1]]) ~ Location, 
        data = adonis_data$unifracU[[2]],
        by = "terms",
        permutations = nperms)
adonis2(as.dist(adonis_data$unifracW[[1]]) ~ Location, 
        data = adonis_data$unifracW[[2]],
        by = "terms",
        permutations = nperms)

# Bray-Curtis with location + more variables
plots$beta_ado$adu_bray_loc = 
  adonis2(as.dist(adonis_data$bray[[1]]) ~ Location+AgeC+AgeM+EntColi+BMI, 
        data = adonis_data$bray[[2]],
        by = "margin",
        permutations = nperms) %>%
  slice(1:(n()-2)) %>%
  arrange(-R2) %>%
  rownames_to_column(var="Variable") %>%
  mutate(Variable = fct_reorder(Variable, R2)) %>%
  mutate(Sig = stars.pval(`Pr(>F)`)) %>%
  ggplot(., aes(x=Variable, y=R2))+
  geom_bar(stat = "identity") +
  geom_text(aes(label = Sig, y = R2*1.1), size=6, vjust=0.5, hjust=0) +
  ylim(c(0, 0.035))+
  coord_flip()+
  scale_x_discrete(labels = c("Location"="Environment", "AgeC"="TsincePreg.", "AgeM"="Age"))+
  labs(x="")+
  ggtitle("Bray-Curtis")+
  theme(plot.title = element_text(hjust = 0.5))

# UniFracU with location + more variables
plots$beta_ado$adu_uniU_loc =
  adonis2(as.dist(adonis_data$unifracU[[1]]) ~ Location+AgeC+AgeM+EntColi+BMI, 
        data = adonis_data$unifracU[[2]],
        by = "margin",
        permutations = nperms) %>%
  slice(1:(n()-2)) %>%
  arrange(-R2) %>%
  rownames_to_column(var="Variable") %>%
  mutate(Variable = fct_reorder(Variable, R2)) %>%
  mutate(Sig = stars.pval(`Pr(>F)`)) %>%
  ggplot(., aes(x=Variable, y=R2))+
  geom_bar(stat = "identity") +
  geom_text(aes(label = Sig, y = R2*1.1), size=6, vjust=0.5, hjust=0) +
  ylim(c(0, 0.045))+
  coord_flip()+
  scale_x_discrete(labels = c("Location"="Environment", "AgeC"="TsincePreg.", "AgeM"="Age"))+
  labs(x="")+
  ggtitle("UniFracU")+
  theme(plot.title = element_text(hjust = 0.5))

# UniFracW with location + more variables
plots$beta_ado$adu_uniW_loc =
  adonis2(as.dist(adonis_data$unifracW[[1]]) ~ Location+AgeC+AgeM+EntColi+BMI, 
        data = adonis_data$unifracW[[2]],
        by = "margin",
        permutations = nperms) %>%
  slice(1:(n()-2)) %>%
  arrange(-R2) %>%
  rownames_to_column(var="Variable") %>%
  mutate(Variable = fct_reorder(Variable, R2)) %>%
  mutate(Sig = stars.pval(`Pr(>F)`)) %>%
  ggplot(., aes(x=Variable, y=R2))+
  geom_bar(stat = "identity") +
  geom_text(aes(label = Sig, y = R2*1.1), size=6, vjust=0.5, hjust=0) +
  ylim(c(0, 0.045))+
  coord_flip()+
  scale_x_discrete(labels = c("Location"="Environment", "AgeC"="TsincePreg.", "AgeM"="Age"))+
  labs(x="")+
  ggtitle("UniFracW")+
  theme(plot.title = element_text(hjust = 0.5))

##### BETA DIVERSITY: mds adults #####

point_size = 1

ncomps = 2
mds$adults = list()
mds$adults$bray = cmdscale(beta$bray[samples_to_keep, samples_to_keep], eig = TRUE, k = ncomps)
mds$adults$unifracU = cmdscale(beta$unifracU[samples_to_keep, samples_to_keep], eig = TRUE, k = ncomps)
mds$adults$unifracW = cmdscale(beta$unifracW[samples_to_keep, samples_to_keep], eig = TRUE, k = ncomps)

mds_points$adults = list()
for (type in names(mds$adults)) {
  mds_points$adults[[type]] = as.data.frame(mds$adults[[type]]$points)
  colnames(mds_points$adults[[type]]) = paste0("PC", 1:ncol(mds_points$adults[[type]]))
  mds_points$adults[[type]] = merge(mds_points$adults[[type]], meta, by = 0)
  rownames(mds_points$adults[[type]]) = mds_points$adults[[type]]$Row.names
  mds_points$adults[[type]] = mds_points$adults[[type]][,-c(1)]
}

plots$beta$ent_coli_pcoa = ellipse_and_centroid(mds_points$adults$bray, "PC1", "PC2", "EntColi", F, 
                                                point_size = point_size, point_alpha = 0.5)+
  scale_color_manual(values=c("#888888", "purple"))

##### BETA DIVERSITY: within group similarities #####

plots$beta_extra = list()

# Nicer metric names for plotting
metric_names = list()
metric_names$bray = "Bray-Curtis"
metric_names$unifracU = "UniFracU"
metric_names$unifracW = "UniFracW"

for (metric in c("bray", "unifracU", "unifracW")) {
  dists = in_vs_out_group_dists(beta[[metric]], meta, "Location")
  # Are individuals closer to each other in rural as compared to urban?
  p1 = dists %>%
    filter(type == "in_group") %>%
    ggplot(., aes(x=Location.x, y=value))+
    geom_boxplot()+
    stat_compare_means()+
    ggtitle("In group distances, all individuals")
  # Are mothers closer to each other in rural as compared to urban?
  p2 = dists %>%
    filter(type == "in_group") %>%
    filter(AdultChild.x == "Adult", AdultChild.y == "Adult") %>%
    ggplot(., aes(x=Location.x, y=value))+
    geom_boxplot()+
    stat_compare_means()+
    ggtitle("In group distances, mothers")
  # Are infants closer to each other in rural as compared to urban?
  p3 = dists %>%
    filter(type == "in_group") %>%
    filter(AdultChild.x == "Child", AdultChild.y == "Child") %>%
    ggplot(., aes(x=Location.x, y=value))+
    geom_boxplot()+
    stat_compare_means()+
    ggtitle("In group distances, infants")
  
  plots$beta_extra[[paste0("in_group_", metric)]] = ggarrange(p1, p2, p3, ncol = 3)
  
  # Are infants closer to adults in rural as compared to urban?
  dists2 = dists %>%
    filter(type == "in_group") %>%
    filter(AdultChild.x != AdultChild.y) %>%
    filter(TimePoint.x == TimePoint.y) %>%
    mutate(LocationTime = interaction(Location.x, TimePoint.x))
  
  dist_comp = my_compare_means(dists2, "value", "Location.x", "TimePoint.x") %>%
    mutate(group1 = interaction("Rural", TimePoint.x))%>%
    mutate(group2 = interaction("Urban", TimePoint.x)) %>%
    mutate(p.format = formatC(p, 2))
  print(dist_comp)
  
  plots$beta_extra[[paste0("adu_child_dist_", metric)]] = dists2 %>%
    mutate(grp4 = interaction(Location.x, TimePoint.x)) %>%
      ggplot(., aes(x=grp4, y=value))+
      geom_violin(aes(color=grp4, fill=grp4))+
      geom_boxplot(width=0.1, outlier.shape = NA, aes(color=grp4))+
      scale_color_manual(values = cols$adults)+
      scale_fill_manual(values = cols$infants)+
      stat_pvalue_manual(dist_comp, y.position = 1.1*max(dists2$value), size = pval_size, label="p.format")+
      ggtitle("Adult-infant distances")+
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(0, 1.2)+
      theme(legend.position="none")+
      labs(y=metric_names[[metric]], x=NULL) +
      scale_x_discrete(labels=c(
        "Rural\nT1", "Urban\nT1",
        "Rural\nT2", "Urban\nT2"))
}

cat("Rural infants are ahead in 'maturation', is it because they're older?\n")
meta %>%
  filter(AdultChild == "Child") %>%
  group_by(Group3) %>%
  summarize(mean(AgeC))
cat("They are only 9 days older on average, is this difference significant?\n")
meta %>%
  filter(AdultChild == "Child") %>%
  t.test(AgeC ~ Location, .)

##### BETA DIVERSITY: related vs unrelated #####

dodge_w = 0.85
dodge = position_dodge(width = dodge_w)

for (metric in c("bray", "unifracU", "unifracW")) {
  # Distances between mothers and children at same timepoint
  # Groups if same index
  dists = in_vs_out_group_dists(beta[[metric]], meta, "Index", list(
    "AdultChild.x != AdultChild.y",
    "TimePoint.x == TimePoint.y"
  ))
  
  dist_comp = my_compare_means(dists, "value", "type", "TimePoint") %>%
    mutate(group1 = "in_group")%>%
    mutate(group2 = "out_group") %>%
    mutate(xmin= c(1-dodge_w/4, 2-dodge_w/4)) %>%
    mutate(xmax= c(1+dodge_w/4, 2+dodge_w/4)) %>%
    mutate(p.format = formatC(p, 2))
  
  # Are related mother-child pairs more similar than unrelated pairs?
  plots$beta_extra[[paste0("relatedness_effect_", metric)]] =  
    ggplot(dists, aes(x=TimePoint, y=value, color=type))+
      geom_violin(position=dodge, aes(fill=type), alpha=0.5)+
      geom_boxplot(position=dodge, width = 0.1, outlier.shape = NA)+
      stat_pvalue_manual(dist_comp, y.position = 1.1*max(dists$value), size = pval_size, label="p.format")+
      scale_color_discrete(name = "Relationship", labels = c("Related", "Not related"))+
      scale_fill_discrete(name = "Relationship", labels = c("Related", "Not related"))+
      labs(y=metric_names[[metric]])+
      ylim(c(0,1.2))
  
  # Are infants closer to their mothers in rural as compared to urban?
  plots$beta_extra[[paste0("mother_own_child_dists_", metric)]] =
    dists %>%
    filter(type == "in_group") %>%
    ggplot(., aes(x=Location.x, y=value))+
      geom_boxplot()+
      labs(y="Mother-child distances")+
      stat_compare_means()+
      ylim(c(0,1.2))
}

##### BETA DIVERSITY: adonis infants #####

samples_to_keep = rownames(filter(meta, AdultChild == "Child"))

adonis_data = list()
adonis_data$bray = prep_for_adonis(beta$bray, meta, samples_to_keep, c("Gender", "AgeC", "AgeM", "FirstMeal", "SinceDivers.", "BMI"))
adonis_data$unifracU = prep_for_adonis(beta$unifracU, meta, samples_to_keep, c("Gender", "AgeC", "AgeM", "FirstMeal", "SinceDivers.", "BMI"))
adonis_data$unifracW = prep_for_adonis(beta$unifracW, meta, samples_to_keep, c("Gender", "AgeC", "AgeM", "FirstMeal", "SinceDivers.", "BMI"))

nperms = 2000
set.seed(1337)

# Just location
adonis2(as.dist(adonis_data$bray[[1]]) ~ Location, 
        data = adonis_data$bray[[2]],
        by = "terms",
        permutations = nperms)
adonis2(as.dist(adonis_data$unifracU[[1]]) ~ Location, 
        data = adonis_data$unifracU[[2]],
        by = "terms",
        permutations = nperms)
adonis2(as.dist(adonis_data$unifracW[[1]]) ~ Location, 
        data = adonis_data$unifracW[[2]],
        by = "terms",
        permutations = nperms)

# Bray-Curtis with location + more variables
plots$beta_ado$chi_bray_loc = 
  adonis2(as.dist(adonis_data$bray[[1]]) ~ Location+AgeC+Gender+BMI, 
          data = adonis_data$bray[[2]],
          by = "margin",
          permutations = nperms) %>%
  slice(1:(n()-2)) %>%
  arrange(-R2) %>%
  rownames_to_column(var="Variable") %>%
  mutate(Variable = fct_reorder(Variable, R2)) %>%
  mutate(Sig = stars.pval(`Pr(>F)`)) %>%
  ggplot(., aes(x=Variable, y=R2))+
  geom_bar(stat = "identity") +
  geom_text(aes(label = Sig, y = R2*1.1), size=6, vjust=0.5, hjust=0) +
  ylim(c(0, 0.10))+
  coord_flip()+
  scale_x_discrete(labels = c("Location"="Environment", "AgeC" = "Age"))+
  labs(x="")+
  ggtitle("Bray-Curtis")+
  theme(plot.title = element_text(hjust = 0.5))

# UniFracU with location + more variables
plots$beta_ado$chi_uniU_loc =
  adonis2(as.dist(adonis_data$unifracU[[1]]) ~ Location+AgeC+Gender+BMI, 
          data = adonis_data$unifracU[[2]],
          by = "margin",
          permutations = nperms) %>%
  slice(1:(n()-2)) %>%
  arrange(-R2) %>%
  rownames_to_column(var="Variable") %>%
  mutate(Variable = fct_reorder(Variable, R2)) %>%
  mutate(Sig = stars.pval(`Pr(>F)`)) %>%
  ggplot(., aes(x=Variable, y=R2))+
  geom_bar(stat = "identity") +
  geom_text(aes(label = Sig, y = R2*1.1), size=6, vjust=0.5, hjust=0) +
  ylim(c(0, 0.16))+
  coord_flip()+
  scale_x_discrete(labels = c("Location"="Environment", "AgeC" = "Age"))+
  labs(x="")+
  ggtitle("UniFracU")+
  theme(plot.title = element_text(hjust = 0.5))

# UniFracW with location + more variables
plots$beta_ado$chi_uniW_loc =
  adonis2(as.dist(adonis_data$unifracW[[1]]) ~ Location+AgeC+Gender+BMI, 
          data = adonis_data$unifracW[[2]],
          by = "margin",
          permutations = nperms) %>%
  slice(1:(n()-2)) %>%
  arrange(-R2) %>%
  rownames_to_column(var="Variable") %>%
  mutate(Variable = fct_reorder(Variable, R2)) %>%
  mutate(Sig = stars.pval(`Pr(>F)`)) %>%
  ggplot(., aes(x=Variable, y=R2))+
  geom_bar(stat = "identity") +
  geom_text(aes(label = Sig, y = R2*1.1), size=6, vjust=0.5, hjust=0) +
  ylim(c(0, 0.16))+
  coord_flip()+
  scale_x_discrete(labels = c("Location"="Environment", "AgeC" = "Age"))+
  labs(x="")+
  ggtitle("UniFracW")+
  theme(plot.title = element_text(hjust = 0.5))

##### BETA DIVERSITY: adonis infants T2 only #####

samples_to_keep = rownames(filter(meta, Group4 == "Child.2"))

adonis_data = list()
adonis_data$bray = prep_for_adonis(beta$bray, meta, samples_to_keep, c("Gender", "AgeC", "AgeM", "FirstMeal", "SinceDivers.", "BMI", "IngredN"))
adonis_data$unifracU = prep_for_adonis(beta$unifracU, meta, samples_to_keep, c("Gender", "AgeC", "AgeM", "FirstMeal", "SinceDivers.", "BMI", "IngredN"))
adonis_data$unifracW = prep_for_adonis(beta$unifracW, meta, samples_to_keep, c("Gender", "AgeC", "AgeM", "FirstMeal", "SinceDivers.", "BMI", "IngredN"))

nperms = 1000
set.seed(1337)

# Just location
adonis2(as.dist(adonis_data$bray[[1]]) ~ Location, 
        data = adonis_data$bray[[2]],
        by = "terms",
        permutations = nperms)
adonis2(as.dist(adonis_data$unifracU[[1]]) ~ Location, 
        data = adonis_data$unifracU[[2]],
        by = "terms",
        permutations = nperms)
adonis2(as.dist(adonis_data$unifracW[[1]]) ~ Location, 
        data = adonis_data$unifracW[[2]],
        by = "terms",
        permutations = nperms)

# Bray-Curtis with Location
plots$beta_ado$chi2_bray_loc =
  adonis2(as.dist(adonis_data$bray[[1]]) ~ Location+AgeC+Gender+BMI+SharingMealAvg+IngredN,
          data = adonis_data$bray[[2]],
          by = "margin",
          permutations = nperms) %>%
  slice(1:(n()-2)) %>%
  arrange(-R2) %>%
  rownames_to_column(var="Variable") %>%
  mutate(Variable = fct_reorder(Variable, R2)) %>%
  mutate(Sig = stars.pval(`Pr(>F)`)) %>%
  ggplot(., aes(x=Variable, y=R2))+
  geom_bar(stat = "identity") +
  geom_text(aes(label = Sig, y = R2*1.1), size=6, vjust=0.5, hjust=0) +
  ylim(c(0, 0.14))+
  coord_flip()+
  scale_x_discrete(labels = c("SubLocation"="Area"))+
  labs(x="")+
  ggtitle("Bray-Curtis")+
  theme(plot.title = element_text(hjust = 0.5))

# UniFracU with Location
plots$beta_ado$chi2_uniU_loc =
  adonis2(as.dist(adonis_data$unifracU[[1]]) ~ Location+AgeC+Gender+BMI+SharingMealAvg+IngredN, 
          data = adonis_data$unifracU[[2]],
          by = "margin",
          permutations = nperms) %>%
  slice(1:(n()-2)) %>%
  arrange(-R2) %>%
  rownames_to_column(var="Variable") %>%
  mutate(Variable = fct_reorder(Variable, R2)) %>%
  mutate(Sig = stars.pval(`Pr(>F)`)) %>%
  ggplot(., aes(x=Variable, y=R2))+
  geom_bar(stat = "identity") +
  geom_text(aes(label = Sig, y = R2*1.1), size=6, vjust=0.5, hjust=0) +
  ylim(c(0, 0.14))+
  coord_flip()+
  scale_x_discrete(labels = c("SubLocation"="Area"))+
  labs(x="")+
  ggtitle("UniFracU")+
  theme(plot.title = element_text(hjust = 0.5))

# UniFracW with Location
plots$beta_ado$chi2_uniW_loc =
  adonis2(as.dist(adonis_data$unifracW[[1]]) ~ Location+AgeC+Gender+BMI+SharingMealAvg+IngredN, 
          data = adonis_data$unifracW[[2]],
          by = "margin",
          permutations = nperms) %>%
  slice(1:(n()-2)) %>%
  arrange(-R2) %>%
  rownames_to_column(var="Variable") %>%
  mutate(Variable = fct_reorder(Variable, R2)) %>%
  mutate(Sig = stars.pval(`Pr(>F)`)) %>%
  ggplot(., aes(x=Variable, y=R2))+
  geom_bar(stat = "identity") +
  geom_text(aes(label = Sig, y = R2*1.1), size=6, vjust=0.5, hjust=0) +
  ylim(c(0, 0.14))+
  coord_flip()+
  scale_x_discrete(labels = c("SubLocation"="Area"))+
  labs(x="")+
  ggtitle("UniFracW")+
  theme(plot.title = element_text(hjust = 0.5))

##### BETA DIVERSITY: mds infants #####

point_size = 1

mds$all$bray = cmdscale(beta$bray, eig = TRUE, k = ncomps)
mds$all$unifracU = cmdscale(beta$unifracU, eig = TRUE, k = ncomps)
mds$all$unifracW = cmdscale(beta$unifracW, eig = TRUE, k = ncomps)

for (type in c("bray", "unifracU", "unifracW")) {
  mds_points$all[[type]] = as.data.frame(mds$all[[type]]$points)
  colnames(mds_points$all[[type]]) = paste0("PC", 1:ncol(mds_points$all[[type]]))
  mds_points$all[[type]] = merge(mds_points$all[[type]], meta, by = 0)
  rownames(mds_points$all[[type]]) = mds_points$all[[type]]$Row.names
  mds_points$all[[type]] = mds_points$all[[type]][,-c(1)]
}

# Bray
ggplot(mds_points$all$bray, aes(x=PC1, y=PC2, color=Group4))+
  geom_point(size=point_size)+
  stat_ellipse(level = 0.68)
# UnifracU
ggplot(mds_points$all$unifracU, aes(x=PC1, y=PC2, color=Group4))+
  geom_point(size=point_size)+
  stat_ellipse(level = 0.68)
# UnifracW
ggplot(mds_points$all$unifracW, aes(x=PC1, y=PC2, color=Group4))+
  geom_point(size=point_size)+
  stat_ellipse(level = 0.68)

plots$beta$adu_chi_pcoa = ellipse_and_centroid(mds_points$all$bray, "PC1", "PC2", "Group4", F, 
                                               point_size = point_size, point_alpha = 0.5)+
  scale_color_manual(values=cols$times,
                     labels=c("Mothers T1", "Infants T1", "Mothers T2", "Infants T2"),
                     name=NULL)+
  theme(legend.position = "left")

# Infants only MDS

samples_to_keep = rownames(filter(meta, AdultChild == "Child"))

mds$infants = list()
mds$infants$bray = cmdscale(beta$bray[samples_to_keep, samples_to_keep], eig = TRUE, k = ncomps)
mds$infants$unifracU = cmdscale(beta$unifracU[samples_to_keep, samples_to_keep], eig = TRUE, k = ncomps)
mds$infants$unifracW = cmdscale(beta$unifracW[samples_to_keep, samples_to_keep], eig = TRUE, k = ncomps)

mds_points$infants = list()
for (type in names(mds$infants )) {
  mds_points$infants[[type]] = as.data.frame(mds$infants[[type]]$points)
  colnames(mds_points$infants[[type]]) = paste0("PC", 1:ncol(mds_points$infants[[type]]))
  mds_points$infants[[type]] = merge(mds_points$infants[[type]], meta, by = 0)
  rownames(mds_points$infants[[type]]) = mds_points$infants[[type]]$Row.names
  mds_points$infants[[type]] = mds_points$infants[[type]][,-c(1)]
}

plots$beta$chi_pcoa = ellipse_and_centroid(mds_points$infants$bray, "PC1", "PC2", "Group3", F, 
                                           point_size = point_size, point_alpha = 0.5)+
  scale_color_manual(values=cols$infants, 
                     labels=c("Rural\ninfants T1", "Urban\ninfants T1", "Rural\ninfants T2", "Urban\ninfants T2"),
                     name=NULL)

##### DIFFERENTIAL ABUNDANCE: consider sparsity #####

sparsity_hists = list()
for (level in c("asv", "spe", "gen", "phy")) {
  sparsity_hists[[level]] = ab[[level]] %>%
    summarize_all(function(x) sum(x == 0)) %>%
    t() %>%
    as.data.frame() %>%
    ggplot(., aes(x=V1))+
    geom_histogram(binwidth = 5) +
    ggtitle(level) +
    xlab("N of samples with zero abundance")
}

plots$extra$sparsity = ggarrange(plotlist = sparsity_hists, ncol = 1, nrow = 4)

taxa$gen_defined = taxa$gen %>%
  filter(Genus != "__") %>%
  filter(!grepl("uncultured", String))

nonzeros_per_group = ab$gen %>%
  dplyr::select(rownames(taxa$gen_defined)) %>%
  merge(meta, by="row.names") %>%
  group_by(Group3) %>%
  summarize(across(starts_with("gen", ignore.case=F), function(x) sum(x != 0))) %>%
  as.data.frame()
rownames(nonzeros_per_group) = nonzeros_per_group$Group3
nonzeros_per_group = nonzeros_per_group[,-c(1)]

taxa$gen_less_sparse = taxa$gen_defined[colnames(nonzeros_per_group)[which(colSums(nonzeros_per_group > 5) > 0)], ]
ab$gen_less_sparse = ab$gen[, rownames(taxa$gen_less_sparse)]

plots$extra$sparsity_kept_counts =  data.frame("X1" = rowSums(ab$gen_less_sparse) / rowSums(ab$gen)) %>%
  ggplot(., aes(x=X1))+
    geom_histogram(aes(y=cumsum(..count..)/nrow(ab$gen)))

print(paste("Before filtering:", nrow(taxa$gen), "| After filtering:", nrow(taxa$gen_less_sparse)))
print(paste("Removed taxa account for: ", mean(100 - 100*rowSums(ab$gen_less_sparse) / rowSums(ab$gen)), "on average"))

##### DIFFERENTIAL ABUNDANCE: maaslin adults #####

plots$maaslin = list()

dir.create(paste0(outpath, "/maaslin"))

ab_mothers = ab$gen_less_sparse %>%
  filter(rownames(.) %in% rownames(subset(meta, AdultChild == "Adult")))
colnames(ab_mothers) = str_replace_all(taxa$gen_less_sparse[colnames(ab_mothers), "Genus"], "g__", "")
meta$Index = as.factor(meta$Index)

maaslin_mothers = Maaslin2(ab_mothers, meta, paste0(outpath, "/maaslin/mothers"),
         fixed_effects = c("Location", "AgeC", "EntColi"),
         random_effects = c("Index", "TimePoint"),
         max_significance = 0.1,
         min_prevalence = 0,
         analysis_method = "LM",  # Default: LM
         normalization = "TSS",   # Default: TSS
         transform = "LOG")      # Default: LOG

res = data.frame(maaslin_mothers$results) %>%
  filter(qval < 0.1) %>%
  mutate(feature = fct_reorder(feature, qval, min)) %>%
  mutate(logq = -sign(coef)*log(qval))

plots$maaslin$adults = ggplot(res, aes(x=feature, y=value, fill=logq, size=abs(coef)))+
  geom_point(pch=21) +
  scale_fill_gradient2(low="blue", mid="light grey", high="red", limits = c(-max(abs(res$logq)), max(abs(res$logq)))) +
  scale_size(breaks = c(0.25, 0.5, 0.75, 1))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.direction="horizontal",
        legend.position = "top",
        legend.key.height = unit(2, "mm"))+
  labs(x=NULL, y=NULL)+
  scale_y_discrete(labels=c("TRUE"="EntColi+", "AgeC" = "TsincePreg."))

##### DIFFERENTIAL ABUNDANCE: maaslin infants #####

ab_infants = ab$gen_less_sparse %>%
  filter(rownames(.) %in% rownames(subset(meta, AdultChild == "Child")))
colnames(ab_infants) = str_replace_all(taxa$gen_less_sparse[colnames(ab_infants), "Genus"], "g__", "")
meta["UrbanAgeC"] = (meta$Location == "Urban") * meta$AgeC

maaslin_infants = Maaslin2(ab_infants, meta, paste0(outpath, "/maaslin/infants"),
                           fixed_effects = c("Location", "AgeC", "UrbanAgeC"),
                           random_effects = c("Index", "TimePoint"),
                           max_significance = 0.1,
                           min_prevalence = 0,
                           analysis_method = "LM",  # Default: LM
                           normalization = "TSS",   # Default: TSS
                           transform = "LOG")      # Default: LOG

# It seems the differences in microbiome between urban and rural and in how their
# microbiome develops over time depends on the acquisition of rare species mostly
# because alpha diversity grows faster in rural children but very few taxa
# have abundance associated with rural/urban or rural/urban*age

res = data.frame(maaslin_infants$results) %>%
  filter(qval < 0.1) %>%
  mutate(feature = fct_reorder(feature, qval, min)) %>%
  mutate(logq = -sign(coef)*log(qval))

plots$maaslin$infants =
  ggplot(res, aes(x=feature, y=value, fill=logq, size=abs(coef)))+
  geom_point(pch=21) +
  scale_fill_gradient2(low="blue", mid="light grey", high="red", limits = c(-max(abs(res$logq)), max(abs(res$logq)))) +
  scale_size(breaks = c(0.2, 0.38, 0.58, 0.7))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.direction="horizontal",
        legend.position = "top",
        legend.key.height = unit(2, "mm"))+
  labs(x=NULL, y=NULL)+
  scale_y_discrete(labels = c("AgeC" = "Age", "UrbanAgeC" = "Urban*Age"))

##### PLOTTING #####

w = 174 # mm
h = 230

dir.create(paste0(outpath, "/extras"))

# Figure S1
p1 = (plots$qc$coverage | plots$qc$raref_gen)
p2 = (plots$qc$raref_spe | plots$qc$raref_asv)
p3 = (plots$beta$all_pcoa | guide_area() + plot_layout(guides = 'collect')) + plot_layout(widths = c(2, 1))
figS1 = p1 / p2 / p3 + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = "bold", size=14))
ggsave(paste0(outpath, "/figS1.pdf"), plot=figS1,
       width = w, height = h*0.7, units = "mm", dpi = 300)

# Figure 1 phylum plots
ggsave(paste0(outpath, "/fig1.pdf"), plot = plots$phy_barplots$combined, 
       width = w, height = h*0.35, units = "mm", dpi = 300)

# All alpha div metrics, comparing environment
figS2 = ggarrange(plotlist = plots$alpha[c(
  "adults_chao", "infants_chao", 
  "adults_shan", "infants_shan",
  "adults_PD", "infants_PD")], 
  nrow=3, ncol=2, labels = "AUTO")
ggsave(paste0(outpath, "/figS2.pdf"), plot = figS2,
       width = w, height = 0.8*h, units = "mm", dpi = 300)

# All alpha div metrics comparing timepoints
alpha_div_all_metrics_time = ggarrange(plotlist = plots$alpha[c(
  "adults_chao_time", "infants_chao_time", 
  "adults_shan_time", "infants_shan_time",
  "adults_PD_time", "infants_PD_time")], 
  nrow=3, ncol=2, labels = "AUTO")
ggsave(paste0(outpath, "/extras/alpha_div_all_metrics_time.pdf"), plot = alpha_div_all_metrics_time,
       width = w, height = 0.8*h, units = "mm", dpi = 300)

# Time variables vs diversity
alpha_time_vars =
  ggarrange(plots$alpha_reg$time_since_preg, plots$alpha_reg$time_child_age, plots$alpha_reg$time_since_divers, 
            common.legend = T, ncol=3)
ggsave(paste0(outpath, "/extras/alpha_time_vars.pdf"), plot = alpha_time_vars,
       width = w, height = h*0.25, units = "mm", dpi = 300)

# Ability of beta div metrics to discriminate between various groups
beta_grouping_hists = ggarrange(
  annotate_figure(plots$beta$hists_AdultChild, top = text_grob("Group = Adult/Child", face = "bold", size = 14)),
  annotate_figure(plots$beta$hists_individuals_adults, top = text_grob("Group = individuals (adults)", face = "bold", size = 14)),
  annotate_figure(plots$beta$hists_individuals_children, top = text_grob("Group = individuals (infants)", face = "bold", size = 14)),
  ncol=3, common.legend = T)
ggsave(paste0(outpath, "/extras/beta_grouping_hists.pdf"), plot = beta_grouping_hists,
       width = w*1.5, height = h, units = "mm", dpi = 300)

# All adonis barplots
figS3 = ggarrange(plots$beta_ado$adu_bray_loc, plots$beta_ado$adu_uniU_loc, plots$beta_ado$adu_uniW_loc, 
                  plots$beta_ado$chi_bray_loc, plots$beta_ado$chi_uniU_loc, plots$beta_ado$chi_uniW_loc, 
                  ncol=3, nrow=2, labels = "AUTO")
ggsave(paste0(outpath, "/figS3.pdf"), plot = figS3,
       width = w, height = h*0.4, units = "mm", dpi = 300)


in_group_plots = ggarrange(plots$beta_extra$in_group_bray, plots$beta_extra$in_group_unifracU, plots$beta_extra$in_group_unifracW,
                           nrow = 3)
ggsave(paste0(outpath, "/extras/in_group_plots.pdf"), plot = in_group_plots,
       width = w, height = h*0.8, units = "mm", dpi = 300)

# Mother-Child distance violin plots
figS4 = 
  (plots$beta_extra$adu_child_dist_bray | plots$beta_extra$relatedness_effect_bray) /
  (plots$beta_extra$adu_child_dist_unifracU | plots$beta_extra$relatedness_effect_unifracU) /
  (plots$beta_extra$adu_child_dist_unifracW |plots$beta_extra$relatedness_effect_unifracW) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = "bold", size=14))
ggsave(paste0(outpath, "/figS4.pdf"), plot = figS4,
       width = w, height = h*0.8, units = "mm", dpi = 300)

com_offsetX = -5
com_offsetY = 0
# Main figure mothers
alignedT = align_plots(plots$alpha$adults_chao, plots$beta_ado$adu_bray_loc, align = "h", axis = "t")
alignedB = align_plots(plots$alpha_reg$ent_coli, plots$beta$ent_coli_pcoa, align = "h", axis = "b")
alignedL = align_plots(alignedT[[1]], plots$alpha_reg$time_since_preg, alignedB[[1]], plots$maaslin$adults, 
                       align = "v", axis = "l")
alignedR = align_plots(alignedT[[2]], alignedB[[2]], alignedL[[4]], 
                       align = "v", axis = "r")
left = plot_grid(plotlist = alignedL[1:3], nrow=3, align="v", axis="lr", rel_heights = c(1.4,1,1.2), 
                 labels = c("A", "B", "C"), hjust = com_offsetX-4, vjust = com_offsetY+c(1.5,0,0))
right = plot_grid(plotlist = alignedR[1:2], nrow=2, align="v", axis="lr", labels= c("D", "E"),
                  hjust = com_offsetX, vjust = com_offsetY+c(1.5,0))
top = plot_grid(left, right, ncol=2)
full = plot_grid(top, alignedR[[3]], nrow = 2, rel_heights = c(1.9, 1), labels=c("", "F"),
                 hjust = com_offsetX-6, vjust = com_offsetY)
ggsave(paste0(outpath, "/fig2.pdf"), plot = full, 
       height = 0.8*h, width = w, units="mm", dpi = 300)

# Main figure infants
alignedT = align_plots(plots$alpha$infants_chao, plots$beta_ado$chi_bray_loc, align = "h", axis = "tb")
alignedB = align_plots(plots$beta$adu_chi_pcoa, plots$beta$chi_pcoa, align = "h", axis = "b")
alignedL = align_plots(alignedT[[1]], alignedB[[1]], plots$maaslin$infants, align = "v", axis = "l")
alignedR = align_plots(alignedT[[2]], alignedB[[2]], alignedL[[3]], align = "v", axis = "r")
left = plot_grid(plotlist = alignedL[1:2], nrow=2, align="v", axis="lr", rel_heights = c(1,1), 
                 labels = c("A", "C"), hjust = com_offsetX-4, vjust = com_offsetY+c(1.5,0))
right = plot_grid(plotlist = alignedR[1:2], nrow=2, align="v", axis="lr", labels= c("B", "D"),
                  hjust = com_offsetX, vjust = com_offsetY+c(1.5,0))
top = plot_grid(left, right, ncol=2)
size_legend = get_legend(plots$maaslin$infants+guides(fill="none"))
full = plot_grid(top, alignedR[[3]], nrow = 2, rel_heights = c(1.5, 1), labels=c("", "E"),
                 hjust = com_offsetX-6, vjust = com_offsetY)
ggsave(paste0(outpath, "/fig3.pdf"), plot = full, 
       height = 0.8*h, width = w, units="mm", dpi = 300)
