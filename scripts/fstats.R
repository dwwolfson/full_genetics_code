# F stats

library(vcfR)
library(tidyverse)
library(StAMPP)
library(gplots)
library(ade4)
library(adegenet)
library(dartR)
library(hierfstat)
library(ggpubr)

# these can been 
vcfR<-read.vcfR("data/stacks_filtered_pre_LD_prune/populations.snps.vcf")

#convert to genlight
gen<-vcfR2genlight(vcfR)

# read in sample info
pops<-read.table(file="data/popmaps/comp_popmap_states_flyways", sep='\t', col.names=c("ID", "state", "flyway"))

#assign populations (a StaMPP requirement)
gen@pop<-as.factor(pops$flyway)

#generate flyway Fst matrix
flyway_fst<-stamppFst(gen, nclusters=4)

# generate flyway Nei distance measures
flyway_nei<-stamppNeisD(gen)
saveRDS(flyway_nei, "output/Nei/flyways_nei.Rda")
flyway_nei<-readRDS("output/Nei/flyways_nei.Rda")

# save out the Fst information to file
saveRDS(flyway_fst, "output/Fst_output/flyways.Rda")
flyway_fst<-readRDS("output/Fst_output/flyways.Rda")

flyway_nj<-nj(flyway_fst$Fsts)
flyway_nj_nei<-nj(flyway_nei)

# Fst based
plot.phylo(flyway_nj, edge.width = 2, font=3, tip.color = "blue", cex=1.7, adj=0)
plot.phylo(flyway_nj, "cladogram")
plot.phylo(flyway_nj, "unrooted", edge.width = 2, font=3, tip.color = "blue", cex=1.7, adj=0)
plot.phylo(flyway_nj, "radial")
plot.phylo(flyway_nj, "tidy")

# Nei based
plot.phylo(flyway_nj_nei, edge.width = 2, font=3, tip.color = "blue", cex=1.7, adj=0)

# Individual-based phylo stuff
gen@pop<-as.factor(pops$ID)
ind_nei<-stamppNeisD(gen)
ind_nei_nj<-nj(ind_nei)

ind_phylo_plot<-plot.phylo(ind_nei_nj)
col_df<-data.frame(ID=ind_nei_nj$tip.label, flyway=pops$flyway)
col_df<-col_df %>% 
  mutate(color=ifelse(grepl("IP", flyway), "green",
         ifelse(grepl("HP", flyway), "orange",
         ifelse(grepl("RMP", flyway), "purple", 
         ifelse(grepl("PCP", flyway), "blue", 'flag')))))

plot.phylo(ind_nei_nj, edge.width = 1, font=0.1, tip.color = col_df$color, cex=1)
plot.phylo(ind_nei_nj, "fan", edge.width = 1, font=0.1, tip.color = col_df$color, cex=1)

# export phylip format
stamppPhylip(ind_nei, "output/Nei/phylip_individuals/NeisD.txt")
stamppPhylip(ind_nei, "output/Nei/phylip_individuals/NeisD.phy")

#extract the pairwise matrix
m<-flyway_fst$Fsts
#fill in upper triangle of the matrix
m[upper.tri(m)] <- t(m)[upper.tri(m)]

#melt to tidy format for ggplotting
heat <- reshape2::melt(m)

#plot as heatmap with exact values labeling each cell
fst_plot<-ggplot(data = heat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  geom_text(data=heat,aes(label=round(value, 3)))+
  theme_minimal()+
  scale_fill_gradient2(low = "white", high = "red", space = "Lab", name="Fst") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face='bold', size=14),
        axis.text.y = element_text(angle = 45, hjust = 1, face='bold', size=14))+
  labs(x="", y="")

ggsave("figures/Fst_pops.tiff",
       compression="lzw", dpi=300, bg="white")
saveRDS(fst_plot, "output/fst_ggplot.Rdata")
     
# Alternatively, add p-values to the other off-diagonals
m1<-m
m1[upper.tri(m1)]<-0  #all p-values were 0
heat1<-reshape2::melt(m1)

fst_plot1<-ggplot(data = heat1, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  geom_text(data=heat1,aes(label=round(value, 3)))+
  theme_minimal()+
  scale_fill_gradient2(low = "white", high = "red", space = "Lab", name="Fst") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face='bold', size=14),
        axis.text.y = element_text(angle = 45, hjust = 1, face='bold', size=14))+
  labs(x="", y="")

ggsave("figures/Fst_pvals.png",bg="white",
       dpi=300)
################################################################


# generate within-IP Fst info

# assign populations to states in the genlight object
gen@pop<-as.factor(pops$state)

# run stamppFst
state_fst<-stamppFst(gen, nclusters=4)

# filter out non-IP 
m<-state_fst$Fsts
m[upper.tri(m)] <- t(m)[upper.tri(m)]
m1<-reshape2::melt(m)
m1<-m1 %>% 
  filter(Var2 %in%c('MN', 'OH', 'WI', 'IA', 'MI', 'AR', 'NE')) %>% 
  filter(Var1 %in%c('MN', 'OH', 'WI', 'IA', 'MI', 'AR', 'NE'))

fst_states<-ggplot(data = m1, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  geom_text(data=m1,aes(label=round(value, 3)))+
  theme_minimal()+
  scale_fill_gradient2(low = "white", high = "red", space = "Lab", name="Fst") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face='bold', size=14),
        axis.text.y = element_text(angle = 45, hjust = 1, face='bold', size=14))+
  labs(x="", y="")

ggsave("figures/Fst_states.png", bg="white", dpi=300)



#################################################
 # estimate Fis inbreeding coefficients

# switch back to flyway assignment for pops
gen@pop<-as.factor(pops$flyway)

# switch from genlight to genind object
flyway_genind<-gl2gi(gen)

flyway_Fis<-inbreeding(flyway_genind)

#write to file
saveRDS(flyway_Fis, "output/Fis_output.Rda")

Fbar<-sapply(flyway_Fis, mean)
# not what I'd expect

##
# test for HWE (but this is global, so not sure this makes sense to filter on)
obj.test<-pegas::hw.test(flyway_genind)
hw_df<-as.data.frame(obj.test)
hw_df$locus<-rownames(hw_df)


# write out to file
saveRDS(obj.test, "output/pegas_HW/hw.test.output.Rda")

# now run summary with adegenet on the genind object
adegenet_summary_flyway<-summary(flyway_genind)
# I think this isnt' what I want because it gives a stat for each locus?

# write out to file
saveRDS(adegenet_summary_flyway, "output/Fst_output/adegenet_summary.Rda")
adegenet_summary_flyway<-readRDS("output/Fst_output/adegenet_summary.Rda")


# try hierfstat approach as well to see if numbers agree
fstats_hierf<-basic.stats(flyway_genind, digits=3)

# write to file
saveRDS(fstats_hierf, "output/hierfstat/fstats.Rda")

# Fis
Fis_flyway<-apply(fstats_hierf$Fis, MARGIN=2, FUN=mean, na.rm=T) %>% 
  round(digits=2)
fis_inb<-as.data.frame(fstats_hierf$Fis)
fis_inb<-reshape2::melt(fis_inb)
fis_inb %>% 
  ggplot(aes(x=variable, y=value))+geom_boxplot()



# Average observed hetero per site
Ho_flyway<-apply(fstats_hierf$Ho, MARGIN=2, FUN=mean, na.rm=T) %>% 
  round(digits=2)
obs_het<-as.data.frame(basic_stat$Ho)
obs_het<-reshape2::melt(obs_het)
obs_het %>% 
  ggplot(aes(x=variable, y=value))+geom_boxplot()+ggtitle("Observed Heterozygosity")

# Expected hetero
He_flyway<-apply(fstats_hierf$Hs, MARGIN=2, FUN=mean, na.rm=T) %>% 
  round(digits=2)
exp_het<-as.data.frame(basic_stat$Hs)
exp_het<-reshape2::melt(exp_het)
exp_het %>% 
  ggplot(aes(x=variable, y=value))+geom_boxplot()+ggtitle("Expected Heterozygosity")

# Combine together
obs_het$metric="Observed Heterozygosity"
exp_het$metric="Expected Heterozygosity"

hets<-rbind.data.frame(obs_het, exp_het)
hets %>% 
  ggplot(aes(x=variable, y=value, color=metric))+
  geom_boxplot()+
  theme_bw()
  



hetero_plot<-hets %>% 
  ggplot(aes(x=variable, y=value, fill=metric))+
  geom_boxplot()+
  theme_pubr()+
  geom_pwc(method="wilcox_test",
           label="p.adj.format",
           p.adjust.method="bonferroni")+
  labs(x="Groups", y="Heterozygosity per Locus\n", fill="Metric")

ggsave( "figures/heterozygosity.tiff",
       dpi=300, compression="lzw")
ggsave( "figures/heterozygosity.png",
       dpi=300)

# If we wanted to compare across populations
my_comparisons<-list(c("HP", "IP"), 
                     c("HP", "PCP"),
                     c("HP", "RMP"),
                     c("IP", "PCP"),
                     c("IP", "RMP"),
                     c("PCP", "RMP"))
p<-hets %>% 
  ggplot(aes(x=variable, y=value, color=metric))+
  geom_boxplot()+
  theme_pubr()+
  stat_compare_means(comparisons=my_comparisons)

ggadjust_pvalue(p, p.adjust.method="bonferroni")

# I think the below is if you want to consider all the panels at the same time
# bxp<-ggboxplot(
#   hets, x='metric', y='value', color='metric',
#   facet.by = 'variable'
# )
# bxp<-bxp+geom_pwc(method="wilcox_test")
# ggadjust_pvalue(
#   bxp, p.adjust.method="bonferroni",
#   label = "{p.adj.format}{p.adj.signif}", hide.ns = TRUE
# )


# Fis
Fis_df<-as.data.frame(basic_stat$Fis)
Fis_df<-reshape2::melt(Fis_df)
Fis_df %>% 
  ggplot(aes(x=variable, y=value))+geom_boxplot()+ggtitle("FisInbreeding")

# dartR genetic diversity metrics
gen<-gl.compliance.check(gen)

dartr_diversity<-gl.report.heterozygosity(
  gen,
  method = "pop",
  plot.out = TRUE,
  save2tmp = T,
  verbose = 5
)

# save to file
write_csv(dartr_diversity, "output/dartR/diversity_metrics.csv")


