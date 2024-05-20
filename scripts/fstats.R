# F stats

library(vcfR)
library(tidyverse)
library(StAMPP)
library(gplots)
library(ade4)
library(adegenet)

# these can been 
vcfR<-read.vcfR("data/stacks_filtered_pre_LD_prune/populations.snps.vcf")

#convert to genlight
gen<-vcfR2genlight(vcfR)

# read in sample info
pops<-read.table(file="data/popmaps/comp_popmap_states_flyways", sep='\t', col.names=c("ID", "state", "flyway"))

#assign populations (a StaMPP requirement)
gen@pop<-as.factor(pops$flyway)

#generate flyway Fst matrix
#flyway_fst<-stamppFst(gen, nclusters=4)
flyway_fst<-loadRDS("output/Fst_output/flyways.Rda")

flyway_nj<-nj(flyway_fst$Fsts)

plot.phylo(flyway_nj, edge.width = 2, font=3, tip.color = "blue", cex=1.7, adj=0)
plot.phylo(flyway_nj, "cladogram")
plot.phylo(flyway_nj, "unrooted")
plot.phylo(flyway_nj, "radial")
plot.phylo(flyway_nj, "tidy")


# save out the Fst information to file
saveRDS(flyway_fst, "output/Fst_output/flyways.Rda")

# generate within-IP Fst info
# filter out non-IP 
# assign populations to states in the genlight object
# run stamppFst

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
     