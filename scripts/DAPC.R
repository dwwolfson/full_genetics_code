# Calculate DAPC

library(adegenet)
library(tidyverse)
library(vcfR)
library(ggpubr)

# import the vcf that has been filtered in stacks and also snpfiltr (although only minimally in snpfiltr)
my_vcf<-read.vcfR("output/SNPfiltR/distance_thin.snps.vcf.gz")

# convert vcf file to genlight format
gen<-vcfR2genlight(my_vcf)

pops<-read.table(file="data/popmaps/comp_popmap_states_flyways", sep='\t', col.names=c("ID", "state", "flyway"))

# combine flyway info with indidual sample IDs
df<-data.frame(ID=indNames(gen))
df<-left_join(df, pops)

# Flyway PCA
strata(gen)<-df
setPop(gen)<-~flyway # set population name to flyway

# k-means to find clusters
grp<-find.clusters(gen, max.n.clust = 9) #this took a long time to run
#saved 200 PCs and a cluster number of 4

#write output to file
saveRDS(grp, "output/find.cluster_grp.Rda")

my_dapc<-dapc(gen, grp$grp)  #this took a long time to run
# saved 200 PCs and 3 discriminant functions
# save to file
saveRDS(my_dapc, "output/DAPC/dapc_4_clusters.Rda")
 

# amount of variation explained by the 3 discriminant functions
percent= my_dapc$eig/sum(my_dapc$eig)*100
barplot(percent, ylab="Percent of genetic variance explained by eigenvectors", 
        names.arg=round(percent,2))


dapc_df<-as.data.frame(my_dapc$ind.coord)
dapc_df$ID<-row.names(dapc_df) 
dapc_df<-left_join(dapc_df, pops)

ggplot(dapc_df, aes(x=LD1, y=LD2, fill=flyway))+ geom_point(size=2, pch=21)+
  labs(x="DPC1 (60.8%)",y="DPC2 (23.9%)", fill="Groups")+
  theme_pubr()+
  geom_hline(yintercept=0, linetype="dashed", col="black", linewidth=1)+
  geom_vline(xintercept=0, linetype="dashed", col="black", linewidth=1)+
  theme(axis.text.x = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
                legend.title = element_text(size=16),
                legend.text = element_text(size=16))


ggsave("figures/dapc_flyways.png",
       dpi=300)


ggsave("figures/dapc_flyways.tiff",
       dpi=300, compression="lzw")

####################
 # run again and use the lowest BIC for 2 clusters instead
grp1<-find.clusters(gen, max.n.clust = 2) 
# 200 pc's retained and a cluster size of 2
 
# write to file
saveRDS(grp1, "output/find.cluster_grp_2_groups.Rda")

my_dapc_2_clust<-dapc(gen, grp1$grp) 
# retained 150 PCs
# only 1 eigenvalue
scatter(my_dapc_2_clust)
assignplot

table(df$flyway, grp1$grp)
table(df$flyway, grp$grp)

#write to file
saveRDS(my_dapc_2_clust, "output/dapc_2_clust.Rda")
