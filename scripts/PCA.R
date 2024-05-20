# Calculate PCA

library(adegenet)
library(tidyverse)
library(vcfR)
library(factoextra)
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

# make allele frequency matrix
x.swans<-tab(gen, freq=T, NA.method="mean")

# This is PCA of the covariance matrix because it is centered and not scaled
pca.swans<-dudi.pca(df = x.swans, center = T, scale = F, scannf = FALSE, nf = 2)

pca_plot<-fviz_pca_ind(pca.swans, label="none", habillage=df$flyway, addEllipses = T, ellipse.level=0.95)+
  theme_pubr()+
  labs(x="PC1 (4.1%)", y="PC2 (1.8%)")+
  ggtitle("")

ggsave("figures/pca.tiff",
       dpi=300, compression="lzw")
