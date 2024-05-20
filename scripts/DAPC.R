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
grp<-find.clusters(gen, max.n.clust = 9)
#write output to file
saveRDS(grp, "output/find.cluster_grp.Rda")

 dapc(gen, grp$grp)
 
 
 dapc.sub<-dapc(d1, grp$grp)
 # 150 pc's saved and 3 discriminant functions saved
 
 pdf(file="figures/dapc_5k_loci_150pcs.pdf", width = 8, height = 8)
 scatter(dapc.sub, scree.da=F, bg="white", pch=20, cell=0,
         cstar=0, solid=0.4, clab=0, cex=3, leg=T,
         txt.leg = c("High Plains Flock",
                     "Interior Population",
                     "Rocky Mountain Population",
                     "Pacific Population"))
 myInset <- function(){
   temp <- dapc.sub$pca.eig
   temp <- 100* cumsum(temp)/sum(temp)
   plot(temp, col=rep(c("black","lightgrey"),
                      c(dapc.sub$n.pca,1000)), ylim=c(0,100),
        xlab="PCA axis", ylab="Cumulated variance (%)",
        cex=1, pch=20, type="h", lwd=2)
 }
 add.scatter(myInset(), posi="topleft",
             ratio=.2,
             bg=transp("white"))
 dev.off()