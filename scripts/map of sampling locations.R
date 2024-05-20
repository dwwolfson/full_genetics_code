# deraad style map

library(tidyverse)
#remotes::install_github("wmgeolab/rgeoboundaries")
library(rgeoboundaries)
library(elevatr)
library(raster)
library(RColorBrewer)
library(gridExtra)
library(sf)
library(mapmixture)
library(tidyterra)
library(ggspatial)
library(rnaturalearth)

# read in location data
df<-read_csv("data/master_genetic_samples.csv")

df$Latitude<-round(df$Latitude, 3)
df$Longitude<-round(df$Longitude, 3)

state_dfs<-split(df, df$Flyway)
sample_df<-data.frame(NULL)
for (i in names(state_dfs)){
  samps<-state_dfs[[i]] %>% dplyr::group_by(Latitude, Longitude) %>% dplyr::summarize(count=dplyr::n())
  dat<-cbind(rep(i, times=nrow(samps)), samps)
  sample_df<-as.data.frame(rbind(sample_df, dat))
}

colnames(sample_df)[1]<-'Flyway'


# get worldwide data for making make map
#pac<-map_data("world")

# pull in raster to use as basemap
earth<-terra::rast("data/spatial_data/NE1_50M_SR_W/NE1_50M_SR_W.tif")

boundary<-c(xmin=-77, xmax=-155, ymin=34, ymax=64) |>  transform_bbox(bbox=_, 4326)
base_crop<-crop(earth, boundary)


# import boundary lines
states<-ne_states(country = "United States of America", returnclass = "sf")
provinces<-ne_states(country="Canada", returnclass = "sf")

# add in lakes
lakes <- rnaturalearth::ne_download(type = 'lakes', 
                                    scale=10,
                                    category = 'physical') %>% 
  sf::st_as_sf(lakes110, crs = 4263) %>% 
  select(name) %>% 
  filter(name == 'Lake Huron' 
         | name == 'Lake Ontario'
         | name == 'Lake Michigan'
         | name == 'Lake Erie'
         | name == 'Lake Superior'
         | name == 'Great Bear Lake'
         | name == 'Great Slave Lake'
         | name == 'Lake Winnipeg')


# make comprehensive map

# should write out the flyway labels instead of abbreviations
sample_df<-sample_df %>% 
  mutate(Flyway=recode(.$Flyway, "IP" = "Interior Population",
                       "PCP" = "Pacific Coast Population",
                       "RMP" = "Rocky Mountain Population"))


sampling_map<-ggplot()+
  layer_spatial(base_crop, alpha=0.7)+
  geom_sf(data=states, fill=NA)+
  geom_sf(data=provinces, fill=NA)+
  geom_sf(data=lakes, color="black", fill="lightblue")+
  #geom_polygon(data=pac, aes(x=long, y=lat, group=group), fill=NA, col="black", cex=.1)+ # we aren't using this anymore
  theme_classic()+
  scale_size_binned(breaks=c(1,2,3,4,8,12,14))+
  geom_point(data=sample_df, aes(x=Longitude, y=Latitude, fill=Flyway, size=count),
             alpha=0.8, pch=21, colour="black")+
  scale_color_manual(values=brewer.pal(3, "Set2"))+
  guides(fill=guide_legend(override.aes = list(size=4)))+
  # geom_rect(aes(xmin = -103, xmax = -78, ymin = 34, ymax = 51), #IP
  #           lwd=1.5, color=brewer.pal(3, "Set2")[1], fill=NA)+ 
  # geom_rect(aes(xmin = -155, xmax = -118, ymin = 45, ymax = 65), #PCP
  #           lwd=1.5, color=brewer.pal(3, "Set2")[2], fill=NA)+  
  # geom_rect(aes(xmin = -121, xmax = -106, ymin = 40, ymax = 58), #RMP
  #           lwd=1.5, color=brewer.pal(3, "Set2")[3], fill=NA)+
  coord_sf(
  xlim = c(boundary[["xmin"]], boundary[["xmax"]]),
  ylim = c(boundary[["ymin"]], boundary[["ymax"]]), 
  expand=F)+
  scale_x_continuous(label=I)+
  scale_y_continuous(label=I)+
  labs(x="\nLongitude", y="Latitude\n", size="Number of\nsamples")+
  theme(axis.text.x = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=16, face="bold"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = c(0.15, 0.3))

ggsave("figures/samples_map.tiff", sampling_map,
       dpi=300, compression="lzw")
  