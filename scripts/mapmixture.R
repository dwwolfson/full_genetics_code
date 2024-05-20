# devtools::install_github("Tom-Jenkins/mapmixture")

# Load package
library(mapmixture)
library(here)
library(tidyverse)
library(rnaturalearth)
library(terra)
library(sf)
library(tidyterra)
library(gridExtra)

# # Read in admixture file format 1
# file <- system.file("extdata", "admixture1.csv", package = "mapmixture")
# admixture1 <- read.csv(file)
# 
# # Read in coordinates file
# file <- system.file("extdata", "coordinates.csv", package = "mapmixture")
# coordinates <- read.csv(file)
# 
# # Run mapmixture
# map1 <- mapmixture(admixture1, coordinates, crs = 3035)
#  map1

#df<-read_tsv(here("data/clumpak_output_3_4_2024/K=4/CLUMPP.files/ClumppIndFile.output"))

# @author Roman Lustrik {roman.lustrik@@biolitika.si}
readClumpp <- function(x) {
  x <- readLines(x)
  
  # extract individuals
  ind <- sapply(unname(sapply(x, function(x) strsplit(x = x, split = "\\("))), "[[", 1)
  ind <- strsplit(trimws(ind), "  ")
  ind <- sapply(ind, "[[", 1)
  
  # extract population
  # http://stackoverflow.com/questions/28267400/extract-a-string-of-words-between-two-specific-words-in-r
  pop <- sub(".*) *(.*?) *:.*", "\\1", x)
  
  # extract assignment values
  pis <- do.call("rbind", strsplit(unlist(lapply(strsplit(x, ":  "), "[[", 2)), " "))
  pis <- as.data.frame(apply(pis, MARGIN = 2, as.numeric))
  names(pis) <- paste("q", 1:ncol(pis), sep = "")
  
  cbind(ind, pop, pis)
}

df<-readClumpp(here("data/clumpak_output_3_4_2024/K=4/CLUMPP.files/ClumppIndFile.output"))

# write out structure output and add coordinates/sample ID
#write_csv(df, here("output/clumpp_output/k=4_3_4_2024/clump4.csv"))

# add in info from popmap used during structure runs
comp_popmap <- read_delim("D:/Data_Documents/!Swan_stuff/Genetics/popmaps/comp_popmap_states_flyways", 
                               delim = "\t", escape_double = FALSE, 
                               col_names = FALSE, trim_ws = TRUE)

names(comp_popmap)<-c("sample_ID", "state", "flyway")
df<-cbind(df, comp_popmap)

# code the arkansas swans as ON-MB (ontario-manitoba)
df<-df %>% 
  mutate(state=recode(state, AR="ON-MB"))


# there is only 1 ID sample, switch it to WY (n=17)
df<-df %>% 
  mutate(state=recode(state, ID="WY"))

my_admix_db<-data.frame(site=df$state, Ind=df$sample_ID, 
                        Cluster1=df$q1, Cluster2=df$q2, Cluster3=df$q3, Cluster4=df$q4)

my_coords<-read_csv(here("data/mapmixture_coordinates.csv"))

my_map<-mapmixture(my_admix_db, my_coords, crs=4326,
           cluster_cols = c("red", "green", "blue", "orange"),
           pie_size=3)+
  # Adjust theme options
  theme(
    legend.position = "top",
    plot.margin = margin(l = 10, r = 10),
  )+
  # Adjust the size of the legend keys
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1)))



# Use a raster as the basemap
earth<-terra::rast(here("data/spatial_data/NE1_50M_SR_W/NE1_50M_SR_W.tif"))

boundary<-c(xmin=-75, xmax=-160, ymin=35, ymax=70) |>  transform_bbox(bbox=_, 4326)
base_crop<-crop(earth, boundary)

rast_mapmix<-mapmixture(my_admix_db, my_coords, crs=4326, basemap=base_crop,
                        cluster_cols = c("red", "green", "blue", "orange"),
                        pie_size=3, 
                        arrow_size=4,
                        scalebar_size = 2,
                        axis_title_size = 16,
                        axis_text_size = 14,
                        boundary=c(xmin=-75, xmax=-160, ymin=35, ymax=70))+
  # Adjust theme options
  theme(
    legend.position = "top",
    plot.margin = margin(l = 10, r = 10),
  )+
  # Adjust the size of the legend keys
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1)))

# Add in structure plot

# more descriptive names
my_admix_db<-my_admix_db %>% 
  mutate(site=recode(site, "MN"="Minnesota",
                     "OH"="Ohio", 
                     "WI"="Wisconsin", 
                     "MI"="Michigan", 
                     "IA"="Iowa", 
                     "ON-MB"="Ontario", 
                     "NE"="Nebraska",
                     "WA"="Washington", 
                     "WY"="Wyoming",
                     "AL"="Alberta",
                     "AK"="Alaska"))


# Traditional structure barplot
structure_barplot <- structure_plot(
  admixture_df = my_admix_db,
  type = "structure",
  cluster_cols =  c("red", "green", "blue", "orange"),
  site_dividers = TRUE,
  divider_width = 0.8,
  site_order = c(
    c("Alaska", "Washington",
      "Alberta", "Wyoming",
      "Nebraska",
      "Minnesota", "Ontario", "Iowa", 
      "Wisconsin", "Michigan", "Ohio")
  ),
  labels = "site",
  flip_axis = F,
  site_ticks_size = -0.05,
  site_labels_y = -0.35,
  site_labels_size = 5
)+
  # Adjust theme options
  theme(
    #axis.text.x=element_text(angle=-25, vjust=0.5),
    axis.title.y = element_text(size = 12, hjust = 1),
    axis.text.y = element_text(size = 8),
  )

# Arrange plots
 grid.arrange(rast_mapmix, structure_barplot, nrow = 2, heights = c(4,1))

############################################################################################

 # parameters
 
 crs=4326
 
 #crs=9311
 # EPSG:9311
 # NAD27 / US National Atlas Equal Area
 # https://spatialreference.org/ref/epsg/9311/
 # 
 # Type: PROJECTED_CRS
 # WGS84 Bounds: 167.65, 15.56, -65.69, 74.71
 # Scope: Statistical analysis.
 # Area: United States (USA) - onshore and offshore.
 # Projection method name: Lambert Azimuthal Equal Area (Spherical)
 # Axes: Easting, Northing (X,Y). Directions: east, north. UoM: metre.
 # Base CRS: EPSG:4267
 
 boundary<-c(xmin=-75, xmax=-160, ymin=35, ymax=70) |>  transform_bbox(bbox=_, crs)
 
 # more descriptive names
 my_coords<-my_coords %>% 
   rename(.,site=state)
 
 my_coords<-my_coords %>% 
   mutate(site=recode(site, "MN"="Minnesota",
                      "OH"="Ohio", 
                      "WI"="Wisconsin", 
                      "MI"="Michigan", 
                      "IA"="Iowa", 
                      "ON-MB"="Ontario", 
                      "NE"="Nebraska",
                      "WA"="Washington", 
                      "WY"="Wyoming",
                      "AL"="Alberta",
                      "AK"="Alaska"))

 # Run mapmixture helper functions to prepare admixture and coordinates data
 admixture_df <- standardise_data(my_admix_db, type = "admixture") |> transform_admix_data(data = _)
 coords_df <- standardise_data(my_coords, type = "coordinates")
 admix_coords <- merge_coords_data(coords_df, admixture_df) |> transform_df_coords(df = _, crs = crs)



# import boundary lines
states<-ne_states(country = "United States of America", returnclass = "sf")
provinces<-ne_states(country="Canada", returnclass = "sf")


# pie_df<-left_join(my_admix_db, my_coords, by=join_by(site==state)) %>% 
#   select(-Ind) %>% 
#   select(site, Lat, Lon, Cluster1, Cluster2, Cluster3, Cluster4) %>% 
#   rename(lat=Lat, lon=Lon)

base_lines<-ggplot()+
  geom_spatraster_rgb(data=base_crop)+
  geom_sf(data=states, fill=NA)+
  geom_sf(data=provinces, fill=NA)+
  coord_sf(
    xlim = c(boundary[["xmin"]], boundary[["xmax"]]),
    ylim = c(boundary[["ymin"]], boundary[["ymax"]]), 
    expand=F
  )+
  add_pie_charts(admix_coords,
                 admix_columns = 4:ncol(admix_coords),
                 lat_column = "lat",
                 lon_column="lon",
                 pie_colours = c("red", "green", "blue", "orange"),
                 border=0.3,
                 opacity=1,
                 pie_size=2)


# Arrange plots
mapmix1<-grid.arrange(base_lines, structure_barplot, nrow = 2, heights = c(4,1))
ggsave(here("figures/mapmixture_base_boundaries_structure.tiff"), mapmix1,
       dpi=300, compression="lzw")


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

base_lines_lakes<-ggplot()+
  geom_spatraster_rgb(data=base_crop)+
  geom_sf(data=states, fill=NA)+
  geom_sf(data=provinces, fill=NA)+
  geom_sf(data=lakes, color="black", fill="lightblue")+
  coord_sf(
    xlim = c(boundary[["xmin"]], boundary[["xmax"]]),
    ylim = c(boundary[["ymin"]], boundary[["ymax"]]), 
    expand=F
  )+
  add_pie_charts(admix_coords,
                 admix_columns = 4:ncol(admix_coords),
                 lat_column = "lat",
                 lon_column="lon",
                 pie_colours = c("red", "green", "blue", "orange"),
                 border=0.3,
                 opacity=1,
                 pie_size=3)

# Arrange plots
mapmix2<-grid.arrange(base_lines_lakes, structure_barplot, nrow = 2, heights = c(4,1))
ggsave(here("figures/mapmixture_base_boundaries_structure_lakes.tiff"), mapmix2,
       dpi=300, compression="lzw")

