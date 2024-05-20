# Join Burke museum data to know date range of sample collection

library(tidyverse)
df<-read_csv("data/Burke_samples_received.csv")
dat<-read_csv("data/idigbio_uwbm_output.csv")

df<-df %>% 
  left_join(., dat,
            by=c("UWBM"="dwc:catalogNumber")) %>% 
  select(UWBM:sex, 'dwc:coordinateUncertaintyInMeters', 
         'dwc:county', 'idigbio:eventDate',
         'dwc:stateProvince',
         'idigbio:geoPoint')

write_csv(df, "data/joined_Burke_sample_data.csv")
