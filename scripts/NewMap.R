##Original script by Eric Anderson

library(RColorBrewer)
library(raster)
library(tidyverse)
library(rgdal)
library(ggspatial)  # use GitHub version

####Map plotting
range <- readOGR("../data/shapefile/",layer="YWAR") # Shape file from NatureServe
breeding <- subset(range,SEASONAL==2)
breeding <- crop(breeding,extent(-170,-55,30,70))
wintering <- subset(range, SEASONAL == 3)

# Load some map layers and things
nat_earth <- stack("../data/maps/HYP_LR_SR_W_DR.tif")
ne_coast <- readOGR("../data/maps/ne_10m_coastline",
                    "ne_10m_coastline")
state_prov <- readOGR("../data/maps/ne_10m_admin_1_states_provinces_lines",
                      "ne_10m_admin_1_states_provinces_lines")

# and immediately drop all but canada and the US:
state_prov <- state_prov[state_prov$adm0_name %in% c("United States of America", "Canada"),]
country_bound <- readOGR("../data/maps/ne_10m_admin_0_boundary_lines_land",
                         "ne_10m_admin_0_boundary_lines_land")


# project and crop the natural earth raster
ne_crop <- crop(nat_earth, extent(-220, -20, 0, 90))  # way big so we can clip a rectangle out of it

lamproj <- "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#ne_crop_lam <- projectRaster(ne_crop, crs =  CRS(lamproj))
ne_crop_lam <- stack("../data/maps/map.tif")

# this puts the base map down as an image on the bottom
lambert_base <- ggplot() + 
  geom_spraster_rgb(ne_crop_lam, interpolate = TRUE) +
  coord_fixed()


# clip country and province bounds to put them in there
#library(stplanr)
country_bounds_clipped <- stplanr::gclip(shp = country_bound, bb = ne_crop)
state_prov_clipped <- stplanr::gclip(shp = state_prov, bb = ne_crop)


with_bounds <- lambert_base + 
  geom_spatial(country_bounds_clipped, crsto = as.CRS(ne_crop_lam), color = 'gray30', size = 0.2) +
  geom_spatial(state_prov_clipped, crsto = as.CRS(ne_crop_lam), color = 'gray30', size = 0.1) 

# get the kriged Q values as a raster stack These come from 01.4..
setwd("~/Documents/MigratoryBirds/RAD/YWAR/fluidigm/08.06.17")
stack <- readRDS("processed/K5breeding_locprior.rds")

# and the last thing we have to do is plot the rasters over the top with some transparency.  
# I have Clearly I have written my own function to map base color and the value in the 
# raster to three RGB values and alpha's for them.  
source("~/Documents/scripts/wifl-popgen-funcs.R")

# that function of mine is expecting values between 0 and 1 inclusive,
# so, I need to set any values in the rasters over 1 to 1.
stack[stack > 1.0] <- 1.0
stack[stack < 0.0] <- 0.0
crs(stack) <- as.CRS(country_bound)  # just give it lat-long
ws_proj <- projectRaster(stack, crs = CRS(lamproj))  # project it to lambert


my_cols <- c('#ff7f00','#377eb8','#e41a1c','#FF69B4','#984ea3','#ffff33') # orange, blue, red, pink,purple

rgba_list <- lapply(1:nlayers(ws_proj), function(i) {
  r <- ws_proj[[i]]
  mi <- cellStats(r, min)  # prepare to stretch everything so that the min value becomes zero and the max remains as it is
  ma <- cellStats(r, max)
  r2 <- (r - mi) * (ma / (ma - mi))
  rgba_brick(r = r2, color = my_cols[i])
})
names(rgba_list) <- names(ws_proj)



with_breeding_groups <- with_bounds +
  geom_spraster_rgb(rgba_list[[1]], interpolate = TRUE) +
  geom_spraster_rgb(rgba_list[[2]], interpolate = TRUE) +
  geom_spraster_rgb(rgba_list[[3]], interpolate = TRUE) +
  geom_spraster_rgb(rgba_list[[4]], interpolate = TRUE) +
  geom_spraster_rgb(rgba_list[[5]], interpolate = TRUE) +
theme_void()
with_breeding_groups


# now, let's put the wintering birds on there in color:
winter_birds <- readRDS("processed/Unknown_assignments_K2_12.19.18.rds") %>%
  mutate(projx = xyTransform(x = longitude, y = latitude, from = as.CRS(country_bound), to = as.CRS(ne_crop_lam))[,1],
         projy = xyTransform(x = longitude, y = latitude, from = as.CRS(country_bound), to = as.CRS(ne_crop_lam))[,2]) %>%
  filter(Stage == "Wintering") %>%
  filter(PofZ > 0.8)

###Wrong color assignments - here's a janky way to deal with it
wint_cols <- my_cols


# let's summarize each lat-long point with the number of birds from different groups
# and put them into different layers so that we can print the smaller ones on top
nums_at_lat_longs <- winter_birds %>% 
  group_by(latitude, longitude, projx, projy, collection) %>% 
  tally() %>%
  ungroup() %>%
  arrange(latitude, longitude, desc(n)) %>%
  group_by(latitude, longitude) %>%
  mutate(plot_order = 1:n()) %>%
  split(., .$plot_order)



#for (k in 1:5) {
with_unk <- with_breeding_groups +
  geom_spatial(wintering, crsto = as.CRS(ne_crop_lam), fill = "white", alpha = 0.35) + # gotta add the range here.
  geom_point(data = winter_birds, mapping = aes(x = projx, y = projy, fill = collection), size = 1, shape=21, stroke=0.1,
             position = position_jitter(width = 8e04, height = 8e04)) +
  scale_fill_manual(values = wint_cols)
#with_unk
filename = paste("figures/YWAR_map_pop",k,".pdf",sep='')
ggsave(with_unk,filename="YWAR_map_12.20.19.pdf",height=10,width=10)
#}
dev.copy(jpeg,"YWAR_map_12.20.19.jpg")
