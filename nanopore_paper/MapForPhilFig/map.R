rm(list=ls())

# Load libraries
library(dplyr)
library(ggplot2)
library(rgdal)
library(ggrepel)
library(rgeos)
library(gridExtra)
library(gtable)

# Load functions
source("R/shp_to_df.R")

# Load data
locs <- read.csv("data/Sequencing_locations.csv", stringsAsFactors = FALSE)

# Load shapefiles
tz_outline <- readOGR("data/gis", "TZ_Outline", stringsAsFactors = FALSE)
tz_regs <- readOGR("data/gis", "TZ_Region_2012_pop", stringsAsFactors = FALSE)
ke_outline <- readOGR("data/gis", "KenyaOutline", stringsAsFactors = FALSE)
ke_counties <- readOGR("data/gis", "KenyaCounties", stringsAsFactors = FALSE)
ph_outline <- readOGR("data/gis", "PHL_Outline", stringsAsFactors = FALSE)
ph_province <- readOGR("data/gis", "PHL_province", stringsAsFactors = FALSE)

# Collect outline centroids
tz_centroid <- coordinates(gCentroid(tz_outline, byid=TRUE))
ke_centroid <- coordinates(gCentroid(ke_outline, byid=TRUE))
ph_centroid <- coordinates(gCentroid(ph_outline, byid=TRUE))

# Transform to dataframe for plotting
tz_outline_df <- shp_to_df(tz_outline)
tz_regs_df <- shp_to_df(tz_regs)
ke_outline_df <- shp_to_df(ke_outline)
ke_counties_df <- shp_to_df(ke_counties)
ph_outline_df <- shp_to_df(ph_outline)
ph_province_df <- shp_to_df(ph_province)

# Isolate locations for plotting
africa_points <- locs[which(locs$Country!="Philippines"),]
philippines_points <- locs[which(locs$Country=="Philippines"),]

# Plot African sites
map_af <- ggplot() +
  geom_polygon(data=tz_regs, aes(x=long, y=lat, group=group), fill="beige", color="dimgrey") +
  geom_polygon(data=tz_outline, aes(x=long, y=lat, group=group), color="#323232", fill=NA) +
  geom_polygon(data=ke_counties, aes(x=long, y=lat, group=group), fill="beige", color="dimgrey") +
  geom_polygon(data=ke_outline, aes(x=long, y=lat, group=group), color="#323232", fill=NA) +
  geom_point(data=africa_points, aes(x=long, y=lat, size=n_samples), color="dodgerblue") +
  scale_size_area(name="Number of samples", breaks=c(5,10,15,50), limits=c(0,55)) +
  geom_label_repel(data=africa_points[which(africa_points$Location %in% c("Mugumu District Field Office")),],
                   aes(x=long, y=lat, label=Location), label.size = NA,
                   nudge_y=0.5, nudge_x=-4, vjust="left", hjust="left") +
  geom_label_repel(data=africa_points[which(africa_points$Location %in% c("NMAIST")),],
                   aes(x=long, y=lat, label=Location), label.size = NA,
                   nudge_y=0.5, nudge_x=-2, vjust="left", hjust="left") +
  geom_label_repel(data=africa_points[which(africa_points$Location %in% c("Makeuni District Field Site", "TVLA", "UNITID")),],
                   aes(x=long, y=lat, label=Location), label.size = NA,
                   nudge_y=0.5, nudge_x=1, vjust="left", hjust="left") +
  annotate("label", x=tz_centroid[1], y=tz_centroid[2], label="TANZANIA", size=6,
           color="#323232", fill="beige", label.size = NA) +
  annotate("label", x=ke_centroid[1], y=ke_centroid[2], label="KENYA", size=6,
           color="#323232", fill="beige", label.size = NA) +
  theme_void() +
  theme(legend.position = "bottom") +
  coord_fixed()

# Plot Philippine site
map_ph <- ggplot() +
  geom_polygon(data=ph_province, aes(x=long, y=lat, group=group), fill="beige", color="dimgrey") +
  geom_polygon(data=ph_outline, aes(x=long, y=lat, group=group), color="#323232", fill=NA) +
  geom_point(data=philippines_points, aes(x=long, y=lat, size=n_samples), color="dodgerblue") +
  scale_size_area(name="Number of samples", breaks=c(5,10,15,50), limits=c(0,55)) +
  geom_label_repel(data=philippines_points,
                   aes(x=long, y=lat, label=Location), label.size = NA,
                   nudge_y=0.5, nudge_x=-2, vjust="left", hjust="left") +
  annotate("label", x=ph_centroid[1], y=ph_centroid[2], label="PHILIPPINES", size=6,
           color="#323232", fill="beige", label.size = NA) +
  theme_void() +
  theme(legend.position = "none") +
  coord_fixed()

# Extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

map_legend <- g_legend(map_af)
#tiff("figs/seq_map.tiff", res = 600)
pdf("figs/seq_map.pdf", height=8, width=11)
completed_map <- grid.arrange(arrangeGrob(map_af + theme(legend.position="none"),
                               map_ph, nrow=1),
                   map_legend, nrow=2,heights=c(10, 1))
dev.off()
