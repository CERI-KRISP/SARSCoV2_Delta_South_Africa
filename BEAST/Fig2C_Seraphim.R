# install.packages("devtools"); library(devtools)
# install_github("sdellicour/seraphim/unix_OS") # for Unix systems
## install_github("sdellicour/seraphim/windows") # for Windows systems
# install.packages("diagram")

library(diagram)
require(seraphim)
library(rgdal)
library(sf)

loadstuff<- TRUE
if(loadstuff){

      # 1. Extracting the spatio-temporal information contained in posterior trees

      treefile<- "clusterA_11637_n286STEXP100ML.trees"

      localTreesDirectory = "Tree_extractions_new"
      allTrees = scan(file=treefile, what="", sep="\n", quiet=T)
      burnIn = 0
      randomSampling = FALSE
      nberOfTreesToSample = 100
      mostRecentSamplingDatum = 2021.621918
      coordinateAttributeName = "Location"

      treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName)


      # 2. Extracting the spatio-temporal information embedded in the MCC tree

      treefile<- "clusterA_11637_n286STEXP100ML.tree"

      source("mccExtractions.r")
      mcc_tre = readAnnotatedNexus(treefile)
      mcc_tab = mccTreeExtraction(mcc_tre, mostRecentSamplingDatum)
      write.csv(mcc_tab, "clusterA_MCC.csv", row.names=F, quote=F)

      # 3. Estimating the HPD region for each time slice

      nberOfExtractionFiles = nberOfTreesToSample
      prob = 0.95; precision = 0.025
      startDatum = min(mcc_tab[,"startYear"])

      polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))


      # 4.1 spatial boundaries

      reloadborders<-TRUE
      if(reloadborders){
          template_raster = raster("ZA5_Template_studyarea.asc")
          #template_raster = raster("zaf_pd_2020_1km.tif")
          borders = crop(getData("GADM", country="ZAF", level=1), extent(template_raster))
          #borders <- borders[borders$NAME_1 == "KwaZulu-Natal", ]
          #borders = getData("GADM", country="ZAF", level=1)
          save(borders,file="ZA_borders.Rdata")
      }
      load("ZA_borders.Rdata") #borders

}

#map south Africa
southAfrica <- st_read("gpkg/southafrica.gpkg")
districts <- st_read("gpkg/southafrica_district.gpkg")
muni <- st_read("gpkg/southafrica_seats.gpkg")


OrRd <- c(
  "#fff7ec",
  "#fee8c8",
  "#fdd49e",
  "#fdbb84",
  "#fc8d59",
  "#ef6548",
  "#d7301f",
  "#990000"
)

muni$polycolor <- NA
muni$polycolor[which(muni$Density > 0 & muni$Density <= 4)] <- OrRd[1]
muni$polycolor[which(muni$Density > 4 & muni$Density <= 8)] <- OrRd[2]
muni$polycolor[which(muni$Density > 8 & muni$Density <= 16)] <- OrRd[3]
muni$polycolor[which(muni$Density > 16 & muni$Density <= 32)] <- OrRd[4]
muni$polycolor[which(muni$Density > 32 & muni$Density <= 64)] <- OrRd[5]
muni$polycolor[which(muni$Density > 64 & muni$Density <= 128)] <- OrRd[6]
muni$polycolor[which(muni$Density > 128 & muni$Density <= 256)] <- OrRd[7]
muni$polycolor[which(muni$Density > 256)] <- OrRd[8]


# 4.2 Defining the different colour scales

minYear = min(mcc_tab[,"startYear"]); maxYear = max(mcc_tab[,"endYear"])
endYears_indices = (((mcc_tab[,"endYear"]-minYear)/(maxYear-minYear))*100)+1

##colors have to go to max time slice - AT LEAST
n_number_colours_needed<- max(round(endYears_indices))
n_repeats_discrete<- 10
c1<- rev((brewer.pal(4,"BuPu")))
c2<- (brewer.pal(9,"Reds"))
colours<- rev(rep(c(c1,c2), each=n_repeats_discrete))
colour_scale<- colorRampPalette(colours)(n_number_colours_needed)

endYears_colours = colour_scale[round(endYears_indices)]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons))
{
  date = as.numeric(names(polygons[[i]]))
  polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
  polygons_colours[i] = paste0(colour_scale[polygon_index],"55")
}

# 5. Generating the dispersal history plot

pdf('map_clusterA21.pdf',width=6, height=6.3,bg="white")

ptsize<- 0.8
pitjit<- 0.1
par(mar=c(0,0,0,0), oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
#plot(template_raster, col="grey85", box=F, axes=F, colNA="white", legend=F)
#plot(borders, add=T, lwd=0.1, border="gray10")
plot(southAfrica$geom, col="white", box=F, axes=F, colNA="grey90", legend=F)
for(i in 1:nrow(muni)){
  plot(muni$geom[i], axes=F, col=muni$polycolor[i], add=T, border=NA)
}
plot(districts$geom, add=T, lwd=0.1, border="gray10")
#plot(borders, lwd=0.1, border="gray10")
for (i in length(polygons):1)
{
  plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
}
for (i in 1:dim(mcc_tab)[1])
{
  # curvedarrow(cbind(mcc_tab[i,"startLon"],mcc_tab[i,"startLat"]), cbind(mcc_tab[i,"endLon"],mcc_tab[i,"endLat"]), arr.length=0,
              # arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
  curvedarrow(cbind(mcc_tab[i,"startLon"],mcc_tab[i,"startLat"]), cbind(mcc_tab[i,"endLon"],mcc_tab[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=2.5*1.1, lty=1, lcol="grey22", arr.col=NA, arr.pos=FALSE, curve=0.3, dr=NA, endhead=F)
  curvedarrow(cbind(mcc_tab[i,"startLon"],mcc_tab[i,"startLat"]), cbind(mcc_tab[i,"endLon"],mcc_tab[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=2*1.1, lty=1, lcol=endYears_colours[i], arr.col=NA, arr.pos=FALSE, curve=0.3, dr=NA, endhead=F)
}
for (i in dim(mcc_tab)[1]:1)
{
  xs<- mcc_tab[i,"startLon"]
  ys<- mcc_tab[i,"startLat"]
  xe<- jitter(mcc_tab[i,"endLon"],pitjit)
  ye<- jitter(mcc_tab[i,"endLat"],pitjit)
  if (i == 1)
  {
    points(xs, ys, pch=16, col=colour_scale[1], cex=ptsize)
    points(xs, ys, pch=1, col="gray10", cex=ptsize)
  }
  points(xe, ye, pch=16, col=endYears_colours[i], cex=ptsize)
  points(xe, ye, pch=1, col="gray10", cex=ptsize)
}

xrange<- c(xmin(template_raster), xmax(template_raster))
yrange<- c(ymin(template_raster), ymax(template_raster))
#rect(xrange[1], yrange[1], xrange[2], yrange[2], xpd=T, lwd=0.2)
#axis(1, c(ceiling(xmin(template_raster)), floor(xmax(template_raster))), pos=ymin(template_raster), mgp=c(0,0.2,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
#axis(2, c(ceiling(ymin(template_raster)), floor(ymax(template_raster))), pos=xmin(template_raster), mgp=c(0,0.5,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc_tab[,"startYear"]); rast[2] = max(mcc_tab[,"endYear"])
plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.40,0.80,0.14,0.155),
     legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.5, col.axis="gray30", line=0, mgp=c(0,-0.02,0)))

legend("bottomleft", 
       legend = c("0-4",
                  "4-8",
                  "8-16",
                  "16-32",
                  "32-64",
                  "64-128",
                  "128 - 256",
                  "> 256"), 
       fill = OrRd, 
       bty = "n",
       cex = .6,
       border = "black")

a<-dev.off()
