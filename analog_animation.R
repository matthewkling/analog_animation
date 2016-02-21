library(ecoclim)
library(raster)
library(dplyr)
library(tidyr)
library(FNN)
library(ggplot2)
library(animation)
library(tsne)
library(colormap)

select <- dplyr::select

prepRasters <- function(outfile, agg=8){
      # function loads, crops, aggregates, log-transforms, and re-saves climateNA raster stacks for two time periods
      
      d <- parseMetadata("F:/ClimateNA/1km_Hamman_raw/27_biovars")
      d <- d[!grepl("MAR", d$path),] # this variable doesn't cover mexico, the rest do
      b <- stack(d$path[grepl("6190", d$path)])
      #f <- stack(d$path[grepl("2020", d$path)])
      f <- stack(d$path[grepl("2080", d$path)])
      
      ext <- c(-2723106, 2522840, 3426221, 7149212) # for 27 bio
      b <- crop(b, ext)
      f <- crop(f, ext)
      
      if(agg>1){
            b <- aggregate(b, agg)
            f <- aggregate(f, agg)
      }
      
      # log transform ppt, and round up logs of <1 mm ppt
      logtransform <- c(24, 22, 20, 13, 6, 4, 1) # for 27 biovars
      notransform <- setdiff(1:26, logtransform)
      ltrans <- function(x){
            x[x<1] <- 1
            log10(x)
      }
      b <- stack(b[[notransform]], calc(b[[logtransform]], ltrans))
      f <- stack(f[[notransform]], calc(f[[logtransform]], ltrans))
      
      b <- mask(b, f)
      f <- mask(f, b)
      
      d <- list(b, f)
      saveRDS(d, outfile)
      return(d)
}


prepData <- function(clim, ext=NULL, clip=NULL){
      # function converts a list of two raster stacks into formatted data frames
      
      if(!is.null(ext)){
            ext <- extent(ext)
            clim <- lapply(clim, crop, y=ext)
      }
      if(!is.null(clip)){
            clim <- lapply(clim, mask, mask=clip)
      }
      
      clim <- lapply(clim, function(x) as.data.frame(rasterToPoints(x)))
      clim <- lapply(clim, function(x){
            names(x) <- c("x", "y",  paste0("v", 1:(ncol(x)-2)))
            return(x)})
      b <- clim[[1]]; f <- clim[[2]]
      b$year <- "baseline"
      f$year <- "future"
      v <- na.omit(rbind(b, f))
      v <- v %>% mutate_each(funs(scale), -x, -y, -year)
      vm <- as.matrix(dplyr::select(v, -x, -y, -year))
      
      y <- v %>%
            gather(variable, value, -year, -x, -y) %>%
            spread(year, value) %>%
            mutate(delta=future-baseline)
      
      return(list(vm, y))
}


cluster <- function(data, k=25){
      # function returns kmeans cluster centers of data matrix
      set.seed(123)
      trn <- data[sample(nrow(data), 10000),] # have to reduce size for kmeans to work
      trn <- kmeans(trn, k, iter.max=50)$centers
      return(trn)
}

embed <- function(data){
      # function uses NMDS to embed n-dimensional data into 3-d
      dst <- dist(data)
      set.seed(123)
      fit <- MASS::isoMDS(dst, k=3)
      embedded <- fit$points
      return(embedded)
}

sm <- function(y, clusters, embedded, outfile){
      # function saves a small multiples menu of color palettes to show options for RGB mapping
      yi <- y %>%
            na.omit() %>%
            mutate(value = baseline + delta * .5) %>%
            select(x, y, variable, value) %>%
            spread(variable, value)
      yim <- yi %>%
            select(-x, -y) %>%
            as.matrix()
      nn <- get.knnx(clusters, yim, k=1)
      for(o in 1:6){
            for(i in 1:8){
                  yj <- yi
                  colors <- colors3d(embedded, trans="ecdf", order=o, inversion=i)
                  yj$color <- colors[as.vector(nn$nn.index)]
                  yj$order <- o
                  yj$inversion <- i
                  if(o==1 & i==1) sm <- yj else(sm <- rbind(sm, yj))
            }
      }
      xstep <- (max(range(sm$x))-min(range(sm$x))) * 1.25
      ystep <- (max(range(sm$y))-min(range(sm$y))) * 1.25
      sm$xx <- sm$x + sm$inversion * xstep
      sm$yy <- sm$y + sm$order * ystep
      p <- ggplot(sm, aes(xx, yy)) +
            geom_raster(fill=sm$color) +
            ggmap::theme_nothing() +
            coord_fixed(ratio=1)
      ggsave(outfile, p, width=20, height=15, units="in")
}

animate <- function(y, clusters, embedded, scheme, steps=100, interval=.1, outfile, tag){
      # function generates an animation of climate analog migration
      
      colors <- colors3d(embedded, trans="ecdf", order=scheme[2], inversion=scheme[1])
      
      oopts = if (.Platform$OS.type == "windows") {
            ani.options(ffmpeg = "C:/Program Files/ffmpeg-20160219-git-98a0053-win64-static/bin/ffmpeg.exe")
      }
      
      source("E:/analog_migration/code/saveVid.R")
      
      saveVideo(expr={
            
            # cluster stetup
            cpus <- 7
            cl <- makeCluster(cpus)
            registerDoParallel(cl)
            
            td <- tempdir()
            
            r <- foreach(i = 0:steps,
                         .packages=c("ggplot2", "dplyr", "tidyr", "FNN")) %dopar% {
                               
                               # interpolate values
                               yi <- y %>%
                                     na.omit() %>%
                                     mutate(value = baseline + delta * i/steps) %>%
                                     select(x, y, variable, value) %>%
                                     spread(variable, value)
                               
                               # assign color of NMDS training values to entire dataset based on nearest neighbors
                               nn <- get.knnx(clusters, as.matrix(select(yi, -x, -y)), k=1)
                               yi$color <- colors[as.vector(nn$nn.index)]
                               
                               # plots
                               #barlat <- min(yi$y) - (max(yi$y) - min(yi$y)) /20
                               p <- ggplot(yi, aes(x, y)) +
                                     geom_raster(fill=yi$color) +
                                     #geom_hline(yintercept=barlat, size=8, alpha=.3) +
                                     #geom_point(x=min(yi$x) + (max(yi$x)-min(yi$x)) * i/steps,
                                     #           y=barlat, size=10, shape=15, alpha=.3) +
                                     #annotate(geom="text", x=c(min(yi$x), max(yi$x)), y=barlat, 
                                     #         label=c("1975", "2080"),
                                     #         vjust=.4, color="white", size=8) +
                                     ggmap::theme_nothing() +
                                     scale_x_continuous(expand=c(0,0)) +
                                     scale_y_continuous(expand=c(0,0)) +
                                     coord_fixed(ratio=1)
                               #p <- ggplot(data.frame(x=1,y=1), aes(x,y)) + geom_text(label=i)
                               
                               png(paste0(td, "/Rplot", i+2, ".png"), width=1500, height=938, units="px")
                               plot(p)
                               dev.off()
                         }
            stopCluster(cl)
            
      },
      tmpdr = tempdir(),
      video.name = paste0(basename(outfile), ".mp4"),
      interval = .05, ani.width=1500, ani.height=938,
      other.opts = "-pix_fmt yuv420p -b 4000k" # these parameters set the pixel color format and the bitrate
      )
      
      ani.options(oopts)
}



end2end <- function(dataset="climateNA26", agg=8, 
                    ext=NULL, clip=NULL,
                    nclust=25, 
                    x=NULL, y=NULL,
                    steps=100, interval=.1,
                    outdir="E:/analog_migration/output", tag){
      
      # raster data
      raster_file <- paste0(outdir, "/data_", dataset, "_agg", agg, ".rds")
      if(file.exists(raster_file)){
            d <- readRDS(raster_file)
      } else{
            d <- prepRasters(raster_file, agg)
      }
      
      # embedded data
      #plot(d[[1]][[1]]); ext=drawExtent()
      d <- prepData(d, ext=ext, clip=clip)
      clusters <- cluster(d[[1]], k=nclust)
      embedded <- embed(clusters)
      
      # graphics
      if(is.null(x) | is.null(y)){
            sm(d[[2]], clusters, embedded, paste0(outdir, "/sm_", dataset, "_agg", agg,"_n", nclust, "_", tag, ".png"))
            x <- readline("Enter column number of selected palette:")
            y <- readline("Enter row number of selected palette:")
      }
      scheme <- c(as.integer(x),7-as.integer(y))
      outfile <- paste0(outdir, "/animation_", dataset, "_agg", agg,"_n", nclust, "_", tag)
      animate(d[[2]], clusters, embedded, scheme, steps=steps, interval=interval, outfile=outfile)
}

###### execute

#CA <- rgdal::readOGR("E:/flow/velocity_wind_alignment/data/shapefiles", "cb_2014_us_state_500k")
#CA <- CA[CA$STUSPS=="CA",]
#prj <- crs(raster("E:/analog_migration/ClimateNA_Reference/ClimateNA_DEM.asc"))
#CA <- spTransform(CA, prj)
#ext=extent(-2625008, -1557861, 4614009, 6645031)


#end2end(agg=2, nclust=50, tag="CA", clip=CA)
end2end(agg=8, nclust=50, tag="USA", x=5, y=5, steps=100)



