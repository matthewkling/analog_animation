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

prepRasters <- function(agg=8){
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
      #logtransform <- c(1, 2, 5) # for 7/27 biovars
      logtransform <- c(24, 22, 20, 13, 6, 4, 1) # for 27 biovars
      notransform <- setdiff(1:27, logtransform)
      ltrans <- function(x){
            x[x<1] <- 1
            log10(x)
      }
      b <- stack(b[[notransform]], calc(b[[logtransform]], ltrans))
      f <- stack(f[[notransform]], calc(f[[logtransform]], ltrans))
      
      b <- mask(b, f)
      f <- mask(f, b)
      
      d <- list(b, f)
      saveRDS(d, paste0("E:/analog_migration/animation/data_b26_agg", agg, ".rds"))
      return(d)
}


prepData <- function(clim, ext=NULL){
      if(!is.null(ext)){
            ext <- extent(ext)
            clim <- lapply(clim, crop, y=ext)
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


cluster <- function(vm, k=25){
      set.seed(123)
      trn <- vm[sample(nrow(vm), 10000),] # have to reduce size for kmeans to work
      trn <- kmeans(trn, k)$centers
      return(trn)
}

embed <- function(trn){
      dst <- dist(trn)
      set.seed(123)
      fit <- MASS::isoMDS(dst, k=3)
      embedded <- fit$points
      return(embedded)
}

sm <- function(y, clusters, embedded, outfile){
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

animate <- function(y, clusters, embedded, scheme, steps=100, interval=.1, 
                    outfile, outdir){
      
      colors <- colors3d(embedded, trans="ecdf", order=scheme[2], inversion=scheme[1])
      
      saveGIF({
            for(i in 0:steps){
                  
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
                  p <- ggplot(yi, aes(x, y)) +
                        geom_raster(fill=yi$color) +
                        geom_hline(yintercept=min(yi$y), size=8, alpha=.3) +
                        geom_point(x=min(yi$x) + (max(yi$x)-min(yi$x)) * i/steps,
                                   y=min(yi$y), size=10, shape=15, alpha=.3) +
                        annotate(geom="text", x=c(min(yi$x), max(yi$x)), y=min(yi$y), label=c("1975", "2080"),
                                 vjust=.4, color="white", size=8) +
                        ggmap::theme_nothing() +
                        coord_fixed(ratio=1)
                  plot(p)
                  
            }
            
      },
      movie.name = basename(outfile),
      interval = interval,
      ani.width = 1500,
      ani.height = 1125,
      outdir = dirname(outfile)
      )
      
      # move file from temp dir to intended destination, which buggy saveGIF fails at
      file.rename(paste0(tempdir(), "/", basename(outfile)),
                  paste0(dirname(outfile), "/", basename(outfile)))
}


###### execute

agg <- 1
nclust <- 25
steps <- 100
interval <- .1
outfile <- "E:/analog_migration/charts/CA_agg1.gif"

d <- prepRasters(agg)
d <- readRDS(paste0("E:/analog_migration/animation/data_b26_agg", agg, ".rds"))
#plot(d[[1]][[1]]); ext=drawExtent()

d <- prepData(d, ext=ext)
clusters <- cluster(d[[1]])
embedded <- embed(clusters)

sm(d[[2]], clusters, embedded, paste0("E:/analog_migration/charts/sm_", "test2", ".png"))
scheme <- c(4,2)

animate(d[[2]], clusters, embedded, scheme, steps=steps, interval=interval, outfile=outfile)


