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

raster_recalculate <- F
if(raster_recalculate){
      
      d <- parseMetadata("F:/ClimateNA/1km_Hamman_raw/27_biovars")
      d <- d[!grepl("MAR", d$path),] # this variable doesn't cover mexico, the rest do
      b <- stack(d$path[grepl("6190", d$path)])
      #f <- stack(d$path[grepl("2020", d$path)])
      f <- stack(d$path[grepl("2080", d$path)])
      
      ext <- c(-2723106, 2522840, 3426221, 7149212) # for 27 bio
      b <- crop(b, ext)
      f <- crop(f, ext)
      
      b <- aggregate(b, 8)
      f <- aggregate(f, 8)
      
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
      
      saveRDS(list(b, f), "E:/analog_migration/animation/data_b26.rds")
}
d <- readRDS("E:/analog_migration/animation/data_b26.rds"); b <- d[[1]]; f <- d[[2]]


# prep data
bv <- as.data.frame(rasterToPoints(b))
names(bv) <- c("x", "y",  paste0("v", 1:(ncol(bv)-2)))
fv <- as.data.frame(rasterToPoints(f))
names(fv) <- c("x", "y",  paste0("v", 1:(ncol(fv)-2)))
bv$year <- "baseline"
fv$year <- "future"
v <- na.omit(rbind(bv, fv))
v <- v %>% mutate_each(funs(scale), -x, -y, -year)
vm <- as.matrix(dplyr::select(v, -x, -y, -year))

y <- v %>%
      gather(variable, value, -year, -x, -y) %>%
      spread(year, value) %>%
      mutate(delta=future-baseline)


# training samples
n <- 1000
# sampled from geographic space
set.seed(123)
#train <- vm[sample(nrow(vm), n),]
# sampled from climate space -- random
#train <- apply(vm, 2, function(x) c(min(x), runif(n-2, min(x), max(x)), max(x)))
# sampled from climate space -- uniform
#train <- apply(vm, 2, function(x) sample(seq(min(x), max(x), c(max(x)-min(x))/(n-1)), n, replace=F))
# training points are kmeans cluster centers
set.seed(123)
train <- vm[sample(nrow(vm), 10000),] # have to reduce size for kmeans to work
train <- kmeans(train, n)$centers

# TSNE
#set.seed(123)
#embedded <- tsne(train, k=3)

#NMDS
dst <- dist(train)
set.seed(123)
fit <- MASS::isoMDS(dst, k=3)
embedded <- fit$points

# compare distributions of embedded data
#library(rgl)
#plot3d(embedded, size=10, col=as.integer(factor(v$year[train])))
#plot3d(embedded, size=10, col=colors)

method <- c("n1000_kmeans_nmds_ecdf")

######## small multiples plot to choose colors #########


sm <- function(){
      steps <- 3
      i <- 2
      yi <- y %>%
            na.omit() %>%
            mutate(value = baseline + delta * i/steps) %>%
            select(x, y, variable, value) %>%
            spread(variable, value)
      yim <- yi %>%
            select(-x, -y) %>%
            as.matrix()
      nn <- get.knnx(train, yim, k=1)
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
      ggsave(paste0("E:/analog_migration/animation/sm_", method, ".png"), p, width=20, height=15, units="in")
}

sm()

stop()


######## animation #########

colors <- colors3d(embedded, trans="ecdf", order=2, inversion=2)

steps <- 200

saveGIF({
      for(i in 0:steps){
            
            yi <- yi <- y %>%
                  na.omit() %>%
                  mutate(value = baseline + delta * i/steps) %>%
                  select(x, y, variable, value) %>%
                  spread(variable, value)
            
            yim <- yi %>%
                  select(-x, -y) %>%
                  as.matrix()
            
            # assign color of NMDS training values to entire dataset based on nearest neighbors
            nn <- get.knnx(train, yim, k=1)
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
            print(p)
            
      }
      
},
movie.name = paste0("analog_animation_", method, ".gif"),
interval = .03,
ani.width = 1500,
ani.height = 1125,
outdir = "E:/analog_migration"
)


