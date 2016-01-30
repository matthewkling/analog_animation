
library(ecoclim)
library(raster)
library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(FNN)
library(tsne)
library(animation)


d <- parseMetadata("E:/ephemera/ephemera/data/monthlies", pattern=".grd", drops="_2012")
s <- lapply(d$path, function(x)stack(x))
b <- stack(s[c(1,3)])
f <- stack(s[c(2,4)])

d <- parseMetadata("F:/ClimateNA/1km_Hamman_raw/27_biovars")
b <- stack(d$path[grepl("6190", d$path)])
b <- stack(d$path[grepl("2020", d$path)])
f <- stack(d$path[grepl("2080", d$path)])


# crop, for testing purposes
#plot(b[[2]])
#ext <- drawExtent()
ext <- c(-126.6319, -67.10868, 24.56991, 49.24116) # for 8monthly and 19bio
ext <- c(-2723106, 2522840, 3426221, 7149212) # for 27 bio
b <- crop(b, ext)
f <- crop(f, ext)

b <- aggregate(b, 4) #only aggregate the 1km data
f <- aggregate(f, 4)

# log transform ppt, and round up logs of <1 mm ppt
logtransform <- c(1:12) # for 19 biovars
logtransform <- c(24, 22, 20, 13, 6, 4, 1) # for 27 biovars
notransform <- setdiff(1:27, logtransform) # retry with all 27--this was accidentally 1:19
ltrans <- function(x){
      x[x<1] <- 1
      log10(x)
}
b <- stack(b[[notransform]], calc(b[[logtransform]], ltrans))
f <- stack(f[[notransform]], calc(f[[logtransform]], ltrans))

b <- mask(b, f)
f <- mask(f, b)

# prep data
bv <- as.data.frame(rasterToPoints(b))
names(bv) <- c("x", "y",  paste0("v", 1:(ncol(bv)-2)))
fv <- as.data.frame(rasterToPoints(f))
names(fv) <- c("x", "y",  paste0("v", 1:(ncol(fv)-2)))
bv$year <- "baseline"
fv$year <- "future"
v <- na.omit(rbind(bv, fv))
v <- v %>% mutate_each(funs(scale), -x, -y, -year)
vm <- as.matrix(select(v, -x, -y, -year))

train <- sample(nrow(vm), 1000) # 1000 # should this be random, or stratified in cartesian e-space to avoid oversampling dense areas?



nn <- get.knnx(vm[train,], vm, k=1)


######## this is totally messed up -- go back to nd_3d_animation and reconstruct correct methods
####### one problem: baseline and future points are totally segregated by TSNE


# 3 dimension reduction methods
for(drm in c("tsne", "nmds", "pca")){

      if(drm=="nmds"){
            dst <- dist(vm[train,])
            fit <- isoMDS(dst, k=3)
            embedded <- fit$points
      }
      if(drm=="tsne") embedded <- tsne(vm[train,], k=3)
      if(drm=="pca") embedded <- prcomp(vm[train,])$x[,1:3]

      #rgl::plot3d(embedded,
      #            col=c("black", "red")[as.integer(factor(v$year[train]))],
      #            size=10)

      # assign embedded values to entire dataset based on nearest neighbors
      vnn <- embedded[as.vector(nn$nn.index),]
      vnn <- apply(vnn, 2, scales::rescale)
      colnames(vnn) <- c("V1", "V2", "V3")

      # two rescaling methods
      for(rm in c("ecdf", "rescale")){

            # animation
            steps <- 100
            saveGIF({
                  for(i in 0:steps){

                        d <- cbind(v[, c("x", "y", "year")], as.data.frame(vnn)) %>%
                              gather(variable, value, V1:V3) %>%
                              spread(year, value) %>%
                              mutate(value = baseline + (future - baseline) * i / steps) %>%
                              select(-baseline, -future) %>%
                              spread(variable, value)
                        # would it be bettter to interpolate in hi-dim space than in embedded space?

                        ym <- as.matrix(select(d, V1:V3))


                        #######if(rm=="ecdf") vnnr <- apply(vnn, 2, function(x)ecdf()(x))

                        # color options
                        rotations <- expand.grid(r=1:3,
                                                 g=1:3,
                                                 b=1:3,
                                                 ir=c(T, F),
                                                 ig=c(T, F),
                                                 ib=c(T, F)) %>%
                              rowwise() %>%
                              filter(length(unique(c(r, g, b)))==3)

                        getColors <- function(j){
                              df <- cbind(vnn[train,rotations$r[j]], vnn[train,rotations$g[j]], vnn[train,rotations$b[j]])
                              for(invert in 1:3) if(as.character(rotations[j, invert+3])) df[,invert] <- 1- df[,invert]
                              cols <- rgb(df[,1], df[,2], df[,3])
                              return(cols)
                        }

                        cols <- getColors(1)

                        # assign embedded values to entire dataset based on nearest neighbors
                        ynn <- get.knnx(vnn[train,], ym, k=1) ##################################### calculate nearest neighbor in training dataset
                        cols <- cols[as.vector(ynn$nn.index)]

                        p <- ggplot(d, aes(x, y)) +
                              geom_raster(fill=cols) +
                              ggmap::theme_nothing() +
                              coord_fixed(ratio=1.0) # 1.5
                        print(p)

                  }

            },
            movie.name = paste0("animation_", drm, "_", rm, ".gif"),
            interval = .05,
            ani.width = 2000,
            ani.height = 1500,
            outdir = "E:/analog_migration"
            )

            rgl::plot3d(embedded,
                        col=getColors(1),
                        size=10)


      }
}


