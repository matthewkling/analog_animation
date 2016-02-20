library(animation)





saveVid <- function (expr, video.name = "animation.mp4", img.name = "Rplot", 
                     ffmpeg = ani.options("ffmpeg"), 
                     other.opts = if (grepl("[.]mp4$", video.name)) "-pix_fmt yuv420p", 
                     tmpdr, ...){
      
      require(doParallel)
      require(devtools)
      
      oopt = ani.options(...)
      on.exit(ani.options(oopt))
      owd = setwd(tmpdr)
      on.exit(setwd(owd), add = TRUE)
      if (is.null(ffmpeg)) {
            ffmpeg = "ffmpeg"
      }
      if (!grepl("^[\"']", ffmpeg)) 
            ffmpeg = shQuote(ffmpeg)
      version = try(system(paste(ffmpeg, "-version"), intern = TRUE))
      if (inherits(version, "try-error")) {
            warning("The command \"", ffmpeg, "\" is not available in your system. Please install FFmpeg or avconv first: ", 
                    ifelse(.Platform$OS.type == "windows", "http://ffmpeg.arrozcru.org/autobuilds/", 
                           "http://ffmpeg.org/download.html"))
            return()
      }
      ani.dev = ani.options("ani.dev")
      file.ext = ani.options("ani.type")
      interval = ani.options("interval")
      if (is.character(ani.dev)) 
            ani.dev = get(ani.dev)
      num = ifelse(file.ext == "pdf", "", "%d")
      unlink(paste(img.name, "*.", file.ext, sep = ""))
      img.fmt = paste(img.name, num, ".", file.ext, sep = "")
      img.fmt = file.path(tmpdr, img.fmt)
      ani.options(img.fmt = img.fmt)
      if ((use.dev <- ani.options("use.dev"))) 
            ani.dev(img.fmt, width = ani.options("ani.width"), height = ani.options("ani.height"))
      in_dir(owd, expr)
      #in_dir(tmpdr, expr)
      if (use.dev) 
            dev.off()
      ffmpeg = paste(ffmpeg, "-y", "-r", 1/ani.options("interval"), 
                     "-i", basename(img.fmt), other.opts, basename(video.name))
      message("Executing: ", ffmpeg)
      cmd = system(ffmpeg)
      
      auto_browse = function(output){
            if(!ani.options('autobrowse')) return()
            if (.Platform$OS.type == 'windows') {
                  try(shell.exec(output))
            } else if (Sys.info()['sysname'] == 'Darwin') {
                  system(paste('open ', shQuote(output)))
            } else system(paste('xdg-open ', shQuote(output)))
      }
      
      if (cmd == 0) {
            setwd(owd)
            file.copy(file.path(tmpdr, basename(video.name)), 
                      video.name, overwrite = TRUE)
            message("\n\nVideo has been created at: ", output.path <- normalizePath(video.name))
            auto_browse(output.path)
      }
      invisible(cmd)
}