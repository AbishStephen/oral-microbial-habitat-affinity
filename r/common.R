# Shared data access and vector/raster export; no cluster paths or Python required.
main_script <- function() {
 a <- commandArgs(trailingOnly=FALSE); f <- a[grepl("^--file=",a)]
 if(length(f)) return(normalizePath(gsub("~+~"," ",sub("^--file=","",f[1]),fixed=TRUE),mustWork=TRUE))
 files<-unlist(lapply(sys.frames(),function(x)if(!is.null(x$ofile))as.character(x$ofile)else NULL))
 if(length(files))return(normalizePath(tail(files,1),mustWork=TRUE))
 if(file.exists("r/build_figures.R"))return(normalizePath("r/build_figures.R"))
 stop("Open the project folder or run Rscript r/build_figures.R")
}
REPO <- normalizePath(file.path(dirname(main_script()),".."),mustWork=TRUE)
.local_lib <- file.path(REPO,".R_library")
if(dir.exists(.local_lib)) .libPaths(c(.local_lib,.libPaths()))
required <- c("ggplot2","patchwork","jsonlite","ape","svglite","ragg")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Install required R packages: ",paste(missing,collapse=", "))
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
source(file.path(REPO,"r/style.R"))
update_geom_defaults("text",list(family=STYLE$font_family))
OUT <- file.path(REPO,"results/R")
read_tsv <- function(path) read.delim(file.path(REPO,path),check.names=FALSE,
 stringsAsFactors=FALSE,quote="",comment.char="",na.strings=c("NA",""))
read_json <- function(path) {
 text<-paste(readLines(file.path(REPO,path),warn=FALSE),collapse="\n")
 # Some Python JSON exports encode missing annotations as unquoted NaN.
 text<-gsub("([,:]\\s*)NaN(?=\\s*[,}\\]])","\\1null",text,perl=TRUE)
 jsonlite::fromJSON(text,simplifyVector=FALSE)
}
records <- function(x) do.call(rbind,lapply(x,function(r) as.data.frame(lapply(r,function(v)if(is.null(v)) NA else v),stringsAsFactors=FALSE)))
write_source <- function(d,name) {
 dir.create(file.path(OUT,"source_data"),recursive=TRUE,showWarnings=FALSE)
 disk<-d
 for(k in names(disk)) {
  if(is.double(disk[[k]]))disk[[k]]<-ifelse(is.na(disk[[k]]),NA,trimws(formatC(disk[[k]],format="g",digits=17)))
  else if(is.character(disk[[k]])||is.factor(disk[[k]]))disk[[k]]<-gsub("[\r\n\t]"," ",as.character(disk[[k]]))
 }
 write.table(disk,file.path(OUT,"source_data",name),sep="\t",row.names=FALSE,quote=FALSE,na="")
}
habitat_scale <- function() scale_colour_manual(values=c(plaque=STYLE$plaque,mucosal=STYLE$mucosal),name="HMT affinity")
base_theme <- function(size=STYLE$base_size) theme_minimal(base_size=size,base_family=STYLE$font_family)+
 theme(panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),
 panel.grid.major.x=element_line(colour=STYLE$grid,linewidth=.25),
 text=element_text(colour=STYLE$text),axis.text=element_text(colour=STYLE$muted,size=STYLE$axis_size),
 axis.title=element_text(size=9),plot.title=element_text(face="bold",size=STYLE$panel_title_size),
 plot.subtitle=element_text(size=8.5,colour=STYLE$muted,margin=margin(b=7)),
 plot.caption=element_text(size=STYLE$caption_size,colour=STYLE$muted,hjust=0),
 plot.tag=element_text(face="bold",size=16),legend.position="bottom",legend.title=element_text(size=8.5),
 legend.text=element_text(size=8),plot.margin=margin(10,14,8,10),strip.text=element_text(face="bold",size=9))
blank_plot <- function() ggplot()+theme_void(base_family=STYLE$font_family)
wrap <- function(x,width=40) vapply(x,function(s) paste(strwrap(s,width),collapse="\n"),character(1))
annotation <- function(title,subtitle,caption) plot_annotation(title=title,subtitle=subtitle,caption=caption,
 theme=theme(plot.title=element_text(size=STYLE$title_size,face="bold",family=STYLE$font_family,colour=STYLE$text),
 plot.subtitle=element_text(size=10,colour=STYLE$muted),plot.caption=element_text(size=STYLE$caption_size,hjust=0,
 colour=STYLE$muted),plot.margin=margin(16,18,12,16)))
write_plot <- function(p,name,width,height) {
 dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
 grDevices::cairo_pdf(file.path(OUT,paste0(name,".pdf")),width=width,height=height,family=STYLE$font_family)
 print(p);invisible(dev.off())
 svglite::svglite(file.path(OUT,paste0(name,".svg")),width=width,height=height,
  system_fonts=list(sans=STYLE$font_family));print(p);invisible(dev.off())
 ragg::agg_png(file.path(OUT,paste0(name,".png")),width=width,height=height,units="in",res=STYLE$dpi,background="white")
 print(p);invisible(dev.off())
 message("Rendered ",name," (PDF, SVG, PNG)")
}
# Polygons retain true gene positions. Clipping bars identify truncated genes.
gene_polygons <- function(genes,window=c(-3500,3500),height=.18) {
 pieces <- list(); k <- 0
 for(i in seq_len(nrow(genes))) {
  r <- genes[i,];lo <- max(window[1],r$xstart);hi <- min(window[2],r$xend)
  if(hi<=lo) next
  head <- min(180,(hi-lo)*.38);y <- r$y
  if(r$display_strand=="+") {
   x <- c(lo,hi-head,hi-head,hi,hi-head,hi-head,lo)
   yy <- y+c(-height*.6,-height*.6,-height,0,height,height*.6,height*.6)
  } else {
   x <- c(hi,lo+head,lo+head,lo,lo+head,lo+head,hi)
   yy <- y+c(-height*.6,-height*.6,-height,0,height,height*.6,height*.6)
  }
  k <- k+1;pieces[[k]] <- data.frame(x=x,y=yy,id=i,role=r$role,track=r$track,stringsAsFactors=FALSE)
 }
 do.call(rbind,pieces)
}
