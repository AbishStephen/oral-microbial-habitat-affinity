#!/usr/bin/env Rscript
script<-sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))][1])
script<-gsub("~+~"," ",script,fixed=TRUE)
source(file.path(dirname(normalizePath(script)),"common.R"))
args<-commandArgs(TRUE);chosen<-"all"
if("--figure"%in%args)chosen<-args[match("--figure",args)+1]
if("--output"%in%args)OUT<-normalizePath(args[match("--output",args)+1],mustWork=FALSE)
if(!chosen%in%c("all",as.character(1:5),"supplement"))stop("--figure must be 1,2,3,4,5,supplement or all")
for(i in 1:5)if(chosen%in%c("all",as.character(i))) {
 source(file.path(REPO,"r",paste0("figure",i,".R")))
 # common.R is re-sourced by an individual entry point; restore an explicit output.
 if("--output"%in%args)OUT<-normalizePath(args[match("--output",args)+1],mustWork=FALSE)
 sizes<-list(c(18.2,13.3),c(19.2,13.8),c(19.2,13),c(17.1,11.3),c(17,11.5))
 write_plot(get(paste0("build_figure",i))(),paste0("Figure",i),sizes[[i]][1],sizes[[i]][2])
}
if(chosen%in%c("all","supplement")){
 source(file.path(REPO,"r/supplement.R"))
 if("--output"%in%args)OUT<-normalizePath(args[match("--output",args)+1],mustWork=FALSE)
 write_plot(build_supplement(),"Supplementary_candidate_domain_neighborhoods",11.5,14)
}
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
versions<-data.frame(package=required,version=vapply(required,function(x)as.character(packageVersion(x)),character(1)))
write.table(versions,file.path(OUT,"R_package_versions.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
# sessionInfo can contain local library paths; record only portable version/platform fields.
jsonlite::write_json(list(R=as.character(getRversion()),platform=R.version$platform,packages=as.list(setNames(versions$version,versions$package))),file.path(OUT,"R_environment.json"),pretty=TRUE,auto_unbox=TRUE)

if(chosen=="all") {
 exported<-read.delim(file.path(OUT,"source_data/Supplementary_gene_coordinates.tsv"),check.names=FALSE,quote="",comment.char="")
 stopifnot(nrow(exported)==206,length(unique(exported$locus_id))==30)
 predictions<-read.delim(file.path(OUT,"source_data/Figure5_heldout_predictions.tsv"),check.names=FALSE,quote="")
 stopifnot(nrow(predictions)==47,length(unique(predictions$Genus))==13)
 message("Export round-trip checks passed: 206 gene records and 47 held-out predictions")
}
