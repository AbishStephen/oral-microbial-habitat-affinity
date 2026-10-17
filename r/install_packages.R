#!/usr/bin/env Rscript
# Restore the tested dependency versions into a project-local library.
a<-commandArgs(FALSE);f<-sub("^--file=","",a[grepl("^--file=",a)][1])
f<-gsub("~+~"," ",f,fixed=TRUE)
repo<-normalizePath(file.path(dirname(normalizePath(f)),".."))
lib<-file.path(repo,".R_library");dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"))
if(!requireNamespace("renv",quietly=TRUE))install.packages("renv",lib=lib)
renv::restore(project=repo,library=lib,lockfile=file.path(repo,"renv.lock"),prompt=FALSE)
message("Locked R figure dependencies restored.")
