#!/usr/bin/env Rscript
script<-sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))][1])
script<-gsub("~+~"," ",script,fixed=TRUE)
source(file.path(dirname(normalizePath(script)),"common.R"))
checks<-list()
m<-read_tsv("figure_sources/Figure1/HOMD47_HMT_manifest.tsv")
checks$cohort<-nrow(m)==47&&length(unique(m$Genus))==13&&sum(m$n_genomes)==1121&&sum(m$Ecology=="p")==24&&sum(m$Ecology=="m")==23
primary<-read_tsv("data/frozen/project/study13/01_genome_manifests/EXPANDED_HMT_MANIFEST.tsv")
checks$full_precision_ecology<-isTRUE(all.equal(m$Median_Y,primary$Median_Y[match(m$HMT,primary$HMT)],tolerance=1e-14,check.attributes=FALSE))
focal<-records(read_json("figure_sources/Figure2/Figure2_plotting_data.json")$focal)
checks$exploratory_statistics<-nrow(focal)==14&&all(focal$raw_perm_p_consistency<.05)&&all(focal$adj_perm_p_consistency<.05)&&all(focal$raw_maxT_all_feature_FWER_p_consistency>=.05)&&all(focal$adj_maxT_all_feature_FWER_p_consistency>=.05)
pred<-records(read_json("figure_sources/Figure5/Figure5_model_predictions_and_sister_pairs.json")$pred)
expected<-read_tsv("data/frozen/expected/nested_genusbalanced_sourcebalanced_HMT_predictions.tsv")
prediction_difference<-max(abs(pred$Y_pred-expected$Y_pred[match(pred$HMT,expected$HMT)]))
checks$full_precision_predictions<-prediction_difference<1e-12
macro_auc<-function(d) mean(vapply(split(d,d$Genus),function(z){p<-z$Y_pred[z$Ecology=="plaque"];m<-z$Y_pred[z$Ecology=="mucosal"];mean(outer(p,m,">")+0.5*outer(p,m,"=="))},numeric(1)))
checks$prediction_auc<-abs(macro_auc(pred)-0.8205128205128205)<1e-12
g<-read_tsv("data/frozen/neighborhoods/locus_genes.tsv");l<-read_tsv("data/frozen/neighborhoods/selected_loci.tsv")
checks$supplement_counts<-nrow(g)==206&&nrow(l)==30&&length(unique(l$Pfam))==10&&length(unique(l$HMT))==16&&length(unique(l$GCA))==25
checks$coordinate_lengths<-all(g$display_end-g$display_start==g$end-g$start+1)
f<-as.character(g$focal)%in%c("True","TRUE");checks$focal_orientation<-all(g$display_strand[f]=="+")&&all(g$display_start[f]==0)
source(file.path(REPO,"r/figure4.R"));norm<-normalise_main_loci(read_json("figure_sources/Figure4/Figure4_real_GFF_loci.json"))
checks$main_coordinate_lengths<-all(abs((norm$genes$xend-norm$genes$xstart)-(norm$genes$end-norm$genes$start))<1e-12)
checks$main_reference_loci<-nrow(norm$loci)==8&&sum(norm$loci$flip)==1
if(!all(unlist(checks))) stop("Failed scientific-input checks: ",paste(names(checks)[!unlist(checks)],collapse=", "))
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
jsonlite::write_json(list(passed=TRUE,n_checks=length(checks),maximum_prediction_difference=prediction_difference,numeric_tolerance=1e-12,checks=checks),file.path(OUT,"R_input_checks.json"),auto_unbox=TRUE,pretty=TRUE)
message("Passed ",length(checks)," R scientific-input checks")
