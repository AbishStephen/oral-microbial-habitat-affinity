.args <- commandArgs(FALSE)
.file <- .args[grepl("^--file=",.args)]
if(length(.file)) .rdir <- dirname(normalizePath(gsub("~+~"," ",sub("^--file=","",.file[1]),fixed=TRUE))) else {
 .src <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
 .rdir <- if(!is.null(.src))dirname(normalizePath(.src)) else file.path(getwd(),"r")
}
source(file.path(.rdir,"common.R"))
build_figure3 <- function() {
 d<-read_tsv("figure_sources/Figure3/Figure3_A_all_13_domain_cluster_counts.tsv")
 ord<-c("PF01638","PF04072","PF13676","PF12680","PF01745","PF18075","PF13413","PF22166","PF07366","PF13577","PF10604","PF01337","PF14470")
 d$Pfam<-factor(d$Pfam,levels=rev(ord));stopifnot(nrow(d)==13)
 pa<-ggplot(d,aes(y=Pfam))+geom_segment(aes(x=max_genera_in_a_single_cluster,xend=n_genera,yend=Pfam),linewidth=1.1,colour="#ADC5CD")+
  geom_point(aes(x=n_genera,colour="Domain observed"),size=3)+geom_point(aes(x=max_genera_in_a_single_cluster,colour="Broadest sequence cluster"),shape=18,size=3.2)+
  scale_colour_manual(values=c("Domain observed"=STYLE$mucosal,"Broadest sequence cluster"=STYLE$plaque),name=NULL)+
  scale_y_discrete(labels=function(x)unname(DOMAIN_NAMES[x]))+scale_x_continuous(limits=c(0,13.4),breaks=c(0,3,6,9,13))+
  base_theme()+labs(x="Number of genera",y=NULL,title="Domains encompass diverse protein families",subtitle="",tag="A")
 pcnt<-ggplot(d,aes(x=n_distinct_MMseqs_X5_clusters,y=Pfam))+geom_col(fill="#A2B7C7",width=.5)+
  geom_text(aes(label=n_distinct_MMseqs_X5_clusters),hjust=-.25,size=2.7,colour=STYLE$text)+
  scale_x_continuous(limits=c(0,86),breaks=c(0,40,80))+base_theme()+
  theme(axis.text.y=element_blank(),axis.title.y=element_blank(),legend.position="none")+
  labs(x="Sequence clusters",title="",subtitle="MMseqs X5")
 memberships<-read_tsv("figure_sources/Figure3/Figure3_B_selected_cluster_genus_memberships.tsv")
 selected<-memberships[(memberships$Pfam=="PF01638"&memberships$rank<=6)|(memberships$Pfam=="PF04072"&memberships$rank<=4)|(memberships$Pfam=="PF13676"&memberships$rank<=3),]
 selected<-selected[order(match(selected$Pfam,c("PF01638","PF04072","PF13676")),selected$rank),]
 selected$row<-paste(unname(DOMAIN_NAMES[selected$Pfam]),sprintf("%02d",selected$rank));selected$row<-factor(selected$row,levels=rev(selected$row))
 pts<-do.call(rbind,lapply(seq_len(nrow(selected)),function(i){gs<-strsplit(selected$genera[i],";",fixed=TRUE)[[1]];data.frame(row=selected$row[i],Genus=gs,x=match(gs,GENERA),Phylum=unname(PHYLA[gs]))}))
 spans<-do.call(rbind,lapply(seq_len(nrow(selected)),function(i){xs<-match(strsplit(selected$genera[i],";",fixed=TRUE)[[1]],GENERA);data.frame(row=selected$row[i],x=min(xs),xend=max(xs))}))
 pts$row<-factor(pts$row,levels=levels(selected$row));spans$row<-factor(spans$row,levels=levels(selected$row))
 pb<-ggplot()+geom_segment(data=spans,aes(x=x,xend=xend,y=row,yend=row),colour="#C0D0D5",linewidth=.9)+
  geom_point(data=pts,aes(x=x,y=row,colour=Phylum),size=2.6)+
  geom_text(data=selected,aes(x=13.6,y=row,label=paste0(n_proteins," proteins")),hjust=0,size=2.6,colour=STYLE$muted)+
  scale_colour_manual(values=PHYLUM_COLOURS)+scale_x_continuous(breaks=1:13,labels=GENERA,limits=c(.5,17))+
  base_theme()+theme(axis.text.x=element_text(angle=45,hjust=1,size=7.2),axis.text.y=element_text(size=8),legend.position="none")+
  labs(x=NULL,y=NULL,title="Selected protein clusters across genera",subtitle="Rows are actual sequence clusters; coloured points show observed membership",tag="B")
 ar<-read_tsv("figure_sources/Figure3/Figure3_C_domain_architecture_proportions.tsv")
 ar$label<-paste(unname(DOMAIN_NAMES[ar$Pfam]),ar$main_coannotation,sep=" | ");ar$label<-factor(ar$label,levels=rev(ar$label))
 pc<-ggplot(ar,aes(y=label))+geom_col(aes(x=100),fill=STYLE$grid,width=.5)+
  geom_col(aes(x=percent),fill=STYLE$mucosal,width=.5)+
  geom_text(aes(x=103,label=sprintf("%.1f%% (%s/%s)",percent,format(n_primary,big.mark=","),format(n_proteins,big.mark=","))),hjust=0,size=2.8)+
  scale_x_continuous(limits=c(0,150),breaks=c(0,25,50,75,100))+base_theme()+
  theme(axis.text.y=element_text(size=8.5))+labs(x="Proteins with the indicated co-annotation (%)",y=NULL,
   title="Dominant protein-domain combinations",subtitle="",tag="C")
 om<-read_tsv("figure_sources/Figure3/Figure3_D_omnibus_repertoire_tests.tsv")
 om$layer<-factor(om$layer,levels=rev(c("PFAM","COG","EC","KMOD","KPATH","X5")))
 pd<-ggplot(om,aes(x=perm_p,y=layer))+geom_col(width=.5,fill="#B4C7CF")+
  geom_point(colour=STYLE$mucosal,size=2.3)+geom_vline(xintercept=.05,linetype=2,colour=STYLE$plaque)+
  geom_text(aes(x=1.03,label=sprintf("P = %.3f",perm_p)),hjust=0,size=3)+
  scale_x_continuous(limits=c(0,1.27),breaks=c(0,.25,.5,.75,1))+
  scale_y_discrete(labels=c(PFAM="Pfam",COG="COG",EC="EC",KMOD="KEGG modules",KPATH="KEGG pathways",X5="Sequence clusters"))+
  base_theme()+labs(x="Within-genus permutation P value",y=NULL,title="Genome-wide directional coherence",subtitle="No significant omnibus effect in the six annotation layers",tag="D")
 write_source(d,"Figure3_domain_breadth.tsv");write_source(selected,"Figure3_selected_sequence_clusters.tsv");write_source(ar,"Figure3_domain_coannotations.tsv");write_source(om,"Figure3_omnibus_P_values.tsv")
 a<-wrap_elements(full=(pa|pcnt)+plot_layout(widths=c(3.6,1)))
 ((a|pb)/(pc|pd))+plot_layout(widths=c(1.18,1),heights=c(1.4,1))+
  annotation("","",
   "Sequence clusters are not experimentally established orthologues. Genomes are nested within HMTs; neither domain co-annotation nor omnibus nonsignificance establishes biochemical mechanism.")
}
if(!interactive() && basename(main_script())=="figure3.R") write_plot(build_figure3(),"Figure3",19.2,13)
