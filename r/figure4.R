.args <- commandArgs(FALSE)
.file <- .args[grepl("^--file=",.args)]
if(length(.file)) .rdir <- dirname(normalizePath(gsub("~+~"," ",sub("^--file=","",.file[1]),fixed=TRUE))) else {
 .src <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
 .rdir <- if(!is.null(.src))dirname(normalizePath(.src)) else file.path(getwd(),"r")
}
source(file.path(.rdir,"common.R"))
normalise_main_loci <- function(rows) {
 genes<-list();loci<-list()
 for(i in seq_along(rows)) {
  r<-rows[[i]];g<-records(r$genes);center<-(r$anchor_start+r$anchor_end)/2
  flip<-r$type=="HxlR/SnoaL2"&&r$strand1=="+"&&r$strand2=="-"
  g$xstart<-g$start-center;g$xend<-g$end-center;g$display_strand<-g$strand
  if(flip){g$xstart<-center-g$end;g$xend<-center-g$start;g$display_strand<-ifelse(g$strand=="+","-","+")}
  g$role<-ifelse(g$gene==r$protein1,"gene1",ifelse(g$gene==r$protein2,"gene2",ifelse(grepl("recombinase|transposase",g$product,ignore.case=TRUE),"mobile","other")))
  g$locus_index<-i;g$track<-as.character(i);g$y<-0
  genes[[i]]<-g;loci[[i]]<-data.frame(locus_index=i,Genus=r$Genus,GCA=r$GCA,HMT=r$HMT,distance_bp=r$distance_bp,type=r$type,flip=flip)
 }
 list(genes=do.call(rbind,genes),loci=do.call(rbind,loci))
}
locus_panel <- function(g,l,core=TRUE) {
 n<-nrow(l);l$y<-rev(seq_len(n));g$y<-l$y[match(g$locus_index,l$locus_index)]
 window<-if(core)c(-3500,3500) else c(-2900,2900);polys<-gene_polygons(g,window,.17)
 focal<-g[g$role%in%c("gene1","gene2"),];focal$label<-ifelse(focal$role=="gene1","HxlR","SnoaL_2")
 if(!core) focal$label<-ifelse(focal$role=="gene1",ifelse(focal$locus_index==7,"HxlR","SnoaL_2"),ifelse(focal$locus_index==7,"WHD_ArnR1","Polyketide_cyc2"))
 focal$mid<-pmax(window[1]+150,pmin(window[2]-150,(focal$xstart+focal$xend)/2));focal$label_y<-focal$y+ifelse(focal$role=="gene1",.27,-.28)
 l$left<-paste(l$Genus,l$HMT,sep=" | ");l$right<-paste0(format(l$distance_bp,big.mark=",")," bp")
 p<-ggplot()+geom_segment(data=l,aes(x=window[1],xend=window[2],y=y,yend=y),colour=STYLE$grid,linewidth=.4)+
  geom_polygon(data=polys,aes(x=x,y=y,group=id,fill=role),colour="white",linewidth=.15)+
  geom_text(data=focal,aes(x=mid,y=label_y,label=label),size=3,fontface="bold",colour=STYLE$text)+
  geom_text(data=l,aes(x=window[1]-170,y=y+.05,label=left),hjust=1,size=3,fontface="italic",colour=STYLE$text)+
  geom_text(data=l,aes(x=window[1]-170,y=y-.18,label=GCA),hjust=1,size=2.4,colour=STYLE$muted)+
  geom_text(data=l,aes(x=window[2]+170,y=y,label=right),hjust=0,size=2.7,colour=STYLE$muted)+
  scale_fill_manual(values=c(gene1=STYLE$plaque,gene2=STYLE$partner,other=STYLE$other_gene,mobile=STYLE$mobile),labels=c(gene1="Focal domain-bearing CDS",gene2="Linked domain-bearing CDS",other="Other CDS",mobile="Recombinase / transposase"),name=NULL)+
  scale_x_continuous(limits=c(window[1]-2400,window[2]+1500),breaks=seq(-3000,3000,1000),labels=function(x)x/1000)+
  scale_y_continuous(limits=c(.35,n+.65))+base_theme()+
  theme(axis.text.y=element_blank(),axis.title.y=element_blank(),panel.grid=element_blank(),legend.position=if(core)"bottom" else "none")+
  labs(x="Position relative to locus midpoint (kb)")
 p
}
build_figure4 <- function() {
 rows<-read_json("figure_sources/Figure4/Figure4_real_GFF_loci.json");n<-normalise_main_loci(rows)
 stopifnot(length(rows)==8,nrow(n$loci)==8)
 pa<-locus_panel(n$genes[n$genes$locus_index<=6,],n$loci[1:6,],TRUE)+
  labs(title="HxlR-SnoaL_2 genomic neighbourhoods",subtitle="Six assemblies representing five genera; actual CDS positions and strands",tag="A")
 dist<-read_tsv("figure_sources/Figure4/Figure4_B_all67_distances.tsv");stopifnot(nrow(dist)==67,sum(dist$bp<=500&dist$orientation=="divergent")==66)
 orderg<-c("Actinomyces","Porphyromonas","Prevotella","Streptococcus","Neisseria");dist$genus<-factor(dist$genus,levels=rev(orderg));dist$Phylum<-unname(PHYLA[as.character(dist$genus)])
 dist$offset<-ave(seq_len(nrow(dist)),dist$genus,FUN=function(x)((seq_along(x)*7)%%13-6)*.014)
 dist$y<-as.numeric(dist$genus)+dist$offset
 pb<-ggplot(dist,aes(x=bp,y=y,colour=Phylum))+annotate("rect",xmin=38,xmax=500,ymin=.5,ymax=5.5,fill=STYLE$pale)+
  geom_vline(xintercept=500,linetype=2,colour=STYLE$mobile)+geom_point(size=2,alpha=.8)+
  scale_colour_manual(values=PHYLUM_COLOURS)+scale_x_log10(limits=c(38,3200),breaks=c(50,100,200,500,1000,2000))+
  scale_y_continuous(breaks=1:5,labels=rev(orderg))+base_theme()+theme(legend.position="none")+
  labs(x="Intergenic separation (bp; log scale)",y=NULL,title="Spacing across all linked genomes",subtitle="66 of 67 arrangements are divergent within 500 bp",tag="B")
 pc<-locus_panel(n$genes[n$genes$locus_index>6,],n$loci[7:8,],FALSE)+
  labs(title="Alternative linked-domain arrangements",subtitle="Actual gene coordinates in Fusobacterium and Actinomyces",tag="C")
 types<-read_tsv("figure_sources/Figure4/Figure4_D_8cluster_pair_types.tsv");types$pair<-factor(seq_len(nrow(types)),levels=8:1);types$breadth<-ifelse(types$n_genera>1,"Two genera","One genus")
 pd<-ggplot(types,aes(x=n_genomes,y=pair,fill=breadth))+geom_col(width=.6)+
  geom_text(aes(label=n_genomes),hjust=-.4,size=3)+scale_fill_manual(values=c("One genus"=STYLE$other_gene,"Two genera"=STYLE$mucosal),name=NULL)+
  scale_x_continuous(limits=c(0,30),breaks=c(0,10,20,30))+base_theme()+
  labs(x="Genomes carrying the exact cluster pair",y="Pair type",title="Eight protein-cluster pair types",subtitle="Only one pair type spans more than one genus",tag="D")
 write_source(n$genes,"Figure4_normalised_gene_coordinates.tsv");write_source(n$loci,"Figure4_reference_loci.tsv");write_source(dist,"Figure4_linkage_distances.tsv");write_source(types,"Figure4_cluster_pair_types.tsv")
 ((pa|pb)/(pc|pd))+plot_layout(widths=c(2.05,1),heights=c(1.5,1))+
  annotation("","",
   "Gene adjacency and orientation do not establish shared regulation, identical protein function or independent evolutionary gains. Linkage enrichment with plaque affinity: two-sided P = 0.217.")
}
if(!interactive() && basename(main_script())=="figure4.R") write_plot(build_figure4(),"Figure4",17.1,11.3)
