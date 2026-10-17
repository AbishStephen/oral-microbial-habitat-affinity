.args <- commandArgs(FALSE)
.file <- .args[grepl("^--file=",.args)]
if(length(.file)) .rdir <- dirname(normalizePath(gsub("~+~"," ",sub("^--file=","",.file[1]),fixed=TRUE))) else {
 .src <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
 .rdir <- if(!is.null(.src))dirname(normalizePath(.src)) else file.path(getwd(),"r")
}
source(file.path(.rdir,"common.R"))
build_figure1 <- function() {
 m <- read_tsv("figure_sources/Figure1/HOMD47_HMT_manifest.tsv")
 m$Ecology <- ifelse(m$Ecology=="p","plaque","mucosal")
 m$Ecology_tier <- ifelse(m$Ecology_tier=="H","High","Moderate")
 tr <- ape::read.tree(file.path(REPO,"figure_sources/Figure1/HOMD47_ribosomal_tree.newick"))
 hmt <- sub("\\|.*","",tr$tip.label);m <- m[match(hmt,m$HMT),]
 stopifnot(nrow(m)==47,length(unique(m$Genus))==13,sum(m$n_genomes)==1121,!anyNA(m$HMT))
 m$y <- rev(seq_len(nrow(m)));m$Phylum <- unname(PHYLA[m$Genus])
 m$taxon <- paste(m$Genus,m$Species," | ",m$HMT,sep=" ")
 depth <- ape::node.depth.edgelength(tr);ys <- numeric(length(depth));ys[seq_len(47)] <- m$y
 calc_y <- function(node) {if(node<=47) return(ys[node]);kids <- tr$edge[tr$edge[,1]==node,2];ys[node]<<-mean(vapply(kids,calc_y,numeric(1)));ys[node]}
 calc_y(48)
 e <- data.frame(x=depth[tr$edge[,1]],xend=depth[tr$edge[,2]],y=ys[tr$edge[,2]],yend=ys[tr$edge[,2]])
 ve <- do.call(rbind,lapply(unique(tr$edge[,1]),function(n) {kids<-tr$edge[tr$edge[,1]==n,2];data.frame(x=depth[n],xend=depth[n],y=min(ys[kids]),yend=max(ys[kids]))}))
 maxd <- max(depth);m$x <- depth[seq_len(47)];labelx <- maxd*1.03
 pt <- ggplot()+geom_segment(data=rbind(e,ve),aes(x=x,xend=xend,y=y,yend=yend),colour="#8CA2AC",linewidth=.3)+
  geom_point(data=m,aes(x=x,y=y,colour=Ecology),size=1.3)+
  geom_text(data=m,aes(x=labelx,y=y,label=taxon),hjust=0,size=2.5,family=STYLE$font_family)+
  habitat_scale()+scale_x_continuous(limits=c(0,maxd*2.08),breaks=c(0,.5,1),name="Substitutions per site")+
  scale_y_continuous(limits=c(.3,47.7),expand=c(0,0))+base_theme()+
  theme(axis.text.y=element_blank(),axis.title.y=element_blank(),panel.grid=element_blank(),legend.position="none")+
  labs(title="Phylogeny and ecological affinity",subtitle="47 taxa across six bacterial phyla",tag="A")
 pa <- ggplot(m,aes(y=y,x=Median_Y,colour=Ecology))+geom_vline(xintercept=0,colour=STYLE$muted,linewidth=.35)+
  geom_segment(aes(x=0,xend=Median_Y,yend=y),linewidth=.45,alpha=.6)+
  geom_point(aes(shape=Ecology_tier),size=2.1,stroke=.65,fill="white")+
  geom_text(aes(x=13,label=n_genomes),colour=STYLE$muted,hjust=1,size=2.4)+
  habitat_scale()+scale_shape_manual(values=c(High=16,Moderate=1),name="Confidence")+
  scale_x_continuous(limits=c(-10,13.5),breaks=c(-8,-4,0,4,8),name="Plaque-to-mucosal affinity (log2)")+
  scale_y_continuous(limits=c(.3,47.7),expand=c(0,0))+base_theme()+
  theme(axis.text.y=element_blank(),axis.title.y=element_blank())+
  labs(title="Habitat score",subtitle="Right-hand counts: genomes per HMT")
 counts <- as.data.frame(table(m$Genus,m$Ecology),stringsAsFactors=FALSE);names(counts)<-c("Genus","Ecology","n")
 go <- unique(m$Genus);counts$Genus<-factor(counts$Genus,levels=rev(go))
 pc <- ggplot(counts,aes(x=n,y=Genus,fill=Ecology))+geom_col(width=.66)+
  scale_fill_manual(values=c(plaque=STYLE$plaque,mucosal=STYLE$mucosal),name="HMT affinity")+
  scale_x_continuous(breaks=c(0,2,4,6),expand=expansion(mult=c(0,.08)))+base_theme()+
  theme(axis.text.y=element_text(face="italic",size=9))+
  labs(x="Number of HMTs",y=NULL,title="Matched genera",subtitle="Both habitats in every genus",tag="B")
 write_source(m,"Figure1_plotted_HMTs.tsv")
 (pt|pa|pc)+plot_layout(widths=c(2.1,1.05,.9))+
  annotation("","",
   "")
}
if(!interactive() && basename(main_script())=="figure1.R") write_plot(build_figure1(),"Figure1",18.2,13.3)
