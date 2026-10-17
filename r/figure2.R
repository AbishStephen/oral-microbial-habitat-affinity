.args <- commandArgs(FALSE)
.file <- .args[grepl("^--file=",.args)]
if(length(.file)) .rdir <- dirname(normalizePath(gsub("~+~"," ",sub("^--file=","",.file[1]),fixed=TRUE))) else {
 .src <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
 .rdir <- if(!is.null(.src))dirname(normalizePath(.src)) else file.path(getwd(),"r")
}
source(file.path(.rdir,"common.R"))
build_figure2 <- function() {
 j<-read_json("figure_sources/Figure2/Figure2_plotting_data.json");stopifnot(j$n_features==13917,j$n_focal==14,j$n_genera==13)
 boxes<-data.frame(x=c(.0,.34,.68,.0,.34,.68),y=c(.62,.62,.62,.12,.12,.12),
  label=c("Six abundance summaries\nEcological classification","47 HMTs / 13 genera\n1,121 reference genomes","13,917 nonredundant\nfunctional profiles","14 exploratory annotations\n>=5 genera; >=3 phyla","Raw and genome-adjusted\nnominal consistency tests","4,096 within-genus\nlabel permutations"))
 links<-data.frame(x=c(.29,.63,.825,.68,.34),xend=c(.34,.68,.825,.63,.29),y=c(.79,.79,.62,.29,.29),yend=c(.79,.79,.46,.29,.29))
 pa<-ggplot(boxes)+geom_rect(aes(xmin=x,xmax=x+.29,ymin=y,ymax=y+.34),fill=STYLE$pale,colour=STYLE$grid,linewidth=.5)+
  geom_text(aes(x=x+.145,y=y+.17,label=label),size=3.05,lineheight=1.1,colour=STYLE$text,family=STYLE$font_family)+
  geom_segment(data=links,aes(x=x,y=y,xend=xend,yend=yend),arrow=grid::arrow(length=grid::unit(.1,"inches")),colour=STYLE$muted,linewidth=.45)+
  annotate("text",x=0,y=-.02,label="",hjust=0,size=3,colour=STYLE$muted)+
  coord_cartesian(xlim=c(-.01,1.01),ylim=c(-.08,1.05),clip="off")+theme_void(base_family=STYLE$font_family)+
  base_theme()+theme(panel.grid=element_blank(),axis.text=element_blank(),axis.title=element_blank())+
  labs(title="",subtitle="",tag="A")
 focal<-records(j$focal)
 focal<-focal[order(match(focal$direction,c("plaque","mucosal")),-focal$genus_support),]
 short<-focal$annotation;short[is.na(short)|short==""]<-"reaction"
 short[focal$feature=="COG:COG1216"]<-"GT2"
 short[focal$feature=="KMOD:M00545"]<-"cinnamate module"
 sid<-sub("^RHEA:RHEA:","RHEA:",sub("^PFAM:","",sub("^KMOD:","",focal$feature)))
 focal$label<-paste(sid,short)
 focal$label<-factor(focal$label,levels=rev(focal$label))
 pb<-ggplot(focal,aes(y=label))+geom_col(aes(x=genus_support,fill=direction),width=.57)+
  geom_point(aes(x=phyla_support),shape=21,fill="white",colour=STYLE$text,size=2.2,stroke=.6)+
  scale_fill_manual(values=c(plaque=STYLE$plaque,mucosal=STYLE$mucosal),name="Direction")+
  scale_x_continuous(limits=c(0,13),breaks=c(0,3,6,9,13))+base_theme()+
  theme(axis.text.y=element_text(size=8),legend.position="none")+
  labs(x="Support breadth: bars = genera; circles = phyla",y=NULL,title="",subtitle="",tag="B")
 atlas<-records(j$atlas);z<-do.call(rbind,lapply(j$raw41,unlist));colnames(z)<-unlist(j$genera)
 heat<-data.frame(feature=rep(atlas$id,each=13),annotation=rep(atlas$short,each=13),theme=rep(atlas$group,each=13),
  Genus=rep(colnames(z),times=41),effect=as.vector(t(z)))
 labs_y<-paste(sub(".*:","",atlas$id),atlas$short,sep=" | ")
 heat$feature<-factor(heat$feature,levels=rev(atlas$id));heat$Genus<-factor(heat$Genus,levels=GENERA)
 heat$theme<-factor(heat$theme,levels=unique(atlas$group))
 pc<-ggplot(heat,aes(x=Genus,y=feature,fill=effect))+geom_tile(colour="white",linewidth=.15)+
  scale_fill_gradient2(low=STYLE$mucosal,mid="white",high=STYLE$plaque,limits=c(-1,1),name="Plaque - mucosal\nwithin-genus difference",guide=guide_colourbar(display="rectangles"))+
  scale_y_discrete(labels=setNames(labs_y,atlas$id))+
  facet_grid(theme~.,scales="free_y",space="free_y",switch="y",labeller=labeller(theme=function(x)wrap(x,17)))+
  base_theme()+theme(axis.text.x=element_text(angle=40,hjust=1,size=7.5),axis.text.y=element_text(size=7.8),
   strip.text.y.left=element_text(angle=0,size=7.5),strip.placement="outside",panel.spacing.y=grid::unit(.12,"inches"),panel.grid=element_blank())+
  labs(x=NULL,y=NULL,title="Functional diversity across the full cohort",subtitle="",tag="C")
 key<-c("PFAM:PF01638","PFAM:PF04072","PFAM:PF13577","COG:COG1216")
 rz<-do.call(rbind,lapply(j$raw14,unlist));colnames(rz)<-unlist(j$genera)
 dd<-do.call(rbind,lapply(seq_along(key),function(i){fi<-match(key[i],focal$feature);k<-match(key[i],vapply(j$focal,function(r)r$feature,character(1)));data.frame(feature=key[i],annotation=focal$annotation[fi],Genus=colnames(rz),effect=rz[k,],offset=((seq_len(13)*7)%%13-6)/50)}))
 dd$annotation<-factor(dd$annotation,levels=unique(dd$annotation));dd$Phylum<-unname(PHYLA[dd$Genus])
 means<-aggregate(effect~annotation,dd,mean)
 pd<-ggplot(dd,aes(x=effect,y=offset))+geom_vline(xintercept=0,colour=STYLE$grid)+
  geom_point(aes(colour=Phylum),size=2.3,alpha=.9)+geom_point(data=means,aes(x=effect,y=0),shape=23,size=3.1,fill=STYLE$partner,colour=STYLE$text)+
  facet_wrap(~annotation,ncol=1)+scale_colour_manual(values=PHYLUM_COLOURS)+
  scale_x_continuous(limits=c(-1,1),breaks=c(-1,0,1))+scale_y_continuous(limits=c(-.25,.25))+
  base_theme()+theme(axis.text.y=element_blank(),axis.title.y=element_blank(),panel.grid.major.y=element_blank(),legend.text=element_text(size=7))+
  labs(x="Within-genus plaque - mucosal difference",title="",subtitle="",tag="D")
 write_source(heat,"Figure2_heatmap_effects.tsv");write_source(focal,"Figure2_focal_annotations.tsv");write_source(dd,"Figure2_candidate_genus_effects.tsv")
 ((wrap_elements(full=pa)|pb)/(pc|pd))+plot_layout(widths=c(2.2,1),heights=c(1,2.35))+
  annotation("","",
   "")
}
if(!interactive() && basename(main_script())=="figure2.R") write_plot(build_figure2(),"Figure2",19.2,13.8)
