.args <- commandArgs(FALSE)
.file <- .args[grepl("^--file=",.args)]
if(length(.file)) .rdir <- dirname(normalizePath(gsub("~+~"," ",sub("^--file=","",.file[1]),fixed=TRUE))) else {
 .src <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
 .rdir <- if(!is.null(.src))dirname(normalizePath(.src)) else file.path(getwd(),"r")
}
source(file.path(.rdir,"common.R"))
build_supplement <- function() {
 loci<-read_tsv("data/frozen/neighborhoods/selected_loci.tsv");g<-read_tsv("data/frozen/neighborhoods/locus_genes.tsv")
 pf<-c("PF01638","PF04072","PF13577","PF13676","PF22166","PF13413","PF18075","PF01745","PF01337","PF14470")
 stopifnot(nrow(loci)==30,nrow(g)==206,length(unique(loci$HMT))==16,length(unique(loci$GCA))==25,setequal(loci$Pfam,pf))
 loci<-loci[order(match(loci$Pfam,pf),loci$display_order),];loci$y<-rep(3:1,10)
 g$y<-loci$y[match(g$locus_id,loci$locus_id)];g$track<-g$locus_id
 g$xstart<-g$display_start;g$xend<-g$display_end;g$focal<-as.character(g$focal)%in%c("True","TRUE")
 prefs<-list(PF01638=c("quinone","disulfide reductase","nadh oxidase","nitroreductase"),PF04072="malate dehydrogenase",PF13577=c("multidrug resistance","mdtd"),PF13676=c("defense","restriction","nadase"),PF22166=c("quinone","methyltransferase","epimerase"),PF13413=c("rodz","cell shape"),PF18075="ftse",PF01745="miab",PF01337="ribonuclease",PF14470=c("phage","sh oct","shoct"))
 text<-tolower(paste(ifelse(is.na(g$product),"",g$product),ifelse(is.na(g$gene),"",g$gene)))
 preferred<-vapply(seq_len(nrow(g)),function(i) !g$focal[i]&&abs(g$CDS_offset[i])==1&&any(vapply(prefs[[g$Pfam[i]]],function(k)grepl(k,text[i],fixed=TRUE),logical(1))),logical(1))
 g$role<-ifelse(g$focal,"Domain-bearing CDS",ifelse(preferred,"Prioritised neighbour",ifelse(grepl("hypothetical",text),"Uncharacterised CDS","Other annotated CDS")))
 panels<-paste0(LETTERS[1:10],"  ",unname(DOMAIN_NAMES[pf])," | ",pf,"\n",ifelse(seq_along(pf)<=7,"Plaque-oriented","Mucosa-oriented")," in exploratory screen")
 loci$panel<-factor(panels[match(loci$Pfam,pf)],levels=panels);g$panel<-factor(panels[match(g$Pfam,pf)],levels=panels)
 polys<-gene_polygons(g,c(-3500,5500),.14);polys$panel<-g$panel[polys$id]
 loci$taxon<-paste(loci$Genus,loci$Species);long<-nchar(loci$taxon)>38;loci$taxon[long]<-paste0(substr(loci$Genus[long],1,1),". ",loci$Species[long]);loci$label<-paste(loci$taxon,loci$HMT)
 loci$habitat<-ifelse(loci$Ecology=="plaque","P","M")
 g$label<-ifelse(is.na(g$gene),"",sub("_[0-9]+$","",g$gene));g$label[g$focal]<-unname(DOMAIN_NAMES[g$Pfam[g$focal]])
 maps<-c(ribonuclease="RNase",quinone="quinone red.","disulfide reductase"="disulfide red.","nadh oxidase"="NADH ox.",methyltransferase="methyltransf.","malate dehydrogenase"="malate DH",multidrug="MFS transporter")
 for(k in names(maps)){i<-which(!g$focal&g$label==""&grepl(k,text,fixed=TRUE));g$label[i]<-maps[k]}
 label_rows<-g[(g$focal|preferred|abs(g$CDS_offset)==1)&g$label!="",];label_rows$mid<-(label_rows$xstart+label_rows$xend)/2
 label_rows<-label_rows[label_rows$mid > -3400 & label_rows$mid < 5400,];label_rows$label_y<-label_rows$y+ifelse(label_rows$focal,.23,-.23)
 clips<-rbind(data.frame(x=-3500,y=g$y[g$xstart< -3500],panel=g$panel[g$xstart< -3500]),data.frame(x=5500,y=g$y[g$xend>5500],panel=g$panel[g$xend>5500]))
 p<-ggplot()+geom_segment(data=loci,aes(x=-3500,xend=5500,y=y,yend=y),colour=STYLE$grid,linewidth=.3)+
  geom_polygon(data=polys,aes(x=x,y=y,group=id,fill=role),colour="white",linewidth=.15)+
  geom_segment(data=clips,aes(x=x,xend=x,y=y-.1,yend=y+.1),linewidth=.5,colour=STYLE$muted)+
  geom_text(data=loci,aes(x=-3500,y=y+.47,label=label),hjust=0,fontface="italic",size=2.5,colour=STYLE$text)+
  geom_text(data=loci,aes(x=5500,y=y+.47,label=habitat,colour=Ecology),hjust=1,size=2.5,fontface="bold")+
  geom_text(data=label_rows,aes(x=mid,y=label_y,label=label),size=STYLE$gene_label_size,colour=STYLE$text)+
  scale_fill_manual(values=c("Domain-bearing CDS"=STYLE$text,"Prioritised neighbour"=STYLE$partner,"Other annotated CDS"=STYLE$other_gene,"Uncharacterised CDS"=STYLE$unknown_gene),name=NULL)+
  habitat_scale()+guides(colour="none",fill=guide_legend(nrow=2,byrow=TRUE))+
  scale_x_continuous(limits=c(-3500,5500),breaks=c(-2000,0,2000,4000),labels=c(-2,0,2,4),expand=expansion(mult=.025))+
  scale_y_continuous(limits=c(.55,3.72),expand=c(0,0))+
  facet_wrap(~panel,ncol=2)+base_theme(9)+
  theme(panel.grid=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(size=7),
   strip.text=element_text(hjust=0,size=8.7,lineheight=1.1),panel.spacing.y=grid::unit(.2,"inches"),panel.spacing.x=grid::unit(.2,"inches"),
   legend.text=element_text(size=8),plot.title=element_text(size=15,face="bold"))+
  labs(x="Position relative to focal CDS start (kb)",y=NULL,
   title="Gene arrangements around ten candidate protein domains",
   subtitle="30 representative loci | 16 HMTs | 25 assemblies | measured CDS coordinates and strands",
   caption=NULL)
 write_source(g,"Supplementary_gene_coordinates.tsv");write_source(loci,"Supplementary_selected_loci.tsv")
 p
}
if(!interactive() && basename(main_script())=="supplement.R") write_plot(build_supplement(),"Supplementary_candidate_domain_neighborhoods",11.5,14)
