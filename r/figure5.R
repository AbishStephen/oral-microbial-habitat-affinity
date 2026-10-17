.args <- commandArgs(FALSE)
.file <- .args[grepl("^--file=",.args)]
if(length(.file)) .rdir <- dirname(normalizePath(gsub("~+~"," ",sub("^--file=","",.file[1]),fixed=TRUE))) else {
 .src <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
 .rdir <- if(!is.null(.src))dirname(normalizePath(.src)) else file.path(getwd(),"r")
}
source(file.path(.rdir,"common.R"))
build_figure5 <- function() {
 j<-read_json("figure_sources/Figure5/Figure5_model_predictions_and_sister_pairs.json");pred<-records(j$pred);auc<-records(j$auc)
 stopifnot(nrow(pred)==47,length(unique(pred$Genus))==13,nrow(auc)==13)
 boxes<-data.frame(x=c(0,.205,.41,.615,.82),label=c("Hold out one genus\n13 outer folds","Rank features\ntraining genera only","Tune model\ninner genus folds","Fit training data\npredict held-out HMTs","Evaluate transfer\nfull-refitting nulls"))
 pa<-ggplot(boxes)+geom_rect(aes(xmin=x,xmax=x+.17,ymin=.2,ymax=.8),fill=STYLE$pale,colour=STYLE$grid)+
  geom_text(aes(x=x+.085,y=.5,label=label),size=3.2,family=STYLE$font_family)+
  geom_segment(data=data.frame(x=c(.17,.375,.58,.785)),aes(x=x,xend=x+.035,y=.5,yend=.5),arrow=grid::arrow(length=grid::unit(.09,"inches")),colour=STYLE$muted)+
  xlim(0,1)+ylim(0,1)+theme_void(base_family=STYLE$font_family)+base_theme()+
  theme(panel.grid=element_blank(),axis.text=element_blank(),axis.title=element_blank())+
  labs(title="Genus-held-out prediction design",subtitle="Held-out genera are excluded from feature ranking, tuning, scaling and fitting",tag="A")
 pred$x<-match(pred$Genus,GENERA)
 pred$offset<-ave(pred$x,paste(pred$Genus,pred$Ecology),FUN=function(x)if(length(x)==1)0 else seq(-.19,.19,length.out=length(x)))
 means<-aggregate(Y_pred~Genus+Ecology,pred,mean);ys<-tapply(pred$Y_pred,pred$Genus,range)
 pb<-ggplot(pred,aes(x=x+offset,y=Y_pred,colour=Ecology))+geom_hline(yintercept=0,linetype=2,colour=STYLE$muted,linewidth=.35)+
  geom_point(size=2.7,alpha=.95)+habitat_scale()+scale_x_continuous(breaks=1:13,labels=GENERA,limits=c(.4,13.6))+
  base_theme()+theme(axis.text.x=element_text(angle=30,hjust=1,size=9))+
  labs(x=NULL,y="Held-out predicted ecological affinity (a.u.)",title="Predictions for all 47 held-out HMTs",subtitle="Source-family-balanced outcome; each point is one taxon",tag="B")
 auc$Genus<-factor(auc$Genus,levels=rev(GENERA))
 aa<-rbind(data.frame(Genus=auc$Genus,AUC=auc$sourcebalanced_genusbalanced,model="Unadjusted"),data.frame(Genus=auc$Genus,AUC=auc$sourcebalanced_genusbalanced_covadj,model="Genome-adjusted"))
 pc<-ggplot(auc,aes(y=Genus))+geom_vline(xintercept=.5,linetype=2,colour=STYLE$muted)+
  geom_segment(aes(x=sourcebalanced_genusbalanced,xend=sourcebalanced_genusbalanced_covadj,yend=Genus),linewidth=.8,colour="#ADC5CD")+
  geom_point(data=aa,aes(x=AUC,colour=model,shape=model),size=2.3)+
  scale_colour_manual(values=c(Unadjusted=STYLE$plaque,"Genome-adjusted"=STYLE$mucosal),name=NULL)+
  scale_shape_manual(values=c(Unadjusted=16,"Genome-adjusted"=18),name=NULL)+
  scale_x_continuous(limits=c(-.05,1.05),breaks=c(0,.25,.5,.75,1))+base_theme()+
  theme(axis.text.y=element_text(face="italic",size=8.5))+
  labs(x="Within-genus paired AUC",y=NULL,title="Transfer varies across genera",subtitle="Source-balanced models; dashed line: chance",tag="C")
 scores<-read_tsv("figure_sources/Figure5/model_performance_and_uncertainty.tsv")
 models<-c("Genus-balanced nested","Genus-balanced nested, genome covariate-residualized")
 take<-unlist(lapply(c("original","source-family-balanced"),function(e) vapply(models,function(mm)which(scores$ecology_outcome==e&scores$model==mm),integer(1))))
 scores<-scores[take,];scores$label<-factor(c("Original","Original, adjusted","Source balanced","Source balanced, adjusted"),levels=rev(c("Original","Original, adjusted","Source balanced","Source balanced, adjusted")))
 scores$adjustment<-ifelse(grepl("residualized",scores$model),"Genome-adjusted","Unadjusted")
 pd<-ggplot(scores,aes(y=label,colour=adjustment))+geom_vline(xintercept=.5,linetype=2,colour=STYLE$muted)+
  geom_segment(aes(x=genus_bootstrap_95_low,xend=genus_bootstrap_95_high,yend=label),linewidth=1)+
  geom_point(aes(x=macro_within_genus_AUC),size=3)+
  geom_text(aes(x=1.04,label=sprintf("P = %.3f",permutation_empirical_p)),hjust=0,size=3,colour=STYLE$text)+
  scale_colour_manual(values=c(Unadjusted=STYLE$plaque,"Genome-adjusted"=STYLE$mucosal))+
  scale_x_continuous(limits=c(.35,1.28),breaks=c(.5,.7,.9,1))+base_theme()+theme(legend.position="none")+
  labs(x="Macro AUC with genus-bootstrap 95% interval",y=NULL,title="Prediction and covariate sensitivity",subtitle="Full-refitting permutation P values",tag="D")
 write_source(pred,"Figure5_heldout_predictions.tsv");write_source(auc,"Figure5_genus_AUCs.tsv");write_source(scores,"Figure5_model_performance.tsv")
 (pa/pb/(pc|pd))+plot_layout(heights=c(.65,1.2,1.25))+
  annotation("","",
   "Unadjusted nulls: 200 complete refits per outcome. Genome-adjusted nulls: 100 refits. HMTs, rather than individual genomes, are the prediction units.")
}
if(!interactive() && basename(main_script())=="figure5.R") write_plot(build_figure5(),"Figure5",17,11.5)
