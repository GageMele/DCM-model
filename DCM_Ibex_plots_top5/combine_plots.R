# Cross-model comparison figures from the per-case outputs of plot_case.R.
# Usage: Rscript combine_plots.R <outdir>
outdir<-commandArgs(TRUE)[1];suppressPackageStartupMessages(library(ggplot2))
tags<-list.dirs(outdir,full.names=FALSE,recursive=FALSE)
tags<-tags[file.exists(file.path(outdir,tags,"validation_rows.csv"))]
if(!length(tags))stop("No finished cases in ",outdir)
rows<-do.call(rbind,lapply(tags,function(t)read.csv(file.path(outdir,t,"validation_rows.csv"))))
rows$model<-factor(rows$model,levels=unique(rows$model))
rows$split<-factor(rows$split,levels=c("seen dates","unseen dates"))
cmp<-file.path(outdir,"compare");dir.create(cmp,showWarnings=FALSE)
ratio_axis<-function(x)pmax(x+1,.01)

# Summary table of the plotted held-out predictions.
s<-do.call(rbind,lapply(split(rows,list(rows$model,rows$split),drop=TRUE),function(r)data.frame(
  model=r$model[1],split=r$split[1],rows=nrow(r),casts=length(unique(r$cast_key)),
  rmse=sqrt(mean((r$chl_norm-r$mean)^2)),mae=mean(abs(r$chl_norm-r$mean)),
  cov50=mean(r$chl_norm>=r$q25&r$chl_norm<=r$q75),cov90=mean(r$chl_norm>=r$q05&r$chl_norm<=r$q95),
  width90=median(r$q95-r$q05),pit_below05=mean(r$pit<.05),pit_above95=mean(r$pit>.95))))
write.csv(s,file.path(cmp,"heldout_summary.csv"),row.names=FALSE);print(s,digits=3,row.names=FALSE)

# PIT histograms, model x split.
p<-ggplot(rows,aes(pit))+geom_histogram(breaks=seq(0,1,.05),fill="grey60",colour="white")+
  facet_grid(model~split,scales="free_y")+theme_bw(base_size=9)+
  labs(x="PIT = P(replicate < observed)",y="Rows",title="Calibration of held-out predictions",
    subtitle="Flat = calibrated; U shape = predictive intervals too narrow; hump = too wide")
ggsave(file.path(cmp,"pit_by_model.png"),p,width=9,height=2+1.6*nlevels(rows$model),dpi=130,bg="white")

# 90% interval coverage by depth (target 0.90).
rows$dbin<-floor(rows$depth/10)*10+5
cv<-aggregate(cbind(cov90=chl_norm>=q05&chl_norm<=q95,cov50=chl_norm>=q25&chl_norm<=q75)~model+split+dbin,rows,mean)
p<-ggplot(cv,aes(cov90,dbin,colour=model))+geom_vline(xintercept=.9,linetype=2)+geom_path()+geom_point(size=1)+
  scale_y_reverse()+facet_wrap(~split)+theme_bw(base_size=10)+
  labs(x="Fraction of held-out rows inside the 90% predictive interval",y="Depth (m)",colour=NULL,
    title="Predictive interval coverage by depth (dashed = nominal 0.90)")
ggsave(file.path(cmp,"coverage90_by_depth.png"),p,width=10,height=6,dpi=130,bg="white")

# Pooled predictive density vs observed, model x split.
dens<-do.call(rbind,lapply(tags,function(t){
  r<-read.csv(file.path(outdir,t,"validation_rows.csv"),nrows=1);z<-readRDS(file.path(outdir,t,"validation_yrep40.rds"))
  nd<-min(25,nrow(z$yrep));data.frame(model=r$model,split=r$split,draw=rep(seq_len(nd),each=ncol(z$yrep)),value=as.vector(t(z$yrep[seq_len(nd),,drop=FALSE])))}))
dens$model<-factor(dens$model,levels=levels(rows$model));dens$split<-factor(dens$split,levels=levels(rows$split))
p<-ggplot()+geom_density(data=dens,aes(ratio_axis(value),group=draw),colour="steelblue",alpha=.2,linewidth=.2)+
  geom_density(data=rows,aes(ratio_axis(chl_norm)),colour="black",linewidth=.7)+scale_x_log10()+
  facet_grid(model~split)+theme_bw(base_size=9)+
  labs(x="Chl / surface Chl (log scale)",y="Density",title="Held-out chl distribution: observed (black) vs 25 posterior predictive replicates (blue)")
ggsave(file.path(cmp,"ppc_density_by_model.png"),p,width=9,height=2+1.6*nlevels(rows$model),dpi=130,bg="white")

# Example held-out casts with every model's posterior mean, one strong-DCM cast per region.
for(sp in levels(droplevels(rows$split))) {
  r<-rows[rows$split==sp,,drop=FALSE];one<-r[r$model==levels(droplevels(r$model))[1],]
  pk<-tapply(one$chl_norm,one$cast_key,max);rg<-tapply(as.character(one$Region),one$cast_key,`[`,1)
  pick<-unlist(lapply(split(names(pk),rg[names(pk)]),function(k){p<-sort(pk[k]);names(p)[ceiling(.75*length(p))]}))
  ex<-r[r$cast_key %in% pick,];ex$panel<-paste0(ex$Region,"\n",substr(ex$date,1,10)," MLD ",round(ex$MLD_temp_MIX),"m")
  obs<-ex[ex$model==ex$model[1],]
  p<-ggplot(ex,aes(y=depth))+geom_point(data=obs,aes(x=chl_norm),size=.4,colour="grey30")+
    geom_path(aes(x=mean,colour=model),linewidth=.6)+scale_y_reverse()+facet_wrap(~panel,scales="free_x",ncol=4)+
    theme_bw(base_size=9)+theme(legend.position="bottom")+
    labs(x="chl_norm",y="Depth (m)",colour=NULL,title=paste("Held-out casts,",sp,": observed (points) and each model's posterior mean"))
  ggsave(file.path(cmp,paste0("profiles_",sub(" ","_",sp),".png")),p,width=12,height=6.5,dpi=130,bg="white")
  # Across-cast median and 90th percentile per region: observed vs each model.
  qp<-do.call(rbind,lapply(tags,function(t){f<-file.path(outdir,t,"validation_quantile_profiles.csv")
    m<-read.csv(file.path(outdir,t,"validation_rows.csv"),nrows=1);if(m$split!=sp)return(NULL);x<-read.csv(f);x$model<-m$model;x}))
  if(!is.null(qp)) {
    qp$model<-factor(qp$model,levels=levels(rows$model));ob<-qp[qp$model==qp$model[1],]
    p<-ggplot(qp,aes(y=depth))+geom_path(aes(x=pred_q50,colour=model))+geom_path(aes(x=pred_q90,colour=model),linetype=2)+
      geom_path(data=ob,aes(x=obs_q50),colour="black",linewidth=.9)+geom_path(data=ob,aes(x=obs_q90),colour="black",linewidth=.9,linetype=2)+
      scale_y_reverse()+facet_wrap(~Region,scales="free_x",ncol=4)+theme_bw(base_size=9)+theme(legend.position="bottom")+
      labs(x="chl_norm",y="Depth (m)",colour=NULL,title=paste("Across-cast median (solid) and 90th percentile (dashed),",sp),
        subtitle="Black = observed held-out casts; colours = posterior predictive")
    ggsave(file.path(cmp,paste0("quantile_profiles_",sub(" ","_",sp),".png")),p,width=12,height=6.5,dpi=130,bg="white")
  }
}
cat("Comparison figures in",cmp,"\n")
