# Posterior plots for one fitted case. Read-only: loads fit.rds and the frozen data,
# writes PNGs and CSVs into <outdir>/<tag>/. Nothing in the experiment is changed.
# Usage: Rscript plot_case.R <experiment_dir> <case_id> <outdir> <tag> <model_name>
args<-commandArgs(TRUE)
experiment<-args[1];case_id<-as.integer(args[2]);outdir<-args[3];tag<-args[4];model_name<-args[5]
Sys.setenv(EXPERIMENT_DIR=experiment)
source(file.path(experiment,"common.R"));source(file.path(experiment,"models.R"));source(file.path(experiment,"scoring.R"))
suppressPackageStartupMessages({library(ggplot2)})
cases<-read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
case<-cases[cases$case_id==case_id,,drop=FALSE];stopifnot(nrow(case)==1)
run<-run_directory(experiment,case);out<-file.path(outdir,tag);dir.create(out,recursive=TRUE,showWarnings=FALSE)
split_label<-if(case$split=="primary")"seen dates" else "unseen dates"
cat("Plotting",tag,"=",case$label,"from",run,"\n")
fit<-readRDS(file.path(run,"fit.rds"))
ndraw<-as.integer(Sys.getenv("PLOT_DRAWS","300"))
total<-posterior::ndraws(posterior::as_draws_array(fit));set.seed(5150);ids<-sort(sample.int(total,min(ndraw,total)))
sets<-evaluation_sets(experiment,case)
sets<-list(training=sets[[1]],validation=sets$validation)
title_base<-paste0(model_name," (",split_label,")")
ratio_axis<-function(x)pmax(x+1,.01)   # chl / surface chl, floored for log axes

for(ds in names(sets)) {
  d<-sets[[ds]];d<-d[order(d$cast_key,d$depth),,drop=FALSE]
  pred<-d[setdiff(names(d),c("chl_norm","ratio"))]
  ds_label<-if(ds=="training")"training casts (in-sample fit)" else paste("held-out casts,",split_label)
  # Draw unseen date/ENSO groups once per date group, as the scoring does.
  mu<-y<-matrix(NA_real_,length(ids),nrow(d))
  for(ii in split(seq_len(nrow(d)),group_key(d))) {
    set.seed(9001+ii[1])
    mu[,ii]<-brms::posterior_epred(fit,newdata=pred[ii,,drop=FALSE],draw_ids=ids,re_formula=NULL,
      allow_new_levels=TRUE,sample_new_levels="gaussian")
    set.seed(9001+ii[1])
    y[,ii]<-brms::posterior_predict(fit,newdata=pred[ii,,drop=FALSE],draw_ids=ids,re_formula=NULL,
      allow_new_levels=TRUE,sample_new_levels="gaussian")
  }
  if(case$family=="lognormal"){mu<-mu-1;y<-y-1}
  q<-apply(y,2,quantile,c(.05,.25,.5,.75,.95))
  rows<-data.frame(model=model_name,tag=tag,split=split_label,dataset=ds,
    d[c("cast_key","Region","Ocean","date","depth","MLD_temp_MIX","chl_norm")],
    mean=colMeans(mu),mean_lo=apply(mu,2,quantile,.05),mean_hi=apply(mu,2,quantile,.95),
    q05=q[1,],q25=q[2,],q50=q[3,],q75=q[4,],q95=q[5,],pit=colMeans(y<matrix(d$chl_norm,nrow(y),ncol(y),byrow=TRUE)))
  write.csv(rows,file.path(out,paste0(ds,"_rows.csv")),row.names=FALSE)
  saveRDS(list(rows=rows[c("Region","depth","chl_norm")],yrep=y[seq_len(min(40,nrow(y))),,drop=FALSE]),file.path(out,paste0(ds,"_yrep40.rds")))

  # 1. Example profiles: per region, a typical cast and a strong-DCM cast.
  peak<-tapply(d$chl_norm,as.character(d$cast_key),max,na.rm=TRUE)
  reg<-tapply(as.character(d$Region),as.character(d$cast_key),`[`,1)
  pick<-unlist(lapply(split(names(peak),reg[names(peak)]),function(k){
    p<-sort(peak[k]);unique(c(names(p)[ceiling(length(p)/2)],names(p)[ceiling(.9*length(p))]))}))
  ex<-rows[rows$cast_key %in% pick,,drop=FALSE]
  ex$panel<-paste0(ex$Region,"\n",substr(ex$date,1,10)," MLD ",round(ex$MLD_temp_MIX),"m")
  p<-ggplot(ex,aes(y=depth))+
    geom_ribbon(aes(xmin=q05,xmax=q95),orientation="y",fill="steelblue",alpha=.18)+
    geom_ribbon(aes(xmin=q25,xmax=q75),orientation="y",fill="steelblue",alpha=.3)+
    geom_point(aes(x=chl_norm),size=.4,colour="grey20")+
    geom_path(aes(x=mean),colour="firebrick",linewidth=.7)+
    scale_y_reverse()+facet_wrap(~panel,scales="free_x",ncol=4)+theme_bw(base_size=9)+
    labs(x="chl_norm = (Chl - surface)/surface",y="Depth (m)",title=title_base,
      subtitle=paste0(ds_label,". Points = observed; red = posterior mean; bands = 50% and 90% posterior predictive"))
  ggsave(file.path(out,paste0(ds,"_profiles.png")),p,width=12,height=2.6*ceiling(length(pick)/4),dpi=130,bg="white",limitsize=FALSE)

  # 2. Posterior predictive density vs observed, by region (log ratio axis).
  nd<-min(40,nrow(y));dens<-data.frame(Region=rep(d$Region,nd),draw=rep(seq_len(nd),each=nrow(d)),value=as.vector(t(y[seq_len(nd),,drop=FALSE])))
  p<-ggplot()+geom_density(data=dens,aes(ratio_axis(value),group=draw),colour="steelblue",alpha=.15,linewidth=.2)+
    geom_density(data=d,aes(ratio_axis(chl_norm)),colour="black",linewidth=.8)+
    scale_x_log10()+facet_wrap(~Region,scales="free_y",ncol=4)+theme_bw(base_size=9)+
    labs(x="Chl / surface Chl (log scale)",y="Density",title=title_base,
      subtitle=paste0(ds_label,". Black = observed distribution; blue = ",nd," posterior predictive replicates"))
  ggsave(file.path(out,paste0(ds,"_ppc_density_region.png")),p,width=12,height=6,dpi=130,bg="white")

  # 3. Same, by depth band.
  band<-cut(d$depth,c(9,30,50,70,90,110,130,151),labels=c("10-30","30-50","50-70","70-90","90-110","110-130","130-150"))
  dens$band<-rep(band,nd);dd<-data.frame(band=band,chl_norm=d$chl_norm)
  p<-ggplot()+geom_density(data=dens,aes(ratio_axis(value),group=draw),colour="steelblue",alpha=.15,linewidth=.2)+
    geom_density(data=dd,aes(ratio_axis(chl_norm)),colour="black",linewidth=.8)+
    scale_x_log10()+facet_wrap(~band,scales="free_y",ncol=4)+theme_bw(base_size=9)+
    labs(x="Chl / surface Chl (log scale)",y="Density",title=title_base,subtitle=paste0(ds_label,", by depth band (m)"))
  ggsave(file.path(out,paste0(ds,"_ppc_density_depth.png")),p,width=12,height=6,dpi=130,bg="white")

  # 4. Observed vs posterior predictive quantiles across casts, by region and 5 m band.
  key<-interaction(d$Region,floor(d$depth/5)*5,drop=TRUE)
  qp<-do.call(rbind,lapply(split(seq_len(nrow(d)),key),function(ii){
    if(length(ii)<5)return(NULL);s<-y[,ii,drop=FALSE];md<-apply(s,1,median);q9<-apply(s,1,quantile,.9)
    data.frame(Region=d$Region[ii[1]],depth=floor(d$depth[ii[1]]/5)*5+2.5,n=length(ii),
      obs_q10=quantile(d$chl_norm[ii],.1),obs_q50=median(d$chl_norm[ii]),obs_q90=quantile(d$chl_norm[ii],.9),
      pred_q10=quantile(s,.1),pred_q50=quantile(s,.5),pred_q90=quantile(s,.9),
      pred_q50_lo=quantile(md,.05),pred_q50_hi=quantile(md,.95),pred_q90_lo=quantile(q9,.05),pred_q90_hi=quantile(q9,.95))}))
  write.csv(qp,file.path(out,paste0(ds,"_quantile_profiles.csv")),row.names=FALSE)
  p<-ggplot(qp,aes(y=depth))+
    geom_ribbon(aes(xmin=pred_q50_lo,xmax=pred_q50_hi),orientation="y",fill="steelblue",alpha=.3)+
    geom_ribbon(aes(xmin=pred_q90_lo,xmax=pred_q90_hi),orientation="y",fill="darkorange",alpha=.3)+
    geom_path(aes(x=pred_q10),colour="grey50",linetype=2)+geom_path(aes(x=pred_q50),colour="steelblue")+geom_path(aes(x=pred_q90),colour="darkorange")+
    geom_path(aes(x=obs_q10),colour="black",linetype=2)+geom_path(aes(x=obs_q50),colour="black",linewidth=.8)+geom_path(aes(x=obs_q90),colour="black",linewidth=.8)+
    scale_y_reverse()+facet_wrap(~Region,scales="free_x",ncol=4)+theme_bw(base_size=9)+
    labs(x="chl_norm",y="Depth (m)",title=title_base,
      subtitle=paste0(ds_label,". Across-cast 10/50/90% quantiles per 5 m: black = observed, colour = posterior predictive (bands = 90% uncertainty)"))
  ggsave(file.path(out,paste0(ds,"_quantile_profiles.png")),p,width=12,height=6.5,dpi=130,bg="white")

  # 5. Calibration: PIT histogram (flat = calibrated) and observed vs predicted.
  p<-ggplot(rows,aes(pit))+geom_histogram(breaks=seq(0,1,.05),fill="grey60",colour="white")+
    geom_hline(yintercept=nrow(rows)/20,colour="firebrick",linetype=2)+facet_wrap(~Region,ncol=4,scales="free_y")+theme_bw(base_size=9)+
    labs(x="PIT = P(replicate < observed)",y="Rows",title=title_base,
      subtitle=paste0(ds_label,". Flat = calibrated; U shape = intervals too narrow; hump = too wide"))
  ggsave(file.path(out,paste0(ds,"_pit.png")),p,width=12,height=6,dpi=130,bg="white")
  p<-ggplot(rows,aes(ratio_axis(mean),ratio_axis(chl_norm)))+geom_point(alpha=.15,size=.4)+
    geom_abline(colour="firebrick")+scale_x_log10()+scale_y_log10()+facet_wrap(~Region,ncol=4)+theme_bw(base_size=9)+
    labs(x="Posterior mean Chl/surface",y="Observed Chl/surface",title=title_base,subtitle=ds_label)
  ggsave(file.path(out,paste0(ds,"_obs_vs_pred.png")),p,width=12,height=6,dpi=130,bg="white")
  cat(" ",ds,nrow(d),"rows,",length(unique(d$cast_key)),"casts; cov90",round(mean(rows$chl_norm>=rows$q05&rows$chl_norm<=rows$q95),3),"\n")
  rm(mu,y);invisible(gc())
}

# 6. Parameter posteriors: population-level coefficients and group/smooth SDs
# (the per-depth-bin sigma coefficients of the 1 m model are summarised, not listed).
dr<-posterior::as_draws_df(fit);v<-grep("^(b_|bsp_|bs_|sd_|sds_|sigma|shape)",posterior::variables(dr),value=TRUE)
v<-v[!grepl("^b_sigma_depth_(m|bin)",v)]
if(length(v)) {
  s<-posterior::summarise_draws(posterior::subset_draws(dr,variable=v),median=median,
    ~quantile(.x,c(.05,.25,.75,.95)),rhat=posterior::rhat,ess_bulk=posterior::ess_bulk)
  names(s)<-c("variable","median","q05","q25","q75","q95","rhat","ess_bulk")
  write.csv(s,file.path(out,"parameter_summary.csv"),row.names=FALSE)
  s$type<-sub("_.*","",s$variable);s$variable<-factor(s$variable,levels=rev(s$variable))
  p<-ggplot(s,aes(y=variable))+geom_vline(xintercept=0,colour="grey70")+
    geom_linerange(aes(xmin=q05,xmax=q95),linewidth=.4)+geom_linerange(aes(xmin=q25,xmax=q75),linewidth=1.2,colour="steelblue")+
    geom_point(aes(x=median),size=1)+facet_wrap(~type,scales="free",ncol=3)+theme_bw(base_size=8)+
    labs(x="Posterior (median, 50% and 90% intervals)",y=NULL,title=paste(title_base,"- parameter posteriors"),
      subtitle="b = population coefficients; sd/sds = random-effect and smooth SDs; sigma/shape = residual scale")
  ggsave(file.path(out,"parameters.png"),p,width=13,height=max(5,.16*length(v)+2),dpi=130,bg="white",limitsize=FALSE)
}
writeLines(c(case$label,run),file.path(out,"SOURCE.txt"))
cat("Done",tag,"\n")
