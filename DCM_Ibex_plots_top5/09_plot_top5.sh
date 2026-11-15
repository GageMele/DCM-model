#!/bin/bash
# 09_plot_top5.sh: posterior plots for the 5 best model setups (read-only on all experiments).
# Submits one small Slurm job per fitted case plus one job that builds the comparison figures.
# Nothing in any experiment directory is changed; outputs go to plots/top5_<stamp>/.
set -euo pipefail
root=/ibex/user/meleg/DCM_MODEL; cd "$root"
n150=$(cat logs/latest_next150_experiment.txt); b1m=$(cat logs/latest_bins1m_experiment.txt)
out=$root/plots/top5_$(date -u +%Y%m%dT%H%M%SZ); mkdir -p "$out/logs"
cat > "$out/plot_case.R" <<'EOF_PLOT'
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
EOF_PLOT
cat > "$out/combine_plots.R" <<'EOF_PLOT'
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
EOF_PLOT
# experiment <tab> case <tab> output tag <tab> model name
printf '%s\t%s\t%s\t%s\n' \
  "$n150" 1  LN_half_seen        "LN half-shift" \
  "$n150" 2  LN_half_unseen      "LN half-shift" \
  "$n150" 3  LN_sigdepth_seen    "LN half + sigma by depth" \
  "$n150" 4  LN_sigdepth_unseen  "LN half + sigma by depth" \
  "$n150" 11 LN_distvary_seen    "LN half + distance by depth" \
  "$n150" 12 LN_distvary_unseen  "LN half + distance by depth" \
  "$n150" 13 LN_peak_seen        "LN peak profile" \
  "$b1m"  7  LN_hist1m_seen      "LN 1 m histogram" \
  "$b1m"  8  LN_hist1m_unseen    "LN 1 m histogram" > "$out/models.tsv"
echo "Checking fits"
ok=1
while IFS=$'\t' read -r exp id tag name; do
  f=$(printf '%s/runs/case_%03d/fit.rds' "$exp" "$id")
  if [ -e "$f" ]; then printf '  ok   %-20s %s (%s)\n' "$tag" "$f" "$(du -shL "$f" | cut -f1)"; else echo "  MISSING $f"; ok=0; fi
done < "$out/models.tsv"
[ $ok = 1 ] || { echo "A fit is missing: nothing submitted. Paste this output to Claude."; exit 2; }
n=$(wc -l < "$out/models.tsv")
cat > "$out/plot_one.slurm" <<EOF_SL
#!/bin/bash
source "$n150/ibex_environment.sh"
line=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$out/models.tsv")
IFS=\$'\t' read -r exp id tag name <<< "\$line"
Rscript --vanilla "$out/plot_case.R" "\$exp" "\$id" "$out" "\$tag" "\$name"
EOF_SL
cat > "$out/combine.slurm" <<EOF_SL
#!/bin/bash
source "$n150/ibex_environment.sh"
Rscript --vanilla "$out/combine_plots.R" "$out"
EOF_SL
arr=$(sbatch --parsable --job-name=dcm-plot-top5 --partition=batch --array=1-$n --cpus-per-task=2 --mem=64G \
  --time=12:00:00 -o "$out/logs/case_%a.out" -e "$out/logs/case_%a.err" "$out/plot_one.slurm")
cmb=$(sbatch --parsable --job-name=dcm-plot-compare --partition=batch --dependency=afterany:$arr --cpus-per-task=1 --mem=16G \
  --time=02:00:00 -o "$out/logs/combine.out" -e "$out/logs/combine.err" "$out/combine.slurm")
echo "$out" > logs/latest_plots_top5.txt
echo; echo "Submitted plot array $arr ($n cases) and comparison job $cmb"
echo "Output folder: $out"
