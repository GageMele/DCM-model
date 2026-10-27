# Read-only comparisons from existing CSV outputs. Base R only; no network.
out<-"original_model_review"
historical<-read.csv(file.path(out,"historical_population_profiles_10_150.csv"),check.names=FALSE)
observed<-read.csv(file.path(out,"historical_observed_median_profiles.csv"))
regions<-levels(factor(historical$Region))
# Connect curves in depth order. Sorting a geom_line by chl_norm distorts peaks.
png(file.path(out,"historical_profiles.png"),width=1560,height=1350,res=130)
par(mfrow=c(3,3),mar=c(4,4,3,1),oma=c(3,0,4,0))
for(region in regions) {
  h<-historical[historical$Region==region,];h<-h[order(h$depth),]
  o<-observed[observed$Region==region,];o<-o[order(o$depth),]
  plot(NA,xlim=range(c(h$lower95,h$upper95,o$chl_norm)),ylim=c(150,10),
    xlab="chl_norm",ylab="Depth (m)",main=region)
  polygon(c(h$lower95,rev(h$upper95)),c(h$depth,rev(h$depth)),col="#D3EDE5",border=NA)
  polygon(c(h$lower50,rev(h$upper50)),c(h$depth,rev(h$depth)),col="#7BC8C4",border=NA)
  lines(h$mean,h$depth,col="#005F73",lwd=2)
  lines(o$chl_norm,o$depth,col="#222222",lty=2,lwd=1.5)
  box()
}
plot.new();legend("center",legend=c("Expected population mean","Observed training median","50% posterior interval","95% posterior interval"),
  col=c("#005F73","#222222","#7BC8C4","#D3EDE5"),lty=c(1,2,NA,NA),lwd=c(2,1.5,NA,NA),pch=c(NA,NA,15,15),pt.cex=2,bty="n")
mtext("Historical expected population profiles over 10-150 m",outer=TRUE,side=3,line=1.7,cex=1.3)
mtext("Recorded regional median distance; 50/95% mean-curve uncertainty; group effects zero.",outer=TRUE,side=3,line=.3,cex=.8)
mtext("Dashed training medians mix observed predictor and group values; this is a descriptive reconstruction check.",outer=TRUE,side=1,line=1.1,cex=.8)
dev.off()
diagnostics<-read.csv(file.path(out,"historical_parameter_diagnostics.csv"))
beta<-diagnostics$mean[diagnostics$variable=="b_log10_distance"]
root<-"dcm_snapshot_20261016_090217/experiments"
specs<-list(
  list(label="Current original, 10% / 150 m",file=file.path(root,"baseline150_20261015T110518Z_1736107/runs/case_01/checks/example_population_profiles.csv")),
  list(label="Current additive MLD",file=file.path(root,"equation150_20261015T104335Z_1647064/runs/case_06/checks/example_population_profiles.csv")),
  list(label="Current MLD-varying",file=file.path(root,"equation150_20261015T104335Z_1647064/runs/case_12/checks/example_population_profiles.csv")))
current<-do.call(rbind,lapply(specs,function(s){d<-read.csv(s$file);d$model<-s$label;d}))
regions<-levels(factor(historical$Region));colours<-c("#222222","#005F73","#EE9B00","#9B2226")
comparison_rows<-list()
png(file.path(out,"historical_vs_current_profiles.png"),width=1560,height=1350,res=130)
par(mfrow=c(3,3),mar=c(4,4,3,1),oma=c(4,0,4,0))
for(region in regions) {
  c<-current[current$Region==region,];h<-historical[historical$Region==region,]
  reference_distance<-c$log10_distance[c$model==specs[[1]]$label][1]
  stopifnot(max(abs(c$log10_distance-reference_distance))<1e-12)
  # Historical distance enters linearly, so this posterior-mean adjustment is exact.
  h$adjusted_mean<-h$mean+beta*(reference_distance-h$log10_distance)
  comparison_rows[[length(comparison_rows)+1L]]<-data.frame(Region=region,depth=h$depth,
    model="Historical posterior mean",summary="posterior_mean",log10_distance=reference_distance,
    MLD_temp_MIX=NA_real_,estimate=h$adjusted_mean)
  comparison_rows[[length(comparison_rows)+1L]]<-data.frame(Region=region,depth=c$depth,
    model=c$model,summary="posterior_median",log10_distance=c$log10_distance,
    MLD_temp_MIX=c$MLD_temp_MIX,estimate=c$median)
  xrange<-range(c(h$adjusted_mean,c$median),finite=TRUE)
  plot(h$adjusted_mean,h$depth,type="l",ylim=c(150,10),xlim=xrange,xlab="chl_norm",ylab="Depth (m)",main=region,lwd=2)
  for(i in seq_along(specs)) {
    z<-c[c$model==specs[[i]]$label,];lines(z$median,z$depth,col=colours[i+1],lwd=1.8,lty=i+1)
  }
}
plot.new();legend("center",legend=c("Historical posterior mean",vapply(specs,`[[`,character(1),"label")),col=colours,lty=1:4,lwd=2,bty="n",cex=.9)
mtext("Historical and current population profiles over 10-150 m",outer=TRUE,side=3,line=1.7,cex=1.3)
mtext("Matched regional distance; current lines are saved posterior medians. MLD uses recorded example cast values.",outer=TRUE,side=3,line=.3,cex=.8)
mtext("Descriptive curve comparison: training datasets and depth domains differ; this is not held-out accuracy.",outer=TRUE,side=1,line=1.3,cex=.8)
dev.off()
write.csv(do.call(rbind,comparison_rows),file.path(out,"historical_vs_current_reference_profiles.csv"),row.names=FALSE)

# Reproduce probability-versus-depth plots using correctly labelled exclusive bins.
cuts<-c(.5,1,2,3,4,5);bins<-c("<=1.5x","1.5-2x","2-3x","3-4x","4-5x","5-6x",">6x")
palette<-c("#dddddd","#20255c","#4575b4","#74add1","#fdae61","#d73027","#a50026")
png(file.path(out,"historical_expected_enhancement_probabilities.png"),width=1560,height=1350,res=130)
par(mfrow=c(3,3),mar=c(4,4,3,1),oma=c(3,0,4,0))
for(region in regions) {
  h<-historical[historical$Region==region,]
  p<-as.matrix(h[paste0("prob_mu_gt_",cuts)])
  p<-cbind(1-p[,1],p[,-ncol(p),drop=FALSE]-p[,-1,drop=FALSE],p[,ncol(p)])
  stopifnot(max(abs(rowSums(p)-1))<1e-10,min(p)>=-1e-10)
  cs<-t(apply(p,1,cumsum));lo<-cbind(0,cs[,-ncol(cs)])
  plot(NA,xlim=c(0,1),ylim=c(150,10),xlab="Posterior probability",ylab="Depth (m)",main=region)
  for(i in seq_along(bins))polygon(c(lo[,i],rev(cs[,i])),c(h$depth,rev(h$depth)),col=palette[i],border=NA)
  box()
}
plot.new();legend("center",legend=bins,fill=palette,bty="n",title="Expected chl / surface chl",cex=.95)
mtext("Historical fit: probabilities for the expected population enhancement",outer=TRUE,side=3,line=1.7,cex=1.2)
mtext("Exclusive bins sum to 100%; P(chl_norm > 0.5) means P(enhancement > 1.5x), not P(1.5x to 2x).",outer=TRUE,side=3,line=.3,cex=.8)
dev.off()
cat("PASS: depth-ordered profile curves; matched-distance mean adjustment; exclusive probability bins sum to one; three offline plots.\n")
