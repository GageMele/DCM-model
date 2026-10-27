# Static scientific figures from saved summaries; no model fitting or prediction.
root <- 'DCM_review_20261018'
read <- function(name)read.csv(file.path(root,name),stringsAsFactors=FALSE,fileEncoding='UTF-8-BOM')
d <- read('all_registered_runs.csv')
d$case_id <- suppressWarnings(as.numeric(d$case_id))
d$fit_completed <- tolower(as.character(d$fit_completed))=='true'
d$prior_only <- tolower(as.character(d$prior_only))=='true'
region <- read('all_validation_conditional_by_region.csv')
t <- d[d$campaign=='targeted150' & d$fit_completed & !d$prior_only,,drop=FALSE]
colors <- c('#005F73','#0A9396','#CA6702','#AE2012','#6A4C93','#607744','#777777')
number <- function(x)suppressWarnings(as.numeric(x))
pass <- function(x)tolower(as.character(x))=='true'
save_png <- function(name,width,height,draw) {
  if(Sys.getenv('DCM_PROFILES_ONLY')=='1' && as.integer(substr(name,1,2))<6L)return(invisible(NULL))
  png(file.path(root,paste0(name,'.png')),width=width,height=height,res=145,type='cairo')
  on.exit(dev.off());draw()
}
dot <- function(x,labels,title,reference=NULL,log=FALSE,cex=.75,col=NULL) {
  x<-number(x);y<-rev(seq_along(x));ok<-is.finite(x)
  if(log)ok<-ok & x>0
  if(!any(ok)){plot.new();title(main=title);return(invisible(NULL))}
  plot(x[ok],y[ok],xlab='',ylab='',yaxt='n',ylim=c(.5,length(x)+.5),pch=19,
    col=if(is.null(col))'#005F73'else col[ok],main=title,log=if(log)'x'else'',cex=.8)
  abline(h=y,col='#eeeeee');points(x[ok],y[ok],pch=19,col=if(is.null(col))'#005F73'else col[ok],cex=.8)
  axis(2,at=y,labels=labels,las=2,cex.axis=cex,tick=FALSE)
  if(!is.null(reference))abline(v=reference,lty=2,col='#AE2012')
}

save_png('01_all_posterior_runs',2000,1450,function(){
  p<-d[d$fit_completed & !d$prior_only,,drop=FALSE]
  camps<-unique(p$campaign);cols<-setNames(colors[seq_along(camps)],camps)
  par(mfrow=c(2,3),mar=c(4,4.5,3,1),oma=c(2,0,4,0))
  sets<-list(c('fit_hours','validation_conditional_rmse','Fit hours','Validation RMSE'),
    c('validation_conditional_rmse','validation_conditional_interval_score90','Validation RMSE','90% interval score'),
    c('validation_conditional_coverage90','validation_conditional_interval_width90','90% predictive coverage','90% interval width'),
    c('max_rhat','min_bulk_ess','Maximum R-hat','Minimum bulk ESS'),
    c('k_depth','validation_conditional_shape_rmse','Regional depth k','Validation shape RMSE'),
    c('fit_hours','tree_hit_fraction_review','Fit hours','Post-warmup tree-hit fraction'))
  for(s in sets){
    x<-number(p[[s[1]]]);y<-number(p[[s[2]]]);ok<-is.finite(x)&is.finite(y)
    primary<-p$comparison_cohort %in% c('primary_md5_verified','primary_membership_only')
    if(s[2] %in% c('validation_conditional_rmse','validation_conditional_interval_score90','validation_conditional_interval_width90','validation_conditional_shape_rmse'))ok<-ok&primary
    plot(x[ok],y[ok],xlab=s[3],ylab=s[4],pch=ifelse(pass(p$harmonized_convergence_screen[ok]),19,4),col=cols[p$campaign[ok]],cex=.8)
    if(s[1]=='max_rhat'){abline(v=1.01,lty=2);abline(h=400,lty=2)}
    if(s[1]=='validation_conditional_coverage90')abline(v=.9,lty=2)
  }
  mtext('All completed posterior fits: every campaign retained',outer=TRUE,line=2,font=2,cex=1.4)
  mtext('Circles: convergence screen passed; crosses: flagged. Date-blocked scores are excluded from the primary-score panels.',outer=TRUE,line=.5,cex=.8)
  par(fig=c(0,1,0,.055),new=TRUE,mar=c(0,0,0,0));plot.new();legend('center',legend=camps,col=cols,pch=19,horiz=TRUE,bty='n',cex=.8)
})

save_png('02_targeted_all_cases',2400,2600,function(){
  x<-t;labels<-sprintf('T%02d %s',x$case_id,x$block)
  co<-ifelse(pass(x$harmonized_convergence_screen),'#005F73','#AE2012')
  par(mfrow=c(1,6),mar=c(4,9,5,1),oma=c(2,0,4,0))
  for(s in list(c('validation_conditional_rmse','Validation RMSE'),c('validation_conditional_shape_rmse','Shape RMSE'),
    c('validation_conditional_coverage90','90% coverage'),c('validation_conditional_interval_score90','90% interval score'),
    c('max_rhat','Maximum R-hat'),c('fit_hours','Fit hours')))
    dot(x[[s[1]]],labels,s[2],reference=if(s[1]=='max_rhat')1.01 else if(s[1]=='validation_conditional_coverage90').9 else NULL,col=co,cex=.61)
  mtext('Every completed targeted posterior case in case order',outer=TRUE,line=2,font=2,cex=1.3)
  mtext('Red: diagnostic flag. T51/T52 use a separate 117-cast date-blocked cohort; their raw errors are not a paired comparison with the primary cohort.',outer=TRUE,line=.5,cex=.8)
})

save_png('03_structure_and_equations',2400,1800,function(){
  x<-t[t$case_id<=22,,drop=FALSE]
  labels<-paste0('T',sprintf('%02d',x$case_id),' ',x$equation,' / ',x$variance,ifelse(x$cast=='none','',paste0(' / cast ',x$cast)),ifelse(x$family=='student4',' / t4',''))
  co<-ifelse(pass(x$harmonized_convergence_screen),'#005F73','#AE2012')
  par(mfrow=c(2,3),mar=c(4,18,3,1),oma=c(2,0,3,0))
  for(s in list(c('validation_conditional_rmse','Validation RMSE'),c('validation_conditional_shape_rmse','Shape RMSE'),
    c('validation_conditional_peak_depth_mae','Peak-depth MAE (m)'),c('validation_conditional_coverage90','90% coverage'),
    c('validation_conditional_interval_score90','90% interval score'),c('fit_hours','Fit hours')))
    dot(x[[s[1]]],labels,s[2],reference=if(s[1]=='validation_conditional_coverage90').9 else NULL,col=co,cex=.62)
  mtext('Residual structure, mean equation and depth basis',outer=TRUE,line=1,font=2,cex=1.2)
})

save_png('04_sampler_settings',2300,1800,function(){
  x<-t[t$case_id %in% c(2,8,29:46),,drop=FALSE]
  labels<-sprintf('T%02d %s: d%.2f td%s w%s s%s th%s',x$case_id,x$equation,x$adapt_delta,x$treedepth,x$warmup,x$samples,x$threads)
  co<-ifelse(pass(x$harmonized_convergence_screen),'#005F73','#AE2012')
  par(mfrow=c(2,3),mar=c(4,18,3,1),oma=c(2,0,3,0))
  for(s in list(c('validation_conditional_rmse','Validation RMSE'),c('min_bulk_ess','Minimum bulk ESS'),
    c('max_rhat','Maximum R-hat'),c('tree_hit_fraction_review','Tree-hit fraction'),
    c('fit_hours','Fit hours'),c('requested_core_hours','Requested core-hours')))
    dot(x[[s[1]]],labels,s[2],reference=if(s[1]=='max_rhat')1.01 else if(s[1]=='min_bulk_ess')400 else NULL,col=co,cex=.65)
  mtext('Controlled sampler comparisons; timing remains hardware dependent',outer=TRUE,line=1,font=2,cex=1.2)
})

save_png('05_regional_calibration',2000,1500,function(){
  ids<-sprintf('targeted150:%03d',c(1,2,3,7,8,13,17,18,47,49))
  regions<-unique(region$Region)
  par(mfrow=c(1,2),mar=c(8,13,4,2),oma=c(1,0,3,0))
  for(metric in c('rmse','coverage90')){
    mat<-matrix(NA_real_,length(ids),length(regions))
    for(i in seq_along(ids))for(j in seq_along(regions)){
      z<-region[region$run_id==ids[i] & region$Region==regions[j],,drop=FALSE]
      if(nrow(z))mat[i,j]<-number(z[[metric]][1])
    }
    pal<-if(metric=='coverage90')colorRampPalette(c('#AE2012','#E9D8A6','#94D2BD','#005F73'))(100)else colorRampPalette(c('#f7fcfd','#94D2BD','#005F73'))(100)
    limits<-if(metric=='coverage90')c(0,1)else range(mat,na.rm=TRUE)
    image(seq_along(regions),seq_along(ids),t(mat[nrow(mat):1,,drop=FALSE]),axes=FALSE,xlab='',ylab='',col=pal,zlim=limits)
    axis(1,seq_along(regions),regions,las=2,cex.axis=.8)
    axis(2,seq_along(ids),rev(sub('targeted150:','T',ids)),las=2,cex.axis=.9,tick=FALSE)
    for(i in seq_along(ids))for(j in seq_along(regions))if(is.finite(mat[i,j]))text(j,length(ids)-i+1,sprintf(if(metric=='coverage90')'%.0f%%'else'%.2f',mat[i,j]*if(metric=='coverage90')100 else 1),cex=.7,col=if(mat[i,j]>(limits[1]+.65*diff(limits)))'white'else'#111111')
    title(main=if(metric=='rmse')'Regional cast-balanced RMSE'else'90% predictive interval coverage')
  }
  mtext('Regional calibration can differ despite similar overall coverage',outer=TRUE,line=1,font=2,cex=1.2)
})

# Example casts were chosen by predictor centrality, never by observed response or error.
examples<-read('all_validation_conditional_example_profiles.csv')
profiles<-read('all_population_reference_profiles.csv')
profile_figure<-function(tab,name,ids,observed=FALSE){
  z<-tab[tab$run_id %in% ids,,drop=FALSE];regions<-unique(z$Region)
  save_png(name,2100,1900,function(){
    nr<-ceiling(length(regions)/3);par(mfrow=c(nr,3),mar=c(4,4.5,3,1),oma=c(7,0,4,0))
    for(reg in regions){
      zz<-z[z$Region==reg,,drop=FALSE]
      field<-if(observed)'predicted'else'mean'
      xx<-number(zz[[field]]);if(observed)xx<-c(xx,number(zz$chl_norm))
      plot(NA,xlim=range(xx,finite=TRUE),ylim=rev(range(zz$depth)),type='n',xlab='chl_norm',ylab='Depth (m)',main=reg)
      if(observed){a<-zz[zz$run_id==ids[1],,drop=FALSE];a<-a[order(a$depth),];lines(a$chl_norm,a$depth,col='#222222',lwd=1.7)}
      for(i in seq_along(ids)){
        a<-zz[zz$run_id==ids[i],,drop=FALSE];a<-a[order(a$depth),]
        lines(a[[field]],a$depth,col=colors[i],lwd=1.7,lty=i)
      }
    }
    mtext(if(observed)'Same held-out example cast per region: observed and predicted mean profiles'else'Population mean profiles at the same regional predictor settings',outer=TRUE,line=2,font=2,cex=1.2)
    plot.new()
    legend('center',legend=c(if(observed)'Observed',sub('targeted150:','T',ids)),col=c(if(observed)'#222222',colors[seq_along(ids)]),lty=c(if(observed)1,seq_along(ids)),lwd=2,bty='n',cex=1.1)
  })
}
ids<-sprintf('targeted150:%03d',c(1,2,8,13,18))
profile_figure(examples,'06_heldout_profile_examples',ids,TRUE)
profile_figure(profiles,'07_population_profile_comparison',ids,FALSE)

save_png('08_historical_depth_bridge',2000,1700,function(){
  ids<-c('baseline150:001','targeted150:061');z<-profiles[profiles$run_id %in% ids,,drop=FALSE]
  regs<-unique(z$Region);par(mfrow=c(ceiling(length(regs)/3),3),mar=c(4,4.5,3,1),oma=c(6,0,4,0))
  for(reg in regs){
    zz<-z[z$Region==reg,,drop=FALSE];xr<-range(c(zz$lower95,zz$upper95),finite=TRUE)
    plot(NA,xlim=xr,ylim=c(150,10),type='n',xlab='Expected chl_norm',ylab='Depth (m)',main=reg)
    for(i in seq_along(ids)){
      a<-zz[zz$run_id==ids[i],,drop=FALSE];a<-a[order(a$depth),]
      polygon(c(a$lower95,rev(a$upper95)),c(a$depth,rev(a$depth)),col=adjustcolor(colors[i],alpha.f=.15),border=NA)
      lines(a$mean,a$depth,col=colors[i],lty=i,lwd=1.8)
    }
  }
  mtext('Original equation: fit through 150 m versus 250 m',outer=TRUE,line=2,font=2,cex=1.2)
  mtext('Same 108 training casts; original priors, k=15, 4 chains, iter=5000, warmup=1000, delta=.97, tree cap=14. Shading: 95% mean-curve intervals.',outer=TRUE,line=.5,cex=.7)
  plot.new();legend('center',c('Fit through 150 m','Fit through 250 m'),col=colors[1:2],lty=1:2,lwd=2,bty='n',cex=1.1)
})
cat('Created eight scientific figures.\n')
