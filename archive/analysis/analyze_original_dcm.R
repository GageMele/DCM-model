# Diagnose and plot the saved historical fit; never call brm()/update().
library(brms)
library(posterior)
library(ggplot2)
options(error=function(){traceback(25);quit(status=1L)})
out <- "original_model_review"; dir.create(out,showWarnings=FALSE)
fit <- readRDS("original_model_input/my_model_50_FINAL.rds")
cat("rstan available:",requireNamespace("rstan",quietly=TRUE),"\n")
cat("DATE CLASS:",class(fit$data$date),"\n"); print(head(fit$data$date)); print(attributes(fit$data$date))
draws <- tryCatch(posterior::as_draws_array(fit), error=function(e) {
  message("Extracting hydrated saved samples directly: ",conditionMessage(e))
  sim <- attr(fit$fit,"sim")
  n <- nrow(sim$samples[[1]])
  a <- array(NA_real_,c(n,length(sim$samples),ncol(sim$samples[[1]])),
    dimnames=list(NULL,NULL,names(sim$samples[[1]])))
  for(i in seq_along(sim$samples)) a[,i,] <- as.matrix(sim$samples[[i]])
  posterior::as_draws_array(a)
})
diag <- as.data.frame(posterior::summarise_draws(draws,"mean","sd","rhat","ess_bulk","ess_tail","mcse_mean"))
write.csv(diag,file.path(out,"historical_parameter_diagnostics.csv"),row.names=FALSE)
sim <- attr(fit$fit,"sim")
sampler <- lapply(seq_along(sim$samples),function(i) {
  x <- as.data.frame(attr(sim$samples[[i]],"sampler_params"))
  x$chain <- i; x
})
print(lapply(sampler,names))
sampler <- do.call(rbind,sampler)
write.csv(sampler,file.path(out,"historical_sampler_draws.csv"),row.names=FALSE)
cap <- as.numeric(attr(fit$fit,"stan_args")[[1]]$max_depth)
chains <- do.call(rbind,lapply(split(sampler,sampler$chain),function(x) data.frame(chain=x$chain[1],
  draws=nrow(x),divergences=sum(x$divergent__),tree_hits=sum(x$treedepth__>=cap),
  max_depth=max(x$treedepth__),mean_leapfrogs=mean(x$n_leapfrog__),
  ebfmi=mean(diff(x$energy__)^2)/var(x$energy__))))
write.csv(chains,file.path(out,"historical_chain_diagnostics.csv"),row.names=FALSE)
print(chains)
varying <- diag[is.finite(diag$sd)&diag$sd>0,]
summary <- data.frame(observations=nrow(fit$data),min_depth=min(fit$data$depth),max_depth=max(fit$data$depth),
  chains=sim$chains,iter=sim$iter,warmup=sim$warmup,retained=posterior::ndraws(draws),
  max_rhat=max(varying$rhat,na.rm=TRUE),min_bulk_ess=min(varying$ess_bulk,na.rm=TRUE),
  min_tail_ess=min(varying$ess_tail,na.rm=TRUE),divergences=sum(chains$divergences),
  tree_hits=sum(chains$tree_hits),min_ebfmi=min(chains$ebfmi))
write.csv(summary,file.path(out,"historical_summary.csv"),row.names=FALSE); print(summary)
print(diag[diag$variable %in% c("sigma","b_log10_distance")|grepl("^sd_",diag$variable),])
set.seed(20261016); ids <- sort(sample.int(posterior::ndraws(draws),1000))
# The saved posterior is hydrated, but the installed rstan DLL cannot load here.
# Supply the existing draws through an analysis-only S3 adapter; no sampling or
# alteration of the saved object is involved. brms still builds its own designs.
fit$dcm_saved_draws <- draws
if(!requireNamespace("rstan",quietly=TRUE)) {
  # A local structural definition allows read-only @slot access without the
  # unavailable compiled rstan module. Values and parameter names are unchanged.
  slots<-setNames(rep("ANY",length(setdiff(names(attributes(fit$fit)),"class"))),
    setdiff(names(attributes(fit$fit)),"class"))
  methods::setClass("dcm_saved_stanfit",slots=slots,where=.GlobalEnv)
  attr(fit$fit,"class")<-structure("dcm_saved_stanfit",package=".GlobalEnv")
}
class(fit) <- c("dcm_hydrated",class(fit))
extract_saved <- function(x,variable=NULL,regex=FALSE,...) {
  out <- x$dcm_saved_draws
  if(!is.null(variable))out<-posterior::subset_draws(out,variable=variable,regex=regex)
  out
}
registerS3method("as_draws_array","dcm_hydrated",extract_saved,envir=asNamespace("posterior"))
registerS3method("as_draws_list","dcm_hydrated",function(x,...)posterior::as_draws_list(extract_saved(x,...)),envir=asNamespace("posterior"))
registerS3method("as_draws_matrix","dcm_hydrated",function(x,...)posterior::as_draws_matrix(extract_saved(x,...)),envir=asNamespace("posterior"))
registerS3method("ndraws","dcm_hydrated",function(x,...)posterior::ndraws(x$dcm_saved_draws),envir=asNamespace("posterior"))
regions <- levels(fit$data$Region)
representatives <- lapply(regions,function(region) {
  d <- fit$data[fit$data$Region==region,]
  # Use a real recorded predictor value nearest the regional median.
  d <- d[which.min(abs(d$log10_distance-median(d$log10_distance))),,drop=FALSE]
  x <- d[rep(1,141),,drop=FALSE]; x$depth <- 10:150; x
})
grid <- do.call(rbind,representatives); rownames(grid)<-NULL
mu <- brms::posterior_epred(fit,newdata=grid,re_formula=NA,draw_ids=ids)
q <- apply(mu,2,quantile,probs=c(.025,.25,.5,.75,.975))
pred <- data.frame(grid[c("Region","Ocean","depth","log10_distance")],
  mean=colMeans(mu),lower95=q[1,],lower50=q[2,],median=q[3,],upper50=q[4,],upper95=q[5,])
for(t in c(0,.5,1,2,3,4,5)) pred[[paste0("prob_mu_gt_",t)]] <- colMeans(mu>t)
write.csv(pred,file.path(out,"historical_population_profiles_10_150.csv"),row.names=FALSE)
observed <- aggregate(chl_norm~Region+depth,fit$data[fit$data$depth<=150,],median)
write.csv(observed,file.path(out,"historical_observed_median_profiles.csv"),row.names=FALSE)
p <- ggplot(pred,aes(mean,depth))+geom_ribbon(aes(xmin=lower95,xmax=upper95),orientation="y",fill="#94D2BD",alpha=.4)+
  geom_ribbon(aes(xmin=lower50,xmax=upper50),orientation="y",fill="#0A9396",alpha=.4)+geom_path(colour="#005F73",linewidth=.6)+
  geom_path(data=observed[order(observed$Region,observed$depth),],aes(chl_norm,depth),linetype="dashed",colour="black",linewidth=.5)+scale_y_reverse()+
  facet_wrap(~Region,scales="free_x",ncol=3)+theme_bw()+
  labs(x="chl_norm",y="Depth (m)",title="Historical fitted population profiles over the current 10-150 m range",
    subtitle="Solid: expected profile at recorded regional median distance; bands: 50/95% posterior intervals; dashed: training-data median")
ggsave(file.path(out,"historical_profiles.png"),p,width=13,height=11,dpi=130,bg="white")
# A comparable bounded reconstruction sample: whole profiles cannot be identified
# from the saved fit because its model data omit cast IDs and MLD.
set.seed(20261016)
index <- unlist(lapply(split(seq_len(nrow(fit$data)),fit$data$Region),function(ii) sample(ii,min(80,length(ii)))))
d <- fit$data[index,]; m <- brms::posterior_epred(fit,newdata=d,draw_ids=ids)
y <- brms::posterior_predict(fit,newdata=d,draw_ids=ids)
q <- apply(y,2,quantile,c(.05,.95))
ppc <- data.frame(d,predicted=colMeans(m),lower=q[1,],upper=q[2,])
write.csv(ppc,file.path(out,"historical_bounded_training_predictions.csv"),row.names=FALSE)
cat("BOUNDED TRAINING RMSE:",sqrt(mean((ppc$chl_norm-ppc$predicted)^2)),
    "COVERAGE90:",mean(ppc$chl_norm>=ppc$lower&ppc$chl_norm<=ppc$upper),"\n")
cat("Historical analysis finished. No model fitting.\n")
