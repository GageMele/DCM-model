# Static scientific figures from all saved posterior results.
library(ggplot2)
root <- "DCM_all_run_comparison_20261016"
out <- "dcm_comparison_figures"
dir.create(out, showWarnings = FALSE)
r <- read.csv(file.path(root, "02_all_completed_posterior_runs.csv"), check.names = FALSE)
e <- r[r$campaign == "extensive150", ]
e$screen <- factor(e$harmonized_convergence_screen, levels=c("True","False"), labels=c("Passed convergence screen","Diagnostic flags"))
theme_set(theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(), legend.position="bottom"))
save <- function(plot, name, width=12, height=8) {
  ggsave(file.path(out,name),plot,width=width,height=height,dpi=130,bg="white")
}
p <- ggplot(e,aes(fit_hours,validation_conditional_rmse,colour=block,shape=screen))+
  geom_point(size=2.8,alpha=.85)+scale_x_log10()+
  geom_text(data=e[e$case_id %in% c(6,24,30,36,40,74,85,95),],aes(label=case_id),
    nudge_y=.007,show.legend=FALSE,size=3)+
  labs(title="All 84 completed held-out evaluations",subtitle="Every completed posterior case is shown, including diagnostic flags; labels identify selected cases",
       x="Fit time, hours (log scale)",y="Held-out cast-balanced RMSE",colour="Experiment block",shape=NULL)
save(p,"all_validation_and_cost.png")

k <- e[e$block=="equation_k",]
kk <- rbind(data.frame(k[c("case_id","equation","k_depth","screen")],mode="Conditional",rmse=k$validation_conditional_rmse),
            data.frame(k[c("case_id","equation","k_depth","screen")],mode="Population (group effects zero)",rmse=k$validation_population_rmse))
p <- ggplot(kk,aes(k_depth,rmse,colour=mode,group=mode))+geom_line()+geom_point(aes(shape=screen),size=2)+
  facet_wrap(~equation,ncol=3)+scale_x_continuous(breaks=c(5,7,10,15))+
  labs(title="Equation and k comparisons: all completed cases 1-40",subtitle="Region-MLD tensor cases 25-28 remain unfinished in this snapshot",
       x="Depth k",y="Held-out cast-balanced RMSE",colour=NULL,shape=NULL)
save(p,"all_equation_k_comparisons.png",13,10)

g <- read.csv(file.path(root,"all_validation_conditional_by_region.csv"))
g$Region <- factor(g$Region,levels=c("Central RS","North RS","South RS","Equatorial","Hawiian","Mariana","NWHI","Samoa"))
g$coverage_percent <- 100*g$coverage90
highlight <- g[g$case_id %in% c(6,85,95),]
highlight$case <- factor(highlight$case_id,levels=c(6,85,95),labels=c("6: additive 10%","85: additive 50%","95: Student-t both scaled"))
p <- ggplot(g,aes(Region,coverage_percent))+geom_hline(yintercept=90,linetype="dashed",colour="#333333")+
  geom_point(alpha=.18,size=1.6,colour="grey40",position=position_jitter(width=.1,height=0,seed=42))+
  geom_point(data=highlight,aes(colour=case,shape=case),size=3.4)+
  scale_y_continuous(limits=c(0,101),breaks=seq(0,100,20))+
  labs(title="Overall coverage hides substantial regional differences",subtitle="Grey points: every completed posterior evaluation; dashed line: nominal 90% prediction-interval coverage",
       x=NULL,y="Held-out predictive coverage (%)",colour=NULL,shape=NULL)
save(p,"all_regional_coverage.png",12,6)

r$screen <- factor(r$harmonized_convergence_screen,levels=c("True","False"),labels=c("Passed convergence screen","Diagnostic flags"))
p <- ggplot(r,aes(min_bulk_ess,max_rhat,colour=campaign,shape=screen))+geom_point(size=2.5,alpha=.85)+
  geom_vline(xintercept=400,linetype="dashed")+geom_hline(yintercept=1.01,linetype="dashed")+
  scale_x_log10()+labs(title="All 116 completed posterior fits: convergence and precision",subtitle="Dashed thresholds reproduce the common review screen; tree-depth warnings are reported separately",
  x="Minimum bulk ESS (log scale)",y="Maximum R-hat",colour="Campaign",shape=NULL)
save(p,"all_posterior_diagnostics.png",12,7)

# Compare existing population curves for all 18 equation tests.
camp <- "dcm_snapshot_20261016_090217/experiments/equation150_20261015T104335Z_1647064/runs"
profiles <- lapply(1:18,function(i) {
  x <- read.csv(file.path(camp,sprintf("case_%02d",i),"checks","example_population_profiles.csv"))
  x$case_id <- i; x$equation <- c("A: additive + Region-depth","B: MLD-varying + Region-depth","C: MLD-only shared shape")[(i-1)%/%6+1]
  x
})
profiles <- do.call(rbind,profiles)
write.csv(profiles,file.path(root,"all_legacy_equation_population_profiles.csv"),row.names=FALSE)
p <- ggplot(profiles,aes(depth,median,colour=equation,group=case_id))+geom_line(alpha=.6,linewidth=.6)+
  facet_wrap(~Region,scales="free_y",ncol=3)+labs(title="All 18 legacy equation tests: existing population profiles",
  subtitle="Each line is one run at its recorded example predictor settings; these curves are not held-out accuracy tests",
  x="Depth (m)",y="Predicted chl_norm",colour=NULL)
save(p,"all_legacy_equation_profiles.png",13,10)
cat("Saved five figures and all 18 legacy profile tables.\n")
