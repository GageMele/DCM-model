library(brms)
d <- readRDS("validation_extensive150/actual/data/p10_model.rds")
d$depth_z <- (d$depth-80)/70
d$cast_key <- factor(d$cast_key)
f <- bf(chl_norm ~ Region + distance_z + mld_z + s(depth,by=Region,bs="cr",k=7) +
          (1|Ocean:Region:date)+(1|Ocean:Region:enso),sigma~0+Region,nu=4,family=student())
p <- c(set_prior("normal(0,2)",class="b"),set_prior("student_t(3,0,1.5)",class="sd"),
 set_prior("student_t(3,0,1.5)",class="sds"),set_prior("normal(0,1)",class="b",dpar="sigma"))
sd <- make_standata(f,data=d,prior=p)
cat("FIXED NU WORKS:",sd$nu,"\n")
f <- bf(chl_norm ~ Region + distance_z + mld_z + s(depth,by=Region,bs="cr",k=7) +
          (1|Ocean:Region:date)+(1|Ocean:Region:enso)+ar(time=depth,gr=cast_key,p=1,cov=TRUE))
tm <- system.time(sd <- make_standata(f,data=d))
cat("AR data time",tm,"SIZE",format(object.size(sd),units="MB"),"\n")
cat("AR NAMES:",names(sd),"\n")
for(n in grep("[Tt][Ii][Mm]|[Tt][Gg]|lag|AR|ar",names(sd),value=TRUE)) {
 cat(n,"dim",dim(sd[[n]]),"length",length(sd[[n]]),"\n")
 print(utils::head(sd[[n]]))
}
code <- make_stancode(f,data=d)
writeLines(code,"probe_targeted_ar.stan")
