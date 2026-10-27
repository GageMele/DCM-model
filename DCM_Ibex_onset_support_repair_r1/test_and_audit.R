# Base-R tests and real-cohort audits run on Ibex before any job is submitted.
args<-commandArgs(trailingOnly=TRUE);package<-args[1];experiment<-args[2]
for(f in list.files(package,pattern="\\.R$",full.names=TRUE))parse(f)
source(file.path(package,"onset_metrics.R"))
# Compare well-supported real curves with the previous scorer as a regression check.
previous<-new.env(parent=globalenv())
old_file<-file.path(Sys.getenv("ONSET_REPAIR_DIR"),"scoring_backups","onset_metrics.R")
if(!file.exists(old_file))old_file<-file.path(experiment,"onset_metrics.R")
sys.source(old_file,envir=previous)
for(n in 0:8) {
  z<-if(n)seq.int(20,length.out=n)else numeric()
  c<-onset_curve(z,rep(1,n));o<-onset_from_curve(c,.5)
  stopifnot(length(c$depth)==n)
  if(n<7L)stopifnot(o$status=="insufficient_support",is.na(o$depth))
  else stopifnot(o$status=="left_censored")
}
stopifnot(onset_from_curve(onset_curve(10:16,rep(0,7)))$status=="insufficient_support",
  !length(onset_curve(12.4,1)$depth),!length(onset_curve(160:170,rep(1,11))$depth))
z<-10:150;y<-pmax(0,(z-40)/20)
stopifnot(onset_from_curve(onset_curve(z,y),.5)$depth==50,
  onset_from_curve(onset_curve(z,pmax(0,(z-60)/20)),.5)$depth==70,
  onset_from_curve(onset_curve(z,rep(0,141)),.5)$status=="not_reached_in_window",
  onset_from_curve(onset_curve(z,rep(1,141)),.5)$status=="left_censored")
spike<-rep(0,141);spike[50:51]<-2
stopifnot(onset_from_curve(onset_curve(z,spike),.5)$status=="not_reached_in_window",
  onset_from_curve(onset_curve(c(10:30,50:150),c(rep(0,21),rep(1,101))),.5)$status=="gap_censored")
one<-data.frame(cast_key="short",Region="R",Ocean="O",date="D",depth=20,
  MLD_temp_MIX=30,chl_norm=.7,predicted=.6)
a<-onset_profile_metrics(one,matrix(c(.5,.6),ncol=1))
stopifnot(nrow(a)==3L,all(a$observed_status=="insufficient_support"),
  all(a$predicted_status=="insufficient_support"),all(is.na(a$reached_agreement)),
  all(a$posterior_insufficient_support_draws==2L),all(is.na(a$probability_reached_in_window)))
many<-one[rep(1L,141),];many$cast_key<-"long";many$depth<-z
many$chl_norm<-y;many$predicted<-y;b<-onset_profile_metrics(many,rbind(y,y))
stopifnot(all(b$onset_abs_error==0))
s<-onset_scores(rbind(a[a$threshold==.5,],b[b$threshold==.5,]))
stopifnot(s$casts==2L,s$observed_insufficient_support==1L,s$both_detected==1L,
  s$onset_mae_given_both_detected==0,s$reached_comparable==1L)
files<-list.files(file.path(experiment,"data"),pattern="(validation|reserve)_d150\\.rds$",full.names=TRUE)
stopifnot(length(files)==4L)
audit<-list();k<-0L
for(file in files) {
  d<-readRDS(file);parts<-split(seq_len(nrow(d)),as.character(d$cast_key))
  for(ii in parts) {
    x<-d[ii,,drop=FALSE];c<-onset_curve(x$depth,x$chl_norm);o<-onset_from_curve(c,.5)
    old_curve<-try(previous$onset_curve(x$depth,x$chl_norm),silent=TRUE)
    old_failed<-inherits(old_curve,"try-error")
    if(o$status!="insufficient_support") {
      stopifnot(!old_failed,isTRUE(all.equal(old_curve$depth,c$depth)),
                isTRUE(all.equal(old_curve$value,c$value)))
      old_onset<-previous$onset_from_curve(old_curve,.5)
      stopifnot(identical(old_onset$status,o$status),identical(old_onset$depth,o$depth))
    }
    k<-k+1L;audit[[k]]<-data.frame(dataset=sub("\\.rds$","",basename(file)),
      cast_key=as.character(x$cast_key[1]),Region=as.character(x$Region[1]),
      observations=nrow(x),minimum_depth=min(x$depth),maximum_depth=max(x$depth),
      grid_points=length(c$depth),maximum_supported_run=longest_supported_run(is.finite(c$value)),
      previous_scorer_failed=old_failed,onset_status=o$status)
  }
  stopifnot(sum(vapply(parts,length,integer(1)))==nrow(d))
}
out<-do.call(rbind,audit)
write.csv(out,file.path(Sys.getenv("ONSET_REPAIR_DIR"),"actual_cohort_support_audit.csv"),row.names=FALSE)
print(as.data.frame(table(out$dataset,out$onset_status)),row.names=FALSE)
print(out[out$onset_status=="insufficient_support",],row.names=FALSE)
cat("PASS: short/empty/singleton profiles; shifted/absent/censored/gapped onsets; all actual held-out casts retained and audited.\n")
source(file.path(package,"ledger.R"))
temp<-tempfile("onset_repair_ledger_");dir.create(file.path(temp,"scheduler"),recursive=TRUE)
n<-0L;submit<-function(args){n<<-n+1L;"123456"};lookup<-function(name)character()
stopifnot(ensure_job(temp,"check_01","fixture",character(),submit,lookup)=="123456",
  ensure_job(temp,"check_01","fixture",character(),submit,lookup)=="123456",n==1L)
unlink(file.path(temp,"scheduler","check_01.job"))
stopifnot(ensure_job(temp,"check_01","fixture",character(),submit,function(name)"123456")=="123456",n==1L)
cat("PASS: repeated/interrupted evaluation submissions preserve accepted IDs.\n")
