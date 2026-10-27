# Run on Ibex before preparation. This parses all R files and tests data/metric
# logic without fitting Stan or reading chlorophyll outcomes for selection.
a<-commandArgs(trailingOnly=TRUE);package<-normalizePath(a[1],mustWork=TRUE)
for(p in list.files(package,pattern="\\.R$",full.names=TRUE))parse(p)
source(file.path(package,"models.R"));source(file.path(package,"data.R"))
source(file.path(package,"onset_metrics.R"))
Sys.setenv(EXPERIMENT_DIR=package);source(file.path(package,"scoring.R"))
cases<-read.delim(file.path(package,"cases.tsv"),stringsAsFactors=FALSE)
stopifnot(nrow(cases)==282L,sum(cases$mode=="posterior")==270L,
  sum(cases$mode=="prior")==12L,sum(cases$fraction==1 & cases$mode=="posterior")==11L,
  identical(sort(unique(cases$k_depth)),c(7L,9L,11L,13L,15L)),
  all(cases$cpus==cases$chains*cases$threads),all(cases$fraction>=.1),
  all(cases$iter==cases$warmup+cases$samples),!anyDuplicated(cases$label))

# A rigid translation keeps the maximum derivative, whereas width compression
# changes it. These deterministic checks encode the actual scientific question.
z<-seq(10,150,by=1);f<-function(z)4*exp(-((z-80)/18)^2)
slope<-function(z,y)max(diff(y)/diff(z))
base<-slope(z,f(z));shift<-slope(z,f(z-10));sharp<-slope(z,f((z-10)*1.5))
stopifnot(abs(base-shift)<1e-12,sharp>base)
unsupported<-onset_from_curve(onset_curve(10,0),.5)
stopifnot(unsupported$status=="insufficient_support")

# Exact sampling targets, nested mandatory casts, outcome-free repeatability.
meta<-data.frame(cast_key=sprintf("c%04d",1:1000),Region=rep(c("Central RS","South RS","Other"),c(200,200,600)),
  Ocean="test",date="2026-01-01",enso="neutral",log10_distance=seq_len(1000)/1000,MLD_temp_MIX=20+seq_len(1000)/100)
levels<-c("Central RS","South RS","Other")
k20<-sample_keys(meta,meta,200,levels,20261022)
k40<-sample_keys(meta,meta,400,levels,20261022,mandatory=k20)
k20b<-sample_keys(meta,meta,200,levels,20261023)
stopifnot(length(k20)==200L,length(k40)==400L,all(k20%in%k40),
  identical(k20,sample_keys(meta,meta,200,levels,20261022)),!setequal(k20,k20b))
cat("PASS: 282-case design; R syntax; full-data labelling; shift/sharpness distinction; singleton support; whole-cast targets and repeat selections.\n")
