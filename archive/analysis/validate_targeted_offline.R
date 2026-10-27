# Offline, base-R verification of the already-written data splits and package.
# This script deliberately does not load brms, rstan, Rcpp, or any HTTP package.
package<-"DCM_Ibex_targeted_tests"
for(f in list.files(package,pattern="\\.R$",full.names=TRUE))parse(f)
source(file.path(package,"campaign.R"))
x<-campaign_cases();stopifnot(nrow(x)==63,min(as.numeric(sub("p","",x$subset)))>=10)
root<-"validation_targeted150/actual"
for(split in c("primary","date_blocked")) {
  tr<-readRDS(file.path(root,"data",paste0(split,"_p10_d150.rds")))
  va<-readRDS(file.path(root,"data",paste0(split,"_validation_d150.rds")))
  stopifnot(length(unique(tr$cast_key))==108,length(unique(va$cast_key))>=108,
    !any(tr$cast_key %in% va$cast_key),min(tr$depth)==10,max(tr$depth)==150)
  if(split=="date_blocked") {
    key<-function(z)paste(z$Ocean,z$Region,z$date,sep="\r")
    stopifnot(!any(key(tr) %in% key(va)))
  }
  cat(split,"train",length(unique(tr$cast_key)),"validation",length(unique(va$cast_key)),"disjoint casts; valid\n")
}
original<-readRDS(file.path(root,"data","p10.rds"))$model_data
prepared<-readRDS(file.path(root,"data","primary_p10_d150.rds"))
canon<-function(z){z<-z[order(as.character(z$cast_key),z$depth),c("cast_key","depth","chl_norm")];rownames(z)<-NULL;z$cast_key<-as.character(z$cast_key);z}
stopifnot(identical(canon(original),canon(prepared)))
cat("PASS: all R syntax, 63-case matrix, minimum data fractions, cast/date independence and unchanged training responses.\n")
