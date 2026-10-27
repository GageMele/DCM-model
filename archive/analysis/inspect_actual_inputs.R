# Read-only inspection of the exact local inputs uploaded to Ibex.
root <- 'K:/Gage Mele/chla/DCM distribution model Ibex'
ctd <- read.csv(file.path(root,'CTD_NOAA_WHOI_KAUST.csv'),stringsAsFactors=FALSE,check.names=FALSE)
mld <- read.csv(file.path(root,'mld_full.csv'),stringsAsFactors=FALSE,check.names=FALSE)
cat('CTD rows:',nrow(ctd),'\nCTD columns:\n');print(names(ctd))
cat('MLD rows:',nrow(mld),'\nMLD columns:\n');print(names(mld))
bad <- function(x) is.na(x)|!nzchar(trimws(as.character(x)))
for(dn in c('ctd','mld')) {
  d <- get(dn)
  cols <- intersect(c('LocationCode','DeepCTDID','cast_id','Region.x','Region','UTCDate.x','date','Ocean'),names(d))
  cat('\nMissing identifier counts:',dn,'\n')
  print(vapply(d[cols],function(x)sum(bad(x)),integer(1)))
  invalid <- Reduce(`|`,lapply(d[cols],bad))
  cat('Examples of incomplete identifiers:\n');print(head(unique(d[invalid,cols,drop=FALSE]),15))
}
cat('\nCTD regions:\n');print(table(ctd$Region.x,useNA='ifany'))
cat('\nMLD regions:\n');print(table(mld$Region,useNA='ifany'))
ctd$Region <- ctd$Region.x
ctd$cast_id <- ctd$DeepCTDID
ctd <- ctd[ctd$Region!='KAUST',]
# Align the spelling to the CTD label used by the current model.
mld$Region[mld$Region=='Hawaiian'] <- 'Hawiian'
ctd$Ocean <- ifelse(ctd$Region %in% c('North RS','Central RS','South RS'),'RedSea','Pacific')
ctd_ids <- unique(ctd[c('LocationCode','Region','cast_id','Ocean','UTCDate.x')])
cat('\nMissing LocationCode cast records by region:\n');print(table(ctd_ids$Region[bad(ctd_ids$LocationCode)]))
cat('\nMLD RedSea LocationCode values:\n');print(table(mld$LocationCode[mld$Ocean=='RedSea']))
key3 <- function(x) do.call(paste,c(x[c('Region','cast_id','Ocean')],sep='|'))
mk <- key3(mld);ck <- key3(ctd_ids)
cat('\nDuplicated 3-column MLD identifiers:',sum(duplicated(mk)),'\n')
cat('CTD distinct cast records:',nrow(ctd_ids),'\n')
cat('Unmatched 3-column CTD identifiers:',sum(!ck %in% mk),'\n')
cat('Unmatched examples:\n');print(head(ctd_ids[!ck %in% mk,],15))
cat('Matching MLD metadata for missing CTD LocationCodes:\n');print(head(mld[match(ck[bad(ctd_ids$LocationCode)],mk),],15))
date_count <- table(ck)
cat('Cast identifiers spanning several date/location records:',sum(date_count>1),'\n')
cat('MLD non-finite values:',sum(!is.finite(mld$MLD_temp_MIX)),'\n')
