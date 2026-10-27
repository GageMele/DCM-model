# Resumable Slurm submission. Accepted fits are never automatically replaced.
ids<-function(x)unique(grep("^[0-9]+$",trimws(as.character(x)),value=TRUE))
atomic_lines<-function(text,path) {
  temporary<-paste0(path,".tmp.",Sys.getpid());writeLines(text,temporary)
  if(!file.rename(temporary,path)) {
    if(.Platform$OS.type!="windows")stop("Cannot save scheduler record: ",path)
    stopifnot(file.copy(temporary,path,overwrite=TRUE));unlink(temporary)
  }
}
atomic_table<-function(x,path) {
  temporary<-paste0(path,".tmp.",Sys.getpid());write.table(x,temporary,sep="\t",row.names=FALSE,quote=FALSE)
  if(!file.rename(temporary,path)) {
    if(.Platform$OS.type!="windows")stop("Cannot save manifest")
    stopifnot(file.copy(temporary,path,overwrite=TRUE));unlink(temporary)
  }
}
hash_text<-function(text) {
  file<-tempfile();on.exit(unlink(file));writeLines(text,file);substr(unname(tools::md5sum(file)),1,10)
}
run_read<-function(command,args) {
  output<-suppressWarnings(system2(command,vapply(args,shQuote,character(1)),stdout=TRUE,stderr=TRUE))
  if(!is.null(attr(output,"status"))&&attr(output,"status")!=0)stop(command," failed: ",paste(output,collapse="\n"))
  output
}
active_jobs<-function(known)intersect(ids(known),ids(run_read("squeue",c("--me","--noheader","--format=%A"))))
find_named_job<-function(name) {
  live<-ids(run_read("squeue",c("--me","--noheader",paste0("--name=",name),"--format=%A")))
  archived<-run_read("sacct",c("-X","--noheader","--parsable2","--starttime=2026-10-19",
    paste0("--user=",Sys.getenv("USER")),paste0("--name=",name),"--format=JobIDRaw,JobName%100"))
  recovered<-vapply(strsplit(archived,"|",fixed=TRUE),function(z)if(length(z)>=2L && z[2]==name)z[1]else"",character(1))
  ids(c(live,recovered))
}
terminal_job<-function(id) {
  x<-run_read("sacct",c("-X","--noheader","--parsable2",paste0("--jobs=",id),"--format=JobIDRaw,State"))
  rows<-strsplit(x,"|",fixed=TRUE)
  states<-vapply(rows,function(z)if(length(z)>=2L&&z[1]==id)sub("[ +].*$","",z[2])else"",character(1))
  any(states %in% c("COMPLETED","FAILED","CANCELLED","TIMEOUT","OUT_OF_MEMORY","NODE_FAIL","PREEMPTED","BOOT_FAIL","DEADLINE","REVOKED"))
}
submit_job<-function(args) {
  output<-suppressWarnings(system2("sbatch",vapply(args,shQuote,character(1)),stdout=TRUE,stderr=TRUE))
  status<-attr(output,"status");accepted<-grep("^[0-9]+(;[^[:space:]]+)?$",trimws(output),value=TRUE)
  if((is.null(status)||status==0)&&length(accepted)==1L)return(sub(";.*$","",accepted))
  error<-simpleError(paste(output,collapse="\n"))
  if(!length(accepted)&&any(grepl("sbatch: error:",output,fixed=TRUE)))class(error)<-c("definite_rejection",class(error))
  stop(error)
}
ensure_job<-function(experiment,stage,name,args,submit=submit_job,lookup=find_named_job) {
  record<-file.path(experiment,"scheduler",paste0(stage,".job"))
  intent<-file.path(experiment,"scheduler",paste0(stage,".intent"))
  if(file.exists(record)) {
    id<-ids(readLines(record,warn=FALSE));if(length(id)!=1L)stop("Invalid saved job record")
    return(id)
  }
  found<-ids(lookup(name))
  if(length(found)>1L)stop("Multiple jobs share a unique stage name; refusing submission: ",name)
  if(length(found)==1L){atomic_lines(found,record);return(found)}
  if(file.exists(intent))stop("Earlier submission outcome is uncertain for ",name,". No duplicate submitted. Rerun after its Slurm record becomes visible; inspect scheduler/",stage,".intent if no job was accepted.")
  atomic_lines(c(name,format(Sys.time(),tz="UTC",usetz=TRUE)),intent)
  result<-tryCatch(submit(c("--parsable",paste0("--job-name=",name),args)),error=function(e)e)
  if(inherits(result,"error")) {
    if(inherits(result,"definite_rejection"))unlink(intent)
    stop(result)
  }
  stopifnot(length(ids(result))==1L);atomic_lines(result,record);result
}
all_work_jobs<-function(root) {
  paths<-c(Sys.glob(file.path(root,"experiments","*","manifest.tsv")),
    Sys.glob(file.path(root,"benchmarks","*","manifest.tsv")),
    Sys.glob(file.path(root,"experiments","*","existing_evaluations.tsv")))
  recorded<-unlist(lapply(paths,function(path){x<-read.delim(path,stringsAsFactors=FALSE)
    unlist(x[intersect(names(x),c("fit_job","check_job","job_id"))],use.names=FALSE)
  }),use.names=FALSE)
  pilot<-unlist(lapply(Sys.glob(file.path(root,"runs","*","run_config.rds")),function(path){
    tryCatch(readRDS(path)$job_id,error=function(e)character())
  }),use.names=FALSE)
  ids(c(recorded,pilot))
}
