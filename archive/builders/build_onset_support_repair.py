"""Build a checks-only recovery package; do not start R, Stan, SSH or fits."""
from pathlib import Path, PurePosixPath
import base64, hashlib, io, tarfile, textwrap, json
ROOT=Path(__file__).resolve().parent
OLD=ROOT/'DCM_Ibex_onset_tests'
PKG=ROOT/'DCM_Ibex_onset_support_repair_r1'
PKG.mkdir(exist_ok=True)
def put(name,text):
    (PKG/name).write_text(text.strip()+'\n',encoding='utf-8',newline='\n')

original=(OLD/'onset_metrics.R').read_text(encoding='utf-8-sig')
start=original.index('onset_curve <- function(')
end=original.index('onset_from_curve <- function(',start)
new=r'''
# Only onset metrics need resolved vertical support; all observations stay in RMSE.
longest_supported_run <- function(good) {
  if(!length(good))return(0L)
  x<-rle(good);if(any(x$values))max(x$lengths[x$values])else 0L
}
onset_curve <- function(depth,value) {
  if(length(depth)!=length(value))stop("Depth/value lengths differ")
  empty<-function()list(depth=numeric(),value=numeric(),support_ok=FALSE)
  if(!length(depth))return(empty())
  if(any(!is.finite(depth))||any(!is.finite(value)))stop("Nonfinite input to onset scorer")
  z<-order(depth);depth<-depth[z];value<-value[z]
  if(anyDuplicated(depth))stop("Duplicate depths in onset evaluation")
  lo<-ceiling(max(10,min(depth)));hi<-floor(min(150,max(depth)))
  if(lo>hi)return(empty())
  grid<-seq.int(lo,hi)
  # approx needs two points; a singleton is retained and flagged below.
  y<-if(length(depth)==1L)rep(value,length(grid))else approx(depth,value,xout=grid,rule=1)$y
  j<-findInterval(grid,depth);between<-j>0L & j<length(depth)
  bad<-between & depth[pmin(j+1L,length(depth))]-depth[pmax(j,1L)]>3 & !grid %in% depth
  y[bad]<-NA_real_;sm<-rep(NA_real_,length(y));good<-which(is.finite(y))
  if(length(good))for(ii in split(good,cumsum(c(TRUE,diff(good)>1L))))
    sm[ii]<-if(length(ii)>=5L)stats::runmed(y[ii],k=5L,endrule="median")else y[ii]
  list(depth=grid,value=sm,support_ok=length(grid)>=7L && longest_supported_run(is.finite(sm))>=7L)
}
'''
fixed=original[:start]+new+original[end:]
fixed=fixed.replace('  if(!length(eligible))return(list(depth=NA_real_,status="insufficient_support"))',r'''
  if(length(z)<7L || longest_supported_run(is.finite(y))<7L ||
     longest_supported_run(is.finite(y)&z>=15)<persistence)
    return(list(depth=NA_real_,status="insufficient_support",boundary=NA_real_))''')
fixed=fixed.replace('      mld=tab$MLD_temp_MIX[1],observed_status=a$status,predicted_status=b$status,',r'''
      mld=tab$MLD_temp_MIX[1],observed_rows=nrow(tab),
      minimum_depth=min(tab$depth),maximum_depth=max(tab$depth),
      grid_points=length(observed$depth),
      maximum_supported_run=longest_supported_run(is.finite(observed$value)),
      observed_status=a$status,predicted_status=b$status,''')
fixed=fixed.replace('      posterior_gap_censored_draws=sum(statuses=="gap_censored"),',
 '''      posterior_insufficient_support_draws=sum(statuses=="insufficient_support"),
      posterior_gap_censored_draws=sum(statuses=="gap_censored"),''')
fixed=fixed.replace('  data.frame(casts=nrow(x),observed_detected=',
 '''  data.frame(casts=nrow(x),
    observed_insufficient_support=sum(x$observed_status=="insufficient_support"),
    predicted_insufficient_support=sum(x$predicted_status=="insufficient_support"),
    reached_comparable=sum(!is.na(x$reached_agreement)),observed_detected=''',1)
fixed=fixed.replace('  metadata <- unique(tab[c("cast_key","Region","MLD_temp_MIX","log10_distance")])',
 '''  metadata <- unique(tab[c("cast_key","Region","MLD_temp_MIX","log10_distance")])
  eligible <- unique(out$cast_key[out$observed_status!="insufficient_support" &
                                  out$predicted_status!="insufficient_support"])
  metadata <- metadata[as.character(metadata$cast_key) %in% eligible,,drop=FALSE]
  if(!nrow(metadata)) {
    writeLines("No onset-resolved profiles available for this example plot",
      file.path(checks,paste0(tag,"_onset_plot_skipped.txt")))
    return(invisible(aggregate))
  }''',1)
fixed=fixed.replace('split(metadata,metadata$Region)', 'split(metadata,metadata$Region,drop=TRUE)')
fixed=fixed.replace('split(train,train$Region)', 'split(train,train$Region,drop=TRUE)')
put('onset_metrics.R',fixed)
original_pred=(OLD/'predictive.R').read_text(encoding='utf-8-sig')
pred=original_pred.replace('    if(dataset!="training")save_onset_metrics(tab,mu,checks,tag)\n','')
needle='    for(n in names(sc))result[[paste(tag,n,sep="_")]]<-sc[[n]]'
assert pred.count(needle)==1
pred=pred.replace(needle,needle+r'''
    # Save standard scores before onset-specific work; failed stages remain visible.
    write.csv(result,file.path(run,"test_summary.csv"),row.names=FALSE)
    write.csv(do.call(rbind,scores),file.path(checks,"prediction_scores.csv"),row.names=FALSE)
    writeLines(paste(tag,"onset diagnostics"),file.path(run,"CHECK_STAGE.txt"))
    if(dataset!="training")save_onset_metrics(tab,mu,checks,tag)
''')
put('predictive.R',pred)
put('ledger.R',(OLD/'submit.R').read_text(encoding='utf-8-sig').split('queue_cases<-function(')[0])

put('worker.R',r'''
# Reuse the existing check scripts and saved fits; never enter the fitting script.
experiment<-normalizePath(Sys.getenv("EXPERIMENT_DIR"),mustWork=TRUE)
repair<-normalizePath(Sys.getenv("ONSET_REPAIR_DIR"),mustWork=TRUE)
source(file.path(experiment,"common.R"));source(file.path(experiment,"campaign.R"))
case<-read_case(experiment,Sys.getenv("CASE_ID"));stopifnot(case$mode=="posterior")
run<-file.path(experiment,"runs",sprintf("case_%02d",case$case_id))
dir.create(run,recursive=TRUE,showWarnings=FALSE)
stopifnot(startsWith(normalizePath(run),paste0(experiment,"/runs/")))
if(file.exists(file.path(run,"CHECKS_COMPLETE.txt")))quit(status=0L)
donor<-run
if(case$execution=="reuse") {
  c<-read.delim(file.path(experiment,"source_controls.tsv"),stringsAsFactors=FALSE)
  donor<-as.character(c$source_run[c$case_id==case$case_id]);stopifnot(length(donor)==1L)
}
if(!all(file.exists(file.path(donor,c("FIT_COMPLETE.txt","fit.rds"))))) {
  writeLines("Source fitting did not complete; no replacement fit submitted",
    file.path(run,"CHECKS_FAILED.txt"));stop("No completed source fit")
}
# Preserve earlier check outputs after any original checking job has ended.
backup<-file.path(repair,"case_backups",sprintf("case_%02d",case$case_id),
  paste0(format(Sys.time(),"%Y%m%dT%H%M%SZ",tz="UTC"),"_",Sys.getpid()))
dir.create(backup,recursive=TRUE)
checks<-file.path(run,"checks")
if(dir.exists(checks)) {
  if(!startsWith(normalizePath(checks),paste0(experiment,"/runs/")))stop("Unexpected checks path")
  if(!file.rename(checks,file.path(backup,"checks")))stop("Cannot preserve previous checks")
}
for(n in c("test_summary.csv","diagnostic_summary.csv","CHECKS_FAILED.txt","CHECK_STAGE.txt")) {
  p<-file.path(run,n)
  if(file.exists(p) && !file.copy(p,file.path(backup,n)))stop("Cannot preserve ",n)
}
failure<-file.path(run,"CHECKS_FAILED.txt")
if(file.exists(failure) && !file.rename(failure,file.path(backup,"previous_failure_marker.txt")))
  stop("Cannot archive previous failure marker")
writeLines(paste("Existing fit reused by evaluation job",Sys.getenv("SLURM_JOB_ID")),
  file.path(run,"ONSET_SUPPORT_REPAIR.txt"))
tryCatch({
  for(p in c("brms","posterior","mgcv","ggplot2","bayesplot"))
    if(!requireNamespace(p,quietly=TRUE))stop("Missing package: ",p)
  source(file.path(experiment,if(case$execution=="reuse")"reuse_check.R" else "03_check_case.R"))
  if(!file.exists(file.path(run,"CHECKS_COMPLETE.txt")))stop("Evaluation did not finish")
},error=function(e){writeLines(conditionMessage(e),file.path(run,"CHECKS_FAILED.txt"));stop(e)})
''')
put('test_and_audit.R',r'''
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
''')
put('submit_repair.R',r'''
# The only allowed submission targets here are check.slurm and collect.slurm.
a<-commandArgs(trailingOnly=TRUE);root<-normalizePath(a[1]);package<-normalizePath(a[2])
experiment<-normalizePath(a[3]);repair<-normalizePath(Sys.getenv("ONSET_REPAIR_DIR"))
source(file.path(package,"ledger.R"))
Sys.setenv(DCM_ROOT=root,EXPERIMENT_DIR=experiment,ONSET_REPAIR_DIR=repair)
cases<-read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
original<-read.delim(file.path(experiment,"manifest.tsv"),stringsAsFactors=FALSE)
controls<-read.delim(file.path(experiment,"source_controls.tsv"),stringsAsFactors=FALSE)
stopifnot(nrow(cases)==34L,sum(cases$mode=="posterior")==30L,!anyDuplicated(original$case_id))
path<-file.path(repair,"repair_manifest.tsv")
manifest<-if(file.exists(path))read.delim(path,stringsAsFactors=FALSE)else
  data.frame(case_id=integer(),original_fit_job=character(),original_check_job=character(),check_job=character(),run_dir=character())
token<-hash_text(experiment);base<-c(paste0("--chdir=",root),"--export=ALL")
for(i in which(cases$mode=="posterior")) {
  case<-cases[i,];item<-original[original$case_id==case$case_id,,drop=FALSE];stopifnot(nrow(item)==1L)
  run<-as.character(item$run_dir)
  if(file.exists(file.path(run,"CHECKS_COMPLETE.txt"))) {cat("Already checked:",case$case_id,"\n");next}
  donor<-run;fit_id<-as.character(item$fit_job);work<-ids(c(item$fit_job,item$check_job))
  if(case$execution=="reuse") {
    old<-controls[controls$case_id==case$case_id,,drop=FALSE];stopifnot(nrow(old)==1L)
    donor<-as.character(old$source_run);fit_id<-as.character(old$fit_job);work<-ids(c(work,old$fit_job,old$check_job))
  }
  if(!file.exists(file.path(donor,"FIT_COMPLETE.txt")) && terminal_job(fit_id)) {
    cat("Source fit ended without completion; no replacement fit:",case$case_id,"\n");next
  }
  Sys.setenv(CASE_ID=case$case_id)
  for(attempt in 1:3) {
    waiting<-active_jobs(work)
    if(!file.exists(file.path(donor,"FIT_COMPLETE.txt")) && !terminal_job(fit_id))waiting<-ids(c(waiting,fit_id))
    if(!any(file.exists(file.path(run,c("CHECKS_COMPLETE.txt","CHECKS_FAILED.txt")))) &&
       !terminal_job(as.character(item$check_job)))waiting<-ids(c(waiting,item$check_job))
    dep<-if(length(waiting))paste0("--dependency=afterany:",paste(waiting,collapse=":"))else character()
    job<-tryCatch(ensure_job(repair,sprintf("check_%02d",case$case_id),sprintf("dcm-or1-%s-q%02d",token,case$case_id),
      c(base,dep,file.path(repair,"check.slurm"))),error=function(e)e)
    if(!inherits(job,"error"))break
    invalid<-inherits(job,"definite_rejection")&&grepl("Job dependency problem|Invalid job id",conditionMessage(job))
    if(!invalid||attempt==3L)stop(job)
  }
  j<-match(case$case_id,manifest$case_id)
  row<-data.frame(case_id=case$case_id,original_fit_job=as.character(item$fit_job),
    original_check_job=as.character(item$check_job),check_job=job,run_dir=run)
  if(is.na(j))manifest<-rbind(manifest,row)else if(as.character(manifest$check_job[j])!=job)stop("Repair ledger mismatch")
  atomic_table(manifest,path);cat("Registered evaluation only: case",case$case_id,"job",job,"\n")
}
known<-ids(manifest$check_job);signature<-hash_text(sort(known))
for(attempt in 1:3) {
  waiting<-active_jobs(known)
  dep<-if(length(waiting))paste0("--dependency=afterany:",paste(waiting,collapse=":"))else character()
  collector<-tryCatch(ensure_job(repair,paste0("collect_",signature),paste0("dcm-or1-",token,"-results-",signature),
    c(base,dep,file.path(repair,"collect.slurm"))),error=function(e)e)
  if(!inherits(collector,"error"))break
  invalid<-inherits(collector,"definite_rejection")&&grepl("Job dependency problem|Invalid job id",conditionMessage(collector))
  if(!invalid||attempt==3L)stop(collector)
}
atomic_lines(collector,file.path(repair,"collection_job_id.txt"))
cat("Collector:",collector,"; existing fit registrations unchanged; no fitting jobs submitted or cancelled.\n")
''')
slurm=(OLD/'check.slurm').read_text().replace('"$EXPERIMENT_DIR/03_check_case.R"','"$ONSET_REPAIR_DIR/worker.R"')
slurm=slurm.replace('srun --ntasks=1', '''# Prevent concurrent recovery workers for this case.
exec 9>"$ONSET_REPAIR_DIR/case_${CASE_ID}.lock"
flock -n 9 || { echo 'Another recovery check is active.' >&2; exit 1; }
srun --ntasks=1''')
put('check.slurm',slurm)
collect=(OLD/'collect.slurm').read_text()
put('collect.slurm',collect)
before={n:hashlib.sha256((OLD/n).read_bytes()).hexdigest() for n in ('onset_metrics.R','predictive.R')}
after={n:hashlib.sha256((PKG/n).read_bytes()).hexdigest() for n in before}
put('original_SHA256SUMS.txt','\n'.join(f'{v}  {n}' for n,v in before.items()))
put('patched_SHA256SUMS.txt','\n'.join(f'{v}  {n}' for n,v in after.items()))
put('patch_files.sh',r'''
#!/bin/bash
set -euo pipefail
for name in onset_metrics.R predictive.R; do
  old_hash=$(awk -v n="$name" '$2==n {print $1}' "$ONSET_REPAIR_DIR/original_SHA256SUMS.txt")
  new_hash=$(awk -v n="$name" '$2==n {print $1}' "$ONSET_REPAIR_DIR/patched_SHA256SUMS.txt")
  actual=$(sha256sum "$EXPERIMENT_DIR/$name" | cut -d' ' -f1)
  if [[ "$actual" != "$old_hash" && "$actual" != "$new_hash" ]]; then
    echo "Unexpected scoring version: $name; no scoring files changed." >&2; exit 1
  fi
done
mkdir -p "$ONSET_REPAIR_DIR/scoring_backups"
for name in onset_metrics.R predictive.R; do
  new_hash=$(awk -v n="$name" '$2==n {print $1}' "$ONSET_REPAIR_DIR/patched_SHA256SUMS.txt")
  actual=$(sha256sum "$EXPERIMENT_DIR/$name" | cut -d' ' -f1)
  [[ "$actual" == "$new_hash" ]] && continue
  [[ -e "$ONSET_REPAIR_DIR/scoring_backups/$name" ]] || cp -- "$EXPERIMENT_DIR/$name" "$ONSET_REPAIR_DIR/scoring_backups/$name"
  cp -- "$ONSET_REPAIR_DIR/$name" "$EXPERIMENT_DIR/$name.repair_r1.tmp"
  mv -- "$EXPERIMENT_DIR/$name.repair_r1.tmp" "$EXPERIMENT_DIR/$name"
done
(cd "$EXPERIMENT_DIR"; sha256sum --check "$ONSET_REPAIR_DIR/patched_SHA256SUMS.txt")
''')
put('run_repair.sh',r'''
#!/bin/bash
# Recover checking only; repeating resumes the same accepted registrations.
set -euo pipefail
export DCM_ROOT="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
cd "$DCM_ROOT"
package=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export EXPERIMENT_DIR="${ONSET_EXPERIMENT:-$(cat logs/latest_onset150_experiment.txt)}"
export ONSET_REPAIR_DIR="$EXPERIMENT_DIR/repairs/onset_support_r1"
mkdir -p "$ONSET_REPAIR_DIR/scheduler"
exec 9>"$ONSET_REPAIR_DIR/submission.lock"
flock -n 9 || { echo 'Another recovery submission is active.' >&2; exit 1; }
source "$EXPERIMENT_DIR/ibex_environment.sh"
(cd "$package"; sha256sum --check --quiet SHA256SUMS.txt)
for file in "$package/"*; do
  [[ -f "$file" ]] || continue
  name=$(basename "$file")
  if [[ -f "$ONSET_REPAIR_DIR/$name" ]]; then
    cmp --silent "$file" "$ONSET_REPAIR_DIR/$name" || { echo "Different recovery version: $name" >&2; exit 1; }
  else
    cp -- "$file" "$ONSET_REPAIR_DIR/$name"
  fi
done
for file in "$ONSET_REPAIR_DIR/"*.sh "$ONSET_REPAIR_DIR/"*.slurm; do bash -n "$file"; done
# Test new scoring on fixtures and every real held-out cast before patching.
Rscript --vanilla "$ONSET_REPAIR_DIR/test_and_audit.R" "$ONSET_REPAIR_DIR" "$EXPERIMENT_DIR"
bash "$ONSET_REPAIR_DIR/patch_files.sh"
Rscript --vanilla "$ONSET_REPAIR_DIR/submit_repair.R" "$DCM_ROOT" "$ONSET_REPAIR_DIR" "$EXPERIMENT_DIR"
printf '\nRecovery directory: %s\nSupport audit: %s/actual_cohort_support_audit.csv\n' "$ONSET_REPAIR_DIR" "$ONSET_REPAIR_DIR"
''')
put('export_compact.sh',r'''
#!/bin/bash
# Export all comparison tables while leaving plots and large raw files on Ibex.
set -euo pipefail
cd "${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
experiment="${ONSET_EXPERIMENT:-$(cat logs/latest_onset150_experiment.txt)}"
archive="$PWD/logs/dcm_onset_compact_$(date -u +%Y%m%dT%H%M%SZ)_$$.tar.gz"
list=$(mktemp)
trap 'rm -f -- "$list"' EXIT
find "$experiment/results" -maxdepth 1 -type f -name '*.csv' -print0 > "$list"
find experiments benchmarks \
  -type d \( -name stan_csv -o -name repairs -o -name 'scheduler_backup_*' \) -prune -o \
  -type f \( -name cases.tsv -o -name manifest.tsv -o -name source_controls.tsv \
    -o -name existing_evaluations.tsv -o -name test_summary.csv \
    -o -name diagnostic_summary.csv -o -name run_config.txt \
    -o -name CHECKS_FAILED.txt -o -name FIT_FAILED.txt -o -name CHECKS_COMPLETE.txt \) -print0 >> "$list"
find "$experiment/repairs/onset_support_r1" -maxdepth 1 -type f \
  \( -name actual_cohort_support_audit.csv -o -name repair_manifest.tsv -o -name README.md \) -print0 >> "$list"
find "$experiment/runs" -maxdepth 3 -type f -name parameter_diagnostics.csv -print0 >> "$list"
for name in campaign.R base_models.R data_identity_audit.csv data_checksums.csv; do
  [[ ! -f "$experiment/$name" ]] || printf '%s\0' "$experiment/$name" >> "$list"
done
status=0
tar --null --files-from="$list" -czf "$archive" || status=$?
(( status<=1 )) || exit "$status"
printf '\nCompact comparison archive:\n%s\n' "$archive"
ls -lh -- "$archive"
''')
put('README.md','''
# Onset support repair R1

The old scorer aborted a whole posterior check when a profile had fewer than
seven one-metre grid points. Unsupported onset metrics are now explicitly flagged
insufficient_support, including singletons and insufficient support after 15 m.
Every observation remains in ordinary prediction errors; missing onset is neither
imputed nor treated as absence. Summaries count unsupported and comparable casts.
Standard scores are saved before onset-specific calculations.

Only onset_metrics.R and predictive.R in the selected experiment are patched,
with exact source hashes and backups. Data, cases, formulas, priors, fits, sampler
settings and original job registrations are preserved. Recovery workers archive
failed checking outputs and reuse the existing diagnostic/prediction scripts.
They wait for active original jobs, including unfinished case 18, and skip completed
checks and prior-predictive cases. No fitting jobs are submitted or cancelled.

Local checks: design/source identities, R delimiters, permitted patch targets,
safe archive paths, payload equality and SHA256. Actual R parsing and tests run
on Ibex before queueing. Tests cover short/empty/singleton, shifted/absent/censored/
gapped onsets, posterior uncertainty, denominators and submission recovery. All
real validation/reserve casts are audited, with counts and cast IDs exported.
Local R is not retried after the earlier automatic-review block.

Run bash repair_onset_checks.sh after upload. Repeating resumes accepted checking
registrations. A final collector refreshes all registered campaigns. Later use
bash DCM_Ibex_onset_tests/collect_results.sh to refresh completed results manually.
''')
def balanced(path):
    stack=[];quote=None;esc=False;comment=False
    for ch in path.read_text():
        if comment:
            if ch=='\n':comment=False
            continue
        if quote:
            if esc:esc=False
            elif ch=='\\':esc=True
            elif ch==quote:quote=None
            continue
        if ch=='#':comment=True
        elif ch in ('"',"'",'`'):quote=ch
        elif ch in '({[':stack.append(ch)
        elif ch in ')}]':assert stack and stack.pop()=={')':'(',']':'[','}':'{'}[ch],path
    assert not quote and not stack,path
for p in PKG.glob('*.R'):balanced(p)
assert 'stop("Insufficient vertical support")' not in fixed
assert 'observed_insufficient_support' in fixed and 'predicted_insufficient_support' in fixed
assert pred.index('write.csv(result,file.path(run,"test_summary.csv")')<pred.index('save_onset_metrics(tab,mu,checks,tag)')
assert 'fit.slurm' not in (PKG/'submit_repair.R').read_text()
assert '02_fit_case' not in (PKG/'worker.R').read_text()
assert set(before)=={'onset_metrics.R','predictive.R'}
files=sorted(p for p in PKG.iterdir() if p.is_file() and p.name!='SHA256SUMS.txt')
put('SHA256SUMS.txt',''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in files))
files.append(PKG/'SHA256SUMS.txt');buf=io.BytesIO()
with tarfile.open(fileobj=buf,mode='w:gz',format=tarfile.USTAR_FORMAT) as tar:
    for p in files:
        b=p.read_bytes();info=tarfile.TarInfo(PKG.name+'/'+p.name)
        info.size=len(b);info.mode=0o644;info.mtime=1792580400;tar.addfile(info,io.BytesIO(b))
payload=buf.getvalue();digest=hashlib.sha256(payload).hexdigest()
header='''#!/bin/bash
# Repair onset checks using saved fits; no fitting submissions or cancellations.
set -euo pipefail
bash -n "$0"
root="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
cd "$root"
root="$PWD"
package="$root/DCM_Ibex_onset_support_repair_r1"
mkdir -p logs
payload=$(mktemp "$root/logs/onset_repair_payload.XXXXXXXX.tar.gz")
trap 'rm -f -- "$payload"' EXIT
base64 --decode > "$payload" <<'__ONSET_SUPPORT_REPAIR__'
'''
footer='''
__ONSET_SUPPORT_REPAIR__
printf '%s  %s\\n' '__HASH__' "$payload" | sha256sum --check --status
[[ ! -L "$package" ]] || { echo 'Unexpected package symlink.' >&2; exit 1; }
if [[ ! -e "$package" ]]; then
  stage=$(mktemp -d "$root/.onset_support_install.XXXXXXXX")
  tar -xzf "$payload" -C "$stage" --no-same-owner
  mv -- "$stage/DCM_Ibex_onset_support_repair_r1" "$package"
  rmdir -- "$stage"
fi
expected=$(tar -xOzf "$payload" DCM_Ibex_onset_support_repair_r1/SHA256SUMS.txt | sha256sum | cut -d' ' -f1)
actual=$(sha256sum "$package/SHA256SUMS.txt" | cut -d' ' -f1)
[[ "$actual" == "$expected" ]] || { echo 'Different existing recovery package; no overwrite.' >&2; exit 1; }
(cd "$package"; sha256sum --check --quiet SHA256SUMS.txt)
bash "$package/run_repair.sh"
'''.replace('__HASH__',digest)
out=ROOT/'repair_onset_checks.sh'
text=header+'\n'.join(textwrap.wrap(base64.b64encode(payload).decode(),76))+footer
out.write_text(text,encoding='ascii',newline='\n')
with tarfile.open(fileobj=io.BytesIO(payload),mode='r:gz') as tar:
    for m in tar.getmembers():
        p=PurePosixPath(m.name)
        assert m.isfile() and not p.is_absolute() and '..' not in p.parts and p.parts[0]==PKG.name
        assert tar.extractfile(m).read()==(ROOT/m.name).read_bytes()
assert b'\r' not in out.read_bytes()
print(json.dumps(dict(installer=str(out),bytes=out.stat().st_size,payload_sha256=digest,
  changed_scoring_files=list(before),fitting_jobs_submitted=0)))
