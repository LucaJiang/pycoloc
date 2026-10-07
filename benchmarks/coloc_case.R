#!/usr/bin/env Rscript
# Full public calls on the same portable float64 inputs as coloc_case.py.
# Usage: Rscript --vanilla coloc_case.R CASE_JSON STAGE MODE OUTPUT [REPEATS] [BATCH] [CLOCK_SO]
# Set R_LIBS_USER to independently installed pinned coloc/susieR libraries.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) stop("Expected CASE_JSON STAGE MODE OUTPUT [REPEATS=7] [BATCH=10]")
case_file <- normalizePath(args[[1]], mustWork=TRUE)
stage <- args[[2]]
mode <- args[[3]]
output <- args[[4]]
repeats <- if (length(args) >= 5) as.integer(args[[5]]) else 7L
batch <- if (length(args) >= 6) as.integer(args[[6]]) else 10L
stopifnot(stage %in% c("abf", "bf", "susie"), mode %in% c("validation", "timing", "memory"),
          repeats >= 7L, batch >= 1L)
if (mode == "timing") {
  if (length(args) < 7L) stop("Timing requires a compiled monotonic_clock.c library as argument 7")
  dyn.load(normalizePath(args[[7]], mustWork=TRUE))
}
clock_symbol <- Sys.getenv("COLOC_CLOCK_SYMBOL", "coloc_monotonic")
clock <- function() .Call(clock_symbol)
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(coloc))
stopifnot(as.character(packageVersion("coloc")) == "6.0.3",
          as.character(packageVersion("susieR")) == "0.16.6")
spec <- fromJSON(case_file, simplifyVector=FALSE)
read_array <- function(s) {
  path <- file.path(dirname(case_file), s$file)
  shape <- as.integer(unlist(s$shape))
  count <- prod(shape)
  if (file.info(path)$size != count * 8) stop("Incorrect float64 byte count: ", s$file)
  con <- file(path, "rb")
  on.exit(close(con))
  a <- readBin(con, what="double", n=count, size=8, endian="little")
  if (length(shape) == 2L) matrix(a, nrow=shape[1], ncol=shape[2], byrow=TRUE) else a
}
read_dataset <- function(d) {
  lapply(d, function(x) if (is.list(x) && !is.null(x$file)) read_array(x) else unlist(x))
}
read_fit <- function(f) {
  bf <- read_array(f$lbf_variable)
  ids <- unlist(f$variant_ids)
  if (!is.null(f$null_index)) {
    stopifnot(as.integer(f$null_index) == length(f$variant_ids))
    ids <- c(ids, "null")
  }
  colnames(bf) <- ids
  fit <- list(lbf_variable=bf, converged=f$converged, niter=f$niter,
              sets=list(cs=lapply(f$sets$cs, function(x) as.integer(unlist(x)) + 1L),
                        cs_index=as.integer(unlist(f$sets$cs_index)) + 1L))
  class(fit) <- "susie"
  fit
}
if (stage == "abf") {
  datasets <- lapply(spec$datasets, read_dataset)
} else {
  fits <- lapply(spec$fits, read_fit)
}
priors <- spec$priors[[stage]]
pair_options <- c(priors, list(overlap.min=spec$overlap_min,
                              trim_by_posterior=spec$trim_by_posterior))
run_case <- switch(stage,
  abf=function() do.call(coloc::coloc.abf, c(list(dataset1=datasets[[1]], dataset2=datasets[[2]]), priors)),
  bf=function() do.call(coloc::coloc.bf_bf, c(list(bf1=fits[[1]]$lbf_variable, bf2=fits[[2]]$lbf_variable), pair_options)),
  susie=function() do.call(coloc::coloc.susie, c(list(dataset1=fits[[1]], dataset2=fits[[2]]), pair_options)))
status_of <- function(result) {
  if (stage == "susie" && any(vapply(fits, function(f) length(f$sets$cs) == 0, logical(1)))) return("no_cs")
  if (stage != "abf" && !length(intersect(colnames(fits[[1]]$lbf_variable), colnames(fits[[2]]$lbf_variable)))) return("no_overlap")
  if (is.null(result$summary)) return("unavailable")
  if (stage != "abf" && all(is.na(result$summary$PP.H4.abf))) return("insufficient_overlap")
  if (stage != "abf") {
    counts <- if (stage == "susie") vapply(fits, function(f) length(f$sets$cs), integer(1)) else vapply(fits, function(f) nrow(f$lbf_variable), integer(1))
    if (nrow(result$summary) < prod(counts)) return("partially_trimmed")
  }
  "ok"
}
peak_mib <- function() {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  line <- grep("^VmHWM:", readLines("/proc/self/status"), value=TRUE)
  if (length(line) != 1L) return(NA_real_)
  as.numeric(strsplit(trimws(sub("^VmHWM:", "", line)), "[[:space:]]+")[[1]][1]) / 1024
}
dir.create(dirname(output), recursive=TRUE, showWarnings=FALSE)
record <- list(case_id=spec$case_id, stage=stage, backend="coloc", mode=mode,
               status="error", error=NULL, records=list())
# Capture public printing in memory; disk export follows clocks and HWM capture.
# At the default bounded 71 calls this holds only the small printed summaries.
console <- rawConnection(raw(), open="w")
sink(console)
tryCatch({
  result <- NULL
  if (mode == "timing") {
    start <- clock()
    result <- run_case()
    record$warmup_seconds <- clock() - start
    result <- NULL
    for (repeat_id in seq_len(repeats)) {
      start_epoch <- as.numeric(Sys.time())
      start <- clock()
      for (i in seq_len(batch)) {
        result <- NULL
        result <- run_case()
      }
      elapsed <- clock() - start
      end_epoch <- as.numeric(Sys.time())
      record$records[[repeat_id]] <- list("repeat"=repeat_id-1L, time_seconds=elapsed/batch,
                                         batch_calls=batch, peak_rss_mib=NULL,
                                         start_epoch=start_epoch, end_epoch=end_epoch)
    }
  } else {
    result <- run_case()
  }
  peak <- peak_mib() # Before result serialization or session metadata collection.
  record$status <- status_of(result)
  record$peak_rss_mib <- if (mode == "memory") peak else NULL
  record$n_variants <- if (stage == "abf") unname(result$summary[["nsnps"]]) else
    length(intersect(colnames(fits[[1]]$lbf_variable), colnames(fits[[2]]$lbf_variable)))
  record$L <- if (stage == "abf") NULL else lapply(fits, function(f) nrow(f$lbf_variable))
  record$niter <- if (stage == "abf") NULL else lapply(fits, function(f) f$niter)
  record$converged <- if (stage == "abf") NULL else lapply(fits, function(f) f$converged)
  if (mode == "validation") {
    summary <- if (stage == "abf") as.list(result$summary) else as.data.frame(result$summary)
    record$result <- list(status=record$status, summary=summary,
                          results=if (is.null(result$results)) list() else as.data.frame(result$results),
                          priors=result$priors)
  }
}, error=function(e) { record$error <<- conditionMessage(e) })
sink()
writeBin(rawConnectionValue(console), paste0(output, ".console.log"))
close(console)
probe <- Sys.getenv("COLOC_BLAS_PROBE", "")
if (nzchar(probe)) dyn.load(normalizePath(probe, mustWork=TRUE))
blas_runtime <- if (is.loaded("prusie_blas_info")) .Call("prusie_blas_info") else NULL
description <- utils::packageDescription("coloc")
record$environment <- list(R=R.version.string, version=as.character(packageVersion("coloc")),
                           susieR=as.character(packageVersion("susieR")),
                           coloc_path=find.package("coloc"), susieR_path=find.package("susieR"),
                           blas=extSoftVersion()[["BLAS"]], blas_runtime=blas_runtime,
                           installed_remote_sha=description$RemoteSha,
                           expected_reference_commit="8f20f0bc5e60ffc99e4c2f787bd55fd30cfe7c45",
                           installed_rdb_md5=unname(tools::md5sum(file.path(find.package("coloc"), "R/coloc.rdb"))),
                           session=capture.output(sessionInfo()),
                           threads=as.list(Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"))))
write_json(record, output, auto_unbox=TRUE, pretty=TRUE, digits=NA, null="null", na="null")
cat(record$case_id, stage, mode, record$status, "\n")
if (!is.null(record$error)) quit(status=1L)
