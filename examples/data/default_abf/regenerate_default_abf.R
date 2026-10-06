# Regenerate only the two literal default-prior ABF calls from the pinned fixture.
# Usage: Rscript --vanilla export_default_abf.R coloc_test_data.rda OUTPUT_DIRECTORY
# Use the coloc 6.0.3 / susieR 0.16.6 libraries in reference.lock.json.
args <- commandArgs(TRUE)
stopifnot(length(args)==2)
suppressPackageStartupMessages(library(coloc))
library(jsonlite)
stopifnot(as.character(packageVersion('coloc'))=='6.0.3',
          as.character(packageVersion('susieR'))=='0.16.6')
load(args[1])
dir.create(args[2], recursive=TRUE, showWarnings=FALSE)
actual_defaults <- lapply(formals(coloc::coloc.abf)[c('p1','p2','p12')],eval)
stopifnot(actual_defaults$p1==1e-4, actual_defaults$p2==1e-4,
          actual_defaults$p12==1e-5)
cases <- list()
for(pair in list(c('D1','D2'),c('D3','D4'))) {
  d1 <- coloc_test_data[[pair[1]]]; d2 <- coloc_test_data[[pair[2]]]
  stopifnot(length(d1$snp)==500, length(d2$snp)==500,
            identical(d1$snp,d2$snp), !anyDuplicated(d1$snp))
  warnings <- character()
  # Intentionally no prior arguments: these are actual public defaults.
  result <- withCallingHandlers(coloc::coloc.abf(d1,d2), warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w)); invokeRestart('muffleWarning')
  })
  result$summary <- as.list(result$summary)
  cases[[paste0('abf_default_',paste(pair,collapse='_'))]] <- list(
    datasets=I(pair), call=paste0('coloc::coloc.abf(',pair[1],', ',pair[2],')'),
    supplied_prior_arguments=list(), actual_default_priors=actual_defaults,
    result=result, warnings=I(warnings), error=NULL)
}
write_json(list(reference='official coloc 6.0.3',
  commit='8f20f0bc5e60ffc99e4c2f787bd55fd30cfe7c45',
  input_sha256='81ddbbe3af502779da76d57a7e1bcc4d3a8d172fcf79bcbe15c8151915dad2c1',
  R_version=R.version.string, cases=cases),
  file.path(args[2],'r_default_abf.json'),auto_unbox=TRUE,digits=NA,
  pretty=FALSE,na='null',null='null',matrix='rowmajor')
writeLines(capture.output(sessionInfo()),file.path(args[2],'sessionInfo.private.txt'))
cat('PASS two literal default-prior ABF reference calls\n')
