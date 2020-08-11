# https://books.ropensci.org/drake/hpc.html#parallel-computing-within-targets
envir <- new.env(parent = globalenv())

source("code/levellogger_packages.R", local = envir)
source("code/levellogger_functions.R", local = envir)
source("code/levellogger_plan.R", local = envir)

drake_config(
  envir$plan, # defined in R/plan.R
  verbose = 2,
  envir = envir
)