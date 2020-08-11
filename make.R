# https://books.ropensci.org/drake/hpc.html#parallel-computing-within-targets
envir <- new.env(parent = globalenv())

source("code/levellogger_packages.R", local = envir)
source("code/levellogger_functions.R", local = envir)
source("code/levellogger_plan.R", local = envir)

# Run the plan from code/levellogger_plan.R
make(
  envir$plan,
  verbose = 2,
  envir = envir
)
