library(brms)
library(multcomp)
library(drake)
library(purrr)
library(data.table)
library(ingestr)
library(glue)
# library(httr)
# library(lme4)
# library(tydygraphs)
library(emmeans)
library(ggplot2)
library(patchwork)
library(ggforce)
library(fs)
library(stringr)
# library(mcp)
# library(broom)
library(errors)
# library(glmnet)
library(fst)
# library(disk.frame)
library(parallel)
library(ggthemes)

# setup_disk.frame(4)
options(future.globals.maxSize = Inf)
setDTthreads(2)
threads_fst(2)
options(future.supportsMulticore.unstable = TRUE,
        mc.cores = 2)
options(errors.notation = "plus-minus")

pale_pal <- 
  c(green = "#7DA050",
    orange = "#D19648",
    purple = "#966283",
    teal = "#329985",
    blue = "#5982A0",
    red = "#9B5249")

options(ggplot2.discrete.colour = unname(pale_pal),
        ggplot2.discrete.fill = unname(pale_pal))
