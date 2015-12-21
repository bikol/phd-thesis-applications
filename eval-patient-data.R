# ---- init ----

rm(list=ls())

source("config.R")
source("stats.R")
source("similarities.R")
source("mcn-test.R")
source("utils.R")
source("prototypes.R")

library(parallel)
library(reshape2)
library(dplyr)
library(igraph)

EVALUATION.OUTPUT.FILE = 'patients-eval-output.RData'
EVALUATION.OUTPUT.LOCATION = paste(DATASETS.DIR, EVALUATION.OUTPUT.FILE, sep='/')

# ---- read-datasets ----

printDebug("read datasets")

DATA.COLS.NUM = 17

colClass = c("factor", "numeric", "integer", "integer",
             rep("numeric", times=DATA.COLS.NUM))

ds.training = read.csv(PATIENT.TRAINING.LOCATION, colClasses=colClass)
ds.test     = read.csv(PATIENT.TEST.LOCATION,     colClasses=colClass)

# ---- convert datasets into interval format

ds.training = bind_cols(ds.training[,1:4], convert.to.interval.format(ds.training[,5:ncol(ds.training)]))
ds.test = bind_cols(ds.test[,1:4], convert.to.interval.format(ds.test[,5:ncol(ds.test)]))

# ---- define-prototype-build-strategy ----
printDebug("protype build strategy")

PROTOTYPES = list(
    list(m=matrix(gen.0, nrow=2), type=0),
    list(m=matrix(post.0, nrow=2), type=0),
    list(m=matrix(pre.1, nrow=2), type=1),
    list(m=matrix(post.1, nrow=2), type=1),
    list(m=matrix(pre.2, nrow=2), type=1),
    list(m=matrix(post.2, nrow=2), type=1)
)

# ---- main-evaluation-procedure
source('eval.R')

# ---- save-evaluation ----

printDebug("save evaluation")

save.image(EVALUATION.OUTPUT.LOCATION)
