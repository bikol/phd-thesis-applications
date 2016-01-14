# ---- init ----

rm(list=ls())

source("config.R")
source("stats.R")
source("similarities.R")
source("mcn-test.R")
source("utils.R")

library(parallel)
library(reshape2)
library(dplyr)
library(Matrix)
library(igraph)
library(RCurl)


EVALUATION.OUTPUT.FILE = 'dyslexic-eval-output.RData'
EVALUATION.OUTPUT.LOCATION = paste(DATASETS.DIR, EVALUATION.OUTPUT.FILE, sep='/')

# ---- read-datasets ----

printDebug("read datasets")



KEEL.DATABASE.LOCATION = paste0(DATASETS.DIR,"/keel-dyslexic/")

temp = tempfile()
content = getBinaryURL(paste('http://sci2s.ugr.es/keel/dataset/data/lowQuality/dyslexic_12_4-10-fold.zip'),
                       httpheader = c("User-Agent"="R"))
writeBin(content, temp)
unzip(temp,exdir=KEEL.DATABASE.LOCATION)
unlink(temp)
print(paste("File downloaded:", 'dyslexic_12_4-10-fold.zip'))

DATA.COLS.NUM = 12

# ---- define-prototype-build-strategy ----
printDebug("protype build strategy")
build.prototypes = function(ts){
    return(list(
        list(m=matrix(rep(c(0, 0.5), times=DATA.COLS.NUM), nrow=2), type=0),
        list(m=matrix(rep(c(0.5, 1), times=DATA.COLS.NUM), nrow=2), type=1),
        list(m=matrix(rep(c(0.4, 0.6), times=DATA.COLS.NUM), nrow=2), type=0),
        list(m=matrix(rep(c(0.4, 0.6), times=DATA.COLS.NUM), nrow=2), type=1)
    ))
}

# ---- parallel-init ----

printDebug("paralell init")

if (THREADS > 1)
{
    CL = makeCluster(THREADS, outfile="")

    clusterExport(CL, c("SEED", "OBSUCRE.REPEAT"))
    clusterCall(CL, function(){ set.seed(SEED)})

    clusterExport(cl=CL, list('topo_sort', 'graph_from_adjacency_matrix', 'DescIterBinSearch', 'AscIterBinSearch',
                              'KNN.MAX.CASE.BASE.SIZE', 'invPerm', 'PROBE.SIZE',
                              'IVFC', 'IVFC.NAME', 'KNN', 'KNN.NAME', 'build.prototypes',
                              'DATA.COLS.NUM', 'read.keel', 'KEEL.DATABASE.LOCATION',
                              'multiclassToOutcome', 'filter', 'bind_rows','mutate', '%>%',
                              'printDebug', 'DEBUG'))
    usedLapply = function(...){ clusterApplyLB(CL, ...) }
} else {
    usedLapply = lapply
}




# ---- training-statistics-knn ----
    if(!SKIP.KNN) {
        printDebug("training statistics knn")

        outcomes.knns = usedLapply(1:length(KNN), function(j){
            printDebug(paste0("Using ", KNN.NAME[[j]]," classifier"))
            diags = c()

            fold.size = PROBE.SIZE/10

            diags = sapply(1:10, function(fold) {
                # reading training and test set for 10 fold CV
                train.data = read.keel(paste0(KEEL.DATABASE.LOCATION, '/dyslexic_12_4-10-',fold,'tra.dat'))
                test.data = read.keel(paste0(KEEL.DATABASE.LOCATION, '/dyslexic_12_4-10-',fold,'tst.dat'))

                # normalise into [0,1] interval
                train.data[, 1:(DATA.COLS.NUM*2)] = train.data[, 1:(DATA.COLS.NUM*2)]/10
                test.data[, 1:(DATA.COLS.NUM*2)] = test.data[, 1:(DATA.COLS.NUM*2)]/10

                # multiply observations with multiple classes
                train.data = bind_rows(train.data,
                                       filter(train.data, !is.na(Class2)) %>% mutate(Class1=Class2),
                                       filter(train.data, !is.na(Class3)) %>% mutate(Class1=Class3))

                # convert training set into proper format accepted by classifier
                ts = apply(train.data, 1 , function(x){
                    return(list(m=matrix(as.numeric(x[1:(DATA.COLS.NUM*2)]), nrow=2), type=x[DATA.COLS.NUM*2+1]))
                })

                if(length(ts) > KNN.MAX.CASE.BASE.SIZE) {
                    classifier = KNN[[j]](ts[sample(length(ts), min(length(ts), KNN.MAX.CASE.BASE.SIZE))])
                } else {
                    classifier = KNN[[j]](ts)
                }

                # selection of appropriate columns
                tmp = apply(test.data[ , 1:(DATA.COLS.NUM*2)], 1, function(row) {
                    # matrix in format required by aggregation method is created and passed into classifier
                    return(classifier(matrix(row, nrow=2)))
                })

                converted = apply(cbind(tmp, test.data[ ,(DATA.COLS.NUM*2+1):(DATA.COLS.NUM*2+4)]),
                                  1, multiclassToOutcome)
                return(converted)
            })

            return(unlist(diags))
        })

        printDebug("finished knn classification")

        outcomes.knns        = data.frame(Dummy1=NA, Dummy2=0, Dummy3=1, outcomes.knns)
        names(outcomes.knns) = c('Dummy1', 'ObscureLevel', 'ObscureRepeat', KNN.NAME)

        training.stats.knns = melt(calculate.stats(
            aggregate.outcomes(outcomes.knns)
        ),
        id.vars = "ObscureLevel" ,
        variable.name = "Method",
        value.name = "Value") %>%
            rename(Measure=L1)

        training.stats.knns = suppressWarnings( # suppress different factor levels warning
            left_join(training.stats.knns,
                      KNN.BINDED.DESCRIPTION,
                      by="Method")
        )
    }

# ---- training-statistics-ivfc ----
if(!SKIP.TRAINING) {
    if(!SKIP.IVFC) {
        printDebug("training statistics ivfc")

        outcomes.ivfcs = usedLapply(1:length(IVFC), function(j){
            printDebug(paste0("Using ",IVFC.NAME[[j]]," classifier"))
            diags = c()
            # split training data set into separate probes
            for(i in 1:(nrow(ds.training)/PROBE.SIZE)) {
                # printDebug(paste("Data repeat(probe):", i))

                # select only cases from this probe
                all.data = ds.training[((i-1)*PROBE.SIZE+1):(i*PROBE.SIZE), ]
                # randomise input to obrain different fold for each probe
                shuffle = sample(nrow(all.data))
                all.data = all.data[shuffle, ]

                fold.size = PROBE.SIZE/10

                d = sapply(1:10, function(fold) {
                    # build training and test set for 10 fold CV
                    mask = rep(F, PROBE.SIZE)
                    mask[((fold-1)*fold.size+1):(fold*fold.size)] = T
                    test.data = all.data[mask,]
                    train.data = all.data[!mask,]

                    # convert training set into proper format accepted by classifier
                    ts = apply(train.data,1 , function(x){
                        return(list(m=matrix(as.numeric(x[5:(5+DATA.COLS.NUM*2-1)]), nrow=2), type=x[4]))
                    })

                    classifier = IVFC[[j]](build.prototypes(ts))

                    # selection of appropriate columns
                    tmp = apply(test.data[, 5:(5+DATA.COLS.NUM*2-1)], 1, function(row) {
                        # matrix in format required by aggregation method is created and passed into classifier
                        return(classifier(matrix(row, nrow=2)))
                    })

                    return(tmp)
                })
                # undo the shuffling to enable comparison with expected results
                d = c(d)[invPerm(shuffle)]
                diags[((i-1)*PROBE.SIZE+1):(i*PROBE.SIZE)] = d
            }
            converted = apply(cbind(sapply(diags, '[[','type'), ds.training$MalignancyCharacter),
                              1, diagnosisToOutcome)
            return(converted)
        })

        printDebug("finished ivfc classification")

        outcomes.ivfcs        = data.frame(ds.training[, 1:3], outcomes.ivfcs)
        names(outcomes.ivfcs) = c(names(ds.training)[1:3], IVFC.NAME)

        training.stats.ivfcs = melt(calculate.stats(
            aggregate.outcomes(outcomes.ivfcs)
        ),
        id.vars = "ObscureLevel" ,
        variable.name = "Method",
        value.name = "Value") %>%
            rename(Measure=L1)

        training.stats.ivfcs = suppressWarnings( # suppress different factor levels warning
            left_join(training.stats.ivfcs,
                      IVFC.BINDED.DESCRIPTION,
                      by="Method")
        )
    }
}

# ---- training-statistics-bind ----
if(!SKIP.TRAINING) {
    printDebug("training statistics bind")
    training.stats.all = suppressWarnings( # suppress different factor levels warning
        bind_rows( if(!SKIP.KNN) training.stats.knns else data.frame(),
                   if(!SKIP.IVFC) training.stats.ivfcs else data.frame()
        )
    )
}

# ---- training-statistics-performance-calculation ----
if(!SKIP.TRAINING) {
    printDebug("training statistics performance calculation")

    training.stats.all.perf = aggregate(training.stats.all$Value,
                                        list(Method=training.stats.all$Method,
                                             Measure=training.stats.all$Measure),
                                        mean) %>%
        rename(Value=x)

    training.stats.all.perf = suppressWarnings( # suppress different factor levels warning
        left_join(training.stats.all.perf,
                  distinct(select(training.stats.all,
                                  Method, Class, Subclass, Subsubclass)),
                  by="Method")
    )
}


# ---- convert-statistics-performance-to-wide-format ----

printDebug("convert statistics performance to wide format")

training.stats.all.perf.wide = dcast(training.stats.all.perf,
                                    Method + Class + Subclass + Subsubclass ~ Measure,
                                    value.var="Value")

# ---- parallel-shutdown ----

printDebug("parallel shutdown")

if (THREADS > 1)
    stopCluster(CL)

# ---- save-evaluation ----

printDebug("save evaluation")

save.image(EVALUATION.OUTPUT.LOCATION)
