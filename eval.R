# ---- parallel-init ----

printDebug("paralell init")

if (THREADS > 1)
{
    CL = makeCluster(THREADS, outfile="")

    clusterExport(CL, c("SEED", "OBSUCRE.REPEAT"))
    clusterCall(CL, function(){ set.seed(SEED) })

    clusterExport(cl=CL, list('topo_sort', 'graph_from_adjacency_matrix', 'DescIterBinSearch', 'AscIterBinSearch',
                              'KNN.MAX.CASE.BASE.SIZE', 'PROBE.SIZE',
                              'IVFC', 'IVFC.NAME', 'KNN', 'KNN.NAME',
                              'ds.training', 'ds.test', 'DATA.COLS.NUM',
                              'diagnosisToOutcome', 'PROTOTYPES', 'invPerm',
                              'printDebug', 'DEBUG'))
    # usedLapply = function(...){ parLapplyLB(CL, ...) }
    usedLapply = function(...){ clusterApplyLB(CL, ...) }
} else {
    usedLapply = lapply
}


# ---- training-statistics-knn ----
if(!SKIP.TRAINING) {
    if(!SKIP.KNN) {
        printDebug("training statistics knn")

        outcomes.knns = usedLapply(1:length(KNN), function(j){
            printDebug(paste0("Using ", KNN.NAME[[j]]," classifier"))
            startTime = as.integer(Sys.time())
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
                    ts = apply(train.data[, 4:(5+DATA.COLS.NUM*2-1)], 1 , function(x){
                        return(list(m=matrix(as.numeric(x[2:(2+DATA.COLS.NUM*2-1)]), nrow=2), type=as.integer(x[1])))
                    })

                    if(length(ts) > KNN.MAX.CASE.BASE.SIZE) {
                        classifier = KNN[[j]](ts[sample(length(ts), KNN.MAX.CASE.BASE.SIZE)])
                    } else {
                        classifier = KNN[[j]](ts)
                    }

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
            converted = apply(cbind(diags, ds.training$MalignancyCharacter),
                              1, diagnosisToOutcome)
            printDebug(paste0("Finished ", KNN.NAME[[j]]," classifier in ", as.integer(Sys.time()) - startTime, ' secs'))
            return(converted)
        })

        printDebug("finished knn classification")

        outcomes.knns        = data.frame(ds.training[, 1:3], outcomes.knns)
        names(outcomes.knns) = c(names(ds.training)[1:3], KNN.NAME)

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
}

# ---- training-statistics-ivfc ----
if(!SKIP.TRAINING) {
    if(!SKIP.IVFC) {
        printDebug("training statistics ivfc")

        outcomes.ivfcs = usedLapply(1:length(IVFC), function(j){
            printDebug(paste0("Using ",IVFC.NAME[[j]]," classifier"))
            startTime = as.integer(Sys.time())
            diags = c()
            # split training data set into separate probes
            for(i in 1:(nrow(ds.training)/PROBE.SIZE)) {
                # printDebug(paste("Data repeat(probe):", i))

                # select only cases from this probe
                all.data = ds.training[((i-1)*PROBE.SIZE+1):(i*PROBE.SIZE), ]

                classifier = IVFC[[j]](PROTOTYPES)

                d = apply(all.data[, 5:(5+DATA.COLS.NUM*2-1)], 1, function(row) {
                    # matrix in format required by aggregation method is created and passed into aggr
                    return(classifier(matrix(row, nrow=2)))
                })
                diags[((i-1)*PROBE.SIZE+1):(i*PROBE.SIZE)] = d
            }
            converted = apply(cbind(sapply(diags, '[[','type'), ds.training$MalignancyCharacter),
                              1, diagnosisToOutcome)
            printDebug(paste0("Finished ", IVFC.NAME[[j]]," classifier in ", as.integer(Sys.time()) - startTime, ' secs'))
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

# ---- select-optimized-classifiers ----
if(!SKIP.TRAINING) {
    printDebug("select optimized classifiers")

    optimizedNames = getOptimizedClassifiers(training.stats.all.perf, PERFORMANCE.MEASURE)
    training.stats.all.perf = subset(training.stats.all.perf, Method %in% optimizedNames)
}

############ TEST DATASET ##################

# ---- test-combine-obscuration-levels ----

printDebug("test combine obscuration levels")

ds.test$ObscureLevel = 0

# ---- test-knn-similarities ----
if(!SKIP.KNN) {
    printDebug("test knn similarities")

    if(!SKIP.TRAINING){
        evaluated.knns = unique(subset(training.stats.all.perf, Class=="kNN")$Method)
    } else {
        evaluated.knns = KNN.NAME
    }

    if (THREADS > 1)
        clusterExport(CL, c("evaluated.knns"))

    outcomes.knns = usedLapply(1:length(evaluated.knns), function(j){
        printDebug(paste0("Using ", evaluated.knns[j]," classifier"))
        startTime = as.integer(Sys.time())
        i = which(evaluated.knns[j] == KNN.NAME)
        sim = KNN[[i]]

        diags = c()
        TEST.SIZE = nrow(ds.test)

        # split training data set into separate probes
        for(i in 1:(nrow(ds.training)/PROBE.SIZE)) {
            # printDebug(paste("Data repeat(probe):", i))

            # select only cases from this probe
            train.data = ds.training[((i-1)*PROBE.SIZE+1):(i*PROBE.SIZE), ]

            # convert training set into proper format accepted by classifier
            ts = apply(train.data[,4:(5+DATA.COLS.NUM*2-1)], 1, function(x){
                return(list(m=matrix(as.numeric(x[2:(2+DATA.COLS.NUM*2-1)]), nrow=2), type=as.integer(x[1])))
            })

            if(length(ts) > KNN.MAX.CASE.BASE.SIZE) {
                classifier = sim(ts[sample(length(ts), KNN.MAX.CASE.BASE.SIZE)])
            } else {
                classifier = sim(ts)
            }

            d = apply(ds.test[, 5:(5+DATA.COLS.NUM*2-1)], 1, function(row) {
                # matrix in format required by aggregation method is created and passed into aggr
                return(classifier(matrix(row, nrow=2)))
            })
            diags[((i-1)*TEST.SIZE+1):(i*TEST.SIZE)] = d
        }
        converted = apply(cbind(diags, ds.test$MalignancyCharacter),
                          1, diagnosisToOutcome)
        printDebug(paste0("Finished ", evaluated.knns[j]," classifier in ", as.integer(Sys.time()) - startTime, ' secs'))
        return(converted)
    })

    printDebug("finished kNN classification")

    obscureLevels = unique(ds.training[, 2:3])
    obscureLevels = obscureLevels[rep(1:nrow(obscureLevels), each=nrow(ds.test)), ]
    outcomes.knns        = data.frame(ds.test[, 1:1], obscureLevels, outcomes.knns)
    names(outcomes.knns) = c(names(ds.test)[1:3], evaluated.knns)

    test.stats.knns = melt(calculate.stats(
        aggregate.outcomes(outcomes.knns)
    ),
    id.vars = "ObscureLevel" ,
    variable.name = "Method",
    value.name = "Value") %>%
        rename(Measure=L1) #%>%

    test.stats.knns = suppressWarnings( # suppress different factor levels warning
        left_join(test.stats.knns,
                  KNN.BINDED.DESCRIPTION,
                  by="Method")
    )
}

# ---- test-ivfc-similarities ----
if(!SKIP.IVFC) {
    printDebug("test ivfc similarities")

    if(!SKIP.TRAINING){
        evaluated.ivfcs = unique(subset(training.stats.all.perf, Class=="IVFC")$Method)
    } else {
        evaluated.ivfcs = IVFC.NAME
    }

    if (THREADS > 1)
        clusterExport(CL, c("evaluated.ivfcs"))

    outcomes.ivfcs = usedLapply(1:length(evaluated.ivfcs), function(j){
        printDebug(paste0("Using ",evaluated.ivfcs[j]," classifier"))
        startTime = as.integer(Sys.time())
        i = which(evaluated.ivfcs[j] == IVFC.NAME)
        sim = IVFC[[i]]

        diags = c()
        TEST.SIZE = nrow(ds.test)

        # split training data set into separate probes
        for(i in 1:(nrow(ds.training)/PROBE.SIZE)) {
            # printDebug(paste("Data repeat(probe):", i))

            # select only cases from this probe
            train.data = ds.training[((i-1)*PROBE.SIZE+1):(i*PROBE.SIZE), ]

            # convert training set into proper format accepted by classifier
            ts = apply(train.data[,4:(5+DATA.COLS.NUM*2-1)], 1, function(x){
                return(list(m=matrix(as.numeric(x[2:(2+DATA.COLS.NUM*2-1)]), nrow=2), type=as.integer(x[1])))
            })
            classifier = sim(PROTOTYPES)

            d = apply(ds.test[, 5:(5+DATA.COLS.NUM*2-1)], 1, function(row) {
                # matrix in format required by aggregation method is created and passed into aggr
                return(classifier(matrix(row, nrow=2)))
            })
            diags[((i-1)*TEST.SIZE+1):(i*TEST.SIZE)] = d
        }
        converted = apply(cbind(sapply(diags, '[[','type'), ds.test$MalignancyCharacter),
                          1, diagnosisToOutcome)
        printDebug(paste0("Finished ", evaluated.ivfcs[j]," classifier in ", as.integer(Sys.time()) - startTime, ' secs'))
        return(converted)
    })

    printDebug("finished IVFC classification")

    obscureLevels = unique(ds.training[, 2:3])
    obscureLevels = obscureLevels[rep(1:nrow(obscureLevels), each=nrow(ds.test)), ]
    outcomes.ivfcs        = data.frame(ds.test[, 1:1], obscureLevels, outcomes.ivfcs)
    names(outcomes.ivfcs) = c(names(ds.test)[1:3], evaluated.ivfcs)

    test.stats.ivfcs = melt(calculate.stats(
        aggregate.outcomes(outcomes.ivfcs)
    ),
    id.vars = "ObscureLevel" ,
    variable.name = "Method",
    value.name = "Value") %>%
        rename(Measure=L1) #%>%

    test.stats.ivfcs = suppressWarnings( # suppress different factor levels warning
        left_join(test.stats.ivfcs,
                  IVFC.BINDED.DESCRIPTION,
                  by="Method")
    )
}

# ---- test-statistics-bind ----

printDebug("test statistics bind")

test.stats.all = suppressWarnings( # suppress different factor levels warning
    bind_rows( if(!SKIP.KNN) test.stats.knns else data.frame(),
               if(!SKIP.IVFC) test.stats.ivfcs else data.frame()
    )
)

# ---- test-statistics-performance-calculation ----

printDebug("test statistics performance calculation")

test.stats.all.perf = aggregate(test.stats.all$Value,
                                list(Method=test.stats.all$Method,
                                     Measure=test.stats.all$Measure),
                                mean) %>%
    rename(Value=x)

test.stats.all.perf = suppressWarnings( # suppress different factor levels warning
    left_join(test.stats.all.perf,
              distinct(select(test.stats.all,
                              Method, Class, Subclass, Subsubclass)),
              by="Method")
)

# ---- test-statistics-performance-bind-with-training ----
if(!SKIP.TRAINING) {
    printDebug("test statistics performance bind with training")

    binded.stats.all.perf = suppressWarnings( # suppress different factor levels warning
        left_join(select(training.stats.all.perf, Method, Measure, Value),
                  test.stats.all,
                  by=c("Method", "Measure")) %>%
            rename(Value.training=Value.x, Value.test=Value.y)
    )
}

# ---- convert-statistics-performance-to-wide-format ----

printDebug("convert statistics performance to wide format")

if(!SKIP.TRAINING) {
    training.stats.all.perf.wide = dcast(training.stats.all.perf,
                                         Method + Class + Subclass + Subsubclass ~ Measure,
                                         value.var="Value")
}

test.stats.all.perf.wide          = dcast(test.stats.all.perf,
                                          Method + Class + Subclass + Subsubclass ~ Measure,
                                          value.var="Value")

# # ---- aggregators-selection-and-statistical-tests ----
#
#   to be moved into result visualisation script
#
# printDebug("aggregators selection and statistical tests")
#
# selected.aggrs = subset(test.stats.all.wide,
#                         Class=="Aggregation" &
#                             Decisiveness>=0.95 &
#                             Decisiveness<1.0 &
#                             Sensitivity>Specificity &
#                             Sensitivity>=0.90 &
#                             Specificity>0.8)
#
# perf.selected.aggrs = subset(binded.stats.all.perf,
#                              Measure==PERFORMANCE.MEASURE &
#                                  Method %in% selected.aggrs$Method)
#
#
# perf.all = rbind(perf.selected.aggrs,
#                  subset(binded.stats.all.perf,
#                         Measure==PERFORMANCE.MEASURE & Class=="Model" & Subclass=="Uncertaintified"))
#
# outcomes.all = bind_cols(select(outcomes.models.unc, -(PatientId:ObscureRepeat)),
#                          select(outcomes.aggrs,      -(PatientId:ObscureRepeat)))
#
# pvals = sapply(subset(perf.all, Class=="Aggregation")$Method, function(a1) {
#     sapply(perf.all$Method, function(a2) {
#         mcn.test(outcomes.all[a1], outcomes.all[a2])
#     })
# })
#
# pvals[upper.tri(pvals, diag=TRUE)] = NA
#
# pvals = matrix(p.adjust(pvals, method = "BH",
#                         n=sum(!is.na(pvals)|is.nan(pvals))),
#                nrow=nrow(pvals),
#                dimnames=list(rownames(pvals), colnames(pvals)))

# ---- parallel-shutdown ----

printDebug("parallel shutdown")

if (THREADS > 1)
    stopCluster(CL)