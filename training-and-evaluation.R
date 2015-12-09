# ---- init ----

rm(list=ls())

source("config.R")
source("stats.R")
source("methods.R")
source("aggregators.R")
source("aggregators-optimize.R")
source("similarities.R")
source("mcn-test.R")
source("utils.R")

library(parallel)
library(reshape2)
library(dplyr)
library(Matrix)


# ---- read-datasets ----

printDebug("read datasets")

colClass = c("factor", "numeric", "integer", "integer",
             rep("numeric", times=2*length(METHODS)))

ds.training = read.csv(TRAINING.LOCATION, colClasses=colClass)
ds.test     = read.csv(TEST.LOCATION,     colClasses=colClass)

# ---- parallel-init ----

printDebug("paralell init")

if (THREADS > 1)
{
    CL = makeCluster(THREADS, outfile="")

    clusterExport(CL, c("SEED", "OBSUCRE.REPEAT"))
    clusterCall(CL, function(){ set.seed(SEED) })

    clusterExport(cl=CL, list('KNN.MAX.CASE.BASE.SIZE', 'invPerm', 'PROBE.SIZE', 'SIMILARITIES', 'SIMILARITIES.NAME', 'AGGREGATORS', 'AGGREGATORS.NAME', 'METHODS',
                              'METHODS.NAME', 'ds.training', 'ds.test',
                              'diagnosisToOutcome', 'CUTOFF.CRISP', 'W.AUC',
                              'COMMON.PART', 'INTERVAL.INTERSECTION'))

    usedLapply = function(...){ parLapply(CL, ...) }
} else {
    usedLapply = lapply
}


# ---- training-statistics-similarities ----

printDebug("training statistics similarities")

outcomes.sims = usedLapply(1:length(SIMILARITIES), function(i){

    sim = SIMILARITIES[[i]]

    diags = c()

    # split training data set into separate probes
    for(i in 1:(nrow(ds.training)/PROBE.SIZE)) {
        printDebug(paste("Data repeat(probe):", i))

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
                return(list(m=matrix(as.numeric(x[5:(5+length(METHODS)*2-1)]), nrow=2), type=x[4]))
                })

            if(length(ts) > KNN.MAX.CASE.BASE.SIZE) {
                classifier = sim(ts[sample(length(ts), min(length(ts), KNN.MAX.CASE.BASE.SIZE))])
            } else {
                classifier = sim(ts)
            }

            # selection of appropriate columns
            tmp = apply(test.data[, 5:(5+length(METHODS)*2-1)], 1, function(row) {
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
    return(converted)
})

printDebug("finished classification")

outcomes.sims        = data.frame(ds.training[, 1:3], outcomes.sims)
names(outcomes.sims) = c(names(ds.training)[1:3], SIMILARITIES.NAME)

training.stats.sims = melt(calculate.stats(
    aggregate.outcomes(outcomes.sims)
                        ),
                        id.vars = "ObscureLevel" ,
                        variable.name = "Method",
                        value.name = "Value") %>%
                    rename(Measure=L1)

training.stats.sims = suppressWarnings( # suppress different factor levels warning
        left_join(training.stats.sims,
                SIMILARITIES.BINDED.DESCRIPTION,
                by="Method")
        )


# ---- training-statistics-bind ----

printDebug("training statistics bind")

training.stats.all = training.stats.sims

# ---- training-statistics-performance-calculation ----

printDebug("training statistics performance calculation")

training.stats.all.perf = aggregate(training.stats.all$Value,
                                    list(Method=training.stats.all$Method,
                                         Measure=training.stats.all$Measure),
                                    mean) %>%
                          rename(Value=x)

training.stats.all.perf = left_join(training.stats.all.perf,
                                    distinct(select(training.stats.all,
                                              Method, Class, Subclass, Subsubclass)),
                                    by="Method")

# ---- select-optimized-similarities ----

printDebug("select optimized similarities")

# currently skip optimisation
# optimizedSimilaritiesNames = getOptimizedAggregators(training.stats.all.perf, PERFORMANCE.MEASURE)
optimizedSimilaritiesNames = SIMILARITIES.NAME

training.stats.all.perf = subset(training.stats.all.perf, Method %in% c(optimizedSimilaritiesNames))


############################
############################


# ---- test-combine-obscuration-levels ----

printDebug("test combine obscuration levels")

ds.test$ObscureLevel = 0

# ---- test-statistics-similarities ----

printDebug("test statistics similarities")

similarities.from.training = unique(subset(training.stats.all.perf, Class=="Similarity")$Method)

if (THREADS > 1)
    clusterExport(CL, c("similarities.from.training"))

outcomes.sims = usedLapply(1:length(similarities.from.training), function(j){
    i = which(similarities.from.training[j] == SIMILARITIES.NAME)
    sim = SIMILARITIES[[i]]

    diags = c()
    TEST.SIZE = nrow(ds.test)

    # split training data set into separate probes
    for(i in 1:(nrow(ds.training)/PROBE.SIZE)) {
        printDebug(paste("Data repeat(probe):", i))

        # select only cases from this probe
        train.data = ds.training[((i-1)*PROBE.SIZE+1):(i*PROBE.SIZE), ]

        # convert training set into proper format accepted by classifier
        ts = apply(train.data, 1, function(x){
            return(list(m=matrix(as.numeric(x[5:(5+length(METHODS)*2-1)]), nrow=2), type=x[4]))
        })

        if(length(ts) > KNN.MAX.CASE.BASE.SIZE) {
            classifier = sim(ts[sample(length(ts), min(length(ts), KNN.MAX.CASE.BASE.SIZE))])
        } else {
            classifier = sim(ts)
        }

        d = apply(ds.test[, 5:(5+length(METHODS)*2-1)], 1, function(row) {
            # matrix in format required by aggregation method is created and passed into aggr
            return(classifier(matrix(row, nrow=2)))
        })
        diags[((i-1)*TEST.SIZE+1):(i*TEST.SIZE)] = d
    }
    converted = apply(cbind(diags, ds.test$MalignancyCharacter),
                      1, diagnosisToOutcome)
    return(converted)
})

printDebug("finished classification")

obscureLevels = unique(ds.training[, 2:3])
obscureLevels = obscureLevels[rep(1:nrow(obscureLevels),each=nrow(ds.test)), ]
outcomes.sims        = data.frame(ds.test[, 1:1], obscureLevels, outcomes.sims)
names(outcomes.sims) = c(names(ds.test)[1:3], similarities.from.training)

test.stats.sims = melt(calculate.stats(
    aggregate.outcomes(outcomes.sims)
                ),
                id.vars = "ObscureLevel" ,
                variable.name = "Method",
                value.name = "Value") %>%
        rename(Measure=L1) #%>%
        #select(-ObscureLevel)

test.stats.sims = suppressWarnings( # suppress different factor levels warning
    left_join(test.stats.sims,
              SIMILARITIES.BINDED.DESCRIPTION,
              by="Method")
)

# ---- test-statistics-bind ----

printDebug("test statistics bind")

test.stats.all = test.stats.sims

# ---- training-statistics-performance-calculation ----

printDebug("training statistics performance calculation")

test.stats.all.perf = aggregate(test.stats.all$Value,
                                list(Method=test.stats.all$Method,
                                     Measure=test.stats.all$Measure),
                                mean) %>%
                            rename(Value=x)

test.stats.all.perf = left_join(test.stats.all.perf,
                                    distinct(select(test.stats.all,
                                                    Method, Class, Subclass, Subsubclass)),
                                    by="Method")

# ---- test-statistics-performance-bind-with-training ----

printDebug("test statistics performance bind with training")

binded.stats.all.perf = left_join(select(training.stats.all.perf, Method, Measure, Value),
                                  test.stats.all,
                                  by=c("Method", "Measure")) %>%
                        rename(Value.training=Value.x, Value.test=Value.y)

# ---- convert-statistics-performance-to-wide-format ----

printDebug("convert statistics performance to wide format")

training.stats.all.perf.wide = dcast(training.stats.all.perf,
                                     Method + Class + Subclass + Subsubclass ~ Measure,
                                     value.var="Value")

test.stats.all.wide          = dcast(test.stats.all.perf,
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

# ---- save-evaluation ----

printDebug("save evaluation")

save.image(EVALUATION.OUTPUT.LOCATION)
