source('similarities-helpers.R')


#########################
# Similarity based classifiers
#########################
# *m* passed to all aggregation functions as argument is matrix.

# str(m) gives following output:
# > num [1:2, 1:6] 0.676 0.781 0.639 0.833 0.688 ...
# The lower and upper bounds of i-th interval diagnosis can be obtained via:
# m[1,i] and m[2,i]
# Each classifier returns binary diangosis or *NA*.
#########################

GEN.KNN.SIM = function(sim.params, k, nbs.selector, vote.strategy, optimisation = F) {
    force(sim.params);force(k);force(nbs.selector);force(vote.strategy);force(optimisation)
    SIM = GEN.SIM.JACCARD(sim.params, optimisation=optimisation)
    return(function(training.set){
        force(training.set)
        return(function(m){
            force(m)
            sims = lapply(training.set, function(nb){
                return(SIM(list(lower=nb$m[1,], upper=nb$m[2,]), list(lower=m[1,], upper=m[2,])))
            })
            sel = nbs.selector(sims, k)
            nbsTypes = sapply(training.set[sel], "[[", 'type')
            return(vote.strategy(nbsTypes, sims[sel]))
        })
    })
}


GEN.IVFC.SIM = function(sim.params, classes, sim.aggr, summary.strategy, optimisation = F) {
    force(sim.params);force(classes);force(sim.aggr);force(summary.strategy);force(optimisation)

    SIM = GEN.SIM.JACCARD(sim.params, optimisation=optimisation)
    return(function(prototypes){
        force(prototypes)
        return(function(m){
            force(m)
            sims = lapply(prototypes, function(nb){
                return(SIM(list(lower=nb$m[1,], upper=nb$m[2,]), list(lower=m[1,], upper=m[2,])))
            })

            cls = sapply(prototypes, '[[','type')

            toReturn = sapply(classes, function(cl){
                return(sim.aggr(sims[which(cls == cl)]))
            })
            # columns represent classes, two rows lower and upper similarity bound
            colnames(toReturn) = as.character(classes)
            return(list(type = summary.strategy(toReturn), ivfc = toReturn))
        })
    })
}

KNN.JACCARD = apply(cbind(expand.grid(SIM.PARAMS, KS, NBS.SELECTORS, VOTE.STRATEGIES),
                           expand.grid(SIM.PARAMS.NAME, KS, NBS.SELECTORS.NAME, VOTE.STRATEGIES.NAME)),
                     1, function(row){
                         list(GEN.KNN.SIM(row[[1]],row[[2]], row[[3]], row[[4]]),
                              paste('knn_', row[[6]],'_(', row[[5]], ')_(', row[[7]], ')_', row[[8]], sep=''),
                              'Jaccard', 'Interval')
                     })

KNN.JACCARD.1 = apply(cbind(expand.grid(SIM.PARAMS, NBS.SELECTORS),
                           expand.grid(SIM.PARAMS.NAME, NBS.SELECTORS.NAME)),
                     1, function(row){
                         list(GEN.KNN.SIM(row[[1]], 1, row[[2]], VOTE.STRATEGY.ALL),
                              paste('knn_', 1,'_(', row[[3]], ')_(', row[[4]], ')_', 'all', sep=''),
                              'Jaccard', 'Interval')
                     })

KNN.JACCARD.2 = apply(cbind(expand.grid(SIM.PARAMS, NBS.SELECTORS, VOTE.STRATEGIES[2:3]),
                           expand.grid(SIM.PARAMS.NAME, NBS.SELECTORS.NAME, VOTE.STRATEGIES.NAME[2:3])),
                     1, function(row){
                         list(GEN.KNN.SIM(row[[1]], 2, row[[2]], row[[3]]),
                              paste('knn_', 2,'_(', row[[4]], ')_(', row[[5]], ')_', row[[6]], sep=''),
                              'Jaccard', 'Interval')
                     })



# at least 2 classifiers must be defined
# name and class must not contain '-' and '=' signs (must be valid data.frame column name)
KNN.LIST = c(KNN.JACCARD, KNN.JACCARD.1, KNN.JACCARD.2)

if(is.finite(CLASSIFIER.NUMBER.LIMIT)){
    KNN.LIST = KNN.LIST[sample(length(KNN.LIST), CLASSIFIER.NUMBER.LIMIT)]
}

KNN.LIST = sample(KNN.LIST, length(KNN.LIST))
KNN = sapply(KNN.LIST,'[[',1)
KNN.NAME = sapply(KNN.LIST,'[[',2)
KNN.CLASS = sapply(KNN.LIST,'[[',3)
KNN.SUBCLASS = sapply(KNN.LIST,'[[',4)

KNN.BINDED.DESCRIPTION = data.frame(Method=KNN.NAME,
                                            Class="kNN",
                                            Subclass=KNN.CLASS,
                                            Subsubclass=KNN.SUBCLASS)


IVFC.JACCARD = apply(cbind(expand.grid(SIM.PARAMS, INTERVAL.AGGRS, SUMMARIES),
                          expand.grid(SIM.PARAMS.NAME, INTERVAL.AGGRS.NAME, SUMMARIES.NAME)),
                    1, function(row){
                        list(GEN.IVFC.SIM(row[[1]], list(0, 1), row[[2]], row[[3]]),
                             paste('ivfc_(', row[[4]],')_', row[[5]], '_', row[[6]], sep=''),
                             'Jaccard', 'Interval')
                    })

IVFC.LIST = c(IVFC.JACCARD)

if(is.finite(CLASSIFIER.NUMBER.LIMIT)){
    IVFC.LIST = IVFC.LIST[sample(length(IVFC.LIST), CLASSIFIER.NUMBER.LIMIT)]
}

IVFC.LIST = sample(IVFC.LIST, length(IVFC.LIST))
IVFC = sapply(IVFC.LIST,'[[',1)
IVFC.NAME = sapply(IVFC.LIST,'[[',2)
IVFC.CLASS = sapply(IVFC.LIST,'[[',3)
IVFC.SUBCLASS = sapply(IVFC.LIST,'[[',4)

IVFC.BINDED.DESCRIPTION = data.frame(Method=IVFC.NAME,
                                    Class="IVFC",
                                    Subclass=IVFC.CLASS,
                                    Subsubclass=IVFC.SUBCLASS)


getOptimizedClassifiers = function(inputData, performanceMeasure) {

    df = subset(inputData, Measure==performanceMeasure & Method %in% KNN.NAME)
    if(nrow(df)>0){
        if(PERFORMANCE.MEASURE.DESC){
            KNN.OPT = arrange(df, desc(Value))[1:min(nrow(df), CLASSIFIER.OPTIMISATION.NUMBER), 1]
        } else {
            KNN.OPT = arrange(df, Value)[1:min(nrow(df), CLASSIFIER.OPTIMISATION.NUMBER), 1]
        }
    } else {
        KNN.OPT = c()
    }

    df = subset(inputData, Measure==performanceMeasure & Method %in% IVFC.NAME)
    if(nrow(df)>0){
        if(PERFORMANCE.MEASURE.DESC) {
            IVFC.OPT = arrange(df, desc(Value))[1:min(nrow(df), CLASSIFIER.OPTIMISATION.NUMBER), 1]
        } else {
            IVFC.OPT = arrange(df, Value)[1:min(nrow(df), CLASSIFIER.OPTIMISATION.NUMBER), 1]
        }
    } else {
        IVFC.OPT = c()
    }

    SIMILARITIES.OPT = c(KNN.OPT,
                         IVFC.OPT)

    return(SIMILARITIES.OPT)
}
