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
            ord = nbs.selector(sims, k)
            nbsTypes = sapply(training.set[ord], "[[", 'type')
            return(vote.strategy(nbsTypes, sims[ord]))
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

            result = data.frame(ly=sapply(sims, '[[','ly'),
                                uy=sapply(sims, '[[','uy'),
                                class=sapply(prototypes, '[[','type'))

            toReturn = sapply(classes, function(cl){
                toAggregate = filter(result, class == cl)[,1:2]
                # data is passed in the same format as 'm' matrix
                return(sim.aggr(t(toAggregate)))
            })
            # columns represent classes, two rows lower and upper similarity bound
            colnames(toReturn) = as.character(classes)

            return(list(type = summary.strategy(toReturn), ivfc = toReturn))
        })
    })
}

KNN.BASIC = list(
        list(function(ts){return(function(x){return(sample(2,1)-1)})}, "dummy.random", 'basic', "Interval"),
        list(function(ts){return(function(x){round((x[1,1]+x[2,1])/2)})}, "dummy.mean", 'basic', "Interval"),
        list(GEN.KNN.SIM(SIM.PARAMS.EXACT.MINIMUM.ID, 3, NBS.SELECTOR.ORDER.CEN, VOTE.STRATEGY.MAJRORITY, optimisation=F), "min.id.ex", 'basic', "Interval"),
        list(GEN.KNN.SIM(SIM.PARAMS.EXACT.PRODUCT.ID, 3, NBS.SELECTOR.ORDER.CEN, VOTE.STRATEGY.MAJRORITY, optimisation=F), "prod.id.ex", 'basic', "Interval"),
        list(GEN.KNN.SIM(SIM.PARAMS.EXACT.LUK.ID, 3, NBS.SELECTOR.ORDER.CEN, VOTE.STRATEGY.MAJRORITY, optimisation=F), "luk.id.ex", 'basic', "Interval")
    )


# at least 2 classifiers must be defined
# name and class must not contain '-' and '=' signs (must be valid data.frame column name)
KNN.LIST = c(KNN.BASIC
                     )

KNN.LIST = sample(KNN.LIST, length(KNN.LIST))
KNN = sapply(KNN.LIST,'[[',1)
KNN.NAME = sapply(KNN.LIST,'[[',2)
KNN.CLASS = sapply(KNN.LIST,'[[',3)
KNN.SUBCLASS = sapply(KNN.LIST,'[[',4)

KNN.BINDED.DESCRIPTION = data.frame(Method=KNN.NAME,
                                            Class="kNN",
                                            Subclass=KNN.CLASS,
                                            Subsubclass=KNN.SUBCLASS)


IVFC.BASIC = list(
    list(function(ts){return(function(x){return(list(type=sample(2,1)-1, ivfc=NA))})}, "iv.dummy.random", 'basic', "Interval"),
    list(function(ts){return(function(x){return(list(type=round((x[1,1]+x[2,1])/2), ivfc=NA))})}, "iv.dummy.mean", 'basic', "Interval"),
    list(GEN.IVFC.SIM(SIM.PARAMS.EXACT.MINIMUM.ID, list(0, 1), INTERVAL.AGGR.MEAN, SUMMARY.MAX.CEN, optimisation=F), "iv.min.id.ex", 'basic', "Interval")#,
#     list(GEN.KNN.SIM(SIM.PARAMS.EXACT.PRODUCT.ID, 3, optimisation=F), "prod.id.ex", 'basic', "Interval"),
#     list(GEN.KNN.SIM(SIM.PARAMS.EXACT.LUK.ID, 3, optimisation=F), "luk.id.ex", 'basic', "Interval")
)

IVFC.LIST = c(IVFC.BASIC
)

IVFC.LIST = sample(IVFC.LIST, length(IVFC.LIST))
IVFC = sapply(IVFC.LIST,'[[',1)
IVFC.NAME = sapply(IVFC.LIST,'[[',2)
IVFC.CLASS = sapply(IVFC.LIST,'[[',3)
IVFC.SUBCLASS = sapply(IVFC.LIST,'[[',4)

IVFC.BINDED.DESCRIPTION = data.frame(Method=IVFC.NAME,
                                    Class="IVFC",
                                    Subclass=IVFC.CLASS,
                                    Subsubclass=IVFC.SUBCLASS)
