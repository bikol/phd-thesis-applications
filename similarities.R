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


GEN.KNN.SIM = function(sim.params, k, optimisation = F) {
    force(sim.params);force(k)
    SIM = GEN.SIM.JACCARD(sim.params, optimisation=optimisation)
    return(function(training.set){
        force(training.set)
        return(function(m){
            force(m)
            sims = lapply(training.set, function(nb){
                return(SIM(list(lower=nb$m[1,], upper=nb$m[2,]), list(lower=m[1,], upper=m[2,])))
            })

            # napisać porownywanie przedziałów, narazie dziala i sortuje po pierwszej skladowej listy
            cmp = sapply(sims, function(x){return((x$ly+x$uy)/2)})

            ord = order(cmp, decreasing = T)
            sims.sorted =sims[ord]

            #    dodać opcje remisu dla parzystych k i wyniku NA
            nbsTypes = sapply((training.set[ord])[1:k], "[[", 'type')
            return(as.numeric(names(which.max(table(nbsTypes)))[[1]]))
#             return(list(
#                 type=names(which.max(table(nbsTypes)))[[1]],
#                 ivfc=NULL# wygenerować jakoś
#             ))
        })
    })
}


KNN.BASIC = list(
        list(function(ts){return(function(x){return(sample(2,1)-1)})}, "dummy.random", 'basic', "Interval"),
        list(function(ts){return(function(x){round((x[1,1]+x[2,1])/2)})}, "dummy.mean", 'basic', "Interval"),
        list(GEN.KNN.SIM(SIM.PARAMS.EXACT.MINIMUM.ID, 3, optimisation=F), "min.id.ex", 'basic', "Interval"),
        list(GEN.KNN.SIM(SIM.PARAMS.EXACT.PRODUCT.ID, 3, optimisation=F), "prod.id.ex", 'basic', "Interval"),
        list(GEN.KNN.SIM(SIM.PARAMS.EXACT.LUK.ID, 3, optimisation=F), "luk.id.ex", 'basic', "Interval")
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
    list(function(ts){return(function(x){return(sample(2,1)-1)})}, "dummy.random", 'basic', "Interval"),
    list(function(ts){return(function(x){round((x[1,1]+x[2,1])/2)})}, "dummy.mean", 'basic', "Interval"),
    list(GEN.KNN.SIM(SIM.PARAMS.EXACT.MINIMUM.ID, 3, optimisation=F), "min.id.ex", 'basic', "Interval"),
    list(GEN.KNN.SIM(SIM.PARAMS.EXACT.PRODUCT.ID, 3, optimisation=F), "prod.id.ex", 'basic', "Interval"),
    list(GEN.KNN.SIM(SIM.PARAMS.EXACT.LUK.ID, 3, optimisation=F), "luk.id.ex", 'basic', "Interval")
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
