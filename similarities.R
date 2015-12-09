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
#             str(training.set)
            sims = lapply(training.set, function(nb){
#                 str(nb)
                return(SIM(list(lower=nb$m[1,], upper=nb$m[2,]), list(lower=m[1,], upper=m[2,])))
            })

#             str(sims)
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
        list(function(ts){return(function(x){return(sample(2,1)-1)})}, "dummy.knn1", "knn", "Interval"),
#         list(function(x){return(0)}, "dummy-knn2", "knn", "Interval"),
        list(function(ts){return(function(x){round((x[1,1]+x[2,1])/2)})}, "dummy.knn3", "knn", "Interval"),
        list(GEN.KNN.SIM(SIM.PARAMS.EXACT.MINIMUM.ID, 3, optimisation=F), "min.id.ex","knn","Interval"),
        list(GEN.KNN.SIM(SIM.PARAMS.EXACT.PRODUCT.ID, 3, optimisation=F), "prod.id.ex","knn","Interval"),
        list(GEN.KNN.SIM(SIM.PARAMS.EXACT.SW.m0.5.ID, 3, optimisation=F), "sw.id.ex","knn","Interval")
    )


# at least 2 aggrs must be defined
# name and class must not contain '-' and '=' signs (must be valid data.frame column name)
SIMILARITIES.LIST = c(KNN.BASIC
                     )

SIMILARITIES.LIST = sample(SIMILARITIES.LIST, length(SIMILARITIES.LIST))
SIMILARITIES = sapply(SIMILARITIES.LIST,'[[',1)
SIMILARITIES.NAME = sapply(SIMILARITIES.LIST,'[[',2)
SIMILARITIES.CLASS = sapply(SIMILARITIES.LIST,'[[',3)
SIMILARITIES.SUBCLASS = sapply(SIMILARITIES.LIST,'[[',4)

SIMILARITIES.BINDED.DESCRIPTION = data.frame(Method=SIMILARITIES.NAME,
                                            Class="Similarity",
                                            Subclass=SIMILARITIES.CLASS,
                                            Subsubclass=SIMILARITIES.SUBCLASS)