source('config.R')
source('similarities-helpers.R')

library(parallel)
library(Matrix)

ERR = 1e-2
# sim.params = SIM.PARAMS.EXACT.MINIMUM.ID
# sim.params = SIM.PARAMS.EXACT.PRODUCT.ID
# sim.params = SIM.PARAMS.EXACT.SW.m0.5.ID
# sim.params = SIM.PARAMS.EXACT.LUK.ID
# sim.params = SIM.PARAMS.EXACT.SS.2.ID
# sim.params = SIM.PARAMS.EXACT.SS.5.ID
# sim.params = SIM.PARAMS.EXACT.SS.100.ID

ts = apply(ds.test, 1 , function(x){
    return(list(m=matrix(as.numeric(x[5:(5+length(METHODS)*2-1)]), nrow=2), type=x[4]))
})

if (THREADS > 1)
{
    CL = makeCluster(THREADS, outfile="")

    clusterExport(CL, c("SEED", "OBSUCRE.REPEAT"))
    clusterCall(CL, function(){ set.seed(SEED) })

    clusterExport(cl=CL, list('ERR', 'SIM.PARAMS', 'myConstrOptim', 'ts', 'sim.params', 'sim.optim', 'sim', 'invPerm', 'PROBE.SIZE', 'SIMILARITIES', 'SIMILARITIES.NAME', 'AGGREGATORS', 'AGGREGATORS.NAME', 'METHODS',
                              'METHODS.NAME', 'ds.training', 'ds.test',
                              'diagnosisToOutcome', 'CUTOFF.CRISP', 'W.AUC',
                              'COMMON.PART', 'INTERVAL.INTERSECTION'))

    usedLapply = function(...){ parLapply(CL, ...) }
} else {
    usedLapply = lapply
}

for( iii in 1:length(SIM.PARAMS)){

    sim.params = SIM.PARAMS[[iii]]
    print(paste("sim.params",iii))
#     sim.params = SIM.PARAMS.EXACT.PRODUCT.ID

    sim.optim = GEN.SIM.JACCARD(sim.params, optimisation=T)
    sim = GEN.SIM.JACCARD(sim.params, optimisation=F)

    usedLapply(sample(175, 15), function(i){
        # for(i in sample(175)){
        for(j in sample(175, 15)){
            #         i=88;j=8
#             print(paste(i, j))
            #         print(ts[[i]]$m)
            #         print(ts[[j]]$m)
            t.o = system.time( (result.optim = sim.optim(
                list(lower=ts[[i]]$m[1,],
                     upper=ts[[i]]$m[2,]),
                list(lower=ts[[j]]$m[1,], upper=ts[[j]]$m[2,]))))
            #         print(result.optim)
            #         print(paste(result.optim$ly, result.optim$uy))
            #         print("GEKMA")
            t.k = system.time( (result = sim(
                list(lower=ts[[i]]$m[1,],
                     upper=ts[[i]]$m[2,]),
                list(lower=ts[[j]]$m[1,], upper=ts[[j]]$m[2,]))))
#             print(paste(round(sum(t.o)/sum(t.k),2),"x faster"))
            #         print(paste(result$ly, result$uy))
            if(result$ly-result.optim$ly > ERR | result.optim$uy-result$uy > ERR){
                print(paste(i, j))
                print(ts[[i]]$m)
                print(ts[[j]]$m)
                print(paste(result.optim$ly, result.optim$uy))
                print("GEKMA")
                print(paste(result$ly, result$uy))
                warning(paste("Err", i, j))
            }
        }
    })
}
if (THREADS > 1)
    stopCluster(CL)


# classifier = (GEN.KNN.SIM(sim.params, 3, optimisation=T))(ts)
# result = classifier(ts[[1]]$m)
