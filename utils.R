diagnosisToOutcome = function(p){
    # p is two element vector, p[[1]] contains diagnosis from model, p[[2]] containts
    # reference (real) diagnosis
    if(is.na(p[[1]]) && p[[2]]==1){
        return("N1")
    }else if(is.na(p[[1]]) && p[[2]]==0){
        return("N0")
    }else if(p[[1]]==1 && p[[2]]==1){
        return("TP")
    }else if(p[[1]]==1 && p[[2]]==0){
        return("FP")
    }else if(p[[1]]==0 && p[[2]]==1){
        return("FN")
    }else if(p[[1]]==0 && p[[2]]==0){
        return("TN")
    }
}

multiclassToOutcome = function(p){
    # p is two element vector, p[[1]] contains class label returned by classifier,
    # p[2:length(p)] contain proper labels of given instance
    d = p[2:length(p)]
    d = d[!is.na(d)]
    if(is.na(p[[1]])){
        return("N0")
    }else if (p[[1]] %in% d) {
        return("TN")
    }else {
        return("FP")
    }
}

orig.models.outcomes = function(data){
    toReturn = usedLapply(1:length(METHODS),function(i){
        name = METHODS.NAME[[i]]
        column = (5+(i-1)*2)
        diags = apply(data[,column:(column+1)],1,
                      function(row){if(row[1]!=row[2]) NA else CUTOFF.CRISP(c(row[1],row[2]))})
        converted = apply(cbind(diags,data$MalignancyCharacter),1,diagnosisToOutcome)
        return(converted)
    })
    toReturn = data.frame(data[,1:3],toReturn)
    names(toReturn) = c(names(data)[1:3],METHODS.NAME)
    return(data.frame(toReturn))
}

unc.models.outcomes = function(data){
    toReturn = usedLapply(1:length(METHODS),function(i){
        name = METHODS.NAME[[i]]
        column = (5+(i-1)*2)
        diags = apply(data[,column:(column+1)],1,CUTOFF.CRISP)
        converted = apply(cbind(diags,data$MalignancyCharacter),1,diagnosisToOutcome)
        return(converted)
    })
    toReturn = data.frame(data[,1:3],toReturn)
    names(toReturn) = c(names(data)[1:3],METHODS.NAME)
    return(data.frame(toReturn))
}

count.outcomes = function(data){
    result = c(TP=sum(data=='TP',na.rm=TRUE),
               TN=sum(data=='TN',na.rm=TRUE),
               FN=sum(data=='FN',na.rm=TRUE),
               FP=sum(data=='FP',na.rm=TRUE),
               N0=sum(data=='N0',na.rm=TRUE),
               N1=sum(data=='N1',na.rm=TRUE)
    )
    return(result)
}
aggregate.outcomes = function(outcomes){
    aggregate(outcomes[,4:ncol(outcomes)],
                         by=list(ObscureLevel=outcomes$ObscureLevel,
                                 ObscureRepeat=outcomes$ObscureRepeat),
                         count.outcomes,
                         simplify = FALSE)
}

calculate.stats = function(aggregated){
    result = list()
    for(i in 1:length(STATS)){
        stat = STATS[[i]]; force(stat)
        name = STATS.NAME[[i]]
        stats = aggregate(aggregated[ , 3:ncol(aggregated)],
                          by = list(ObscureLevel=aggregated$ObscureLevel),
                          function(data){
                              outcomes = matrix(c(data, recursive=TRUE),
                                                ncol=6, byrow=TRUE,
                                                dimnames = list(NULL, c("TP","TN","FN","FP","N0","N1")))
                              return(stat(outcomes))
                          },
                          simplify = TRUE)
        result[[name]] = stats
    }
    return(result)
}

min.max = function(vect){
    return(pmin(1, pmax(0, vect)))
}

normalize.data = function(dataset){
    dataset = dataset %>%
        mutate(Age = min.max((Age - 25)/50),
               ADimension = min.max((ADimension-25)/175),
               ASolidDimension = min.max(ASolidDimension/75),
               Color = min.max((Color-1)/3),
               SmEchogenicity = min.max(SmEchogenicity/4),
               SmInnerWallThickness = min.max(SmInnerWallThickness/2),
               TumorVolume = min.max(TumorVolume/1000),
               APapDimension = min.max(APapDimension/20),
               # SeptumThickness = min.max(SeptumThickness/7),
               AgeAfterMenopause = min.max(AgeAfterMenopause/20),
               Ca125 = min.max(Ca125/1000),
               IotaQuality = min.max( (IotaQuality-1)/5)
        )
    return(dataset)
}

convert.to.interval.format = function(dataset){
    colNames= c()
    colNames[seq(1, 2*ncol(dataset), by=2)] = paste(names(dataset), 'min', sep='.')
    colNames[seq(2, 2*ncol(dataset), by=2)] = paste(names(dataset), 'max', sep='.')

    dataset = dataset[,rep(1:ncol(dataset), each=2)]

    for(i in seq(1, ncol(dataset), by=2)){
        tmp = dataset[, i]
        tmp[is.na(tmp)] = 0
        dataset[, i] = tmp

        tmp = dataset[, i+1]
        tmp[is.na(tmp)] = 1
        dataset[, i+1] = tmp
    }
    names(dataset) = colNames
    return(dataset)
}

# auxiliary functions

printDebug = function(msg)
{
    if (DEBUG)
        cat(paste(format(Sys.time(), "[%X]"), msg, "\n"))
}

read.keel = function(file){
    skip = 17
    text = gsub('\\?','[0,10]', readChar(file, 1e6))
    text = gsub('\\[|\\]|\\{|\\}', '', text)
    data = scan(text=text, sep = ",",
                what=list(numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), integer(), integer(), integer(), integer()),
                skip=skip, comment.char='@', fill=T, quiet = T)
    df = data.frame(data)
    names(df)[seq(1, 24, by=2)]=paste0('Indicator', 1:12,'.min')
    names(df)[seq(2, 24, by=2)]=paste0('Indicator', 1:12,'.max')
    names(df)[25:28]= paste0("Class",1:4)
    return(df)
}
