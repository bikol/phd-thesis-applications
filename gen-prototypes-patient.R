source('config.R')
source('utils.R')

library('ggplot2')
library('sampling')

# ---- db-read ----

# list of required columns
COLS.ALL = c("OvarianCancerInFamily", "HormoneReplacementTherapy", "Age",
             "ADimension", "PainAtExamination", "Ascites", "PapBloodFlow",
             "Solid", "ASolidDimension", "InnerWall", "Shadow", "Color", "Septum",
             "SmEchogenicity", "Location", "SmInnerWallThickness", "TumorVolume",
             "Pap", "APapDimension", "SeptumThickness", "AgeAfterMenopause",
             "Ca125", "Ri", "UterusRemoved", "IotaQuality")
COLS.SURE = c("HormoneReplacementTherapy", "Age", "PainAtExamination", "AgeAfterMenopause",
              "UterusRemoved")
COLS.OBSC = COLS.ALL[!COLS.ALL %in% COLS.SURE]


db = read.csv(DATABASE.LOCATION, header=TRUE)
db = db %>%
    # select(GivenId, MalignancyCharacter, one_of(COLS.ALL)) %>%
    rename(PatientId=GivenId)

# ---- db-preprocess ----

db[with(db, which(MalignancyCharacter == 2)), "MalignancyCharacter"] = 1

db[with(db, which(is.na(PapBloodFlow)    & Pap == 0)),    "PapBloodFlow"]    = 0
db[with(db, which(is.na(APapDimension)   & Pap == 0)),    "APapDimension"]   = 0
db[with(db, which(is.na(SeptumThickness) & Septum == 0)), "SeptumThickness"] = 0
db[with(db, which(grepl("^Sms", PatientId) & is.na(Ri))), "Ri"] = 1
db[with(db, which(grepl("^Sz",  PatientId) & is.na(Ri) & Color == 1)), "Ri"] = 1

db = filter(db, !is.na(MalignancyCharacter))

db = normalize.data(db)


adnexes = apply(db, 1, function(row){
    if(row[['MalignancyCharacter']]==0){
        return(0)
    }
    if(!is.na(row[['Histopath']])){
        if(as.integer(row[['Histopath']])==30){
            return(3)
        }
        if(as.integer(row[['Histopath']])==31){
            return(4)
        }
    }

    if(as.integer(row[['Figo']]) %in% c(1,2,3)){
        return(1)
    }
    return(2)

})
mnps = apply(db, 1, function(row){
    if(is.na(row[['MenopauseAge']])){
        return(0)
    }else{
        return(1)
    }
})

db = mutate(db, Adnex = adnexes, Mnp=mnps)

str = strata(db, stratanames=c('Adnex', 'Mnp'), size=rep(100,10), method='srswr')
db.str = db[str$ID_unit,]
# db.str=db



# --------

for( col in COLS.ALL){
# col = 'ADimension'
    print(col)
    result.gen  = matrix(nrow=2, ncol=length(unique(disc)), dimnames=list(NULL, 0:4))
    result.pre  = matrix(nrow=2, ncol=length(unique(disc)), dimnames=list(NULL, 0:4))
    result.post = matrix(nrow=2, ncol=length(unique(disc)), dimnames=list(NULL, 0:4))
    for(i in 0:4){
        data = db[disc==i,]
        gen.l = min.max(mean(data[[col]], na.rm=T) - sd(data[[col]], na.rm=T))
        gen.u = min.max(mean(data[[col]], na.rm=T) + sd(data[[col]], na.rm=T))
        result.gen[,i+1]=c(gen.l,gen.u)
        cat(i, round(mean(data[[col]], na.rm=T),2), round(sd(data[[col]], na.rm=T),2),round(gen.l,2), '-',round(gen.u, 2),sep="\t", fill=T)

        data.pre = db[disc==i&is.na(db$MenopauseAge),]
        pre.l = min.max(mean(data.pre[[col]], na.rm=T) - sd(data.pre[[col]], na.rm=T))
        pre.u = min.max(mean(data.pre[[col]], na.rm=T) + sd(data.pre[[col]], na.rm=T))
        result.pre[,i+1]=c(pre.l,pre.u)
        cat(i, round(mean(data.pre[[col]], na.rm=T),2), round(sd(data.pre[[col]], na.rm=T),2), round(pre.l, 2), '-', round(pre.u, 2),sep="\t", fill=T)

        data.post = db[disc==i&!is.na(db$MenopauseAge),]
        post.l = min.max(mean(data.post[[col]], na.rm=T) - sd(data.post[[col]], na.rm=T))
        post.u = min.max(mean(data.post[[col]], na.rm=T) + sd(data.post[[col]], na.rm=T))
        result.post[,i+1]=c(post.l,post.u)
        cat(i, round(mean(data.post[[col]], na.rm=T),2), round(sd(data.post[[col]], na.rm=T),2), round(post.l, 2), '-', round(post.u, 2),sep="\t", fill=T)
    }
    print(table(db[[col]]))
    toBarplot = function(results){
        results[2,] = results[2,]-results[1,]
        return(results)
    }

    df = data.frame(Disc = db.str$Adnex, Value = db.str[[col]], Mnp = factor(db.str$Mnp))

    print( qplot(data=df, x=Disc, y=Value, size=3,  main=col, color=Mnp) + geom_jitter(size=2, position=position_jitter(width=0.3, height=0.05))
           + geom_errorbar(size=0.75, width=0.5, data=data.frame(t(result.gen), Disc=0:4, Mnp=factor(0), Value=(result.gen[1,]+result.gen[2,])/2), mapping=aes(x=Disc-0.25, ymin=c(X1), ymax=c(X2)), color="blue")
           + geom_errorbar(size=0.75, width=0.5, data=data.frame(t(result.pre), Disc=0:4, Mnp=factor(0), Value=(result.pre[1,]+result.pre[2,])/2), mapping=aes(x=Disc, ymin=c(X1), ymax=c(X2)), color="red", shape=16)
           + geom_errorbar(size=0.75, width=0.5, data=data.frame(t(result.post), Disc=0:4, Mnp=factor(0), Value=(result.post[1,]+result.post[2,])/2), mapping=aes(x=Disc+0.25, ymin=c(X1), ymax=c(X2)), color="cyan", shape=16)
    )


}