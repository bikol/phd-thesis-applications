df = transmute(train.data,
               I1 = train.data[, 1]+train.data[,2],
               I2 = train.data[, 3]+train.data[,4],
               I3 = train.data[, 5]+train.data[,6],
               I4 = train.data[, 7]+train.data[,8],
               I5 = train.data[, 9]+train.data[,10],
               I6 = train.data[, 11]+train.data[,12],
               I7 = train.data[, 13]+train.data[,14],
               I8 = train.data[, 15]+train.data[,16],
               I9 = train.data[, 17]+train.data[,18],
               I10 = train.data[, 19]+train.data[,20],
               I11 = train.data[, 21]+train.data[,22],
               I12 = train.data[, 23]+train.data[,24],
               Class1 = Class1,
               Class2 = Class2,
               Class3 = Class3,
               Class4 = Class4
               )

dff = bind_rows(df,
                filter(df, !is.na(Class2)) %>% mutate(Class1=Class2),
                filter(df, !is.na(Class3)) %>% mutate(Class1=Class3))

for(i in 1:12){
    print(paste0("I",i))
    tmp = filter(dff, Class1==0)[ , i]
    print(paste0('Mean 0: ', round(mean(tmp[[1]])/2, 2)))

    tmp = filter(dff, Class1==1)[ , i]
    print(paste0('Mean 1: ', round(mean(tmp[[1]])/2, 2)))

    tmp = filter(dff, Class1==2)[ , i]
    print(paste0('Mean 2: ', round(mean(tmp[[1]])/2, 2)))

    tmp = filter(dff, Class1==4)[ , i]
    print(paste0('Mean 4: ', round(mean(tmp[[1]])/2, 2)))
}