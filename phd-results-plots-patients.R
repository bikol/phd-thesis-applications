rm(list=ls())

source("config.R")
source("results-overview.R")

library(tikzDevice)
options(tikzPdftexWarnUTF = FALSE)


e.patients = new.env()
e.models = new.env()
e.aggr = new.env()
load("datasets/evaluation-output.RData", e.aggr)
load("datasets/patients-eval-output.RData", e.patients)
load("datasets/models-eval-output.RData", e.models)

patients.training.stats.all = convertClasses(get("training.stats.all", envir = e.patients))
patients.training.stats.all.perf = generatePerf(patients.training.stats.all)
patients.test.stats.all = convertClasses(get("test.stats.all", envir = e.patients))
patients.test.stats.all.perf = convertClasses(get("test.stats.all.perf", envir = e.patients))
patients.test.stats.all.perf.wide = convertClasses(get("test.stats.all.perf.wide", envir = e.patients))

# models.training.stats.all = convertClasses(get("training.stats.all", envir = e.models))
# models.training.stats.all.perf = generatePerf(models.training.stats.all)
# models.test.stats.all = convertClasses(get("test.stats.all", envir = e.models))
# models.test.stats.all.perf = convertClasses(get("test.stats.all.perf", envir = e.models))
# models.test.stats.all.perf.wide = convertClasses(get("test.stats.all.perf.wide", envir = e.models))

required.training.stats.all = classesForOEA(subset(get("training.stats.all", envir = e.aggr), (Class=="Model"|Method=="mean_(dec_(owa_min))_cen_0.025")))
required.training.stats.all.perf = classesForOEA(
    subset(get("training.stats.all.perf", envir = e.aggr),
            Class=="Model" | Method == "mean_(dec_(owa_min))_cen_0.025"
            )
    )
required.test.stats.all.perf = classesForOEA(
    subset(get("test.stats.all", envir = e.aggr),
        Class == "Model" | Method == "mean_(dec_(owa_min))_cen_0.025"
        )
    )

# build non perf version of test results
required.test.stats.all = required.test.stats.all.perf[
    rep(1:nrow(required.test.stats.all.perf), each=length(OBSCURE.PERCENTAGES)), ]
required.test.stats.all$ObscureLevel = OBSCURE.PERCENTAGES


ref.colors =c("#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#36c7ab")

################
# tikz(file="4-training-patients-results.tex", sanitize = T, height = 5.5, width = 6,
#      documentDeclaration = '\\documentclass[11pt,a4paper,oldfontcommands]{memoir}')
# agg.sub = subset(patients.training.stats.all.perf,
#                  Measure==PERFORMANCE.MEASURE & (Class=='Similarity')) %>%
#     group_by(Class, Subsubclass)
#
# bestAggregators = NA
# if (PERFORMANCE.MEASURE.DESC) {
#     bestAggregators = top_n(agg.sub, 5, Value)$Method
# } else {
#     bestAggregators = top_n(agg.sub, 5, -Value)$Method
# }
#
# plotStatsTrainingSet2(subset(bind_rows(patients.training.stats.all, required.training.stats.all),
#                             Method %in% bestAggregators | Class=='Model'),
#                      PERFORMANCE.MEASURE)
# dev.off()

###############
# test cost detailed
###############
agg.sub = subset(patients.training.stats.all.perf,
                 Measure==PERFORMANCE.MEASURE & (Class=='Similarity') & (Subclass=='J_SS_m5_id' | Subclass=='J_SS_m2_id' | Subclass=='J_min_id') ) %>%
    group_by(Class, Subsubclass)

bestAggregators = NA
if (PERFORMANCE.MEASURE.DESC) {
    bestAggregators = top_n(agg.sub, 3, Value)$Method
} else {
    bestAggregators = top_n(agg.sub, 3, -Value)$Method
}

stats = subset(bind_rows(patients.test.stats.all, required.test.stats.all),
               Method %in% bestAggregators | Class=='Model')

shapes.ids = c(21,22,23,24,25,20,32)

ymax.val = 1.1*max(subset(stats, Class=="Model" & Measure==PERFORMANCE.MEASURE)$Value)

stats$Method = gsub("orig. ", "", stats$Method)
stats$Method = gsub("unc. ",  "", stats$Method)
# browser()
stats$Method = factor(stats$Method, levels=c("Alc", "ivfc_(J_min_id)_mean_1_max", "ivfc_(J_SS_m5_id)_cen_max",
                    "ivfc_(J_SS_m5_id)_mean_1_max", "knn_5_(J_min_id)_(lat)_maj",
                    "knn_5_(J_SS_m2_id)_(cen)_maj", "knn_5_(J_SS_m5_id)_(cen)_maj",
                    "LR1", "LR2", "RMI", "SM", "Tim", "OEA"))

levels(stats$Method) = c("Alc", "ivfc_min_mean_1_max", "ivfc_SS_m5_cen_max",
                         "ivfc_SS_m5_mean_1_max", "knn_5_min_lat_maj",
                         "knn_5_SS_m2_cen_maj", "knn_5_m5_cen_maj",
                         "LR1", "LR2", "RMI", "SM", "Tim", "OEA")

a = ggplot(data=subset(stats, Class=="Model" & Subclass=="Original" & Measure==PERFORMANCE.MEASURE),
           aes(x=ObscureLevel, y=Value, group=Method, colour=Method, shape=Method, fill=Method)) +
    geom_line() +
    geom_point(size=3, color="black") +
    # ggtitle("modele referencyjne") +
    scale_shape_manual(values=shapes.ids) +
    scale_color_manual(values=ref.colors) +
    scale_fill_manual(values=ref.colors) +
    xlab("poziom niekompletności danych") +
    ylab(ifelse(PERFORMANCE.MEASURE=="Cost matrix", "całkowity koszt", PERFORMANCE.MEASURE)) +
    theme_bw() +
    theme(legend.position="bottom",
          plot.margin=unit(c(0, 0, 0, 0), "npc"),
          legend.margin=unit(0.05, "npc"),
          panel.margin = unit(0.07, "npc"),
          axis.text.x=element_text(size=9),
          strip.text.x=element_text(size=9),
          axis.title.x = element_text(vjust=-0.5),
          axis.title.y = element_text(vjust=1.5),
          title = element_text(vjust=1.2)) +
    guides(color = guide_legend(nrow = 3,byrow=T),fill = guide_legend(nrow = 3,byrow=T),shape = guide_legend(nrow = 3,byrow=T)) +
    coord_cartesian(xlim=LIM.X, ylim=c(50, ymax.val)) +
    scale_x_continuous(breaks=BREAKS.X)

b = ggplot(data=subset(stats, (Class=='Similarity') & Measure==PERFORMANCE.MEASURE),
           aes(x=ObscureLevel, y=Value, group=Method, colour=Method, shape=Method, fill=Method)) +
    geom_line() +
    geom_point(size=3, color="black") +
    # ggtitle("modele referencyjne") +
    scale_shape_manual(values=shapes.ids) +
    scale_color_manual(values=ref.colors) +
    scale_fill_manual(values=ref.colors) +
    xlab("poziom niekompletności danych") +
    ylab(ifelse(PERFORMANCE.MEASURE=="Cost matrix", "całkowity koszt", PERFORMANCE.MEASURE)) +
    theme_bw() +
    theme(legend.position="bottom",

          plot.margin=unit(c(0, 0, 0, 0), "npc"),
          legend.margin=unit(0.05, "npc"),
          panel.margin = unit(0.07, "npc"),
          axis.text.x=element_text(size=9),
          strip.text.x=element_text(size=9),
          axis.title.x = element_text(vjust=-0.5),
          axis.title.y = element_text(vjust=1.5),
          title = element_text(vjust=1.2)) +
    guides(color = guide_legend(nrow = 3,byrow=T,title=""),fill = guide_legend(nrow = 3,byrow=T,title=""),shape = guide_legend(nrow = 3,byrow=T,title="")) +
    coord_cartesian(xlim=LIM.X, ylim=c(50, ymax.val)) +
    scale_x_continuous(breaks=BREAKS.X)
tikz(file="4-test-patients-detailed-ref.tex", sanitize = T, height = 3, width = 3,
     documentDeclaration = '\\documentclass[11pt,a4paper,oldfontcommands]{memoir}')
print(a)
dev.off()
tikz(file="4-test-patients-detailed-sel.tex", sanitize = T, height = 3, width = 3,
     documentDeclaration = '\\documentclass[11pt,a4paper,oldfontcommands]{memoir}')
print(b)
dev.off()

########
# total cost all
########
stats = bind_rows(patients.test.stats.all.perf, required.test.stats.all.perf)

stats$Measure = factor(stats$Measure,
                       levels=c(PERFORMANCE.MEASURE))

orig.models.perf.measure = subset(stats, Measure==PERFORMANCE.MEASURE & Class=="Model" & Subclass=="Original")

orig.models.perf.measure$Method = gsub("orig. ", "", orig.models.perf.measure$Method)
orig.models.perf.measure$Method = factor(orig.models.perf.measure$Method,
                                         levels=rev(c("Alc", "LR1", "LR2", "RMI", "SM", "Tim", "OEA")))

a0 = ggplot(data=orig.models.perf.measure, aes(x=Method, weight=Value, fill=Method),
            environment=environment()) +
    geom_bar(color="black") +
    # scale_fill_manual(values=rev(ref.colors)) +
    scale_fill_grey(start=1, end=0.8) +
    xlab("modele referencyjne") +
    ylab(NULL) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          legend.position="none",
          plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_text(aes(x=Method,y=Value+5, label=round(Value, 1)),
              hjust=0, vjust=0.5, size=3) +
    geom_text(aes(x=Method, y=0.05*max(orig.models.perf.measure$Value), label=Method),
              hjust=0, vjust=0.5, size=4) +
    coord_cartesian() + coord_flip(ylim=c(0, 1.1*max(orig.models.perf.measure$Value)))

aggrs.perf.measure = subset(stats, Measure==PERFORMANCE.MEASURE & Class=="Similarity") %>%
    group_by(Class, Subclass, Subsubclass)

if (PERFORMANCE.MEASURE.DESC) {
    aggrs.perf.measure = slice(aggrs.perf.measure, which.max(Value))
} else {
    aggrs.perf.measure = slice(aggrs.perf.measure, which.min(Value)) }

c = ggplot(data=aggrs.perf.measure,
           aes(x=factor(paste(Subsubclass, Subclass),
                        levels=rev(c("IVFC J_min_id", "IVFC J_SS_m5_id", "IVFC J_SS_m2_id", "IVFC J_SS_m1_id",
                                 "IVFC J_SS_m05_id", "IVFC J_prod_id", "IVFC J_SS_05_id", "IVFC J_luk_id",
                                 "IVFC J_SS_2_id", "IVFC J_SS_5_id",
                                 "kNN J_min_id", "kNN J_SS_m5_id", "kNN J_SS_m2_id", "kNN J_SS_m1_id",
                                "kNN J_SS_m05_id", "kNN J_prod_id", "kNN J_SS_05_id", "kNN J_luk_id",
                                "kNN J_SS_2_id", "kNN J_SS_5_id"
                        ))),
               weight=Value,
               fill=factor(paste(Subsubclass, Subclass),
                           levels=rev(c("IVFC J_min_id", "IVFC J_SS_m5_id", "IVFC J_SS_m2_id", "IVFC J_SS_m1_id",
                                        "IVFC J_SS_m05_id", "IVFC J_prod_id", "IVFC J_SS_05_id", "IVFC J_luk_id",
                                        "IVFC J_SS_2_id", "IVFC J_SS_5_id",
                                        "kNN J_min_id", "kNN J_SS_m5_id", "kNN J_SS_m2_id", "kNN J_SS_m1_id",
                                        "kNN J_SS_m05_id", "kNN J_prod_id", "kNN J_SS_05_id", "kNN J_luk_id",
                                        "kNN J_SS_2_id", "kNN J_SS_5_id"
                           )))), environment=environment()) +
    geom_bar(width = 0.8, position = position_identity(width=0.8), color="black") +
    scale_fill_grey(start=1, end=0.8) +
    xlab("zaproponowane metody dla różnych t-norm") +
    ylab(ifelse(PERFORMANCE.MEASURE=="Cost matrix", "całkowity koszt", PERFORMANCE.MEASURE)) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          legend.position="none",
          plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_text(aes(x=paste(Subsubclass, Subclass), y=Value+5,
                  label=round(Value, 1)),
              hjust=0, vjust=0.5, size=3) +
    geom_text(aes(x=paste(Subsubclass, Subclass), y=0.05*max(orig.models.perf.measure$Value),
                  label=paste(Subsubclass, Subclass)),
                hjust=0, vjust=0.5, size=4) +
    coord_cartesian() +
    coord_flip(ylim=c(0, 1.1*max(orig.models.perf.measure$Value)))

tikz(file="4-test-patients-all.tex", sanitize = T, height = 5.5, width = 6,
     documentDeclaration = '\\documentclass[11pt,a4paper,oldfontcommands]{memoir}')
grid.arrange(a0, c, nrow=2, heights=unit(c(0.23, 0.77), "npc"))
dev.off()

######
# selected
######

selected.aggrs = subset(patients.test.stats.all.perf.wide,
                        Class=="Similarity" &
                            #                             Decisiveness>=0.95 &
                            #                             Decisiveness<1.0 &
                            Sensitivity>Specificity &
                            Sensitivity>=0.90 &
                            Specificity>0.8)

stats = selected.aggrs
secondaryPerformanceMeasures = c("Sensitivity", "Specificity", "Decisiveness", "Accuracy")

stats = melt(stats, id.vars=c("Method", "Class", "Subclass", "Subsubclass")) %>%
    rename(Measure=variable, Value=value)

stats = transform(stats, Subsubclass=as.character(Subsubclass))

PERF.COLORS = c("#D3DDE2","#B2DF8A","#33A02C","#1B6699")

stats$Measure = factor(stats$Measure,
                       levels=c(PERFORMANCE.MEASURE, "Accuracy",
                                "Sensitivity", "Specificity", "Decisiveness"))

stats.other.measures = subset(stats, Measure %in% secondaryPerformanceMeasures)
stats.other.measures$Measure = revalue(as.factor(stats.other.measures$Measure),
                                       sapply(secondaryPerformanceMeasures,
                                              paste0, # hack to add spacing
                                              paste(rep(" ", LEGEND.SPACING),
                                                    collapse="")))

d = ggplot(data=stats.other.measures, aes(x=Method, weight=Value, fill=Measure)) +
    geom_bar(position=position_dodge(width = 0.9), width=0.9) +
    geom_bar(position=position_dodge(width = 0.9), width=0.9, color="black", show_guide=FALSE) +
    scale_fill_manual(values=PERF.COLORS, name="miara skuteczności    ") +
    coord_cartesian(ylim=c(0.0, 1.1)) +
    xlab(NULL) +
    ylab("wartość") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position="top",
          plot.margin=unit(c(0,0,0.01,0), "npc"),
          strip.text.x=element_text(size=9)) +
    geom_text(aes(y=Value+0.01, label=round(Value, 2), ymax=1),
              position=position_dodge(width=0.9), vjust=0, size=2.8) +
    facet_grid(.~Subclass+Subsubclass, scales="free", space="free")

tikz(file="4-test-patients-selected.tex", sanitize = T, height = 3, width = 6.25,
     documentDeclaration = '\\documentclass[11pt,a4paper,oldfontcommands]{memoir}')
print(d)
dev.off()