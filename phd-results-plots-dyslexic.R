rm(list=ls())

source("results-overview.R")

library(tikzDevice)
options(tikzPdftexWarnUTF = FALSE)

env = new.env()
load("datasets/dyslexic-eval-output.RData", env)

filename = "dyslexic"

training.stats.all.perf = convertClasses(get("training.stats.all.perf", envir = env))
training.stats.all.perf.wide = convertClasses(get("training.stats.all.perf.wide", envir = env))

ref.colors =c("#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#36c7ab")


stats = transmute(training.stats.all.perf.wide, Method, Subclass, Value=`Accuracy (all)`) %>%
    group_by(Subclass) %>%
    top_n(1, jitter(Value))


stats = mutate(stats, Class="Metoda Najbliższego Sąsiada", Err=0)

stats = bind_rows(stats, list(Method="crisp", Subclass="crisp", Value=0.343, Class="GFS", Err=0.0),
                  list(Method="extended", Subclass="extended", Value=0.579, Class="GFS", Err=0.137))

stats$Subclass = factor(stats$Subclass, levels=c("J_min_id", "J_SS_m5_id", "J_SS_m2_id", "J_SS_m1_id",
                                                 "J_SS_m05_id", "J_prod_id", "J_SS_05_id", "J_luk_id",
                                                 "J_SS_2_id", "J_SS_5_id", "crisp", "extended"))
levels(stats$Subclass) = c("J_min", "J_SS_m5", "J_SS_m2", "J_SS_m1",
                           "J_SS_m05", "J_prod", "J_SS_05", "J_luk",
                           "J_SS_2", "J_SS_5", "ostra", "rozszerzona")

a = ggplot(data=stats,
           aes(x=Subclass, y=1-Value, group=Subclass, colour=Subclass, shape=Subclass, fill=Subclass)) +
    geom_errorbar(aes(ymin=1-Value, ymax=1-Value+Err),
                  width=0.9,                    # Width of the error bars
                  position=position_dodge(.9), color="black") +
    geom_bar(stat = "identity", color="black") +
    scale_fill_grey(start=1, end=0.8) +
    xlab("") +
    ylab("poziom błędu") +
    theme_bw() +
    theme(legend.position="none",
          plot.margin=unit(c(0.005, 0.005, 0.0, 0.01), "npc"),
          legend.margin=unit(0.05, "npc"),
          panel.margin = unit(0.07, "npc"),
          axis.text.x=element_text(size=9),
          strip.text.x=element_text(size=9),
          axis.title.x = element_text(vjust=-0.5),
          axis.title.y = element_text(vjust=1.5),
          title = element_text(vjust=1.2)) +
    coord_cartesian(ylim=c(0, 0.7)) +
    facet_grid(.~Class, scales="free", space="free")

tikz(file=paste0("4-",filename,"-all.tex"), sanitize = T, height = 3, width = 6.25,
     documentDeclaration = '\\documentclass[11pt,a4paper,oldfontcommands]{memoir}')
print(a)
dev.off()