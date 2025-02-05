library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(multcompView)
library(cowplot)

#Read in 42 pot experimental data
df.plant <- read.csv("R:/GaeumannomycesGenomics/00_experimental_data/pot_experiment_plantlevel.csv")
df.pot <- read.csv("R:/GaeumannomycesGenomics/00_experimental_data/pot_experiment_potlevel.csv")

#Format data
df.plant.sum <- df.plant %>%
  gather(test, value, -Treatment, -Pot, -Bed, -Row, -Col, -Plant)

df.pot.sum <- df.pot %>%
  gather(test, value, -Treatment, -Pot, -Bed, -Row, -Col) %>%
  filter(test!="Dried.Root.Biomass.per.g")

#Recode treatments with strain codes
df.sum <- bind_rows(df.plant.sum, df.pot.sum) %>%
  mutate(
    Treatment.name=factor(recode(Treatment,
                                 A="Control",
                                 B="Gt-19d1",
                                 C="Gt-8d",
                                 E="Gt-4e",
                                 F="Gt-23d",
                                 G="Gt-LH10"),
                          levels=c("Control", "Gt-8d", "Gt-19d1", 
                                   "Gt-23d", "Gt-4e", "Gt-LH10")),
    type=recode(Treatment,
                A="Control",
                B="Type A",
                C="Type A",
                E="Type B",
                F="Type B",
                G="Type B")
  )

#Check for normality of residuals
for (i in unique(df.sum$test)) {
  
  plot(ggqqplot(residuals(lm(value ~ Treatment, data=df.sum[df.sum$test == i,]))) +
         ggtitle(i))
  
}

#Check for homogeneity of variance
df.levene <- df.sum %>%
  group_by(test) %>%
  levene_test(value ~ Treatment)

#Do multiple comparison testing
df.tukey <- df.sum %>%
  filter(test %in% df.levene$test[which(df.levene$p >= 0.05)]) %>%
  group_by(test) %>%
  tukey_hsd(value ~ Treatment)

df.gh <- df.sum %>%
  filter(test %in% df.levene$test[which(df.levene$p < 0.05)]) %>%
  group_by(test) %>%
  games_howell_test(value ~ Treatment)

#Combine results
df.significance <- bind_rows(df.tukey, df.gh)

#Add letters for significance groups
for (i in 1:length(unique(df.significance$test))) {
  
  tmp <- df.significance[df.significance$test == unique(df.significance$test)[i],]
  
  letters <- data.frame(
    test=unique(df.significance$test)[i],
    tukey=multcompLetters(setNames(tmp$p.adj,
                                   paste0(tmp$group1,
                                          "-", tmp$group2)))$Letters,
    Treatment=names(multcompLetters(setNames(tmp$p.adj,
                                             paste0(tmp$group1,
                                                    "-", tmp$group2)))$Letters)
  )
  
  if (length(unique(letters$tukey)) > 1) {
   
    assign(paste0("df.letters.tmp.", i), letters) 
    
  }
  
}

df.significance.letters <- bind_rows(mget(ls(pattern="df.letters.tmp.")))
df.significance.letters$Treatment.name <- 
  df.sum$Treatment.name[match(df.significance.letters$Treatment, df.sum$Treatment)]

#Format strip labels
test.names <- 
  c(Height.cm="Plant height (cm)",
    Ear.length.cm="Mean ear length (cm)",
    Flag.Leaf.Length.cm="Flag leaf length (cm)",
    Number.of.ears="Number of ears",
    Number.of.tillers="Number of tillers",
    Number.of.roots="Number of roots",
    Roots.per.tiller="Number of roots per tiller",
    Root.length.cm="Mean root length (cm)",
    Average.dried.root.biomass.per.plant.g="Mean dried root\nbiomass per plant (g)",
    TAI="Take-all index (TAI)")

#Plot box and violin plots for all characteristics
gg.pots <- ggplot(df.sum, aes(x=Treatment.name, y=value)) +
  facet_wrap(~ factor(
    test,
    levels=c(c("TAI", "Number.of.roots",
               "Roots.per.tiller", "Root.length.cm",
               "Average.dried.root.biomass.per.plant.g", "Height.cm", 
                "Ear.length.cm", "Flag.Leaf.Length.cm",
                "Number.of.ears", "Number.of.tillers"))
  ),
  scales="free_y",
  ncol=5,
  labeller=as_labeller(test.names)) +
  geom_violin(aes(fill=type),
              linewidth=0.3,
              show.legend=FALSE) +
  geom_boxplot(aes(fill=type),
               linewidth=0.3,
               outlier.size=0.7,
               width=0.1) +
  geom_text(data=df.significance.letters,
            aes(x=Treatment.name, y=Inf, label=tukey),
            vjust=1.5,
            size=2) +
  scale_y_continuous(limits=c(0,NA),
                     expand=expansion(mult=c(0, 0.2))) +
  scale_fill_manual(values=c("grey", "#E69F00", "#f5d896")) +
  theme(legend.position="top",
        legend.title=element_blank(),
        legend.key=element_blank(),
        legend.margin=margin(0, 0, 0, 0),
        legend.direction="horizontal",
        axis.title=element_blank(),
        axis.text=element_text(size=5),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        strip.text=element_text(size=5, face="bold"),
        plot.margin=margin(5.5, 5.5, 5.5, 30))

#Write to file
pdf(file=paste0("R:/GaeumannomycesGenomics/00_experimental_data/pot_experiment-", Sys.Date(), ".pdf"),
     height=3, width=7.5)
ggarrange(
  ggdraw(gg.pots) + 
    draw_label("Above-ground", size=8, angle=90, fontface="bold", x=0.03, y=0.3) + 
    draw_label("Below-ground", size=8, angle=90, fontface="bold", x=0.03, y=0.67),
  labels="c")
dev.off()
