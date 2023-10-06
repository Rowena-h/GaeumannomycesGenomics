library(tidyverse)
library(jsonlite)
library(ComplexUpset)
library(ggpubr)
library(tgutil)

#Read in metadata
strains <- read.csv("R:/GaeumannomycesGenomics/strains", sep="\t", header=FALSE)
metadata <- read.csv(paste0("R:/GaeumannomycesGenomics/05_phylogenomics/raxmlng/metadata.csv"))

#Read in BUSCO output files of missing BUSCOs
buscos <- do.call(rbind, lapply(
  Sys.glob("S://012_busco_asco/*/hifi_assm_*_busco_asco/run_ascomycota_odb10/missing_busco_list.tsv"),
  function(x) { 
    read.csv(x, header=FALSE, comment.char="#") %>%
      mutate(strain=sub("/hifi_assm.*", "", sub(".*012_busco_asco/", "", x))) }
))

#Format
buscos.df <- buscos %>%
  rename(busco=V1) %>%
  mutate(new.strain=strains$V4[match(strain, strains$V2)]) %>%
  filter(!is.na(new.strain))

#Create presence absence table of missing BUSCOs in each strain
presence.absence <- buscos.df %>%
  select(-strain) %>%
  mutate(value=1) %>%
  pivot_wider(names_from="new.strain", values_fill=0) %>%
  select(c("busco", "Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1",
           "Gt-19d1", "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10"))

#Plot upset plot
gg.upset <- 
  upset(
    as.data.frame(presence.absence),
    c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1", "Gt-19d1",
      "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10"),
    sort_sets=FALSE,
    matrix=(
      intersection_matrix(
        geom=geom_point(
          size=2
        ) 
      )
      # +
      #   scale_color_manual(
      #     values=c('TRUE'='black', 'FALSE'='white'),
      #     guide="none"
      #   )
    ),
    base_annotations=list(
      'Intersection size'=intersection_size(
        text=list(
          size=2
        )
      )
    ),
    stripes=upset_stripes(
      geom=geom_segment(linewidth=5),
      data=metadata %>%
        filter(own == "Y") %>%
        select(new.strain, clade) %>%
        rename(set=new.strain),
      mapping=aes(color=clade),
      colors=c(
        'Gh'='#F2F2F2',
        'Ga'='#E5E5E5',
        'GtB'='#E5E5E5',
        'GtA'='#F2F2F2'
      )
    ),
    themes=upset_modify_themes(
      list(
        'Intersection size'=theme(
          axis.title.y=element_text(margin=margin(r=-40))
        ),
        intersections_matrix=theme(
          legend.position="none",
          axis.title.x=element_blank()
        )
      )
    )
  )

#Read in BUSCO descriptor data
busco.names <- read.csv("R:/data/busco_database/v5/data/lineages/ascomycota_odb10/links_to_ODB10.txt", sep="\t")

#Filter for BUSCOs missing in all strains
busco.ids.all <- presence.absence %>%
  filter_at(vars(-busco), all_vars(. == 1)) %>%
  mutate(name=sapply(busco.names$Proteasome.subunit.alpha.type[match(busco, busco.names$X262829at4890)],
                     function(x) paste(strwrap(x, 35), collapse = "\n"))) %>%
  select(busco, name)

#Filter for BUSCOs missing in Gh
busco.ids.Gh <- presence.absence %>%
  filter(`Gh-1B17` == 1 & `Gh-2C17` == 1) %>%
  filter_at(vars(-`Gh-1B17`, -`Gh-2C17`, -busco), all_vars(. == 0)) %>%
  mutate(name=sapply(busco.names$Proteasome.subunit.alpha.type[match(busco, busco.names$X262829at4890)],
                     function(x) paste(strwrap(x, 40), collapse = "\n"))) %>%
  select(busco, name)

#Plot tables of missing BUSCOs and descriptors
busco.table.all <- 
  ggtexttable(busco.ids.all, rows=NULL,
              cols=c("BUSCO", "Descriptor"),
              theme=ttheme(base_size=5,
                           padding=unit(c(3, 1), "mm"),
                           tbody.style=tbody_style(size=6, hjust=0, x=0.01))) %>%
  tab_add_title(text="BUSCOs missing in all strains", face = "plain", size=8)

busco.table.Gh <- 
  ggtexttable(busco.ids.Gh, rows=NULL,
              cols=c("BUSCO", "Descriptor"),
              theme=ttheme(base_size=5,
                           padding=unit(c(3, 1), "mm"),
                           tbody.style=tbody_style(size=6, hjust=0, x=0.01))) %>%
  tab_add_title(text="BUSCOs present in Gt/Ga but missing in Gh", face = "plain", size=8)

#Add tables to upset plot
gg.buscos <- gg.upset +
  annotation_custom(ggplotGrob(busco.table),
                    xmin=-9.5, xmax=-9.5,
                    ymin=20.5, ymax=20.5) +
  annotation_custom(ggplotGrob(busco.table.Gh),
                    xmin=15, xmax=15,
                    ymin=24, ymax=24) +
  ggpreview(width=7, height=5)

#Write to file
pdf(file=paste0("R:/GaeumannomycesGenomics/01_assembly/missing_buscos-",
                Sys.Date(), ".pdf"),
    height=5, width=7)
gg.buscos
dev.off()
