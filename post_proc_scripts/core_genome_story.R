library(data.table)
library(ggplot2)
library(hues)
library(vegan)

runs = list.dirs(path = "outputs/", full.names=FALSE)
runs = runs[2:length(runs)]

whole_data = lapply(runs, function(r)
  list(
    data = fread(paste("outputs",r,"otu.stats", sep="/", collapse="/")),
    rarefy = fread(paste("outputs",r,"rarefy.stats", sep="/", collapse="/")),
    counts = fread(paste("outputs",r,"counts.stats", sep="/", collapse="/"))
  )
)
names(whole_data) = runs

for( run in  names(whole_data) ){
  for(frame in names(whole_data[[run]]))
  {
      whole_data[[run]][[frame]]$phylum = sapply(strsplit(whole_data[[run]][[frame]]$taxonomy, ";"), "[",2)
      whole_data[[run]][[frame]]$class = sapply(strsplit(whole_data[[run]][[frame]]$taxonomy, ";"), "[",3)
      whole_data[[run]][[frame]]$order = sapply(strsplit(whole_data[[run]][[frame]]$taxonomy, ";"), "[",4)
      whole_data[[run]][[frame]]$family = sapply(strsplit(whole_data[[run]][[frame]]$taxonomy, ";"), "[",5)
      whole_data[[run]][[frame]]$genus = sapply(strsplit(whole_data[[run]][[frame]]$taxonomy, ";"), "[",6)
      whole_data[[run]][[frame]][,run := run]
      whole_data[[run]][[frame]][,ani := as.numeric(sapply(strsplit(run,"_"), "[", 2))]
      whole_data[[run]][[frame]][,complete := as.numeric(sapply(strsplit(run,"_"), "[", 3))]
  }
  count_phylum = table(factor(whole_data[[run]]$data$phylum))
  count_class = table(factor(whole_data[[run]]$data$class))
  count_order = table(factor(whole_data[[run]]$data$order))
  count_family = table(factor(whole_data[[run]]$data$family))
  count_genus = table(factor(whole_data[[run]]$data$genus))
  for(frame in names(whole_data[[run]]))
  {
      whole_data[[run]][[frame]][, otu_per_phylum := count_phylum[phylum] ]
      whole_data[[run]][[frame]][, otu_per_class := count_class[class] ]
      whole_data[[run]][[frame]][, otu_per_order := count_order[order] ]
      whole_data[[run]][[frame]][, otu_per_family := count_family[family] ]
      whole_data[[run]][[frame]][, otu_per_genus := count_genus[genus] ]
  }
}


ggplot(whole_data$otus_95_70$rarefy[otu_per_genus > 5] , aes(x=rr_genome_count, y=rr_mean, col = taxonomy, fill = taxonomy, by=otu, ymax=rr_max_95, ymin=rr_min_95)) + geom_point(size=0.3)+geom_line(size=0.3)+facet_grid(singletons~genus)+geom_ribbon(alpha =0.05, size=0.1)+ theme(legend.position = "none")
ggplot(whole_data$otus_95_70$data, aes(x=mean_cog_count, y=nosingle_pan/log(nb_genomes), col = phylum, size=nb_genomes)) + geom_point()+scale_color_manual(values=as.vector(iwanthue(length(unique(whole_data$otus_95_70$data$phylum)))))+scale_y_log10()+scale_x_log10()
ggplot(whole_data$otus_95_70$counts[otu_per_family > 3], aes(x=fract, y=cog_count, col = taxonomy, by=otu)) + geom_point()+geom_line()+facet_wrap(~phylum+family, ncol=4)+ theme(legend.position = "none")+scale_y_log10()


COG_CATS = unlist(list(
"A" = "RNA processing and modification",
"B" = "Chromatin Structure and dynamics",
"C" = "Energy production and conversion",
"D" = "Cell cycle control and mitosis",
"E" = "Amino Acid metabolis and transport",
"F" = "Nucleotide metabolism and transport",
"G" = "Carbohydrate metabolism and transport",
"H" = "Coenzyme metabolis",
"I" = "Lipid metabolism",
"J" = "Translation",
"K" = "Transcription",
"L" = "Replication and repair",
"M" = "Cell wall/membrane/envelop biogenesis",
"N" = "Cell motility",
"O" = "Post-translational modification, protein turnover, chaperone functions",
"P" = "Inorganic ion transport and metabolism",
"Q" = "Secondary Structure",
"T" = "Signal Transduction",
"U" = "Intracellular trafficing and secretion",
"V" = "Defense mechanisms",
"W" = "Extracellular structures",
"Y" = "Nuclear structure",
"Z" = "Cytoskeleton",
"R" = "General Functional Prediction only",
"S" = "Function Unknown",
"no_annot" = "no_annot"
))

goperset = fread("COG_cats.csv")
bla = goperset[,Filter(is.numeric, .SD)]
setkey(whole_data$otus_95_70$data, "otu")
 colnames(bla) = COG_CATS[ colnames(bla)]
tokeep = colnames(bla)[colnames(bla) != "no_annot"]
bla = bla[,..tokeep ]
otus = goperset$V1
otus = otus[rowSums(bla) > 0]
bla = bla[rowSums(bla) > 0,]
bla = bla[grepl("_core",  otus)]
otus = otus[grepl("_core",  otus)]
#otus = otus[rowSums(bla) > 50]
#bla = bla[rowSums(bla) > 50]

bla = bla/whole_data$otus_95_70$data[sub("_singletons", "", sub("_nosingle_pan","",sub("_core","", otus)))]$core_size
#bla = bla[!startsWith(otus,"aniOTU_78_")]
#otus = otus[!startsWith(otus,"aniOTU_78_")]
#bla = bla[,colSums(bla > 0 )>100, with = FALSE]
#bla = bla[,colnames(bla) != "Function Unknown", with = FALSE]


mds = metaMDS(bla, maxit = 100, weakties = TRUE, dist = "euclidean", trymax=40, k=2)
dd = data.table(mds$point)
dd$otu = otus
dd[, type := sapply(strsplit(otu,"_"), "[", 3)]
specs = data.table(mds$species)
specs$otu = colnames(bla)
dd[, otu_id := sub("_singletons","", sub("_nosingle_pan", "", sub("_core","", otu))) ]
dd = merge(dd,whole_data$otus_95_70$data, by.x = "otu_id", by.y = "otu")
dd2 = copy(dd)
dd2[, phylum := NULL]
dd2[, family := NULL]
ggplot(dd, aes(x=MDS1, y=MDS2, col=type)) + geom_point()


ggplot(dd[otu_per_phylum > 3], aes(x=MDS1, y=MDS2, col = paste0(phylum,class, sep =";"), label = otu, size=mean_cog_count))+geom_point(data = dd2, col="black", alpha=0.1)+ geom_point()+scale_color_manual(values=as.vector(iwanthue(length(unique(whole_data$otus_95_70$data[otu_per_phylum > 3]$class)))))+facet_wrap(~phylum)

bla = goperset[,Filter(is.numeric, .SD)]
colnames(bla) = COG_CATS[ colnames(bla)]
tokeep = colnames(bla)[colnames(bla) != "no_annot"]
bla = bla[,..tokeep ]
otus = goperset$V1
otus = otus[rowSums(bla) > 0]
bla = bla[rowSums(bla) > 0,]
bla = bla[grepl("_nosingle_pan",  otus)]
otus = otus[grepl("_nosingle_pan",  otus)]

bla = bla/whole_data$otus_95_70$data[sub("_nosingle_pan","", otus)]$nosingle_pan
#prots = whole_data$otus_95_70$data[sub("_nosingle_pan","",otus)][phylum=="Proteobacteria"]$otu

bla = bla[!startsWith(otus,"aniOTU_78_")]
otus = otus[!startsWith(otus,"aniOTU_78_")]
#bla = bla[sub("_nosingle_pan","",otus) %in% prots]
#otus = otus[sub("_nosingle_pan","",otus) %in% prots]


bla = bla[,colSums(bla) > 0 , with = FALSE]
bla = bla[,colnames(bla) != "Function Unknown", with = FALSE]


mds = metaMDS(bla, maxit = 100, weakties = TRUE, dist = "euclidean", trymax=40, k=2)
dd = data.table(mds$point)
dd$otu = otus
dd[, type := sapply(strsplit(otu,"_"), "[", 3)]
specs = data.table(mds$species)
specs$otu = colnames(bla)
dd[, otu_id := sub("_nosingle_pan", "", sub("_core","", otu)) ]
dd = merge(dd,whole_data$otus_95_70$data, by.x = "otu_id", by.y = "otu")
dd2 = copy(dd)
dd2[, phylum := NULL]


ggplot(dd[otu_per_phylum > 3], aes(x=MDS1, y=MDS2, col = paste0(phylum,class, sep =";"), label = otu, size=mean_cog_count))+geom_point(data = dd2, col="black", alpha=0.1)+ geom_point()+scale_color_manual(values=as.vector(iwanthue(length(unique(whole_data$otus_95_70$data[otu_per_phylum > 3]$class)))))+facet_wrap(~phylum)+xlim(-0.1,0.1)+ylim(-0.1,0.1)
