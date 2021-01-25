library(data.table)
library(pheatmap)
library(RColorBrewer)
library('philentropy')
library(ggplot2)


dd = fread("analyses/motupan_species_rarefact_w_ppanggolin_cogs.tsv")
dd[, tool := "mOTUpan"]
dd = melt(dd, measure.vars = c('core_len','aux_len'), variable.name = "set", value.name = "nb_genes" )
names(dd)[1] = 'nb_org'
pp = fread("other_soft/ppanggolin/rarefaction_species/rarefaction.csv")
pp[, tool := "ppanggolin"]
pp[,'persistent+shell' := persistent+shell]

pp = melt(pp, measure.vars = c('cloud','persistent','persistent+shell' ,"shell"), variable.name = "set", value.name = "nb_genes" )

cols = intersect(names(dd), names(pp))
species_rar = dd[,..cols]

species_rar = rbind(pp[,..cols],dd[,..cols])

plot_dat = species_rar[sample(nrow(species_rar))][!set %in% c('aux_len','cloud') ]
plot_dat[, set := factor(set, levels = c('core_len', 'persistent', 'shell', 'persistent+shell'))]
levels(plot_dat$set) = c("mOTUpan core", "PPanGGOLiN persistents", "PPanGGOLiN shell",   "persistent+shell")
ggplot(plot_dat[set != 'persistent+shell'], aes(x=nb_org, y=nb_genes, by=set))+
  geom_point(alpha =0.1, col = "grey")+geom_smooth(se=FALSE, col='red')+
  facet_grid(~set)+theme_minimal()+xlab("number of genomes")+ylab("number of COGs in set")+
  theme(text = element_text(size = 20), panel.spacing = unit(2, "lines"))


ggsave("~/temp/motupan_vs_ppanggolin_s__prochlo.svg", width = 12, height = 8)
ggsave("~/temp/motupan_vs_ppanggolin_s__prochlo.pdf", width = 12, height = 8)
ggsave("~/temp/motupan_vs_ppanggolin_s__prochlo.jpg", width = 12, height = 8)







roary_rar = fread("analyses/roary_rarefaction.csv")
roary_rar[,'core+softcore' := core+softcore]
melted_roary_rar = melt(roary_rar, measure.vars = c('core+softcore','shell','cloud', 'motupan_core', 'motupan_cloud'), variable.name = "set", value.name = "nb_genes" )
ggplot(melted_roary_rar[grep("cloud", set, inv= TRUE)], aes(x=nb_org, y=nb_genes, col=set))+geom_point()+geom_hline(yintercept=1875.4772727272727)+geom_smooth(se=FALSE)


ppangg_rar = fread("analyses/ppanggolin_rarefaction.csv")
ppangg_rar[,'core+softcore' := core+softcore]
melted_ppangg_rar = melt(ppangg_rar, measure.vars = c('core+softcore','shell','cloud', 'motupan_core', 'motupan_cloud'), variable.name = "set", value.name = "nb_genes" )
ggplot(melted_ppangg_rar[grep("cloud", set, inv= TRUE)], aes(x=nb_org, y=nb_genes, col=set))+geom_point()+geom_hline(yintercept=1861.9772727272727)+geom_smooth(se=FALSE)





ppangg = read.csv("analyses/ppanggolin_matrix_species.csv", row.names=1)
ppangg_cogs = read.csv("analyses/ppanggolin_matrix_species_cogs.csv", row.names=1)
#row_dist = distance(ppangg, method = "jaccard",use.row.names = TRUE, as.dist.obj = TRUE)
#col_dist = distance(t(ppangg), method = "jaccard",use.row.names = TRUE, as.dist.obj = TRUE)
ppangg_genomes = read.csv("analyses/ppanggolin_matrix_species_genomes.csv", row.names=1)
cols = colorRampPalette(brewer.pal(n = 7, name = "Greys"))(100)
row.names(ppangg_genomes) = gsub("-",".",row.names(ppangg_genomes))
selects = names(which(rowSums(ppangg) > 1))

ppangg =  ppangg[selects, ]
ppangg_cogs = ppangg_cogs[selects,]

core_md = ppangg_cogs[ppangg_cogs$motupan == "core",]
access_md = ppangg_cogs[ppangg_cogs$motupan != "core",]

core = ppangg[row.names(core_md),]
access = ppangg[row.names(access_md),]

row_dist = distance(access, method = "jaccard",use.row.names = TRUE, as.dist.obj = TRUE)
col_dist = distance(t(access), method = "jaccard",use.row.names = TRUE, as.dist.obj = TRUE)

access = data.table(access)
access$cog = row.names(access_md)
core = data.table(core)
core$cog = row.names(core_md)


melted_access = melt(access, variable.name = "genome", value.name = "presence", measure.vars=setdiff(colnames(access), "cog"))
melted_core = melt(core, variable.name = "genome", value.name = "presence", measure.vars=setdiff(colnames(core), "cog"))
melted_access[, set := "accessory"]
melted_core[, set := "core"]

big_heat = rbind(melted_core, melted_access)
#big_heat[,cog := gsub("-",".", cog)]
cog_order = c(access$cog[hclust(row_dist, method = "ward.D2")$order], core$cog)

big_heat[, cog  := factor(cog, levels = cog_order)]
big_heat[, genome  := factor(genome, levels = colnames(access)[hclust(col_dist, method = "ward.D2")$order])]

ggplot(big_heat, aes(x=genome,y=cog, fill=presence==1))+geom_tile()+scale_fill_manual(values=c("white", "black"))+facet_grid(vars(set), scale="free", space = "free")+theme_minimal()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

pheatmap(ppangg, annotation_row=ppangg_cogs, annotation_col=ppangg_genomes, color=cols, show_rownames = FALSE, show_colnames=FALSE, clustering_distance_rows=row_dist, clustering_distance_cols=col_dist, treeheight_row = 0, treeheight_col= 0)




ppangg = read.csv("analyses/ppanggolin_matrix_goods.csv", row.names=1)
ppangg_cogs = read.csv("analyses/ppanggolin_matrix_species_cogs.csv", row.names=1)
row_dist = distance(ppangg, method = "jaccard",use.row.names = TRUE, as.dist.obj = TRUE)
col_dist = distance(t(ppangg), method = "jaccard",use.row.names = TRUE, as.dist.obj = TRUE)

ppangg_genomes = read.csv("analyses/ppanggolin_matrix_genomes.csv", row.names=1)

cols = colorRampPalette(brewer.pal(n = 7, name = "Greys"))(100)

pheatmap(ppangg, annotation_row=ppangg_cogs, annotation_col=ppangg_genomes, color=cols, show_rownames = FALSE, clustering_distance_rows=row_dist, clustering_distance_cols=col_dist)

     ann_colors = list(
         Time = c("white", "firebrick"),
         CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
         GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
     )
