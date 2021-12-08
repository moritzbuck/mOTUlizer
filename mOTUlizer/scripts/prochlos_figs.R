library(data.table)
library(pheatmap)
library(RColorBrewer)
#library('philentropy')
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(scales)

dd = fread("analyses/motupan_species_rarefact_w_roary_cogs.tsv")
dd[, tool := "mOTUpan"]
dd = melt(dd, measure.vars = c('core_len','aux_len'), variable.name = "set", value.name = "nb_genes" )
names(dd)[1] = 'nb_org'
pp = fread("other_soft/roary/rarefaction_species/rarefaction.csv")
pp[, tool := "roary"]
pp[,'persistent+shell' := persistent+shell]

pp = melt(pp, measure.vars = c('cloud','persistent','persistent+shell' ,"shell"), variable.name = "set", value.name = "nb_genes" )

cols = intersect(names(dd), names(pp))
species_rar = dd[,..cols]

species_rar = rbind(pp[,..cols],dd[,..cols])

plot_dat = species_rar[sample(nrow(species_rar))][!set %in% c('aux_len','cloud') ]
plot_dat[, set := factor(set, levels = c('core_len', 'persistent', 'shell', 'persistent+shell'))]
levels(plot_dat$set) = c("mOTUpan core", "roary persistents", "roary shell",   "persistent+shell")
ggplot(plot_dat[set != 'persistent+shell'], aes(x=nb_org, y=nb_genes, by=set))+
  geom_point(alpha =0.1, col = "grey")+geom_smooth(se=FALSE, col='red')+
  facet_grid(~set)+theme_minimal()+xlab("number of genomes")+ylab("number of COGs in set")+
  theme(text = element_text(size = 20), panel.spacing = unit(2, "lines"))


ggsave("~/temp/motupan_vs_roary_s__prochlo.svg", width = 12, height = 8)
ggsave("~/temp/motupan_vs_roary_s__prochlo.pdf", width = 12, height = 8)
ggsave("~/temp/motupan_vs_roary_s__prochlo.jpg", width = 12, height = 8)







roary_rar = fread("analyses/roary_rarefaction.csv")
roary_rar[,'core+softcore' := core+softcore]
melted_roary_rar = melt(roary_rar, measure.vars = c('core+softcore','shell','cloud', 'motupan_core', 'motupan_cloud'), variable.name = "set", value.name = "nb_genes" )
ggplot(melted_roary_rar[grep("cloud", set, inv= TRUE)], aes(x=nb_org, y=nb_genes, col=set))+geom_point()+geom_hline(yintercept=1875.4772727272727)+geom_smooth(se=FALSE)


roary_rar = fread("analyses/roary_rarefaction.csv")
roary_rar[,'core+softcore' := core+softcore]
melted_roary_rar = melt(roary_rar, measure.vars = c('core+softcore','shell','cloud', 'motupan_core', 'motupan_cloud'), variable.name = "set", value.name = "nb_genes" )
ggplot(melted_roary_rar[grep("cloud", set, inv= TRUE)], aes(x=nb_org, y=nb_genes, col=set))+geom_point()+geom_hline(yintercept=1861.9772727272727)+geom_smooth(se=FALSE)





roary = read.csv("analyses/roary_matrix_species.csv", row.names=1)
roary_cogs = read.csv("analyses/roary_matrix_species_cogs.csv", row.names=1)
#row_dist = distance(roary, method = "jaccard",use.row.names = TRUE, as.dist.obj = TRUE)
#col_dist = distance(t(roary), method = "jaccard",use.row.names = TRUE, as.dist.obj = TRUE)
roary_genomes = read.csv("analyses/roary_matrix_species_genomes.csv", row.names=1)
cols = colorRampPalette(brewer.pal(n = 7, name = "Greys"))(100)
row.names(roary_genomes) = gsub("-",".",row.names(roary_genomes))
selects = names(which(rowSums(roary) > 1))

roary =  roary[selects, ]
roary_cogs = roary_cogs[selects,]

core_md = roary_cogs[roary_cogs$motupan == "core",]
access_md = roary_cogs[roary_cogs$motupan != "core",]

core = roary[row.names(core_md),]
access = roary[row.names(access_md),]

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

pheatmap(roary, annotation_row=roary_cogs, annotation_col=roary_genomes, color=cols, show_rownames = FALSE, show_colnames=FALSE, clustering_distance_rows=row_dist, clustering_distance_cols=col_dist, treeheight_row = 0, treeheight_col= 0)




roary = read.csv("analyses/roary_matrix_goods.csv", row.names=1)
roary_cogs = read.csv("analyses/roary_matrix_species_cogs.csv", row.names=1)
row_dist = distance(roary, method = "jaccard",use.row.names = TRUE, as.dist.obj = TRUE)
col_dist = distance(t(roary), method = "jaccard",use.row.names = TRUE, as.dist.obj = TRUE)

roary_genomes = read.csv("analyses/roary_matrix_genomes.csv", row.names=1)

cols = colorRampPalette(brewer.pal(n = 7, name = "Greys"))(100)

pheatmap(roary, annotation_row=roary_cogs, annotation_col=roary_genomes, color=cols, show_rownames = FALSE, clustering_distance_rows=row_dist, clustering_distance_cols=col_dist)

     ann_colors = list(
         Time = c("white", "firebrick"),
         CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
         GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
     )

roary_species_cores = as.data.table(read.csv("analyses/roary_species_stats.csv", as.is =TRUE))
roary_species_cores[, diff := (motupan-roary)/mean(c(motupan, roary) )]
roary_species_cores[,clade := sapply(strsplit(taxo,";"), "[", 7)]
ggplot(roary_species_cores[sample(nrow(roary_species_cores))], aes(x=roary,y=motupan, size = nb_genomes, label=clade, col = type))+
  geom_abline(size=2, col="grey")+
  theme_minimal()+scale_size_continuous(trans="log10")+
  geom_point(alpha=0.5)+scale_color_brewer(palette="Dark2")+
  xlab("COGs in roary's persistent set")+
  ylab("COGs in mOTUpan's core set")
ggsave("~/temp/morupan_vs_roary.pdf", width = 7, height = 6)




ggplot(roary_species_cores[sample(nrow(roary_species_cores))], aes(x=mean_completeness,y=diff, size = nb_genomes, label=clade, col = type))+
  geom_hline(yintercept=0, size=2, col="grey")+
  theme_minimal()+scale_size_continuous(trans="log10")+
  geom_point(alpha=0.5)+scale_color_brewer(palette="Dark2")+
  xlab("mean CheckM completeness")+
  ylab("normalized residuals")
ggsave("~/temp/normalized_residuals.pdf", width = 7, height = 6)

meted = melt(roary_species_cores, measure.vars=c('motupan','roary'))
ggplot(meted[sample(nrow(meted))], aes(x=mean_completeness, y=value/est_size, col=variable, size=nb_genomes))+
  geom_line(mapping=aes(by=species), col="black", alpha=0.4, size=0.3)+scale_size_continuous(trans="log10")+
  geom_point(mapping = aes(shape=type), alpha=0.5)+geom_smooth(alpha=0.2, size=3)+
  theme_minimal()+
  xlab('mean completeness')+ylab('core-COGs per Mb')#+ylim(0,1100)

ggsave("~/temp/motupan_vs_roary_core_fract.pdf", width = 7, height = 6)





roary_species_cores = as.data.table(read.csv("analyses/species_stats.csv", as.is =TRUE))
#roary_species_cores = roary_species_cores[mean_completeness_roary > 80]

roary_species_cores[, diff := (motupan_w_roary-roary_core)/mean(c(motupan_w_roary, roary_core) )]
roary_species_cores[,clade := sapply(strsplit(taxo,";"), "[", 7)]
ggplot(roary_species_cores[sample(nrow(roary_species_cores))], aes(x=roary_core,y=motupan_w_roary, size = nb_genomes, label=clade, col = type))+
  geom_abline(size=2, col="grey")+
  theme_minimal()+scale_size_continuous(trans="log10")+
  geom_point(alpha=0.5)+scale_color_brewer(palette="Dark2")+
  xlab("COGs in roary's core set")+
  ylab("COGs in mOTUpan's core set")
ggsave("~/temp/motupan_vs_roary.pdf", width = 7, height = 6)

ggplot(roary_species_cores[sample(nrow(roary_species_cores))], aes(x=mean_completeness_roary,y=diff, size = nb_genomes, label=clade, col = type))+
  geom_hline(yintercept=0, size=2, col="grey")+
  theme_minimal()+scale_size_continuous(trans="log10")+
  geom_point(alpha=0.5)+scale_color_brewer(palette="Dark2")+
  xlab("mean CheckM completeness")+
  ylab("normalized residuals")
ggsave("~/temp/roary_normalized_residuals.pdf", width = 7, height = 6)

meted = melt(roary_species_cores, measure.vars=c('motupan_w_roary','roary_core'))
ggplot(meted[sample(nrow(meted))], aes(x=mean_completeness_roary, y=value/mean_est_roary_cogs, col=variable, size=nb_genomes))+
  geom_line(mapping=aes(by=species), col="black", alpha=0.4, size=0.3)+scale_size_continuous(trans="log10")+
  geom_point(mapping = aes(shape=type), alpha=0.5)+geom_smooth(alpha=0.2, size=3)+
  theme_minimal()+
  xlab('mean completeness')+ylab('Fraction of core COGs')#+ylim(0,1100)

ggsave("~/temp/motupan_vs_roary_core_fract.pdf", width = 7, height = 6)
roary_species_cores[,ppan_persist_p_shell := ppanggolin_persist + ppanggolin_shell]

roary_species_cores[, diff := (motupan_w_ppan-ppanggolin_persist)/mean(c(motupan_w_ppan, ppanggolin_persist) )]
ggplot(roary_species_cores[sample(nrow(roary_species_cores))], aes(x=ppanggolin_persist,y=motupan_w_ppan, size = nb_genomes, label=clade, col = type))+
  geom_abline(size=2, col="grey")+
  theme_minimal()+scale_size_continuous(trans="log10")+
  geom_point(alpha=0.5)+scale_color_brewer(palette="Dark2")+
  xlab("COGs in ppanggolin's persistent set")+
  ylab("COGs in mOTUpan's core set")
ggsave("~/temp/motupan_vs_ppan.pdf", width = 7, height = 6)

ggplot(roary_species_cores[sample(nrow(roary_species_cores))], aes(x=mean_completeness_ppan,y=diff, size = nb_genomes, label=clade, col = type))+
  geom_hline(yintercept=0, size=2, col="grey")+
  theme_minimal()+scale_size_continuous(trans="log10")+
  geom_point(alpha=0.5)+scale_color_brewer(palette="Dark2")+
  xlab("mean CheckM completeness")+
  ylab("normalized residuals")
ggsave("~/temp/ppan_normalized_residuals.pdf", width = 7, height = 6)



meted = melt(roary_species_cores, measure.vars=c('motupan_w_ppan','ppanggolin_persist'))
ggplot(meted[sample(nrow(meted))], aes(x=mean_completeness_ppan, y=value/mean_est_roary_cogs, col=variable, size=nb_genomes))+
  geom_line(mapping=aes(by=species), col="black", alpha=0.4, size=0.3)+scale_size_continuous(trans="log10")+
  geom_point(mapping = aes(shape=type), alpha=0.5)+geom_smooth(alpha=0.2, size=3)+
  theme_minimal()+
  xlab('mean completeness')+ylab('Fraction of core COGs')#+ylim(0,1100)

  ggsave("~/temp/motupan_vs_ppanggolin_core_fract.pdf", width = 7, height = 6)

rar_new = fread("analyses/motupan_rarefact_w_ppanggolin_cogs_3.tsv")
tt = rar_new$empirical_fpr > 0 & rar_new$bootstrapped_fpr > 0
rar_new = rar_new[tt]
mm = melt(rar_new[,c("nb_org", "empirical_fpr", "bootstrapped_fpr")], id.vars = "nb_org")

g1 = ggplot(rar_new, aes(x=empirical_fpr, y=bootstrapped_fpr, col=nb_org))+geom_point()+
  geom_density2d()+geom_smooth(method="lm", se = FALSE, col="purple", size=1.5)+
  geom_abline(col="red", size=1.5)+theme_minimal()+scale_y_log10()+scale_x_log10()+scale_color_gradientn(colors = c("#91bfdb", "#ffffbf", "#fc8d59"), trans = "log10", name = "Genome count:")+
  xlab("Empirical fpr")+ylab("Bootstrapped fpr")
g2 = ggplot(mm, aes(x=nb_org, y=value, col=variable))+geom_point(alpha=0.3)+geom_smooth(se = FALSE, size=1.5)+scale_y_log10()+theme_minimal()+
  scale_color_manual(values = c(muted("#1b9e77"), muted("#d95f02")),  labels = c("Empirical", "Bootstrapped"), name = "fpr type:")+
  xlab("Genome count")+ylab("False positive rate")

ggarrange(g1,g2,labels = c("A", "B"), ncol=2)
ggsave("analyses/sup_fig.pdf", width = 12, height = 8)
