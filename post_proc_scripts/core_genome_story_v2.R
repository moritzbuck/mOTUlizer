library(data.table)
library(ggplot2)
library(hues)
library(vegan)
library(colorspace)
library("ggpubr")


col2alpha <- function
(x,
 maxValue=1,
 ...)
{
   ## Purpose is to extract the alpha value for a set of colors defined in hex space,
   ## for those R tools that use the inconsistent method of defining alpha separate from
   ## the RGB color, although most tools should really be using RGBA format instead...
   if (length(x) == 0) {
      return(x);
   }
   xNA <- is.na(x);
   alphaValues <- col2rgb(x, alpha=TRUE)["alpha",]/255*maxValue;
   alphaValues[xNA] <- 0;
   return(alphaValues);
}

col2hcl <- function
(x,
 maxColorValue=255,
 model=c("polarLUV", "polarLAB"),
 ...)
{
   if (!suppressPackageStartupMessages(require(colorspace))) {
      stop("The colorspace package is required.");
   }
   model <- match.arg(model);

   x1 <- col2rgb(x);
   a1 <- col2alpha(x);
   x2 <- colorspace::sRGB(t(x1)[,1:3,drop=FALSE]/maxColorValue);
   x3 <- rbind(t(colorspace::coords(as(x2, model))), "alpha"=a1);
   if (length(names(x)) > 0) {
      colnames(x3) <- names(x);
   }
   x3[is.na(x3)] <- 0;
   x3;
}

gtdb_dat = fread("~/uppmax/people/0023_anoxicencyclo/4500_assembly_analysis/cool_dat.csv")
gtdb_dat = gtdb_dat[!is.na(mean_cogs)]
gtdb_dat = gtdb_dat[mean_cogs < 8000]
gtdb_dat[set == "anoxi", set := "freshwater MAGs"]
gtdb_dat = gtdb_dat[sample(nrow(gtdb_dat))]

count_phylum = table(factor(gtdb_dat$phylum))
count_class = table(factor(gtdb_dat$class))
count_order = table(factor(gtdb_dat$order))
count_family = table(factor(gtdb_dat$family))
count_genus = table(factor(gtdb_dat$genus))
gtdb_dat[, otu_per_phylum := count_phylum[phylum] ]
gtdb_dat[, otu_per_class := count_class[class] ]
gtdb_dat[, otu_per_order := count_order[order] ]
gtdb_dat[, otu_per_family := count_family[family] ]
gtdb_dat[, otu_per_genus := count_genus[genus] ]

selected_phylum = levels(factor(gtdb_dat[otu_per_phylum > 3]$phylum))
selected_classes = levels(factor(gtdb_dat[otu_per_phylum > 3]$class))
tt = gtdb_dat[otu_per_phylum > 3, .(class, phylum)]
class2phylum = tt$phylum
names(class2phylum) = tt$class
class2phylum = class2phylum[!duplicated(names(class2phylum))]
phylum2class = sapply(selected_phylum, function(x) names(class2phylum)[class2phylum == x])

base_cols = seq(0,360, 360/(length(selected_phylum)+1))[1:length(selected_phylum)]
incr = as.numeric(base_cols[2])

base_cols = sample(base_cols)
names(base_cols) =  selected_phylum

phyla_of_fig1 = c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Proteobacteria", "Patescibacteria", "Cyanobacteria")

tt = lapply(selected_phylum, function(x){
    n = length(phylum2class[[x]])
    if(n > 1)
    {
      cols = iwanthue(n, hmin = base_cols[x], hmax = base_cols[x]+(incr/2), lmin = 30)
    }
    else
    {
      cols = iwanthue(2, hmin = base_cols[x], hmax = base_cols[x]+(incr/2), lmin = 30)[1]
    }
    names(cols)  = phylum2class[[x]]
    cols
})
class2col = unlist(tt)
names(class2col) = paste(unlist(class2phylum[names(class2col)]), names(class2col), sep=";")

bb = lapply(tt, function(x) apply(coords(hex2RGB(x)),2,mean))
phylum2col = sapply(bb, function(x) rgb(x[1], x[2], x[3]))
names(phylum2col) = selected_phylum

write.table(class2col, "class2col.csv", quote = FALSE, col.names = FALSE, sep=",")
write.table(phylum2col, "phylum2col.csv", quote = FALSE, col.names = FALSE, sep=",")

anon_dat = copy(gtdb_dat)
anon_dat[,class := NULL]
anon_dat[,phylum := NULL]
anon_dat[,set := NULL]

ggplot(gtdb_dat[phylum %in% selected_phylum], aes(x=core_len, col=paste(phylum,class , sep= ";")))+
  facet_wrap(~set)+
  guides(col=guide_legend(ncol=2))+
  geom_histogram(data = anon_dat, alpha = 0.4, fill = "grey", col=NA)+theme_minimal()+
  geom_freqpoly(size=3)+
  xlab("core size (number of COGs)")+
  scale_color_manual(values = class2col)
ggsave("mOTUlizer/doc/figs/Supplemental_fig_1b_core_sizes_hist.pdf", width = 8, height = 6)

anon_dat = copy(gtdb_dat)
anon_dat[,class := NULL]
anon_dat[,phylum := NULL]

fig2_h = 15
fig2_w = 3

ggplot(gtdb_dat[phylum %in% phyla_of_fig1], aes(x=mean_cogs, y=mean_variable, shape = set, col=paste(phylum,class , sep= ";"), size = nb_genomes))+
  geom_point(data = anon_dat, alpha = 0.4, col = '#dcdcdc')+
  geom_point()+
#  geom_smooth(data = anon_dat, method = "lm", size = 1, alpha = 0.8, se=F)+
  facet_wrap(~phylum, ncol=1)+
  guides(col=guide_legend(ncol=2))+
  theme_minimal()+
  scale_y_log10()+scale_x_log10()+
  theme(legend.position = "none")+
  ylab("variable genome size (mean number number of COGs per genome)")+
  xlab("genome size (mean number of COGs per genome)")+
  scale_color_manual(values = class2col)

ggsave("mOTUlizer/doc/figs/Fig_1a_var_vs_genome_size.pdf", width = fig2_w, height = fig2_h)

ggplot(gtdb_dat[phylum %in% phyla_of_fig1], aes(shape = set, x=mean_cogs, y=mean_variable/mean_cogs, col=paste(phylum,class , sep= ";"), size = nb_genomes))+
  geom_point(data = anon_dat, alpha = 0.2, col = '#dcdcdc')+
  geom_point( alpha = 0.7)+
#  geom_smooth(data = anon_dat, method = "lm", size = 1, alpha = 0.8, se=F)+
  facet_wrap(~phylum, ncol=1)+
  theme_minimal()+
  scale_y_log10()+scale_x_log10()+
  ylab("variable fraction of genome")+
  guides(col=guide_legend(ncol=2))+
  xlab("genome size (mean number of COGs per genome)")+
  theme(legend.position = "none")+
  scale_color_manual(values = class2col)#+ylim(0,1)

ggsave("mOTUlizer/doc/figs/Fig_1b_variable_fract_vs_genome_size.pdf", width = fig2_w, height = fig2_h)

ggplot(gtdb_dat[phylum %in% phyla_of_fig1], aes( x=mean_cogs, y=aux_sinsingle/mean_variable/sqrt(nb_genomes), shape = set, col=paste(phylum,class , sep= ";"), size = nb_genomes))+
  geom_point(data = anon_dat, alpha = 0.4, col = "#dcdcdc")+
  geom_point( alpha = 0.7)+
#  geom_smooth(data = anon_dat, method = "lm", size = 1, alpha = 0.8, se=F)+
  facet_wrap(~phylum, ncol=1)+
  theme_minimal()+
  scale_y_log10()+scale_x_log10()+
  guides(col=guide_legend(ncol=2))+
  ylab("auxiliary genome diversity (COGs per position per genome^0.5)")+
  xlab("genome size (mean number of COGs per genome)")+
  theme(legend.position = "none")+
  scale_color_manual(values = class2col)#+ylim(0,1)

ggsave("mOTUlizer/doc/figs/Fig_1c_pangediv_vs_genome_size.pdf", width = fig2_w, height = fig2_h)

ggplot(gtdb_dat[phylum %in% phyla_of_fig1], aes( x=mean_cogs, y=aux_sinsingle, shape = set, col=paste(phylum,class , sep= ";"), size = nb_genomes))+
  geom_point(data = anon_dat, alpha = 0.4, col = "#dcdcdc")+
  geom_point( alpha = 0.7)+
#  geom_smooth(data = anon_dat, method = "lm", size = 1, alpha = 0.8, se=F)+
  facet_wrap(~phylum, ncol=1)+
  theme_minimal()+
  scale_y_log10()+scale_x_log10()+
  guides(col=guide_legend(ncol=2))+
  ylab("absolute auxiliary genome size")+
  xlab("genome size (mean number of COGs per genome)")+
  theme(legend.position = "none")+
  scale_color_manual(values = class2col)#+ylim(0,1)

ggsave("mOTUlizer/doc/figs/Fig_1d_pangesize_vs_genome_size.pdf", width = fig2_w, height = fig2_h)

p1 = ggplot(gtdb_dat[phylum %in% phyla_of_fig1], aes( x=mean_cogs, y=aux_sinsingle, shape = set, col=paste(phylum,class , sep= ";"), size = nb_genomes))+
    geom_point( alpha = 0.7)+
    theme_minimal()+
    guides(shape = guide_legend(override.aes = list(size = 5)), col = guide_legend(override.aes = list(size = 5), ncol=1))+
    scale_color_manual(values = class2col)#+ylim(0,1)

plot(as_ggplot(get_legend(p1)))
ggsave("mOTUlizer/doc/figs/Fig_1d_legend.pdf", width = fig2_w, height = fig2_h)

ggplot(gtdb_dat[class %in% valids_classes], aes(x=mean_variable*log(nb_genomes), y=aux_sinsingle, col=class, size = nb_genomes))+
  geom_point(data = anon_dat, alpha = 0.4, col = "grey")+
  geom_point( alpha = 0.7)+scale_color_manual(values = classes_cscale)+
  facet_wrap(~phylum, ncol=3)+
  geom_smooth(data = anon_dat, method = "lm", size = 1, alpha = 0.8, se=F, col="red")+
  scale_y_log10()+scale_x_log10()+
  theme_minimal()

ggplot(gtdb_dat[phylum %in% selected_phylum], aes(shape = set, x=mean_cogs, y=mean_variable/mean_cogs, col=paste(phylum,class , sep= ";"), size = nb_genomes))+
  geom_point(data = anon_dat, alpha = 0.2, col = '#dcdcdc')+
  geom_point( alpha = 0.7)+
#  geom_smooth(data = anon_dat, method = "lm", size = 1, alpha = 0.8, se=F)+
  theme_minimal()+
  scale_y_log10()+scale_x_log10()+
  ylab("variable fraction of genome")+
  guides(col=guide_legend(ncol=2))+
  xlab("genome size (mean number of COGs per genome)")+
  theme(legend.position = "none")+
  scale_color_manual(values = class2col)#+ylim(0,1)


ggplot(gtdb_dat[set != 'gtdb'], aes(x=mean_cogs, y=mean_variable/mean_cogs, col=core_gc-aux_sinsingle_gc, size = nb_genomes))+
  geom_point( alpha = 0.7)

defreshMG
