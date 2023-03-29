


dat_tmp <- dat[dat$Image==30,]
dat_tmp$Phenotype = dat_tmp$cellType
dat_tmp$Phenotype[which(dat_tmp$cellType%in%paste0('c',c(1,8)))] = "CD4 T"
dat_tmp$Phenotype[which(dat_tmp$cellType%in%paste0('c',5:6))] = "CD8 T"
dat_tmp$Phenotype[which(dat_tmp$cellType%in%paste0('c',7))] = "Other"
dat_tmp$Phenotype[which(dat_tmp$cellType%in%paste0('c',4))] = "B"
dat_tmp$Phenotype[which(dat_tmp$cellType%in%paste0('c',c(2:3)))] = "Tumour"
dat_tmp$Phenotype[which(dat_tmp$cellType%in%paste0('c',c(9:18)))] = "Other"

data_plot <- dat_tmp[,c('x','y','Phenotype')]

ggplot2::ggplot() +
  ggplot2::geom_point(mapping=ggplot2::aes(x=data_plot[,1],y=data_plot[,2],
                                           col=data_plot[,3]),
                      size=2) +
  ggplot2::theme_bw() +
  ggplot2::xlim(c(min(data_plot[,1]),max(data_plot[,1]))) +
  ggplot2::ylim(c(min(data_plot[,2]),max(data_plot[,2]))) +
  ggplot2::guides(color='none')+#ggplot2::guide_legend(title="Phenotype")) +
  ggplot2::theme(axis.title = ggplot2::element_blank(),
                 axis.text = ggplot2::element_blank()) +
  ggplot2::scale_color_manual(values = c( "red", "blue", "darkcyan", "darkgreen", 'gray'),
                              breaks=c("Tumour", "CD4 T", "CD8 T", "B", "Other"))


phenotypes <- dat_tmp$Phenotype
coord_x <- dat_tmp$x
coord_y <- dat_tmp$y
imagep2 <- SPIAT::format_image_to_spe(format = "general",
                                      phenotypes = phenotypes,
                                      coord_x = coord_x,coord_y = coord_y)
my_colors <- c("red", "blue", "darkcyan", "darkgreen")

SPIAT::plot_cell_categories(spe_object = imagep2,
                            categories_of_interest =
                              c("Tumour", "CD4 T", "CD8 T", "B"),
                            colour_vector = my_colors, feature_colname = "Phenotype")

nrow(dat_tmp)
table(dat_tmp$Phenotype)




data_plot <- dat[which(dat$Person==unique(dat$Person)[1]),]
data_plot$cellType <-
  ifelse(data_plot$Phenotype%in%c("Tumour", "CD4 T", "CD8 T", "B"),
         data_plot$Phenotype,'Other')
data_plot <- data_plot[,c('cellx','celly','cellType')]
plotPP(data_plot,ptSize = 2, dropAxes = T,
       colorGuide = ,
       colors = c( "darkgreen", "blue", "darkcyan",'gray', "red")) +
  scale_color_discrete(c( "darkgreen", "blue", "darkcyan",'gray', "red"),
                       breaks=c('Tumour', 'CD4 T', 'CD8 T','B','OTHER'))


########################



dat_tmp <- dat0[dat0$Image==1,]
dat_tmp$Phenotype = dat_tmp$cellType
dat_tmp$Phenotype[which(dat_tmp$cellType%in%paste0('c',c(1,8)))] = "CD4 T"
dat_tmp$Phenotype[which(dat_tmp$cellType%in%paste0('c',5:6))] = "CD8 T"
dat_tmp$Phenotype[which(dat_tmp$cellType%in%paste0('c',7))] = "Other1"
dat_tmp$Phenotype[which(dat_tmp$cellType%in%paste0('c',4))] = "B"
dat_tmp$Phenotype[which(dat_tmp$cellType%in%paste0('c',c(2:3)))] = "Tumour"
dat_tmp$Phenotype[which(dat_tmp$cellType%in%paste0('c',c(9:18)))] = "Other"

phenotypes <- dat_tmp$Phenotype
coord_x <- dat_tmp$x
coord_y <- dat_tmp$y
imagep2 <- SPIAT::format_image_to_spe(format = "general",
                                      phenotypes = phenotypes,
                                      coord_x = coord_x,coord_y = coord_y)
my_colors <- c("red", "blue", "darkcyan", "darkgreen")

SPIAT::plot_cell_categories(spe_object = imagep2,
                            categories_of_interest =
                              c("Tumour", "CD4 T", "CD8 T", "B"),
                            colour_vector = my_colors, feature_colname = "Phenotype")

nrow(dat_tmp)
table(dat_tmp$Phenotype)
