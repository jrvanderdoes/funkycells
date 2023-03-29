
#datap2 =data1[which(data1$Person==unique(data1$Person)[3]),]
datap2 =data0[which(data0$Person==unique(data0$Person)[5]),]

# I have tested out and used the SPIAT package for plotting:

phenotypes <- datap2$Phenotype
coord_x <- datap2$cellx
coord_y <- datap2$celly
imagep2 <- SPIAT::format_image_to_spe(format = "general",
                                      phenotypes = phenotypes,
                                      coord_x = coord_x,coord_y = coord_y)
my_colors <- c("red", "blue", "darkcyan", "darkgreen")

# pdf(paste("0_2", ".pdf", sep=""))
plot_cell_categories(spe_object = imagep2,
                     categories_of_interest =
                       c("Tumour", "CD4 T", "CD8 T", "B"),
                     colour_vector = my_colors, feature_colname = "Phenotype")
