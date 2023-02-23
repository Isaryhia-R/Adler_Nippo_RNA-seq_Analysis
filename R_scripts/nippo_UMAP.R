library(plotly)
library(Seurat)
library(dplyr)


# Use counts as input since Seurat normalizes

counts_0 = read.delim('data/working_matrices/Nippo_Allsamples_02-13-23_counts_matrix_rm_col_filabv1TPM.txt',header=TRUE,row.names="gene_id", stringsAsFactors = FALSE)
head(counts_0)


meta = read.delim('Metadata_samples_2-20-23.txt',header=TRUE, stringsAsFactors = FALSE)
head(meta)


meta$sample_name = gsub(" ",".",meta$sample_name)
meta$sample_name = gsub("-",".",meta$sample_name)



############ Vitro sets day 0 ONLY ("VT_0d.*", meta$stage8_sampletype)
########### vitro timecourse ("VT.*d", meta$stage8_sampletype)


counts=counts_0 #[,!grepl ("IJ", meta$stage2_age)] ## subset the data matrix 

head (counts) 


head (counts) 

#subset meta

meta_1=meta[match(colnames(counts), meta$sample_name),]  ## subset the metadata to match the data
head(meta_1)

meta_1$sample_name = gsub("_",".",meta_1$sample_name)


# sanity check that subsetted matrix and subsetted metadata match
table(meta_1$sample_name == colnames(counts)) 



##################### Make Seurat object ##################### 
nippo <- CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)

# Fill in some metadata
nippo@meta.data$SampleType = meta_1$stage3_host_sex_genotype
nippo@meta.data$HostType = meta_1$stage0d3_host_genotype

nippo@meta.data$VitroDay = meta_1$stage6_post.infection_day


  
##################### Seurat processing ##################### 
nippo = NormalizeData(nippo, normalization.method = "LogNormalize")
nippo = FindVariableFeatures(nippo, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(nippo)
nippo = ScaleData(nippo, features = all.genes)

nippo = RunPCA(nippo, features = VariableFeatures(object = nippo), npcs = 50)

ElbowPlot(nippo, ndims = 50)

# choose # dims where standard deviation of principal kinda levels out
# can play around with number and see how plot changes
nippo = RunUMAP(nippo, dims = 1:10) 

##################### Plot ##################### 
colors = c("#800000","#DC143C",
           "#FF8C00", "gold", "#EE82EE",
           "#31a354","#9370DB","#0000FF",
           "green2","#1E90FF","#40E0D0",
           "mediumvioletred","palevioletred",
           "seagreen1","peachpuff","rosybrown2",
           "aquamarine4", "antiquewhite4", "chocolate3","sienna", "darkcyan",
           "navy","olivedrab4","greenyellow",
           "sandybrown", "darkseagreen", "khaki",
           "khaki3","deeppink", "lightcyan2", "lightcoral"
           ,"saddlebrown","gray8","peru","thistle1", 
           "powderblue","blueviolet", "magenta2")

##old colors not consistent with pca###  
#("#800000","#DC143C",
#"#FF8C00", "gold", 
#"#31a354","rosybrown2","#0000FF",
#"green2","#1E90FF","#40E0D0",
#"mediumvioletred","palevioletred1",
#"peachpuff","#9370DB",
#"antiquewhite4","sienna", "darkcyan",
#"navy","olivedrab4","greenyellow",
#"sandybrown", "darkseagreen", "khaki",
#"khaki3","deeppink", "lightcyan2", "lightcoral"
#,"saddlebrown","gray8","peru","thistle1", 
#"powderblue","blueviolet", "magenta2" )


DimPlot(nippo, reduction = "umap", pt.size = 3, cols = colors,
      shape.by = "HostType", group.by = "VitroDay", )

#####################  3D UMAP ##################### 
# Need to re-run UMAP with 3 dimensions (3L)
nippo <- RunUMAP(nippo, dims = 1:10, n.components = 3L)

DimPlot(nippo, group.by = "SampleType", cols = colors) + NoAxes() + NoLegend()

# Extract tSNE information from Seurat Object
umap_1 <- nippo[["umap"]]@cell.embeddings[,1]
umap_2 <- nippo[["umap"]]@cell.embeddings[,2]
umap_3 <- nippo[["umap"]]@cell.embeddings[,3]

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = nippo, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "SampleType"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Make font preferences settings for axes
 font.pref <- list(size=12, family="Arial, sans-serif", color="black")
x.list <- list(title = "UMAP_1", titlefont = font.pref)
y.list <- list(title = "UMAP_2", titlefont = font.pref)
z.list <- list(title = "UMAP_3",titlefont = font.pref)


# Plot your data, in this example my Seurat object had 21 clusters (0-20)
fig <- plot_ly(data = plot.data, 
               x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
               color = ~SampleType, 
               colors = colors,
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 4, width=3), # controls size of points
               text=~label, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text")  %>%
layout(scene=list(xaxis = x.list,
                    yaxis = y.list,
                    zaxis = z.list,
                    camera = list(eye = list(x = 0, y = 0, z = 2.5))))#When you visualize your plotly object, hovering your mouse pointer over a point shows cell names
fig

###
###
###
###