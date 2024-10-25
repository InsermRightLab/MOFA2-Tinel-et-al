## MOFA2 downstream analyses:
#--------------------------------------------------------------
## ----------------------------------------------------------------------------
library(ggplot2)
library(MOFA2)
library(msigdbr)
library(fgsea)
library(org.Mm.eg.db)
library(DT)
library(ggpubr)
library(scales)
library(dplyr)
library(ggrepel)
library(reshape2)

## -----------------------------------------------------------------------------
sessionInfo()

## -----------------------------------------------------------------------------

#filepath <- system.file("extdata", "model.hdf5_300822_100inits_Model61", package = "MOFA2")
#model <- load_model("model.hdf5_300822_100inits_Model61")
#outfile = file.path(getwd(),"model.hdf5_300822_100inits_Model61")
#model <- MOFAobject

#We load the saved model 61 obtained by the script "MOFA2_Model_Training_25.10.24.R"
model <- readRDS("model.hdf5_300822_100inits_Model61")



## -----------------------------------------------------------------------------

plot_data_overview(model)

views_names(model) <- c("Blood.RNAseq", "Blood.mRNAs", "Blood.miRNAs", "Biopsy.mRNAs", "Biopsy.miRNAs", "Urine.mRNAs")

cols = c('Blood.RNAseq'='#a50026','Blood.mRNAs'='#d73027',"Blood.miRNAs"="#f46d43"
         ,"Biopsy.mRNAs"="#4575b4","Biopsy.miRNAs"="#74add1",
         "Urine.mRNAs"="#ffffbf")

plot_data_overview(model, colors=cols)

#views_names(MOFAobject) <- c("All Blood miRNAs", "All Biopsy miRNAs")

## -----------------------------------------------------------------------------
Nsamples = sum(get_dimensions(model)[["N"]])


## -----------------------------------------------------------------------------
# Total variance explained per view
head(get_variance_explained(model)$r2_total[[1]])
## -----------------------------------------------------------------------------
plot_variance_explained(model, x="view", y="factor", plot_total = TRUE)[[2]]
#We reproduce this plot in graphpad to put nice colors on it


## -----------------------------------------------------------------------------
# Variance explained for every factor in per view
head(get_variance_explained(model)$r2_per_factor[[1]])

plot_variance_explained(model, x="view", y="factor", x.cex=6,
                        min_r2 = 0,
                        max_r2 = 20)+
 coord_flip()+scale_x_discrete(limits=c("Biopsy.miRNAs","Biopsy.mRNAs", "Blood.miRNAs","Blood.mRNAs", "Blood.RNAseq", "Urine.mRNAs"))


## -----------------------------------------------------------------------------
# Correlation of each factor 
col3 <- colorRampPalette(c("white", "white", "#0F1487")) 

plot_factor_cor(model, col = col3(20), "pearson")
#Correlation plot looks much better (uncorrelated factors) after implementation of the function "scale.views"



#METADATA########
#We complete the metadata file with a lot of different variables such as:

#central Banff scores i ti ah g ptc C4d ci ct v Thrombi (0,1),DSA(O,1),
#fraction of sclerotic glomeruli calculated with total glomeruli and sclerotic glomeruli ,  
#time after transplantation calculated with specific R code,
####We replaced missing value by 0 and the time after transplantation >5000 days by 5000 days

# Tacrolimus (0,1)
#Primary renal disease
#Type of induction
#time of cold ischemia calculated in minutes
#creatinine calculated in umol/L
#blood concentration of erythrovytes, platelets and white blood cells
## -We add metadata to the model----------------------------------------------------------------------------

sample_metadata <-read.csv2(file="Step_1_metadata130922.csv", header=T, na.string="NA", sep=";", dec=",")
samples_metadata(model) <- sample_metadata
head(samples_metadata(model), n=3)

#Correlating the factors according to the patients outcomes provided by the metadata

#Extract the values of each factor for each patient
p <- plot_factor(model, 
                 factors = 1:8,
                 color_by = "sample")
print(p)

IDsPatients<-as.data.frame(p$data$sample)
Factors<-as.data.frame(p$data$factor)
Values<-as.data.frame(p$data$value)

#We turn it into a dataframe and a pivot table
merged_df_factors <- cbind(IDsPatients,Factors, Values)
pivot_table <- dcast(merged_df_factors, p$data$sample ~ p$data$factor, value.var = "p$data$value")



#We merge it with metadata
factors_vs_metadata <- merge(pivot_table,sample_metadata, by.x=c("p$data$sample"), by.y=c("sample"))

#write.csv2(factors_vs_metadata, file="factorsvsmetadata15.06.23.csv")

#We figured out that factor 1 is related to HLA DSA, ti and cg lesions
########################################ACCORDING TO Banffs lesions##############
p <- plot_factor(model, 
                 factors = 1,
                 color_by = "Donor.Specific.Antibodies",
                 dot_size = 1,        # change dot size
                 dodge = TRUE,           # dodge points with different colors
                 legend = TRUE,          # remove legend
                 add_violin = TRUE,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

print(p)

#Already done! 
p <- plot_factor(model, 
                 factors = 1,
                 color_by = "Central...T.cell.mediated.rejection...ti",
                 dot_size = 1,        # change dot size
                 dodge = TRUE,           # dodge points with different colors
                 legend = TRUE,          # remove legend
                 add_violin = TRUE,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

print(p)

p <- plot_factor(model, 
                 factors = 1,
                 color_by = "Central...Antibody.mediated.rejection...cg",
                 dot_size = 1,        # change dot size
                 dodge = TRUE,           # dodge points with different colors
                 legend = TRUE,          # remove legend
                 add_violin = TRUE,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

print(p)

plot_data_heatmap(model,
                  view = "Blood.mRNAs",         # view of interest
                  factor = 1,             # factor of interest
                  features = 11,          # number of features to plot (they are selected by weight)
                  cluster_rows = F, cluster_cols = T,
                  show_rownames = T, show_colnames = F,
                  transpose = F,
                  annotation_samples = c("Donor.Specific.Antibodies", "Phenotypes")
)


plot_data_heatmap(model,
                  view = "Blood.mRNAs",         # view of interest
                  factor = 1,             # factor of interest
                  features = 11,          # number of features to plot (they are selected by weight)
                  cluster_rows = F, cluster_cols = T,
                  show_rownames = T, show_colnames = T,
                  transpose = F,
                  annotation_samples = c("Donor.Specific.Antibodies", "Phenotypes")
)




#We figured out that factor 2 is related to C4d lesions
########################################ACCORDING TO Banffs lesions##############
p <- plot_factor(model, 
                 factors = 2,
                 color_by = "Central...Antibody.mediated.rejection...C4d",
                 dot_size = 1,        # change dot size
                 dodge = TRUE,           # dodge points with different colors
                 legend = TRUE,          # remove legend
                 add_violin = TRUE,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

print(p)



plot_data_heatmap(model,
                  view = "Biopsy.mRNAs",         # view of interest
                  factor = 2,             # factor of interest
                  features = 15,          # number of features to plot (they are selected by weight)
                  cluster_rows = F, cluster_cols = T,
                  show_rownames = T, show_colnames = F,
                  transpose = F,
                  annotation_samples = c("Central...Antibody.mediated.rejection...C4d")
)

plot_data_heatmap(model,
                  view = "Biopsy.mRNAs",         # view of interest
                  factor = 2,             # factor of interest
                  features = 15,          # number of features to plot (they are selected by weight)
                  cluster_rows = F, cluster_cols = T,
                  show_rownames = T, show_colnames = T,
                  transpose = F,
                  annotation_samples = c("Central...Antibody.mediated.rejection...C4d")
)


#Scatter plots showing uniquely the features with weight>0.9 from the main compartment (blood for factor 1)
#C4d


plot_data_scatter(model,
                  view = "Urine.mRNAs",         # view of interest
                  factor = 2,             # factor of interest
                  features = 1,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,          # add linear regression
                  color_by = "Central...Antibody.mediated.rejection...C4d"
)+
  scale_color_manual(values=c("0"="#9CAFD2", "1"="#DACDB8","2"="#DACDB8", "3"="#DACDB8", "N"="#8C898D")) +
  scale_fill_manual(values=c("0"="#9CAFD2", "1"="#DACDB8","2"="#DACDB8", "3"="#DACDB8", "N"="#8C898D"))+
  facet_wrap(~ feature, ncol = 5, scales = "free_x", dir="h") 




#We figured out that factor 4 is related to Induction 

p <- plot_factor(model, 
                 factors = 4,
                 color_by = "Type.of.primary.immunosuppression...Induction",
                 dot_size = 1,        # change dot size
                 dodge = TRUE,           # dodge points with different colors
                 legend = TRUE,          # remove legend
                 add_violin = TRUE,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

print(p)



plot_data_heatmap(model,
                  view = "Blood.RNAseq",         # view of interest
                  factor = 4,             # factor of interest
                  features = 11,          # number of features to plot (they are selected by weight)
                  cluster_rows = F, cluster_cols = T,
                  show_rownames = T, show_colnames = T,
                  transpose = F,
                  annotation_samples = c("Type.of.primary.immunosuppression...Induction")
)




sessionInfo()
