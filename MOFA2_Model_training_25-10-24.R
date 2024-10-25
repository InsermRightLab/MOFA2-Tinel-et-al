#MOFA script: 

#With miRNA databases fixed for decimal issues and removal of outliers (Biopsy N=2, Blood N=1)
#With miRNA databases filtered out for miRNAs without any variability (cv=0)
     #Blood miRNAs N=3: hsa-miR-580-001621; hsa-miR-596-001550; hsa-miR-657-001512
     #Biopsy miRNAs N=2: hsa-miR-33a-002135; hsa-miR-596-001550
#With miRNA and urine mRNAs datasets beeing log transformed
#Without miRNA databases filtered out by cv for the top500, after log transform, with %CV specific formula
#With RNAseq, blood mRNA and biopsy mRNA databases being filtered out by cv: 
     #top 5000 highly variable features is retained
#Collinearity between Blood RNAseq and microarray is 1547/5000 mRNAs, i.e. 30.94%

#Without addition of the RegressCovariates option (only in MOFA1, removed from MOFA2 because not performing well)
#With Initialization and run of n=100 models and selection of the best MOFA model

#SystemRequirements Python (>=3), numpy, pandas, h5py, scipy, argparse, sklearn, mofapy2
#---------------------------------------------------------------------------------------

## -----------------------------------------------------------------------------
library(MOFA2)
#library(MOFA)
library(reticulate)


##-----------------------------------------------------------------------------------------
#Databases load:
DF0<-read.csv2(file="Blood_RNAseq_Step1_filter_VSTnorm_top5k.csv", row.names=1, header=T, na.string="NA", sep=";", dec=",")
DF1<-read.csv2(file="blood_mRNA_Step1_ct_top5k_170621.csv", row.names=1, header=T, na.string="NA", sep=";", dec=",")
DF2<-read.csv2(file="Blood_miRNA_Step1_ct_norm_exp_log_all751_130721.csv", row.names=1, header=T, na.string="NA", sep=";", dec=",")
DF3<-read.csv2(file="biopsy_mRNA_Step1_ct_top5k_170621.csv", row.names=1, header=T, na.string="NA", sep=";", dec=",")
DF4<-read.csv2(file="biopsy_miRNA_Step1_ct_norm_exp_log_all752_130721.csv", row.names=1, header=T, na.string="NA", sep=";", dec=",")
DF5<-read.csv2(file="Urine_mRNA_Step1_norm_exp_log_010721.csv", row.names=1, header=T, na.string="NA", sep=";", dec=",")


##-----------------------------------------------------------------------------------------
#All samples must be in the same order:
#Here we choose to sort samples by ascending order from their BIOS ID

DF0 <- DF0[, c("BIOS.03.0002", "BIOS.03.0004", "BIOS.03.0005", "BIOS.03.0007", "BIOS.03.0008",
"BIOS.03.0010", "BIOS.03.0011", "BIOS.03.0013", "BIOS.03.0014", "BIOS.03.0015", "BIOS.03.0017",
"BIOS.03.0018", "BIOS.03.0019", "BIOS.03.0020", "BIOS.03.0021", "BIOS.03.0023", "BIOS.03.0025",
"BIOS.03.0026", "BIOS.03.0028", "BIOS.03.0029", "BIOS.03.0030", "BIOS.03.0031", "BIOS.03.0032",
"BIOS.03.0033", "BIOS.03.0034", "BIOS.03.0035", "BIOS.03.0038", "BIOS.03.0039", "BIOS.03.0040",
"BIOS.03.0041", "BIOS.03.0043", "BIOS.03.0044", "BIOS.03.0045", "BIOS.03.0047", "BIOS.03.0050",
"BIOS.03.0051", "BIOS.03.0052", "BIOS.03.0054", "BIOS.03.0055", "BIOS.03.0056", "BIOS.03.0058",
"BIOS.06.0002", "BIOS.06.0009", "BIOS.06.0016", "BIOS.06.0025", "BIOS.06.0026", "BIOS.06.0029",
"BIOS.06.0032", "BIOS.06.0034", "BIOS.06.0035", "BIOS.06.0037", "BIOS.06.0050", "BIOS.06.0055",
"BIOS.06.0056", "BIOS.06.0060", "BIOS.06.0065", "BIOS.06.0070", "BIOS.06.0087", "BIOS.06.0088",
"BIOS.06.0095", "BIOS.06.0096", "BIOS.06.0098", "BIOS.06.0100", "BIOS.06.0114", "BIOS.06.0130",
"BIOS.06.0136", "BIOS.06.0166", "BIOS.06.0176", "BIOS.06.0183", "BIOS.06.0187", "BIOS.06.0191",
"BIOS.06.0193", "BIOS.06.0196", "BIOS.06.0197", "BIOS.06.0198", "BIOS.06.0199", "BIOS.06.0201",
"BIOS.06.0211", "BIOS.06.0219", "BIOS.06.0222", "BIOS.06.0224", "BIOS.06.0226", "BIOS.06.0237",
"BIOS.06.0238", "BIOS.06.0240", "BIOS.06.0242", "BIOS.06.0250", "BIOS.06.0263", "BIOS.06.0264",
"BIOS.06.0265", "BIOS.06.0277", "BIOS.06.0290", "BIOS.06.0295", "BIOS.06.0311", "BIOS.09.0001",
"BIOS.09.0002", "BIOS.09.0003", "BIOS.09.0004", "BIOS.09.0005", "BIOS.09.0006", "BIOS.09.0007",
"BIOS.09.0008", "BIOS.09.0011", "BIOS.09.0012", "BIOS.09.0013", "BIOS.09.0014", "BIOS.09.0015",
"BIOS.09.0016", "BIOS.09.0018", "BIOS.09.0019", "BIOS.09.0021", "BIOS.09.0022", "BIOS.09.0023",
"BIOS.09.0026", "BIOS.09.0027", "BIOS.09.0029", "BIOS.09.0031", "BIOS.09.0032", "BIOS.09.0033",
"BIOS.09.0035", "BIOS.09.0036", "BIOS.09.0039", "BIOS.09.0041", "BIOS.13.0001", "BIOS.13.0002",
"BIOS.13.0003", "BIOS.13.0005", "BIOS.13.0007", "BIOS.13.0010", "BIOS.13.0012", "BIOS.13.0015")]


DF1 <- DF1[, c("BIOS.03.0002", "BIOS.03.0004", "BIOS.03.0005", "BIOS.03.0007", "BIOS.03.0008",
               "BIOS.03.0010", "BIOS.03.0011", "BIOS.03.0013", "BIOS.03.0014", "BIOS.03.0015", "BIOS.03.0017",
               "BIOS.03.0018", "BIOS.03.0019", "BIOS.03.0020", "BIOS.03.0021", "BIOS.03.0023", "BIOS.03.0025",
               "BIOS.03.0026", "BIOS.03.0028", "BIOS.03.0029", "BIOS.03.0030", "BIOS.03.0031", "BIOS.03.0032",
               "BIOS.03.0033", "BIOS.03.0034", "BIOS.03.0035", "BIOS.03.0038", "BIOS.03.0039", "BIOS.03.0040",
               "BIOS.03.0041", "BIOS.03.0043", "BIOS.03.0044", "BIOS.03.0045", "BIOS.03.0047", "BIOS.03.0050",
               "BIOS.03.0051", "BIOS.03.0052", "BIOS.03.0054", "BIOS.03.0055", "BIOS.03.0056", "BIOS.03.0058",
               "BIOS.06.0002", "BIOS.06.0009", "BIOS.06.0016", "BIOS.06.0025", "BIOS.06.0026", "BIOS.06.0029",
               "BIOS.06.0032", "BIOS.06.0034", "BIOS.06.0035", "BIOS.06.0037", "BIOS.06.0050", "BIOS.06.0055",
               "BIOS.06.0056", "BIOS.06.0060", "BIOS.06.0065", "BIOS.06.0070", "BIOS.06.0087", "BIOS.06.0088",
               "BIOS.06.0095", "BIOS.06.0096", "BIOS.06.0098", "BIOS.06.0100", "BIOS.06.0114", "BIOS.06.0130",
               "BIOS.06.0136", "BIOS.06.0166", "BIOS.06.0176", "BIOS.06.0183", "BIOS.06.0187", "BIOS.06.0191",
               "BIOS.06.0193", "BIOS.06.0196", "BIOS.06.0197", "BIOS.06.0198", "BIOS.06.0199", "BIOS.06.0201",
               "BIOS.06.0211", "BIOS.06.0219", "BIOS.06.0222", "BIOS.06.0224", "BIOS.06.0226", "BIOS.06.0237",
               "BIOS.06.0238", "BIOS.06.0240", "BIOS.06.0242", "BIOS.06.0250", "BIOS.06.0263", "BIOS.06.0264",
               "BIOS.06.0265", "BIOS.06.0277", "BIOS.06.0290", "BIOS.06.0295", "BIOS.06.0311", "BIOS.09.0001",
               "BIOS.09.0002", "BIOS.09.0003", "BIOS.09.0004", "BIOS.09.0005", "BIOS.09.0006", "BIOS.09.0007",
               "BIOS.09.0008", "BIOS.09.0011", "BIOS.09.0012", "BIOS.09.0013", "BIOS.09.0014", "BIOS.09.0015",
               "BIOS.09.0016", "BIOS.09.0018", "BIOS.09.0019", "BIOS.09.0021", "BIOS.09.0022", "BIOS.09.0023",
               "BIOS.09.0026", "BIOS.09.0027", "BIOS.09.0029", "BIOS.09.0031", "BIOS.09.0032", "BIOS.09.0033",
               "BIOS.09.0035", "BIOS.09.0036", "BIOS.09.0039", "BIOS.09.0041", "BIOS.13.0001", "BIOS.13.0002",
               "BIOS.13.0003", "BIOS.13.0005", "BIOS.13.0007", "BIOS.13.0010", "BIOS.13.0012", "BIOS.13.0015")]


DF2 <- DF2[, c("BIOS.03.0002", "BIOS.03.0004", "BIOS.03.0005", "BIOS.03.0007", "BIOS.03.0008",
               "BIOS.03.0010", "BIOS.03.0011", "BIOS.03.0013", "BIOS.03.0014", "BIOS.03.0015", "BIOS.03.0017",
               "BIOS.03.0018", "BIOS.03.0019", "BIOS.03.0020", "BIOS.03.0021", "BIOS.03.0023", "BIOS.03.0025",
               "BIOS.03.0026", "BIOS.03.0028", "BIOS.03.0029", "BIOS.03.0030", "BIOS.03.0031", "BIOS.03.0032",
               "BIOS.03.0033", "BIOS.03.0034", "BIOS.03.0035", "BIOS.03.0038", "BIOS.03.0039", "BIOS.03.0040",
               "BIOS.03.0041", "BIOS.03.0043", "BIOS.03.0044", "BIOS.03.0045", "BIOS.03.0047", "BIOS.03.0050",
               "BIOS.03.0051", "BIOS.03.0052", "BIOS.03.0054", "BIOS.03.0055", "BIOS.03.0056", "BIOS.03.0058",
               "BIOS.06.0002", "BIOS.06.0009", "BIOS.06.0016", "BIOS.06.0025", "BIOS.06.0026", "BIOS.06.0029",
               "BIOS.06.0032", "BIOS.06.0034", "BIOS.06.0035", "BIOS.06.0037", "BIOS.06.0050", "BIOS.06.0055",
               "BIOS.06.0056", "BIOS.06.0060", "BIOS.06.0065", "BIOS.06.0070", "BIOS.06.0087", "BIOS.06.0088",
               "BIOS.06.0095", "BIOS.06.0096", "BIOS.06.0098", "BIOS.06.0100", "BIOS.06.0114", "BIOS.06.0130",
               "BIOS.06.0136", "BIOS.06.0166", "BIOS.06.0176", "BIOS.06.0183", "BIOS.06.0187", "BIOS.06.0191",
               "BIOS.06.0193", "BIOS.06.0196", "BIOS.06.0197", "BIOS.06.0198", "BIOS.06.0199", "BIOS.06.0201",
               "BIOS.06.0211", "BIOS.06.0219", "BIOS.06.0222", "BIOS.06.0224", "BIOS.06.0226", "BIOS.06.0237",
               "BIOS.06.0238", "BIOS.06.0240", "BIOS.06.0242", "BIOS.06.0250", "BIOS.06.0263", "BIOS.06.0264",
               "BIOS.06.0265", "BIOS.06.0277", "BIOS.06.0290", "BIOS.06.0295", "BIOS.06.0311", "BIOS.09.0001",
               "BIOS.09.0002", "BIOS.09.0003", "BIOS.09.0004", "BIOS.09.0005", "BIOS.09.0006", "BIOS.09.0007",
               "BIOS.09.0008", "BIOS.09.0011", "BIOS.09.0012", "BIOS.09.0013", "BIOS.09.0014", "BIOS.09.0015",
               "BIOS.09.0016", "BIOS.09.0018", "BIOS.09.0019", "BIOS.09.0021", "BIOS.09.0022", "BIOS.09.0023",
               "BIOS.09.0026", "BIOS.09.0027", "BIOS.09.0029", "BIOS.09.0031", "BIOS.09.0032", "BIOS.09.0033",
               "BIOS.09.0035", "BIOS.09.0036", "BIOS.09.0039", "BIOS.09.0041", "BIOS.13.0001", "BIOS.13.0002",
               "BIOS.13.0003", "BIOS.13.0005", "BIOS.13.0007", "BIOS.13.0010", "BIOS.13.0012", "BIOS.13.0015")]


DF3 <- DF3[, c("BIOS.03.0002", "BIOS.03.0004", "BIOS.03.0005", "BIOS.03.0007", "BIOS.03.0008",
               "BIOS.03.0010", "BIOS.03.0011", "BIOS.03.0013", "BIOS.03.0014", "BIOS.03.0015", "BIOS.03.0017",
               "BIOS.03.0018", "BIOS.03.0019", "BIOS.03.0020", "BIOS.03.0021", "BIOS.03.0023", "BIOS.03.0025",
               "BIOS.03.0026", "BIOS.03.0028", "BIOS.03.0029", "BIOS.03.0030", "BIOS.03.0031", "BIOS.03.0032",
               "BIOS.03.0033", "BIOS.03.0034", "BIOS.03.0035", "BIOS.03.0038", "BIOS.03.0039", "BIOS.03.0040",
               "BIOS.03.0041", "BIOS.03.0043", "BIOS.03.0044", "BIOS.03.0045", "BIOS.03.0047", "BIOS.03.0050",
               "BIOS.03.0051", "BIOS.03.0052", "BIOS.03.0054", "BIOS.03.0055", "BIOS.03.0056", "BIOS.03.0058",
               "BIOS.06.0002", "BIOS.06.0009", "BIOS.06.0016", "BIOS.06.0025", "BIOS.06.0026", "BIOS.06.0029",
               "BIOS.06.0032", "BIOS.06.0034", "BIOS.06.0035", "BIOS.06.0037", "BIOS.06.0050", "BIOS.06.0055",
               "BIOS.06.0056", "BIOS.06.0060", "BIOS.06.0065", "BIOS.06.0070", "BIOS.06.0087", "BIOS.06.0088",
               "BIOS.06.0095", "BIOS.06.0096", "BIOS.06.0098", "BIOS.06.0100", "BIOS.06.0114", "BIOS.06.0130",
               "BIOS.06.0136", "BIOS.06.0166", "BIOS.06.0176", "BIOS.06.0183", "BIOS.06.0187", "BIOS.06.0191",
               "BIOS.06.0193", "BIOS.06.0196", "BIOS.06.0197", "BIOS.06.0198", "BIOS.06.0199", "BIOS.06.0201",
               "BIOS.06.0211", "BIOS.06.0219", "BIOS.06.0222", "BIOS.06.0224", "BIOS.06.0226", "BIOS.06.0237",
               "BIOS.06.0238", "BIOS.06.0240", "BIOS.06.0242", "BIOS.06.0250", "BIOS.06.0263", "BIOS.06.0264",
               "BIOS.06.0265", "BIOS.06.0277", "BIOS.06.0290", "BIOS.06.0295", "BIOS.06.0311", "BIOS.09.0001",
               "BIOS.09.0002", "BIOS.09.0003", "BIOS.09.0004", "BIOS.09.0005", "BIOS.09.0006", "BIOS.09.0007",
               "BIOS.09.0008", "BIOS.09.0011", "BIOS.09.0012", "BIOS.09.0013", "BIOS.09.0014", "BIOS.09.0015",
               "BIOS.09.0016", "BIOS.09.0018", "BIOS.09.0019", "BIOS.09.0021", "BIOS.09.0022", "BIOS.09.0023",
               "BIOS.09.0026", "BIOS.09.0027", "BIOS.09.0029", "BIOS.09.0031", "BIOS.09.0032", "BIOS.09.0033",
               "BIOS.09.0035", "BIOS.09.0036", "BIOS.09.0039", "BIOS.09.0041", "BIOS.13.0001", "BIOS.13.0002",
               "BIOS.13.0003", "BIOS.13.0005", "BIOS.13.0007", "BIOS.13.0010", "BIOS.13.0012", "BIOS.13.0015")]


DF4 <- DF4[, c("BIOS.03.0002", "BIOS.03.0004", "BIOS.03.0005", "BIOS.03.0007", "BIOS.03.0008",
               "BIOS.03.0010", "BIOS.03.0011", "BIOS.03.0013", "BIOS.03.0014", "BIOS.03.0015", "BIOS.03.0017",
               "BIOS.03.0018", "BIOS.03.0019", "BIOS.03.0020", "BIOS.03.0021", "BIOS.03.0023", "BIOS.03.0025",
               "BIOS.03.0026", "BIOS.03.0028", "BIOS.03.0029", "BIOS.03.0030", "BIOS.03.0031", "BIOS.03.0032",
               "BIOS.03.0033", "BIOS.03.0034", "BIOS.03.0035", "BIOS.03.0038", "BIOS.03.0039", "BIOS.03.0040",
               "BIOS.03.0041", "BIOS.03.0043", "BIOS.03.0044", "BIOS.03.0045", "BIOS.03.0047", "BIOS.03.0050",
               "BIOS.03.0051", "BIOS.03.0052", "BIOS.03.0054", "BIOS.03.0055", "BIOS.03.0056", "BIOS.03.0058",
               "BIOS.06.0002", "BIOS.06.0009", "BIOS.06.0016", "BIOS.06.0025", "BIOS.06.0026", "BIOS.06.0029",
               "BIOS.06.0032", "BIOS.06.0034", "BIOS.06.0035", "BIOS.06.0037", "BIOS.06.0050", "BIOS.06.0055",
               "BIOS.06.0056", "BIOS.06.0060", "BIOS.06.0065", "BIOS.06.0070", "BIOS.06.0087", "BIOS.06.0088",
               "BIOS.06.0095", "BIOS.06.0096", "BIOS.06.0098", "BIOS.06.0100", "BIOS.06.0114", "BIOS.06.0130",
               "BIOS.06.0136", "BIOS.06.0166", "BIOS.06.0176", "BIOS.06.0183", "BIOS.06.0187", "BIOS.06.0191",
               "BIOS.06.0193", "BIOS.06.0196", "BIOS.06.0197", "BIOS.06.0198", "BIOS.06.0199", "BIOS.06.0201",
               "BIOS.06.0211", "BIOS.06.0219", "BIOS.06.0222", "BIOS.06.0224", "BIOS.06.0226", "BIOS.06.0237",
               "BIOS.06.0238", "BIOS.06.0240", "BIOS.06.0242", "BIOS.06.0250", "BIOS.06.0263", "BIOS.06.0264",
               "BIOS.06.0265", "BIOS.06.0277", "BIOS.06.0290", "BIOS.06.0295", "BIOS.06.0311", "BIOS.09.0001",
               "BIOS.09.0002", "BIOS.09.0003", "BIOS.09.0004", "BIOS.09.0005", "BIOS.09.0006", "BIOS.09.0007",
               "BIOS.09.0008", "BIOS.09.0011", "BIOS.09.0012", "BIOS.09.0013", "BIOS.09.0014", "BIOS.09.0015",
               "BIOS.09.0016", "BIOS.09.0018", "BIOS.09.0019", "BIOS.09.0021", "BIOS.09.0022", "BIOS.09.0023",
               "BIOS.09.0026", "BIOS.09.0027", "BIOS.09.0029", "BIOS.09.0031", "BIOS.09.0032", "BIOS.09.0033",
               "BIOS.09.0035", "BIOS.09.0036", "BIOS.09.0039", "BIOS.09.0041", "BIOS.13.0001", "BIOS.13.0002",
               "BIOS.13.0003", "BIOS.13.0005", "BIOS.13.0007", "BIOS.13.0010", "BIOS.13.0012", "BIOS.13.0015")]


DF5 <- DF5[, c("BIOS.03.0002", "BIOS.03.0004", "BIOS.03.0005", "BIOS.03.0007", "BIOS.03.0008",
               "BIOS.03.0010", "BIOS.03.0011", "BIOS.03.0013", "BIOS.03.0014", "BIOS.03.0015", "BIOS.03.0017",
               "BIOS.03.0018", "BIOS.03.0019", "BIOS.03.0020", "BIOS.03.0021", "BIOS.03.0023", "BIOS.03.0025",
               "BIOS.03.0026", "BIOS.03.0028", "BIOS.03.0029", "BIOS.03.0030", "BIOS.03.0031", "BIOS.03.0032",
               "BIOS.03.0033", "BIOS.03.0034", "BIOS.03.0035", "BIOS.03.0038", "BIOS.03.0039", "BIOS.03.0040",
               "BIOS.03.0041", "BIOS.03.0043", "BIOS.03.0044", "BIOS.03.0045", "BIOS.03.0047", "BIOS.03.0050",
               "BIOS.03.0051", "BIOS.03.0052", "BIOS.03.0054", "BIOS.03.0055", "BIOS.03.0056", "BIOS.03.0058",
               "BIOS.06.0002", "BIOS.06.0009", "BIOS.06.0016", "BIOS.06.0025", "BIOS.06.0026", "BIOS.06.0029",
               "BIOS.06.0032", "BIOS.06.0034", "BIOS.06.0035", "BIOS.06.0037", "BIOS.06.0050", "BIOS.06.0055",
               "BIOS.06.0056", "BIOS.06.0060", "BIOS.06.0065", "BIOS.06.0070", "BIOS.06.0087", "BIOS.06.0088",
               "BIOS.06.0095", "BIOS.06.0096", "BIOS.06.0098", "BIOS.06.0100", "BIOS.06.0114", "BIOS.06.0130",
               "BIOS.06.0136", "BIOS.06.0166", "BIOS.06.0176", "BIOS.06.0183", "BIOS.06.0187", "BIOS.06.0191",
               "BIOS.06.0193", "BIOS.06.0196", "BIOS.06.0197", "BIOS.06.0198", "BIOS.06.0199", "BIOS.06.0201",
               "BIOS.06.0211", "BIOS.06.0219", "BIOS.06.0222", "BIOS.06.0224", "BIOS.06.0226", "BIOS.06.0237",
               "BIOS.06.0238", "BIOS.06.0240", "BIOS.06.0242", "BIOS.06.0250", "BIOS.06.0263", "BIOS.06.0264",
               "BIOS.06.0265", "BIOS.06.0277", "BIOS.06.0290", "BIOS.06.0295", "BIOS.06.0311", "BIOS.09.0001",
               "BIOS.09.0002", "BIOS.09.0003", "BIOS.09.0004", "BIOS.09.0005", "BIOS.09.0006", "BIOS.09.0007",
               "BIOS.09.0008", "BIOS.09.0011", "BIOS.09.0012", "BIOS.09.0013", "BIOS.09.0014", "BIOS.09.0015",
               "BIOS.09.0016", "BIOS.09.0018", "BIOS.09.0019", "BIOS.09.0021", "BIOS.09.0022", "BIOS.09.0023",
               "BIOS.09.0026", "BIOS.09.0027", "BIOS.09.0029", "BIOS.09.0031", "BIOS.09.0032", "BIOS.09.0033",
               "BIOS.09.0035", "BIOS.09.0036", "BIOS.09.0039", "BIOS.09.0041", "BIOS.13.0001", "BIOS.13.0002",
               "BIOS.13.0003", "BIOS.13.0005", "BIOS.13.0007", "BIOS.13.0010", "BIOS.13.0012", "BIOS.13.0015")]


#------------------------------------------------------------------------------
#To create a MOFA object, database need to be a list of matrices: 
#DFgroup <-as.matrix(DFgroup)
DF0 <-as.matrix(DF0)
DF1 <-as.matrix(DF1)
DF2 <-as.matrix(DF2)
DF3 <-as.matrix(DF3)
DF4 <-as.matrix(DF4)
DF5 <-as.matrix(DF5)
DFmofa <- list(DF0, DF1, DF2, DF3, DF4, DF5)


#To check that all data are numeric (mandatory for the creation of a MOFA object):
sapply(DF0, class)
sapply(DF1, class)
sapply(DF2, class)
sapply(DF3, class)
sapply(DF4, class)
sapply(DF5, class)


## To create the MOFA object:
N = ncol(DFmofa[[1]])
MOFAobject <- create_mofa(DFmofa, groups=NULL)
views_names(MOFAobject) <- c("Top 5k Blood RNAseq", "Top 5k Blood mRNAs", "All Blood miRNAs", "Top 5k Biopsy mRNAs", "All Biopsy miRNAs", "All Urine mRNAs")
views_names(MOFAobject) <- c("F", "E", "D", "C", "B", "A")

## -----------------------------------------------------------------------------
library(ggplot2)
#Function to do a tile plot showing the missing value structure of the input data:
plot_data_overview(MOFAobject)
plot_data_overview(MOFAobject, colors=c("dodgerblue2", "turquoise4", "darkturquoise", "darkred", "coral1", "orange"))

## -----------------------------------------------------------------------------
#filepath <- system.file("extdata", "test_data.RData", package = "MOFA2")
#load(filepath)
#head(dt)

## -----------------------------------------------------------------------------
print(MOFAobject)

## -----------------------------------------------------------------------------
data_opts <- get_default_data_options(MOFAobject)
#scale_views: if views have different ranges/variances, it is good practice to scale each view
#to unit variance. Default is FALSE
data_opts$scale_views <- TRUE
head(data_opts)

## -----------------------------------------------------------------------------
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10
head(model_opts)

## -----------------------------------------------------------------------------
train_opts <- get_default_training_options(MOFAobject)
train_opts$drop_factor_threshold <- 0.02 
#0.01 implies that factors explaining less than 1% of variance in each view will be dropped
head(train_opts)



## -----------------------------------------------------------------------------
#To run the model n times:
n_inits <- 100

MOFAlist <- lapply(seq_len(n_inits), function(it) {
        train_opts$seed <- 1234 + it  
        MOFAobject <- prepare_mofa (
                object = MOFAobject,
                data_options = data_opts,
                model_options = model_opts,
                training_options = train_opts 
        )
        run_mofa(MOFAobject, use_basilisk = TRUE)
})






## -----------------------------------------------------------------------------
#Compare different random inits and select the best model
#Different objects of MOFA are compared in terms of the final value of the ELBO statistics. 
#For model selection the model with the highest ELBO value is selected.
#compare_elbo(models, log = FALSE, return_data = FALSE)
# Using an existing trained model on simulated data
#file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#model1 <- load_model(file)
#model2 <- load_model(file)
#compare_elbo(list(model1,model2))

compare_elbo(MOFAlist)
compare_elbo(MOFAlist, log=FALSE)
#...Plotting the absolute value of the ELBO for every model (the smaller the better)
#To return a data.frame with the ELBO values per model:
#Plotting the log2 of the negative of the ELBO (the higher the better)
compare_elbo(MOFAlist, log=TRUE, return_data=TRUE)
compare_elbo(MOFAlist, log=FALSE, return_data=TRUE)

## -----------------------------------------------------------------------------
#Get an overview of how robust the factors are between different model:
#Different MOFA objects are compared in terms of correlation between their factors
#Using an existing trained model on simulated data
#file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#model1 <- load_model(file)
#model2 <- load_model(file)
#compare_factors(list(model1,model2))
compare_factors(MOFAlist)


## -----------------------------------------------------------------------------
#Choose the model with the best ELBO value as is done by select_model:
MOFAobject.trained <- select_model(MOFAlist, plot =FALSE)
#...if plot=TRUE: Plotting the absolute value of the ELBO for every model (the smaller the better)  
MOFAobject.trained
#Extract the value of the ELBO statistics after model training:
get_elbo(MOFAobject.trained)
#[1] -1493193   #Model 61 according to the absolute ELBO value (the smaller)


## -----------------------------------------------------------------------------

#outfile = file.path(getwd(),"model.hdf5_090821_1inits")
#MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)
#NB: Connecting to the mofapy2 package using basilisk. 
#Set "use_basilisk" to FALSE if you prefer to manually set the python binary using 'reticulate'.

saveRDS(MOFAobject.trained, file = "model.hdf5_300822_100inits_Model61")
saveRDS(MOFAlist, file = "MOFAlist_300822_100inits_all")

#--------------------------------------------------------------------------------


col3 <- colorRampPalette(c("blue", "white", "red")) 
col0 <- colorRampPalette(c("cyan", "white", "#00007F")) 
plot_factor_cor(MOFAobject.trained, col = col0(20), "pearson")
#Correlation plot looks much better (uncorrelated factors) after implementation of the function "scale.views"


## -----------------------------------------------------------------------------
sessionInfo()



