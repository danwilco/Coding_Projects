##ASICS metabolomics analysis
##pre-processing, filtering, multivariate analysis
##v2 16/05/22 IJA

#####XCMS processing####
#load required packages
library(xcms)
library(pander)
library(magrittr)
library(SummarizedExperiment)
library(CAMERA)

#set working directory and load files
setwd("D:/ICU keto/21Mar22_RPpos/mzXML") #example folder
mzXML <- dir("Samples", full.names = TRUE, recursive = TRUE)

#load raw data using readMSData method from MSnbase package 
raw_data <- readMSData(files = mzXML, mode = "onDisk")

#get base peak chromatograms
bpis <- chromatogram(raw_data, aggregationFun = "max")
plot(bpis) #visualise base peak chromatograms

#define peakwidth parameters using widest and narrowest peaks 
rtr <- c(8.5*60, 10.5*60)
mzr <- c(293.0981+0.01, 293.0981-0.01) #replace with values from manual inspection

chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr) #plot peaks

#run chromatogram extraction by CentWaveParam method
cwp <- CentWaveParam(peakwidth = c(10, 40),ppm = 20) #use above plot to determine peakwidth
xdata <- findChromPeaks(raw_data, param = cwp)

#align data 
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))

#get base peak chromatograms
bpis_adj <- chromatogram(xdata, aggregationFun = "max")                                 
#does the object have adjusted retention times?
hasAdjustedRtime(xdata) #should return TRUE

##correspondence
#define bandwidth paramters for peak density method
pdp <- PeakDensityParam(sampleGroups = rep(1, length(fileNames(xdata))),
                        minFraction = 0.4, bw = 10)

#perform the correspondence using defined pdp 
xdata <- groupChromPeaks(xdata, param = pdp)

#fill missing peaks using default settings
xdata <- fillChromPeaks(xdata)

#convert XCMSnExp object to xcmsSet for use with CAMERA
xset <- as(xdata, "xcmsSet")
xsa <- annotate(xset, cor_eic_th = 0)
peaklist <- getPeaklist(xsa)
indexNA <- is.na(peaklist)
peaklist[indexNA] <- 0 #replace NAs in data
is.na(peaklist) #sanity check

#write peaklist as csv file for further analysis
write.csv(peaklist, "ICU_rppos.csv")

####Filtering using pmp####
#load required packages
library(pmp)
library(S4Vectors)

#set wd and load files
setwd("D:/ICU keto/21Mar22_RPpos/pmp")
data <- read.csv("data.csv", row.names = 1) #samples in columns and features in rows
dat <- list()
dat$dataMatrix <- data

#replace missing values with NA
dat$dataMatrix[dat$dataMatrix == 0] <- NA

#import metadata and store as list
meta <- read.csv("sample_meta.csv")
dat$sampleMeta <- meta

#make SummarizedExperiment object combining metabolites and metadata
dat <- SummarizedExperiment(assays = list(dat$dataMatrix),
                            colData = DataFrame(dat$sampleMeta))

#check missing samples and find % for reference
sum(is.na(assay(dat)))
sum(is.na(assay(dat)))/length(assay(dat))*100

#filter by QC - peaks must be present in at least 70% of QC samples
dat_qc <- filter_peaks_by_fraction(df = dat,
                                   min_frac = 0.7,
                                   classes = dat$Class,
                                   method = "QC",
                                   qc_label = "QC")
#filter by RSD < 30%
dat_rsd <- filter_peaks_by_rsd(df = dat_qc,
                               max_rsd = 30,
                               classes = dat$Class,
                               qc_label = "QC")

#filter by blank - FC 20 is equivalent to 0.05 blank/QC ratio
rppos_blank <- filter_peaks_by_blank(rppos_filtered,
                                     fold_change = 20,
                                     classes = rppos_filtered$Group,
                                     blank_label = "B",
                                     qc_label = "QC",
                                     remove_peaks = TRUE,
                                     remove_samples = FALSE,
                                     fraction_in_blank = 0)

#data normalisation using probabilistic quotient normalisation
rppos_norm <- pqn_normalisation(df = rppos_blank,
                                classes = rppos_blank$Group,
                                qc_label = "QC")

#missing value imputation
rppos_mv_imputed <- mv_imputation(rppos_norm,
                                  method = "knn")

#data scaling via glog transformation
rppos_glog <- glog_transformation(df = rppos_mv_imputed,
                                  classes = rppos_mv_imputed$Group,
                                  qc_label = "QC")

#write glog file as csv for further analysis
write.csv(rppos_glog, "rppos_glog.csv")

####Multivariate analysis####
#load required packages
library(mixOmics) #for multivariate functions
library(dplyr) #for filtering

#filter data to remove QCs from analysis
glog <- rppos_glog %>%
  filter(Group == "Control"|Group == "Keto")

#example comparison: Keto vs Control on Day 10
glog_d10 <- glog %>%
  filter(Day == "D10") #keep only samples from Day 10

#set up dataframe for PCA and PLS-DA
X_d10 <- glog_d10[,3:ncol(glog_d10)] #X input needs to be numeric only
Y_d10 <- as.factor(glog_d10$Group) #classifying samples based on group (i.e. control vs keto). Y input does not need to be numeric

#preliminary analysis with PCA
pca_d10 <- pca(X_d10, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca_d10) #check most variance is explained in first component

plotIndiv(pca_d10, group = glog_d10$Group, ind.names = FALSE,
          legend = TRUE, ellipse = TRUE, 
          legend.title = "Group", title = "Keto v Control, Day 10") #plots PCA

#PLS-DA
pls_d10 <- plsda(X_d10, Y_d10, ncomp = 10) #PLS-DA function

#tuning parameters for optimisation and validation of model
pls.perf <- perf(pls_d10, validation = "Mfold", folds = 8,
                 progressBar = TRUE, nrepeat = 100)

plot(pls.perf) #view error rates for all components
title(main = "Error rates Day 10 comparison")

#sparse PLS-DA - can be used when error from PLS-DA is high
list.keepX <- c(5:10, seq(20, 100, 10)) #variable selection prior to separation


tune.splsd10 <- tune.splsda(X_d10, Y_d10, ncomp = 6,
                            validation = "Mfold", 
                            folds = 8, progressBar = TRUE, dist = "max.dist",
                            measure = "BER", test.keepX = list.keepX,
                            nrepeat = 100) #validation of model
#ncomp should be component with lowest error from PLS-DA validation

plot(tune.splsd10) #view error rates

select.keepX <- tune.splsd10$choice.keepX[1:5] #replace 5 with number of component with lowest error

#execute final sPLS model
splsd10.final <- splsda(X_d10, Y_d10, ncomp = 5, keepX = select.keepX) #replace ncomp with number of component with lowest error
plotIndiv(splsd10.final, ind.names = FALSE, legend = TRUE,
          title = "Keto vs Control, Day 10", 
          ellipse = TRUE, legend.title = "Group") #plot sPLS-DA

#get VIP scores from separation model
vip_d10 <- as.data.frame(vip(splsd10.final))

vip_d10$Mean <- rowMeans(vip_d10, [c(1:5)]) #get mean VIP score
vip_d10 <- vip_d10 %>%
  filter(Mean > 1) #keep only metabolites with VIP > 1

#export csv file
write.csv(vip_d10, "vip_d10.csv")