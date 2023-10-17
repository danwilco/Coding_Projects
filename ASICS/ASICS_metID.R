##metID ASICS
##16/05/2022 v2
##IJA

#load packages
library(metID)
library(dplyr)

#set working directory to where your files are stored
setwd("C:/Users/mzxia1/Documents/ITU/RP pos/metID d10")

#load databases from folder - databases are publically available
param1 <- identify_metabolites_params(
  ms1.match.ppm = 5, 
  rt.match.tol = 10000,
  polarity = "positive", #switch to negative depending on data
  ce = "all",
  column = "rp", #switch to hilic depending on column
  total.score.tol = 0.5,
  database = "hmdbMS1Database0.0.1"
) #HMDB MS1 level

param2 <- identify_metabolites_params(
  ms1.match.ppm = 5, 
  rt.match.tol = 10000,
  polarity = "positive",
  ce = "all",
  column = "rp",
  total.score.tol = 0.5,
  database = "keggMS1Database_1.0"
) #kegg MS1 level

param3 <- identify_metabolites_params(
  ms1.match.ppm = 5,
  rt.match.tol = 10000,
  polarity = "positive",
  ce = "all",
  column = "rp",
  total.score.tol = 0.5,
  database = "hmdbDatabase0.0.2"
) #HMDB MS2 level

param4 <- identify_metabolites_params(
  ms1.match.ppm = 5,
  rt.match.tol = 10000,
  polarity = "positive",
  ce = "all",
  column = "rp",
  total.score.tol = 0.5,
  database = "massbankDatabase0.0.2"
) #MassBank

param5 <- identify_metabolites_params(
  ms1.match.ppm = 5,
  rt.match.tol = 10000,
  polarity = "positive",
  ce = "all",
  column = "rp",
  total.score.tol = 0.5,
  database = "monaDatabase0.0.2"
) #MoNA

param6 <- identify_metabolites_params(
  ms1.match.ppm = 5,
  rt.match.tol = 10000,
  polarity = "positive",
  ce = "all",
  column = "rp",
  total.score.tol = 0.5,
  database = "orbitrapDatabase0.0.1"
)
#rt match tol set to 10000 to exclude from analysis - if using inhouse data of standards, could include RT information

#run database comparison
result <- identify_metabolite_all(
  ms1.data = "ms1_d10.csv",
  ms2.data = "PooledQC_MSn_3_neg.mzXML",
  parameter = c(param1, param2, param3, param4, param5, param6),
  path = "."
)

#get annotation table
annotation_table1 <- get_identification_table(result[[1]], 
                                              candidate.num = 1,
                                              type = "new")
annotation_table2 <- get_identification_table(result[[2]],
                                              candidate.num = 1,
                                              type = "new")
annotation_table3 <- get_identification_table(result[[3]],
                                              candidate.num = 1,
                                              type = "new")
annotation_table4 <- get_identification_table(result[[4]],
                                              candidate.num = 1,
                                              type = "new")
annotation_table5 <- get_identification_table(result[[5]],
                                              candidate.num = 1,
                                              type = "new")
annotation_table6 <-get_identification_table(result[[6]],
                                             candidate.num = 1,
                                             type = "new")
#if no annotations, environment will say NULL

#filter and combine tables
annotation_table1 <- annotation_table1 %>%
  filter(!is.na(Compound.name)) #removes unidentified compounds, do for each table that has annotations

#if one peak has multiple annotations, we only want the one with highest confidence
annotation_table2 <- 
  annotation_table2 %>% 
  filter(!(name %in% annotation_table1$name))

#combines into one table
annotation_table_all <- rbind(annotation_table1, 
                              annotation_table2) #include as many as have annotations

#export table to csv file
write.csv(annotation_table_all, "annotation_table_rppos.csv")
