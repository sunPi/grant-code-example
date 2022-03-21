#-------------------------------------------------------------------------------
# Data Analysis Pipeline - Complete on patients, explants, cell lines
#-------------------------------------------------------------------------------
# XGBOOST Parameters
# - verbosity: Verbosity of printing messages. Valid values of 0 (silent), 1 (warning), 2 (info), and 3 (debug).
# - booster [default= gbtree ]Which booster to use. Can be gbtree, gblinear or dart; gbtree and dart use tree based models while gblinear uses linear fun
infolder <- "./outputs/r-objects/"
epochs           <- 70
bst              <- "gbtree"
e                <- 0.3
g                <- 0
md               <- 9
mcwt             <- 1
ss               <- 1
csbt             <- 1
verbose          <- 1

Sys.setenv(CUDA="11.1")

#------------------ Requirements, dependencies and libraries
pkgs <- c('pacman','docopt','dplyr','xgboost', 
          'data.table', 'caret', 'ggplot2', 
          'tidyr','stringr', 'readxl', 'here', 
          'ggplot2', 'caret', 'igraph', 'DiagrammeR', 
          'Ckmeans.1d.dp', 'pROC')

#------------------ Gather citations/meta data
inf <- sessionInfo()
citation("pROC")


suppressMessages(if (!require("BiocManager", character.only = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install()
  } else {
  ipkgs <- sapply(pkgs, function(...) require(..., character.only = TRUE))
  if (any(!ipkgs)) {
    BiocManager::install(pkgs[!ipkgs])
    install.packages(pkgs[!ipkgs])
  } else {
    message("\n\nCool! your machine has everything is needed.\n\n")
  }
  })


print("Loading required packages...")
library(pacman)
pacman::p_load(pkgs, install = TRUE, character.only = TRUE)
pacman::p_loaded()



# "PARPi Drug Screening Pipeline - Patients, explants, cell lines
# 
# Usage: explants.R --path=<value> --infolder=<folder> --epochs=<value> --max_depth=<value> --eta=<value> --gamma=<value> --colsample_bytree=<value> --min_child_weight=<value> --subsample=<value> --verbose=<value>
# 
# Options:
#   -h --help                  Show this screen.
#   --path=<value>
#   --infolder=<folder>        Folder where the outputs are placed.
#   --epochs=<value>           Number of training epochs for the model.
#   --max_depth=<value>
#   --eta=<value>
#   --gamma=<value>
#   --colsample_bytree=<value>
#   --min_child_weight=<value>
#   --subsample=<value>
#   --verbose=<value>          If set to 1 prints all messages.
#   --version                  
# "-> doc
# 
# arguments <- docopt(doc, quoted_args = TRUE, help = TRUE)
# print(arguments)

#------------------ Functions
source("./fun/util.R")
source("./fun/global.R")

#------------------ Load dataset and parameters into R environment 
# infolder          <- arguments$infolder
# path              <- arguments$path 
# epochs            <- as.numeric(arguments$epochs)
# md                <- as.numeric(arguments$max_depth)
# e                 <- as.numeric(arguments$eta)           
# gm                <- as.numeric(arguments$gamma)
# csbt              <- as.numeric(arguments$colsample_bytree)
# mcwt              <- as.numeric(arguments$min_child_weight)
# ss                <- as.numeric(arguments$subsample)   
# verbose           <- as.numeric(arguments$verbose)
# --- AUX

paths <- list(datasets   = here('datasets/'),
              outputs    = here('outputs/'),
              models     = here('outputs', 'models/'),
              model_eval = here('outputs', 'models/', 'model_eval/'),
              pdf        = here('outputs', 'pdf/'),
              grant      = here('outputs', 'pdf', 'grant_proposal/'),
              robj       = here('outputs', 'r-objects/')
              )
#saveRDS(object = paths, file = paste0(paths$robj,'dir.RDS'))

#------------------ Recode the data into numeric variables
# Processing Patients Dataset
patients  <- as.data.frame(read_excel("./datasets/DDR Gene List.xlsx", sheet = 3))
rcd <- patients[,-c(1:4)]
rcd[is.na(rcd)] <- 0
meta <- patients[,c(1:4)]
for (i in names(rcd)) {
  rcd[,i] <- recode(rcd[,i], "deletion" = -2, "loss" = -1, "wild-type" = 0, "gain" = 1, "amplification" = 2, .default = 0, .missing = 0)
}
rcd <- apply(rcd, 2, as.integer)
patients <- as.data.frame(cbind(meta,rcd))

# explants <- as.data.frame(read_excel("./datasets/DDR Gene List.xlsx", sheet = 4)) 
# rcd <- explants[,-c(1:4)]
# meta <- explants[,c(1:4)]
# for (i in names(rcd)) {
#   rcd[,i] <- recode(rcd[,i], "Sd" = -2, "Cd" = -1, "wild-type" = 0, "C" = 1, "S" = 2, .default = 0, .missing = NULL)
# }
# explants <- as.data.frame(cbind(meta,rcd))

# Get the explants data 
explants <- readRDS("../../../global-robj/all.explants.RDS")
# Get the cell lines data
cell.lines <- readRDS("./outputs/r-objects/dataframe_v2.RDS")
# All CL w/o Talazoparib
cell.lines$public.exp$all.exp.cl <- cell.lines$public.exp$all.exp.cl[-which(cell.lines$public.exp$all.exp.cl$DrugName == "Talazoparib"),]

#-------------- PARP data analysis
pipe.data <- list(patients       = patients,
                  explants       = explants,
                  cell.lines = rbind(cell.lines$public.exp$all.exp.cl,
                                     cell.lines$public.exp$all.public.cl)
)
write.csv(pipe.data$patients, "./outputs/patient_dataframe.csv", row.names = FALSE)
write.csv(pipe.data$explants, "./outputs/explants_dataframe.csv", row.names = FALSE)
write.csv(pipe.data$cell.lines, "./outputs/cell.lines_dataframe.csv", row.names = FALSE)
unique(pipe.data$cell.lines$SampleName)
#writexl::write_xlsx(pipe.data, "./datasets/pipedata.xlsx", col_names = TRUE)

#------------------ Prepare the gene set variables
geneset <- list(patients = colnames(pipe.data$patients[,5:ncol(pipe.data$patients)]), 
                explants = colnames(pipe.data$explants[,5:ncol(pipe.data$explants)]),
                cell.lines = colnames(pipe.data$cell.lines[,5:ncol(pipe.data$cell.lines)])
                )
saveRDS(geneset$patients, "./outputs/r-objects/pat.geneset")
write.csv(geneset$patients, "./outputs/pat_geneset.csv", row.names = FALSE)
saveRDS(geneset$explants, "./outputs/r-objects/exp+cl.geneset")
write.csv(geneset$explants, "./outputs/exp+cl_geneset.csv", row.names = FALSE)

#------------------ Test normalisation
# pipe.data$gh2ax$ruca.c.count.mean$ZScore <- nor.min.max(pipe.data$gh2ax$ruca.c.count.mean$ZScore)

#------------------ Sample the dataset (70%-30%)
set.seed(376)
qty <- 0.9
#------------------ Set the parameters
params    <- list(booster         = bst,
                  objective        = "reg:squarederror",
                  eta              = e,
                  gamma            = g,
                  min_child_weight = mcwt,
                  max_depth        = md,
                  colsample_bytree = csbt,
                  subsample        = ss)

# Retrain the model using binary-transformed labels

i <- 1
for(z in pipe.data$patients$ZScore){
  pipe.data$patients$ZScore[i] <- ToBinary(z)
  i <- i + 1
}

i <- 1
for(z in pipe.data$explants$ZScore){
  pipe.data$explants$ZScore[i] <- ToBinary(z)
  i <- i + 1
}

i <- 1
for(z in pipe.data$cell.lines$ZScore){
  pipe.data$cell.lines$ZScore[i] <- ToBinary(z)
  i <- i + 1
}
params$objective <- "binary:logistic"
epochs <- 70
params$eta <- 0.3
source("./pipeline/functions.R")
source("./fun/global.R")

i <- 1
complete.response <- list()
for(data in pipe.data){
  print(paste0("Analysing ", names(pipe.data[i])))
  complete.response <- append(complete.response, list(trainSensitivity(x = data,
                                                     qty = qty,
                                                     params = params,
                                                     epochs = epochs,
                                                     geneset = geneset[[i]],
                                                     save_name_dir = paste0(paths$pdf, names(pipe.data[i]),'.response.pdf'), # REPLACE THIS WITH THE CORRECT FILE/ANALYSIS NAME
                                                     verbose = 1))
  )
  i <- i + 1
}

names(complete.response) <- names(pipe.data)
cat("Succesfully trained",length(complete),"models named as",names(pipe.data))
cat("Model summary and SHAP metrics saved under ", paths$pdf)


metrics <- c("accuracy" = complete.response$cell.lines$metrics$accuracy, 
        "preision" = complete.response$cell.lines$metrics$precision, 
        "recall"   = complete.response$cell.lines$metrics$recall, 
        "F1"       = complete.response$cell.lines$metrics$F1) %>% 
  as.data.frame() %>% 
  t() %>%
  as.data.frame()

overall <- complete.response$cell.lines$metrics$confusion.matrix$overall %>% 
  as.data.frame() %>% 
  t() %>%
  as.data.frame()

overall

byclass <- complete.response$cell.lines$metrics$confusion.matrix$byClass %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame()

writexl::write_xlsx(list(metrics,overall,byclass),path = "./cl.metrics.xlsx", col_names = TRUE)

response   <- complete.response$cell.lines$meso$test.label
prediction <- complete.response$cell.lines$preds 

roc_obj <- roc(response = response, prediction)
auc(roc_obj)
## Area under the curve
roc_df <- data.frame(
  TPR=rev(roc_obj$sensitivities), 
  FPR=rev(1 - roc_obj$specificities), 
  labels=roc_obj$response, 
  scores=roc_obj$predictor)


im1 <- complete.response$patients$explainer$imatrix
im2 <- complete.response$cell.lines$explainer$imatrix
im1 <- common_genes(im1, im2)$im1
im2 <- common_genes(im1, im2)$im2
col_names <- c("Patients_HRD_Genes", "Patients_SHAP_Values", "Cell_Lines_HRD_Genes", "Cell_Lines_SHAP_Values") 

im.corr.response <- build_im_corr(im1, im2, col_names)
im.corr.response <- log_trans(im.corr.response)

hrd.response <- explore_2DSHAP(imatrix.corr    = im.corr.response,
                                 col_names     = col_names,
                                 title         = "Patients SHAP Values vs Cell Line SHAP Values on selected HRD Genes",
                                 log           = TRUE,
                                 do_lm         = TRUE,
                                 level         = 0.99,
                                 save_file     = FALSE,
                                 dir_save_name = paste0(paths$grant,"RESPONSE_2D.Patients.vs.Cell.Lines.HRD.matched"),
                                 ftype         = "pdf")

# Test the fitted model

complete$patients$preds
complete$patients$meso$learning.set$test$ZScore

complete$cell.lines$preds
complete$cell.lines$meso$learning.set$test$ZScore

complete$patients$metrics
complete$cell.lines$metrics

complete$patients$meso$train.matrix
complete$cell.lines$meso$train.matrix

complete$patients$meso$test.matrix
complete$cell.lines$meso$test.matrix
