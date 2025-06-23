
library(lme4) ## mixed models
library(refund) ## fpca.face 
library(dplyr) ## organize lapply results
library(mgcv) ## smoothing in step 2
library(mvtnorm) ## joint CI
library(parallel) ## mcapply
library(Rfast)
library(caret)
library(data.table)
library(sanic)
library(SuperGauss)
library(SimCorMultRes)

simParams <- expand.grid(c(50), # 100 30, 
                         c("fastkfold"), #c(TRUE, "kfold", "fastCV"), # "fastBoot"
                         c(TRUE), # TRUE, FALSE
                         c(FALSE),  # FALSE, TRUE, "fastBoot" # FALSE, TRUE, "fastBoot"
                         c(50), # 15, 30, 100
                         c(100), # L
                         c(0.75) # 0.25, 0.75
)
# currently trying: does lowering rho = 0.75 to rho = 0.5 improve coverage
colnames(simParams) <- c("n", "cv", "rho.smooth", "var.type", "Ni", "L", "rho")

cluserInd <- FALSE

# paths
runNum <- 1
save.wd <- "/Users/loewingergc/Desktop/NIMH Research/functional_GEE/fastFGEE/Vignette"
wd <- "/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final/Resubmission/Code/Utility Functions/"
iter <- 1
data_path <-"/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/science_paper" # path to original Matlab files for data
source("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/photometry_sim_fLME_fn_multi.R") # simulations
source(paste0(wd, "fui_fast.R")) # function code
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_fGLMM.R')
source("/Users/loewingergc/Desktop/NIMH Research/Photometry/code/existing_funs.R") # for time align functions
source("/Users/loewingergc/Desktop/NIMH Research/functional_GEE/gee1step-main/R/fun.gee1step_joint.R")
source("/Users/loewingergc/Desktop/NIMH Research/functional_GEE/sims/sim_func_data_Li.R")
source("/Users/loewingergc/Desktop/NIMH Research/functional_GEE/sims/binary_Var_splines_mgcv.R")
wd4 <- "/Users/loewingergc/Desktop/NIMH Research/Causal/msm_hr/Pen_Spline/Sims/Code/"

wd3 <- "/Users/loewingergc/Desktop/NIMH Research/Causal/msm_hr/Pen_Spline/Sims/Code/"
source(paste0(wd3,"sandwich_normal_pen_fGEE.R"))
wd <- "/Users/loewingergc/Desktop/NIMH Research/functional_GEE/gee1step-main/sources/"
#colnames(simParams) <- c("fix_smooth", "resid_var_indep", "resid_var_subj", "n", "beta_mult", "cv", "weight")


# simulation parameters
reps <- 300 # replicates of simulations
boots <- 5000 # for naive bootstrap permutation test
nknots <- 10 # 3 and 45 works well
n_sim <- simParams[runNum, 1] #c(100, 250, 10)[runNum]#simParams[runNum, 4]   # 10, 20, 30, 50, 70, 
N_i <- simParams[runNum, 5]
cv <- simParams[runNum, 2] # c(TRUE, "fastBoot", "fastCV") 
rho.smooth <- simParams[runNum, 3] #c(TRUE, FALSE) 
var.type <- simParams[runNum, 4] #"boot" # TRUE # "boot # FALSE 
bs <- "ps"
comparison_methods <- TRUE # whether to fit other models besides fLME
L <- simParams$L[runNum] #200 # originally 100
p <- 3 # number of covariates
nknots_min <- 10
z11 <- 5
z12 <- 3 # 1.5
sigma_sq <- 10 # for AR1 this needs to be small enough not to dominate AR1 part
n_delta <- "n.5" #10 # variability in sample size
cov.type <- "ar1"
comb_pen <- TRUE
rho <- simParams[runNum, 7]
rho.vec <- rep(rho, L) # ar9 uses this
long_rho <- rho
family <- "binomial"
divFactor <- 3

##################
# simulate data
##################
set.seed(iter)
beta_idx <- 2:3 # this corresponds to slope effects
n_delta <- round(N_i/2)

if(n_delta == "n.5"){
  n_delta <- floor(N_i/2)
}

set.seed(iter)
dat_sim <-  GenerateData(Nsubj=n_sim, 
                         numFunctPoints = L, 
                         min_visit=N_i,
                         max_visit=N_i,
                         numLongiPoints = N_i, 
                         sigma_sq = sigma_sq, 
                         sigma_z11 = z11, sigma_z12 = z12, 
                         sigma_z21 = 2, sigma_z22 = 1,
                         corstr = cov.type,
                         divFactor = divFactor, # divide betas by this for binomial to avoid all 1s
                         divVec = NULL,
                         rho.func = long_rho, # within function correlation
                         family = family,
                         true.cov = FALSE,
                         rho.vec = rho.vec)

Y <- as.data.frame(dat_sim$data$Y)
colnames(Y) <- paste0("Y_", 1:ncol(Y))
L <- ncol(Y)
X <- model.matrix(~X1 + X2, data = data.frame(dat_sim$data$Cov))
d = data.frame(Y = I(Y), 
               ID = dat_sim$data$subjID, 
               dat_sim$data$Cov,
               time = dat_sim$time)

setwd(save.wd)
saveRDS(d, "binary_ar1_data")
