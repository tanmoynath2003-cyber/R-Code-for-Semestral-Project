library(MASS)              
library(clusterGeneration) 
library(LinDA)
library(corrplot)
library(MicroBVS)
library(ggplot2)
library(tidyr)
library(foreach)
library(doSNOW)

#Correlation Plots

qzip_matrix <- function(U, p, Lambda_mat) {
  X <- matrix(0, nrow = nrow(U), ncol = ncol(U))
  mask <- U > p
  U_rescaled <- (U - p) / (1 - p)
  if(any(mask)) {
    X[mask] <- qpois(U_rescaled[mask], Lambda_mat[mask])
  }
  return(X)
}

simulate_with_covariates <- function(n_samples, n_taxa, rho, p_zero = 0.2, n_signif_taxa = 2, eff_list) {
  treatment <- rbinom(n_samples, 1, 0.5)
  age       <- rnorm(n_samples, 35, 10)
  bmi       <- rnorm(n_samples, 24, 4)
  plaque    <- runif(n_samples, 0, 3) 
  
  meta_df <- data.frame(Treatment = factor(treatment), Age = age, BMI = bmi, Plaque = plaque)
  
  beta_0 <- runif(n_taxa, 1, 3) 
  B <- matrix(0, nrow = 4, ncol = n_taxa) 
  
  actual_signif <- min(n_signif_taxa, n_taxa)
  true_idx <- seq_len(actual_signif)
  
  if(actual_signif > 0){
    B[1, true_idx] <- eff_list$cat
    B[2, true_idx] <- eff_list$age
    B[3, true_idx] <- eff_list$bmi
    B[4, true_idx] <- eff_list$plq
  }
  
  Lambda_mat <- matrix(NA, nrow = n_samples, ncol = n_taxa)
  X_design   <- cbind(treatment, scale(age), scale(bmi), scale(plaque))
  
  for(j in 1:n_taxa) {
    Lambda_mat[, j] <- exp(beta_0[j] + X_design %*% B[, j])
  }
  
  if (!isSymmetric(rho)) {
    rho[lower.tri(rho)] <- t(rho)[lower.tri(rho)]
  }
  
  Z <- mvrnorm(n_samples, mu = rep(0, n_taxa), Sigma = rho)
  U <- pnorm(Z)
  
  OTU_table <- qzip_matrix(U, p_zero, Lambda_mat) 
  
  colnames(OTU_table) <- paste0("Taxon_", 1:n_taxa)
  rownames(OTU_table) <- paste0("Sample_", 1:n_samples)
  
  return(list(otu = OTU_table, meta = meta_df, truth = true_idx))
}

set.seed(123)

N <- 100    
D <- 5      
S <- 2      

eff_list <- list(cat = 1.5, age = 0, bmi = 0.7, plq = 0.9)

rho_matrix <- matrix(c(
  1.0,  0.8,  0.0,  0.0,  0.0,  
  0.8,  1.0,  0.0,  0.0,  0.0,  
  0.0,  0.0,  1.0, -0.6,  0.2,  
  0.0,  0.0, -0.6,  1.0,  0.0,  
  0.0,  0.0,  0.2,  0.0,  1.0   
), nrow = 5, ncol = 5, byrow = TRUE)

sim_data <- simulate_with_covariates(
  n_samples = N, 
  n_taxa = D, 
  rho = rho_matrix, 
  n_signif_taxa = S, 
  eff_list = eff_list
)

colnames(rho_matrix) <- rownames(rho_matrix) <- paste0("Taxon_", 1:D)
otu_tab <- as.data.frame(sim_data$otu)
meta_dat <- sim_data$meta

res_linda <- linda(
  otu.tab = t(otu_tab),      
  meta = meta_dat, 
  formula = '~Treatment + Age + BMI + Plaque', 
  alpha = 0.05
)

clr <- function(x) {
  log_x <- log(x + 0.5) 
  return(log_x - mean(log_x))
}

otu_clr <- t(apply(otu_tab, 1, clr))
fit_multi <- lm(otu_clr ~ Treatment + Age + BMI + Plaque, data = meta_dat)
residuals_multi <- residuals(fit_multi)
residual_correlation_multi <- cor(residuals_multi)

par(mfrow = c(1, 2)) 

corrplot(rho_matrix, 
         method = "color", 
         title = "True Sigma Matrix", 
         mar = c(0, 0, 2, 0),
         addCoef.col = "black", 
         tl.col = "black")

corrplot(residual_correlation_multi, 
         method = "color", 
         title = "Multivariate Residual Correlation", 
         mar = c(0, 0, 2, 0),
         addCoef.col = "black")





#LINDA VS MicroBVS

N_CORR_MATRICES <- 10
M_DATASETS      <- 10
TOTAL_RUNS      <- N_CORR_MATRICES * M_DATASETS

N_SAMPLES   <- 50
N_TAXA      <- 50
N_SIGNIF    <- 5
EFFECT_SIZE <- 0.5

MCMC_ITERATIONS <- 10000
MCMC_THIN       <- 10
Bn              <- 0.5

EXPECTED_SAVED  <- MCMC_ITERATIONS / MCMC_THIN
BURNIN_COUNT    <- floor(EXPECTED_SAVED * Bn)

print(paste("Settings: Total Iterations =", MCMC_ITERATIONS, "| Saved =", EXPECTED_SAVED, "| Burn-in =", BURNIN_COUNT))

calc_perf <- function(detected_names, true_names, total_possible) {
  TP <- length(intersect(detected_names, true_names))
  FP <- length(setdiff(detected_names, true_names))
  FN <- length(setdiff(true_names, detected_names))
  
  sens <- ifelse(length(true_names) > 0, TP / length(true_names), 0)
  prec <- ifelse(length(detected_names) > 0, TP / length(detected_names), 0)
  fdr  <- 1 - prec
  f1   <- ifelse((prec + sens) > 0, 2 * (prec * sens) / (prec + sens), 0)
  
  return(c(Sensitivity = sens, FDR = fdr, F1_Score = f1))
}

simulate_data_sequenced <- function(n_samples, n_taxa, rho, n_signif=2, effect_size=1.5, depth=5000) {
  treatment <- rbinom(n_samples, 1, 0.5)
  noise_var <- rnorm(n_samples)
  meta_df   <- data.frame(Treatment = factor(treatment), NoiseVar = noise_var)
  
  beta_0 <- runif(n_taxa, 3, 5) 
  beta_treatment <- rep(0, n_taxa)
  true_taxa_indices <- seq_len(min(n_signif, n_taxa))
  if(length(true_taxa_indices) > 0) beta_treatment[true_taxa_indices] <- effect_size 
  true_taxa_names <- paste0("Taxon_", true_taxa_indices)
  
  Lambda_mat <- matrix(NA, nrow=n_samples, ncol=n_taxa)
  for(j in 1:n_taxa) Lambda_mat[, j] <- exp(beta_0[j] + (beta_treatment[j] * treatment))
  
  Sigma <- rho
  if (!isSymmetric(Sigma)) Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
  Z <- mvrnorm(n_samples, mu=rep(0, n_taxa), Sigma=Sigma)
  U <- pnorm(Z)
  
  raw_counts <- matrix(0, nrow=n_samples, ncol=n_taxa)
  for(j in 1:n_taxa) raw_counts[,j] <- qpois(U[,j], Lambda_mat[,j])
  
  sequenced_counts <- matrix(0, nrow=n_samples, ncol=n_taxa)
  for(i in 1:n_samples) {
    total <- sum(raw_counts[i, ])
    if(total > 0) probs <- raw_counts[i, ] / total else probs <- rep(1/n_taxa, n_taxa)
    sequenced_counts[i, ] <- rmultinom(1, size=depth, prob=probs)
  }
  
  colnames(sequenced_counts) <- paste0("Taxon_", 1:n_taxa)
  rownames(sequenced_counts) <- paste0("Sample_", 1:n_samples)
  
  return(list(otu = sequenced_counts, meta = meta_df, truth = true_taxa_names))
}

num_cores <- parallel::detectCores() - 1 
cl <- makeCluster(num_cores)
registerDoSNOW(cl)

print(paste("Parallel Cluster initiated with", num_cores, "cores."))
pb <- txtProgressBar(max = TOTAL_RUNS, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

start_time <- Sys.time()

results_log <- foreach(mat_id = 1:N_CORR_MATRICES, .combine = rbind, 
                       .packages = c('MASS', 'LinDA', 'MicroBVS', 'clusterGeneration'),
                       .options.snow = opts) %:%
  foreach(data_id = 1:M_DATASETS, .combine = rbind) %dopar% {
    
    unique_seed <- 42 + (mat_id * 10000) + data_id
    set.seed(unique_seed)
    
    set.seed(42 + mat_id)
    rho <- rcorrmatrix(N_TAXA)
    set.seed(unique_seed)
    
    sim <- simulate_data_sequenced(N_SAMPLES, N_TAXA, rho, n_signif = N_SIGNIF, 
                                   effect_size = EFFECT_SIZE) 
    
    otu_tab    <- sim$otu
    meta_dat   <- sim$meta
    true_names <- sim$truth
    
    iter_res <- data.frame()
    
    try({
      t0 <- Sys.time()
      res_linda  <- linda(t(otu_tab), meta_dat, formula = '~Treatment+NoiseVar', alpha = 0.05)
      time_linda <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      target     <- grep("Treatment", names(res_linda$output), value=TRUE)[1]
      linda_detected <- rownames(res_linda$output[[target]])[res_linda$output[[target]]$reject]
      perf       <- calc_perf(linda_detected, true_names, N_TAXA)
      iter_res   <- rbind(iter_res, data.frame(MatrixID=mat_id, DatasetID=data_id, Method="LinDA", Sensitivity=perf["Sensitivity"], FDR=perf["FDR"], F1=perf["F1_Score"], Time=time_linda))
    }, silent=TRUE)
    
    try({
      t0 <- Sys.time()
      z_mat <- as.matrix(otu_tab)
      keep_cols <- colMeans(z_mat > 0) > 0.05
      
      if(sum(keep_cols) >= 2) {
        z_filt <- z_mat[, keep_cols, drop=FALSE]
        x_mat  <- model.matrix(~ Treatment + NoiseVar, data = meta_dat)[, -1, drop = FALSE]
        
        fit <- DMbvs_R(
          z = z_filt, 
          x = x_mat, 
          prior = "BB", 
          iterations = MCMC_ITERATIONS,
          thin = MCMC_THIN,
          seed = unique_seed
        )
        
        n_saved <- dim(fit$zeta)[3]
        manual_burnin <- floor(n_saved * Bn)
        if(manual_burnin >= n_saved) manual_burnin <- n_saved - 10
        
        res_mbvs <- selected_DM(dm_obj = fit, burnin = manual_burnin, plotting = FALSE)
        
        pips <- res_mbvs$mppi
        if(ncol(pips) != ncol(z_filt)) pips <- t(pips)
        rn <- if(nrow(pips) == ncol(x_mat)+1) c("Intercept", colnames(x_mat)) else colnames(x_mat)
        rownames(pips) <- rn
        colnames(pips) <- colnames(z_filt)
        
        target_row    <- grep("Treatment", rownames(pips), value=TRUE)[1]
        mbvs_detected <- names(which(pips[target_row, ] > 0.5))
        
        perf <- calc_perf(mbvs_detected, true_names, N_TAXA)
        
        iter_res <- rbind(iter_res, data.frame(
          MatrixID=mat_id, DatasetID=data_id, Method="MicroBVS", 
          Sensitivity=perf["Sensitivity"], FDR=perf["FDR"], F1=perf["F1_Score"], Time=as.numeric(difftime(Sys.time(), t0, units="secs"))
        ))
      }
    }, silent=TRUE)
    
    return(iter_res)
  }

close(pb)
stopCluster(cl)
print(paste("Total Simulation Time:", round(difftime(Sys.time(), start_time, units="mins"), 2), "minutes"))

if(nrow(results_log) > 0) {
  print("--- Final Performance Summary ---")
  print(aggregate(cbind(Sensitivity, FDR, F1, Time) ~ Method, data = results_log, mean))
  
  results_long <- pivot_longer(results_log, cols = c("Sensitivity", "FDR", "F1"), 
                               names_to = "Metric", values_to = "Score")
  
  p1 <- ggplot(results_long, aes(x = Method, y = Score, fill = Method)) +
    geom_boxplot(alpha=0.7, outlier.shape = 1) +
    facet_wrap(~Metric, scales = "free_y") +
    labs(title="LinDA vs MicroBVS", 
         subtitle = paste0("N=", N_SAMPLES, ", Taxa=", N_TAXA, ", Effect=", EFFECT_SIZE),
         y="Score (0 to 1)") +
    theme_minimal() +
    scale_fill_manual(values = c("LinDA" = "#E69F00", "MicroBVS" = "#56B4E9"))
  
  print(p1)
} else {
  print("Simulation produced no results.")
}
