# BAYESIAN CAUSAL FOREST (BCF)

# This script has two parts, which are run separately:
# Part A: treatment effect estimation over the 3 x 3 main grid (20 seeds).
# Part B: permutation-based explanation fidelity over the same grid (5 seeds).
#
# covariate matrix cleaned and subsampled by the Python
# notebook: lbidd_covariates_clean.csv

#  1: Packages

required_packages <- c("bcf", "ranger", "ggplot2")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(bcf)
  library(ranger)
  library(ggplot2)
})


#  2: Configuration

COVARIATE_FILE <- "lbidd_covariates_clean.csv"

RESULTS_FILE          <- "bcf_estimation_results.csv"
FAILURES_FILE         <- "bcf_estimation_failures.csv"
SUMMARY_FILE          <- "bcf_estimation_summary.csv"
SUMMARY_FILE_ROBUST   <- "bcf_estimation_summary_robust.csv"
OVERLAP_FILE          <- "bcf_overlap_diagnostics.csv"
REPORTING_TABLE_FILE  <- "bcf_reporting_table.csv"
FIGURE_FILE           <- "bcf_npehe_vs_signal.png"

SEEDS       <- 1:20
PRED_SCALES <- c(0.05, 0.30, 1.00)

CONFOUND_CONFIGS <- list(
  list(type = "none",       strength = 0.0),
  list(type = "prognostic", strength = 5.0),
  list(type = "predictive", strength = 5.0)
)

TEST_SIZE        <- 0.20
PROPENSITY_TREES <- 500

BCF_NBURN   <- 1000
BCF_NSIM    <- 1000
BCF_NCHAINS <- 2

SEEDS_EXPLAIN         <- 1:3
N_EXPLAIN             <- 20
# NSIM_SHAP             <- 10
EXPLAIN_RESULTS_FILE  <- "bcf_explanation_results.csv"
EXPLAIN_FAILURES_FILE <- "bcf_explanation_failures.csv"


#  3: Utility functions

sigmoid <- function(x) {
  plogis(x)
}

zscore <- function(x) {
  x_sd <- stats::sd(x)
  if (!is.finite(x_sd) || x_sd < 1e-8) {
    return(rep(0, length(x)))
  }
  (x - mean(x)) / x_sd
}

append_csv <- function(row, file_path) {
  exists_already <- file.exists(file_path)
  utils::write.table(
    row, file = file_path, sep = ",", row.names = FALSE,
    col.names = !exists_already, append = exists_already, quote = TRUE
  )
}

make_run_key <- function(pred_scale, confound_type, confound_strength, seed) {
  paste(pred_scale, confound_type, confound_strength, seed, sep = "_")
}

ensure_mixed_z <- function(z) {
  # predict.bcf() expects a treatment vector a tiny explanation sample may
  # contain one treatment group only, so change one placeholder value when
  # needed. This argument is required by the prediction interface, tau itself
  # is extracted from the returned treatment effect draws
  z <- as.numeric(z)
  if (length(z) >= 2 && length(unique(z)) < 2) {
    z[1] <- 1 - z[1]
  }
  z
}


#  4: Loading Python-cleaned covariates

load_covariates_option_A <- function(path = COVARIATE_FILE) {
  if (!file.exists(path)) {
    stop("Covariate file not found: ", path)
  }
  covariates <- utils::read.csv(path, check.names = FALSE)
  if (nrow(covariates) == 0 || ncol(covariates) == 0) {
    stop("The covariate file is empty.")
  }
  numeric_columns <- vapply(covariates, is.numeric, logical(1))
  if (!all(numeric_columns)) {
    stop(
      "All covariates must be numeric. Non-numeric columns: ",
      paste(names(covariates)[!numeric_columns], collapse = ", ")
    )
  }
  if (anyNA(covariates)) {
    stop("The covariate matrix contains missing values.")
  }
  message(
    "Loaded covariates: ", nrow(covariates), " rows x ",
    ncol(covariates), " columns."
  )
  covariates
}


#  5: Semi-synthetic data-generating process (DGP)

generate_semi_synthetic_dataset <- function(
    covariates,
    pred_scale = 1.0,
    confound_type = "none",
    confound_strength = 0.0,
    noise_sd = 1.0,
    n_relevant_per_group = 5,
    random_seed = NULL
) {
  valid_types <- c("none", "prognostic", "predictive")
  if (!confound_type %in% valid_types) {
    stop("Unknown confound_type: ", confound_type)
  }
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }
  
  X <- as.matrix(covariates)
  storage.mode(X) <- "double"
  n <- nrow(X)
  d <- ncol(X)
  n_required <- 3 * n_relevant_per_group
  if (d < n_required) {
    stop("At least ", n_required, " covariates are required.")
  }
  
  selected_indices <- sample(seq_len(d), size = n_required, replace = FALSE)
  Iprog <- selected_indices[1:n_relevant_per_group]
  I0 <- selected_indices[(n_relevant_per_group + 1):(2 * n_relevant_per_group)]
  I1 <- selected_indices[(2 * n_relevant_per_group + 1):(3 * n_relevant_per_group)]
  
  prog_weights  <- runif(n_relevant_per_group, -1, 1)
  pred0_weights <- runif(n_relevant_per_group, -1, 1)
  pred1_weights <- runif(n_relevant_per_group, -1, 1)
  
  prog_term  <- as.vector(X[, Iprog, drop = FALSE] %*% prog_weights)
  pred0_term <- as.vector(X[, I0,    drop = FALSE] %*% pred0_weights)
  pred1_term <- as.vector(X[, I1,    drop = FALSE] %*% pred1_weights)
  
  Y0 <- prog_term + pred_scale * pred0_term
  Y1 <- prog_term + pred_scale * pred1_term
  true_ITE <- Y1 - Y0
  
  if (confound_type == "none" || confound_strength == 0) {
    remaining_indices <- setdiff(seq_len(d), selected_indices)
    if (length(remaining_indices) > 0) {
      irrelevant_index <- sample(remaining_indices, size = 1)
      gamma <- zscore(X[, irrelevant_index])
    } else {
      gamma <- rep(0, n)
    }
  } else if (confound_type == "prognostic") {
    gamma <- zscore(prog_term)
  } else {
    gamma <- zscore(pred1_term - pred0_term)
  }
  
  propensity <- sigmoid(confound_strength * gamma)
  W <- rbinom(n = n, size = 1, prob = propensity)
  noise <- rnorm(n = n, mean = 0, sd = noise_sd)
  Y <- W * Y1 + (1 - W) * Y0 + noise
  
  df <- as.data.frame(X, check.names = FALSE)
  df$W <- W
  df$Y <- Y
  df$Y0 <- Y0
  df$Y1 <- Y1
  df$true_ITE <- true_ITE
  
  list(
    df = df,
    feature_cols = colnames(X),
    ground_truth = list(
      Iprog = colnames(X)[Iprog],
      I0 = colnames(X)[I0],
      I1 = colnames(X)[I1],
      Ipred = c(colnames(X)[I0], colnames(X)[I1])
    )
  )
}

#  Propensity model 
# fit_propensity_model <- function(X_train, treatment, X_test, seed) {
# if (length(unique(treatment)) < 2L) {
#    stop("Training data contain only one treatment group.")
#  }
  
#  model <- ranger::ranger(
##    x = as.data.frame(X_train),
#    y = factor(treatment, levels = c(0, 1)),
#    probability = TRUE,
#    num.trees = PROPENSITY_TREES,
#    seed = seed
#  )
  
  train_prob <- predict(model, data = as.data.frame(X_train))$predictions[, "1"]
  test_prob <- predict(model, data = as.data.frame(X_test))$predictions[, "1"]
  
  list(
    model = model,
    pihat_train = as.numeric(train_prob),
    pihat_test = as.numeric(test_prob)
  )
}

#  7: BCF posterior prediction helper

predict_bcf_tau <- function(fit, X_new_control, X_new_moderate, pi_new, z_new, tree_dir) {
  X_new_control  <- as.matrix(X_new_control)
  X_new_moderate <- as.matrix(X_new_moderate)
  storage.mode(X_new_control)  <- "double"
  storage.mode(X_new_moderate) <- "double"
  
  if (length(pi_new) != nrow(X_new_control)) stop("pi_new has the wrong length.")
  if (length(z_new)  != nrow(X_new_control)) stop("z_new has the wrong length.")
  
  z_new <- ensure_mixed_z(z_new)
  
  prediction <- predict(
    fit,
    x_predict_control = X_new_control,
    x_predict_moderate = X_new_moderate,
    pi_pred = pi_new,
    z_pred = z_new,
    save_tree_directory = tree_dir
  )
  
  tau_draws <- prediction$tau
  if (is.null(tau_draws)) stop("predict.bcf() did not return tau draws.")
  if (is.matrix(tau_draws)) return(as.numeric(colMeans(tau_draws)))
  as.numeric(tau_draws)
}


# 8.  fit one BCF model

fit_bcf_estimation_only <- function(
    dataset_list,
    test_size = TEST_SIZE,
    random_seed = 42,
    nburn = BCF_NBURN,
    nsim = BCF_NSIM,
    n_chains = BCF_NCHAINS,
    propensity_trees = PROPENSITY_TREES,
    verbose_bcf = FALSE
) {
  df <- dataset_list$df
  feature_cols <- dataset_list$feature_cols
  
  tree_dir <- file.path(
    tempdir(),
    paste0("bcf_trees_seed_", random_seed, "_",
           format(Sys.time(), "%Y%m%d_%H%M%S"), "_", sample.int(1e8, 1))
  )
  dir.create(tree_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(tree_dir, recursive = TRUE, force = TRUE), add = TRUE)
  
  set.seed(random_seed)
  n <- nrow(df)
  test_indices <- sample(seq_len(n), size = floor(test_size * n), replace = FALSE)
  train_indices <- setdiff(seq_len(n), test_indices)
  
  train_df <- df[train_indices, , drop = FALSE]
  test_df  <- df[test_indices, , drop = FALSE]
  
  X_train <- as.matrix(train_df[, feature_cols, drop = FALSE]); storage.mode(X_train) <- "double"
  X_test  <- as.matrix(test_df[, feature_cols, drop = FALSE]);  storage.mode(X_test)  <- "double"
  W_train <- as.numeric(train_df$W)
  W_test  <- as.numeric(test_df$W)
  Y_train <- as.numeric(train_df$Y)
  true_ite_test <- as.numeric(test_df$true_ITE)
  
  if (length(unique(W_train)) < 2) {
    stop("Training data contain only one treatment group.")
  }
  
  cat("Stage 1/4: fitting propensity model...\n"); flush.console()
  propensity_model <- ranger::ranger(
    x = as.data.frame(X_train), y = factor(W_train, levels = c(0, 1)),
    probability = TRUE, num.trees = propensity_trees, seed = random_seed
  )
  pihat_train <- as.numeric(predict(propensity_model, data = as.data.frame(X_train))$predictions[, "1"])
  pihat_test  <- as.numeric(predict(propensity_model, data = as.data.frame(X_test))$predictions[, "1"])
  
  cat("Stage 2/4: fitting BCF MCMC...\n"); flush.console()
  fit <- bcf::bcf(
    y = Y_train, z = W_train, x_control = X_train, x_moderate = X_train,
    pihat = pihat_train, nburn = nburn, nsim = nsim, n_chains = n_chains,
    random_seed = random_seed, include_pi = "control", verbose = verbose_bcf,
    save_tree_directory = tree_dir
  )
  
  cat("Stage 3/4: predicting held-out treatment effects...\n"); flush.console()
  predicted_ite <- predict_bcf_tau(fit, X_test, X_test, pihat_test, W_test, tree_dir)
  
  if (length(predicted_ite) != length(true_ite_test)) {
    stop("Predicted ITE length does not match test-set length.")
  }
  if (any(!is.finite(predicted_ite))) {
    stop("Predicted ITE contains invalid values.")
  }
  
  cat("Stage 4/4: calculating PEHE and nPEHE...\n"); flush.console()
  pehe <- sqrt(mean((predicted_ite - true_ite_test)^2))
  true_ite_sd <- sd(true_ite_test)
  if (!is.finite(true_ite_sd) || true_ite_sd < 1e-8) {
    stop("True ITE standard deviation is invalid.")
  }
  npehe <- pehe / true_ite_sd
  
  cat("BCF estimation run completed.\n"); flush.console()
  
  list(
    pehe = pehe, npehe = npehe, true_ite_sd = true_ite_sd,
    n_train = length(train_indices), n_test = length(test_indices)
  )
}


#  Time one experiment

time_one_bcf_experiment <- function(
    covariates, pred_scale = 1.0, confound_type = "predictive",
    confound_strength = 5.0, random_seed = 1,
    nburn = BCF_NBURN, nsim = BCF_NSIM, n_chains = BCF_NCHAINS
) {
  dataset <- generate_semi_synthetic_dataset(
    covariates = covariates, pred_scale = pred_scale,
    confound_type = confound_type, confound_strength = confound_strength,
    random_seed = random_seed
  )
  start_time <- Sys.time()
  result <- fit_bcf_estimation_only(
    dataset_list = dataset, random_seed = random_seed,
    nburn = nburn, nsim = nsim, n_chains = n_chains
  )
  elapsed_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  estimated_hours <- (elapsed_minutes * length(SEEDS) * length(PRED_SCALES) *
                        length(CONFOUND_CONFIGS)) / 60
  
  cat("\nOne run took ", round(elapsed_minutes, 2), " minutes.\n", sep = "")
  cat("Estimated 180-run runtime: ", round(estimated_hours, 2), " hours.\n", sep = "")
  print(result)
  
  invisible(list(result = result, elapsed_minutes = elapsed_minutes,
                 estimated_hours = estimated_hours))
}


# checkpointed estimation grid

run_bcf_grid_checkpointed <- function(
    covariates, seeds = SEEDS, pred_scales = PRED_SCALES,
    confound_configs = CONFOUND_CONFIGS, results_file = RESULTS_FILE,
    failures_file = FAILURES_FILE, nburn = BCF_NBURN, nsim = BCF_NSIM,
    n_chains = BCF_NCHAINS, propensity_trees = PROPENSITY_TREES
) {
  completed_keys <- character(0)
  if (file.exists(results_file)) {
    existing_results <- utils::read.csv(results_file, stringsAsFactors = FALSE)
    required_columns <- c("pred_scale", "confound_type", "confound_strength", "seed")
    if (!all(required_columns %in% names(existing_results))) {
      stop("Existing results file has incompatible columns.")
    }
    completed_keys <- make_run_key(
      existing_results$pred_scale, existing_results$confound_type,
      existing_results$confound_strength, existing_results$seed
    )
    message("Found ", length(completed_keys), " previously completed runs.")
  }
  
  total_runs <- length(seeds) * length(pred_scales) * length(confound_configs)
  counter <- 0
  
  for (pred_scale in pred_scales) {
    for (config in confound_configs) {
      confound_type <- config$type
      confound_strength <- config$strength
      
      for (seed in seeds) {
        counter <- counter + 1
        current_key <- make_run_key(pred_scale, confound_type, confound_strength, seed)
        
        if (current_key %in% completed_keys) {
          cat(sprintf("[%d/%d] SKIP: %s\n", counter, total_runs, current_key))
          next
        }
        
        cat(sprintf(
          "\n[%d/%d] signal=%.2f, confounding=%s(%.1f), seed=%d\n",
          counter, total_runs, pred_scale, confound_type, confound_strength, seed
        ))
        
        run_start <- Sys.time()
        run_result <- tryCatch({
          dataset <- generate_semi_synthetic_dataset(
            covariates = covariates, pred_scale = pred_scale,
            confound_type = confound_type, confound_strength = confound_strength,
            random_seed = seed
          )
          fitted <- fit_bcf_estimation_only(
            dataset_list = dataset, random_seed = seed, nburn = nburn,
            nsim = nsim, n_chains = n_chains, propensity_trees = propensity_trees
          )
          runtime_minutes <- as.numeric(difftime(Sys.time(), run_start, units = "mins"))
          
          data.frame(
            pred_scale = round(pred_scale, 4), confound_type = confound_type,
            confound_strength = round(confound_strength, 4), seed = seed,
            pehe = round(fitted$pehe, 4), npehe = round(fitted$npehe, 4),
            true_ite_sd = round(fitted$true_ite_sd, 4),
            n_train = fitted$n_train, n_test = fitted$n_test,
            runtime_minutes = round(runtime_minutes, 4), stringsAsFactors = FALSE
          )
        }, error = function(e) {
          runtime_minutes <- as.numeric(difftime(Sys.time(), run_start, units = "mins"))
          failure <- data.frame(
            pred_scale = pred_scale, confound_type = confound_type,
            confound_strength = confound_strength, seed = seed,
            runtime_minutes = round(runtime_minutes, 4),
            error_message = conditionMessage(e), stringsAsFactors = FALSE
          )
          append_csv(failure, failures_file)
          message("Run failed: ", conditionMessage(e))
          NULL
        })
        
        if (!is.null(run_result)) {
          append_csv(run_result, results_file)
          completed_keys <- c(completed_keys, current_key)
          cat(sprintf("Saved: nPEHE=%.4f, runtime=%.2f minutes.\n",
                      run_result$npehe, run_result$runtime_minutes))
        }
        invisible(gc())
      }
    }
  }
  message("BCF grid finished.")
}


#  result

validate_bcf_results <- function(results_file = RESULTS_FILE) {
  if (!file.exists(results_file)) stop("Results file not found: ", results_file)
  results <- utils::read.csv(results_file, stringsAsFactors = FALSE)
  
  keys <- make_run_key(results$pred_scale, results$confound_type,
                       results$confound_strength, results$seed)
  duplicate_count <- sum(duplicated(keys))
  if (duplicate_count > 0) {
    warning(duplicate_count, " duplicate result rows were found.")
  }
  
  counts <- aggregate(
    results$seed,
    by = results[c("pred_scale", "confound_type", "confound_strength")],
    FUN = length
  )
  names(counts)[4] <- "n_completed"
  
  cat("Expected total runs: 180\n")
  cat("Saved rows: ", nrow(results), "\n", sep = "")
  print(counts)
  invisible(counts)
}

summarise_bcf_results <- function(results_file = RESULTS_FILE, summary_file = SUMMARY_FILE) {
  if (!file.exists(results_file)) stop("Results file not found: ", results_file)
  results <- utils::read.csv(results_file, stringsAsFactors = FALSE)
  grouping <- c("pred_scale", "confound_type", "confound_strength")
  metrics <- c("pehe", "npehe", "runtime_minutes")
  
  means <- aggregate(results[metrics], by = results[grouping], FUN = mean, na.rm = TRUE)
  standard_deviations <- aggregate(results[metrics], by = results[grouping], FUN = sd, na.rm = TRUE)
  minima <- aggregate(results[metrics], by = results[grouping], FUN = min, na.rm = TRUE)
  maxima <- aggregate(results[metrics], by = results[grouping], FUN = max, na.rm = TRUE)
  counts <- aggregate(results$seed, by = results[grouping], FUN = length)
  
  names(means)[4:6] <- c("pehe_mean", "npehe_mean", "runtime_minutes_mean")
  names(standard_deviations)[4:6] <- c("pehe_std", "npehe_std", "runtime_minutes_std")
  names(minima)[4:6] <- c("pehe_min", "npehe_min", "runtime_minutes_min")
  names(maxima)[4:6] <- c("pehe_max", "npehe_max", "runtime_minutes_max")
  names(counts)[4] <- "n_completed"
  
  summary_df <- Reduce(function(l, r) merge(l, r, by = grouping, all = TRUE),
                       list(means, standard_deviations, minima, maxima, counts))
  summary_df <- summary_df[order(summary_df$pred_scale, summary_df$confound_type), ]
  
  numeric_columns <- vapply(summary_df, is.numeric, logical(1))
  summary_df[numeric_columns] <- lapply(summary_df[numeric_columns], function(x) round(x, 4))
  
  utils::write.csv(summary_df, summary_file, row.names = FALSE)
  print(summary_df)
  invisible(summary_df)
}

summarise_bcf_results_robust <- function(
    results_file = RESULTS_FILE, summary_file = SUMMARY_FILE_ROBUST
) {
  results <- utils::read.csv(results_file, stringsAsFactors = FALSE)
  grouping <- c("pred_scale", "confound_type", "confound_strength")
  
  robust_summary <- do.call(rbind, lapply(
    split(results, results[grouping]),
    function(g) {
      if (nrow(g) == 0) return(NULL)
      npehe <- g$npehe
      q1 <- stats::quantile(npehe, 0.25, na.rm = TRUE)
      q3 <- stats::quantile(npehe, 0.75, na.rm = TRUE)
      iqr <- q3 - q1
      upper_fence <- q3 + 3 * iqr
      n_extreme <- sum(npehe > upper_fence, na.rm = TRUE)
      
      data.frame(
        pred_scale = g$pred_scale[1], confound_type = g$confound_type[1],
        confound_strength = g$confound_strength[1], n = nrow(g),
        npehe_mean = round(mean(npehe, na.rm = TRUE), 4),
        npehe_median = round(stats::median(npehe, na.rm = TRUE), 4),
        npehe_sd = round(stats::sd(npehe, na.rm = TRUE), 4),
        npehe_iqr = round(iqr, 4), n_extreme_outliers = n_extreme,
        pct_extreme_outliers = round(100 * n_extreme / nrow(g), 1)
      )
    }
  ))
  
  rownames(robust_summary) <- NULL
  robust_summary <- robust_summary[
    order(robust_summary$pred_scale, robust_summary$confound_type), ]
  
  utils::write.csv(robust_summary, summary_file, row.names = FALSE)
  print(robust_summary)
  invisible(robust_summary)
}

diagnose_bcf_overlap <- function(
    covariates, seeds = SEEDS, pred_scales = PRED_SCALES,
    confound_configs = CONFOUND_CONFIGS, propensity_trees = PROPENSITY_TREES,
    output_file = OVERLAP_FILE
) {
  rows <- list()
  counter <- 0
  
  for (pred_scale in pred_scales) {
    for (config in confound_configs) {
      for (seed in seeds) {
        counter <- counter + 1
        
        dataset <- generate_semi_synthetic_dataset(
          covariates = covariates, pred_scale = pred_scale,
          confound_type = config$type, confound_strength = config$strength,
          random_seed = seed
        )
        df <- dataset$df
        feature_cols <- dataset$feature_cols
        
        set.seed(seed)
        n <- nrow(df)
        test_indices <- sample(seq_len(n), size = floor(TEST_SIZE * n))
        train_indices <- setdiff(seq_len(n), test_indices)
        
        X_train <- as.matrix(df[train_indices, feature_cols, drop = FALSE])
        storage.mode(X_train) <- "double"
        W_train <- as.numeric(df$W[train_indices])
        
        propensity_model <- ranger::ranger(
          x = as.data.frame(X_train), y = factor(W_train, levels = c(0, 1)),
          probability = TRUE, num.trees = propensity_trees, seed = seed
        )
        pihat_train <- predict(propensity_model, data = as.data.frame(X_train))$predictions[, "1"]
        
        rows[[counter]] <- data.frame(
          pred_scale = pred_scale, confound_type = config$type,
          confound_strength = config$strength, seed = seed,
          pihat_min = round(min(pihat_train), 4), pihat_max = round(max(pihat_train), 4),
          pct_below_01 = round(100 * mean(pihat_train < 0.01), 2),
          pct_above_99 = round(100 * mean(pihat_train > 0.99), 2)
        )
        
        cat(sprintf(
          "[%d/%d] pred_scale=%.2f confound=%s: pihat range [%.4f, %.4f], extreme=%.1f%%\n",
          counter, length(pred_scales) * length(confound_configs) * length(seeds),
          pred_scale, config$type, min(pihat_train), max(pihat_train),
          100 * mean(pihat_train < 0.01 | pihat_train > 0.99)
        ))
      }
    }
  }
  
  out <- do.call(rbind, rows)
  utils::write.csv(out, output_file, row.names = FALSE)
  invisible(out)
}

create_bcf_reporting_table <- function(summary_file = SUMMARY_FILE, output_file = REPORTING_TABLE_FILE) {
  summary_df <- utils::read.csv(summary_file, stringsAsFactors = FALSE)
  reporting_table <- data.frame(
    predictive_signal = summary_df$pred_scale,
    confounding_mechanism = summary_df$confound_type,
    n_completed = summary_df$n_completed,
    pehe_mean = round(summary_df$pehe_mean, 4), pehe_std = round(summary_df$pehe_std, 4),
    npehe_mean = round(summary_df$npehe_mean, 4), npehe_std = round(summary_df$npehe_std, 4)
  )
  utils::write.csv(reporting_table, output_file, row.names = FALSE)
  print(reporting_table)
  invisible(reporting_table)
}

make_bcf_npehe_figure <- function(summary_file = SUMMARY_FILE, figure_file = FIGURE_FILE) {
  if (!file.exists(summary_file)) stop("Summary file not found: ", summary_file)
  summary_df <- utils::read.csv(summary_file, stringsAsFactors = FALSE)
  summary_df$confound_type <- factor(
    summary_df$confound_type, levels = c("none", "prognostic", "predictive"),
    labels = c("None", "Prognostic", "Predictive")
  )
  
  plot_object <- ggplot(summary_df, aes(x = pred_scale, y = npehe_mean,
                                        colour = confound_type, group = confound_type)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = pmax(0, npehe_mean - npehe_std), ymax = npehe_mean + npehe_std),
                  width = 0.025) +
    scale_x_continuous(breaks = c(0.05, 0.30, 1.00)) +
    labs(title = "BCF treatment-effect estimation accuracy",
         subtitle = "Mean normalised PEHE \u00b1 1 SD across 20 seeds",
         x = "Predictive signal strength", y = "Normalised PEHE",
         colour = "Confounding mechanism") +
    theme_minimal(base_size = 12)
  
  ggsave(filename = figure_file, plot = plot_object, width = 7.5, height = 5.2, dpi = 300)
  print(plot_object)
  invisible(plot_object)
}


# EXPLANATION STAGE

#   Attrpred / Attrprog  


compute_attr_metrics <- function(importance, feature_names, pred_covariates, prog_covariates) {
  importance <- abs(importance)
  names(importance) <- feature_names
  total <- sum(importance)
  if (!is.finite(total) || total < 1e-12) {
    return(list(attrpred = NA, attrprog = NA))
  }
  list(
    attrpred = sum(importance[pred_covariates]) / total,
    attrprog = sum(importance[prog_covariates]) / total
  )
}



# SECTION 12: Part B — permutation importance (replaces SHAP)     
# SHAP-based Shapley sampling was found computationally expensive at
# full-grid scale, because predict.bcf() carries a large fixed cost per call
# that barely depends on row count, and Shapley sampling requires many
# coalition evaluations per explicand. Permutation importance needs only ONE
# baseline batch + ONE permuted batch per feature, combined into a SINGLE
# predict.bcf() call, exploiting that same fixed-cost structure.
#
# pihat_baseline and z_baseline must be PER-ROW VECTORS (length = nrow(
# X_baseline)), one value per explicand not a single averaged scalar.
# Using each explicand's own propensity score and treatment status keeps
# every row internally consistent when only one covariate is permuted
# collapsing to an averaged scalar both breaks the row-count math

compute_bcf_permutation_importance <- function(fit, X_baseline, pihat_baseline, z_baseline, tree_dir) {
  n <- nrow(X_baseline)
  p <- ncol(X_baseline)
  feature_names <- colnames(X_baseline)
  
  if (length(pihat_baseline) != n) stop("pihat_baseline must have length nrow(X_baseline).")
  if (length(z_baseline)     != n) stop("z_baseline must have length nrow(X_baseline).")
  
  blocks <- list(X_baseline)
  for (j in seq_len(p)) {
    X_perm <- X_baseline
    X_perm[, j] <- sample(X_perm[, j])
    blocks[[j + 1]] <- X_perm
  }
  big_matrix <- do.call(rbind, blocks)
  storage.mode(big_matrix) <- "double"
  
  # Repeat the PER-ROW vectors (not a scalar) across all p+1 blocks, so
  # every row of big_matrix has the matching individual's own pihat/z.
  pihat_vec <- rep(pihat_baseline, times = p + 1)   # length (p+1)*n
  z_vec     <- rep(z_baseline,     times = p + 1)   # length (p+1)*n
  z_vec     <- ensure_mixed_z(z_vec)
  
  if (length(pihat_vec) != nrow(big_matrix)) {
    stop("Internal error: pihat_vec length does not match big_matrix rows.")
  }
  
  cat("Permutation importance: ONE predict.bcf() call on", nrow(big_matrix), "rows...\n")
  flush.console()
  t0 <- Sys.time()
  preds <- predict_bcf_tau(fit, big_matrix, big_matrix, pihat_vec, z_vec, tree_dir)
  cat("Batched predict took", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n")
  
  baseline_preds <- preds[1:n]
  importance <- numeric(p)
  for (j in seq_len(p)) {
    perm_preds <- preds[(j * n + 1):((j + 1) * n)]
    importance[j] <- mean(abs(perm_preds - baseline_preds))
  }
  names(importance) <- feature_names
  importance
}



#  fit one BCF model WITH permutation importance 

fit_bcf_and_explain <- function(
    dataset_list, test_size = TEST_SIZE, random_seed = 42,
    nburn = BCF_NBURN, nsim = BCF_NSIM, n_chains = BCF_NCHAINS,
    propensity_trees = PROPENSITY_TREES, n_explain = N_EXPLAIN
) {
  df <- dataset_list$df
  feature_cols <- dataset_list$feature_cols
  ground_truth <- dataset_list$ground_truth
  
  tree_dir <- file.path(tempdir(), paste0("bcf_perm_", random_seed, "_", sample.int(1e8, 1)))
  dir.create(tree_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(tree_dir, recursive = TRUE, force = TRUE), add = TRUE)
  
  set.seed(random_seed)
  n <- nrow(df)
  test_indices <- sample(seq_len(n), size = floor(test_size * n))
  train_indices <- setdiff(seq_len(n), test_indices)
  
  X_train <- as.matrix(df[train_indices, feature_cols, drop = FALSE]); storage.mode(X_train) <- "double"
  X_test  <- as.matrix(df[test_indices,  feature_cols, drop = FALSE]); storage.mode(X_test)  <- "double"
  W_train <- as.numeric(df$W[train_indices]); W_test <- as.numeric(df$W[test_indices])
  Y_train <- as.numeric(df$Y[train_indices])
  true_ite_test <- as.numeric(df$true_ITE[test_indices])
  
  propensity_model <- ranger::ranger(
    x = as.data.frame(X_train), y = factor(W_train, levels = c(0, 1)),
    probability = TRUE, num.trees = propensity_trees, seed = random_seed
  )
  pihat_train <- predict(propensity_model, data = as.data.frame(X_train))$predictions[, "1"]
  pihat_test  <- predict(propensity_model, data = as.data.frame(X_test))$predictions[, "1"]
  
  cat("Fitting BCF MCMC...\n"); flush.console()
  fit <- bcf::bcf(
    y = Y_train, z = W_train, x_control = X_train, x_moderate = X_train,
    pihat = pihat_train, nburn = nburn, nsim = nsim, n_chains = n_chains,
    random_seed = random_seed, include_pi = "control", verbose = FALSE,
    save_tree_directory = tree_dir
  )
  
  cat("Predicting held-out treatment effects...\n"); flush.console()
  predicted_ite <- predict_bcf_tau(fit, X_test, X_test, pihat_test, W_test, tree_dir)
  pehe <- sqrt(mean((predicted_ite - true_ite_test)^2))
  true_ite_sd <- sd(true_ite_test)
  npehe <- pehe / true_ite_sd
  .
  explicand_idx <- sample(seq_len(nrow(X_test)), size = min(n_explain, nrow(X_test)))
  X_explain <- X_test[explicand_idx, , drop = FALSE]
  pihat_explain <- pihat_test[explicand_idx]   # per-row vector, length = n_explain
  z_explain     <- W_test[explicand_idx]        # per-row vector, length = n_explain
  
  cat("Computing permutation importance on", nrow(X_explain), "explicands...\n"); flush.console()
  importance <- compute_bcf_permutation_importance(fit, X_explain, pihat_explain, z_explain, tree_dir)
  
  attr_metrics <- compute_attr_metrics(
    importance, feature_cols, ground_truth$Ipred, ground_truth$Iprog
  )
  
  list(
    pehe = pehe, npehe = npehe,
    attrpred = attr_metrics$attrpred, attrprog = attr_metrics$attrprog,
    n_explain = length(explicand_idx)
  )
}


#  checkpointed explanation grid  


run_bcf_explanation_grid_checkpointed <- function(
    covariates, seeds = SEEDS_EXPLAIN, pred_scales = PRED_SCALES,
    confound_configs = CONFOUND_CONFIGS, results_file = EXPLAIN_RESULTS_FILE,
    failures_file = EXPLAIN_FAILURES_FILE, nburn = BCF_NBURN, nsim = BCF_NSIM,
    n_chains = BCF_NCHAINS, propensity_trees = PROPENSITY_TREES,
    n_explain = N_EXPLAIN
) {
  completed_keys <- character(0)
  if (file.exists(results_file)) {
    existing_results <- utils::read.csv(results_file, stringsAsFactors = FALSE)
    completed_keys <- make_run_key(
      existing_results$pred_scale, existing_results$confound_type,
      existing_results$confound_strength, existing_results$seed
    )
    message("Found ", length(completed_keys), " previously completed explanation runs.")
  }
  
  total_runs <- length(seeds) * length(pred_scales) * length(confound_configs)
  counter <- 0
  
  for (pred_scale in pred_scales) {
    for (config in confound_configs) {
      confound_type <- config$type
      confound_strength <- config$strength
      
      for (seed in seeds) {
        counter <- counter + 1
        current_key <- make_run_key(pred_scale, confound_type, confound_strength, seed)
        
        if (current_key %in% completed_keys) {
          cat(sprintf("[%d/%d] SKIP: %s\n", counter, total_runs, current_key))
          next
        }
        
        cat(sprintf(
          "\n[%d/%d] (explanation) signal=%.2f, confounding=%s(%.1f), seed=%d\n",
          counter, total_runs, pred_scale, confound_type, confound_strength, seed
        ))
        
        run_start <- Sys.time()
        run_result <- tryCatch({
          dataset <- generate_semi_synthetic_dataset(
            covariates = covariates, pred_scale = pred_scale,
            confound_type = confound_type, confound_strength = confound_strength,
            random_seed = seed
          )
          fitted <- fit_bcf_and_explain(
            dataset_list = dataset, random_seed = seed, nburn = nburn,
            nsim = nsim, n_chains = n_chains, propensity_trees = propensity_trees,
            n_explain = n_explain
          )
          runtime_minutes <- as.numeric(difftime(Sys.time(), run_start, units = "mins"))
          
          data.frame(
            pred_scale = round(pred_scale, 4), confound_type = confound_type,
            confound_strength = round(confound_strength, 4), seed = seed,
            pehe = round(fitted$pehe, 4), npehe = round(fitted$npehe, 4),
            attrpred = round(fitted$attrpred, 4),
            attrprog = round(fitted$attrprog, 4),
            n_explain = fitted$n_explain,
            runtime_minutes = round(runtime_minutes, 4), stringsAsFactors = FALSE
          )
        }, error = function(e) {
          runtime_minutes <- as.numeric(difftime(Sys.time(), run_start, units = "mins"))
          failure <- data.frame(
            pred_scale = pred_scale, confound_type = confound_type,
            confound_strength = confound_strength, seed = seed,
            runtime_minutes = round(runtime_minutes, 4),
            error_message = conditionMessage(e), stringsAsFactors = FALSE
          )
          append_csv(failure, failures_file)
          message("Explanation run failed: ", conditionMessage(e))
          NULL
        })
        
        if (!is.null(run_result)) {
          append_csv(run_result, results_file)
          completed_keys <- c(completed_keys, current_key)
          cat(sprintf("Saved: nPEHE=%.4f, Attrpred=%.4f, Attrprog=%.4f, runtime=%.2f min.\n",
                      run_result$npehe, run_result$attrpred,
                      run_result$attrprog, run_result$runtime_minutes))
        }
        invisible(gc())
      }
    }
  }
  message("BCF explanation grid finished.")
}


# summary and figures


summarise_bcf_explanation_results <- function(
    results_file = EXPLAIN_RESULTS_FILE, summary_file = "bcf_explanation_summary.csv"
) {
  if (!file.exists(results_file)) stop("Explanation results file not found: ", results_file)
  results <- utils::read.csv(results_file, stringsAsFactors = FALSE)
  grouping <- c("pred_scale", "confound_type", "confound_strength")
  metrics <- c("npehe", "attrpred", "attrprog")
  
  means   <- aggregate(results[metrics], by = results[grouping], FUN = mean,   na.rm = TRUE)
  medians <- aggregate(results[metrics], by = results[grouping], FUN = median, na.rm = TRUE)
  sds     <- aggregate(results[metrics], by = results[grouping], FUN = sd,     na.rm = TRUE)
  counts  <- aggregate(results$seed,     by = results[grouping], FUN = length)
  
  names(means)[4:6]   <- paste0(metrics, "_mean")
  names(medians)[4:6] <- paste0(metrics, "_median")
  names(sds)[4:6]     <- paste0(metrics, "_std")
  names(counts)[4]    <- "n_completed"
  
  summary_df <- Reduce(function(l, r) merge(l, r, by = grouping, all = TRUE),
                       list(means, medians, sds, counts))
  summary_df <- summary_df[order(summary_df$pred_scale, summary_df$confound_type), ]
  
  numeric_columns <- vapply(summary_df, is.numeric, logical(1))
  summary_df[numeric_columns] <- lapply(summary_df[numeric_columns], function(x) round(x, 4))
  
  utils::write.csv(summary_df, summary_file, row.names = FALSE)
  print(summary_df)
  invisible(summary_df)
}


# =============================================================================
# SECTION 15: Part B — figures
# =============================================================================
# Figure 1: BCF Attrpred vs predictive signal strength, by confounding type
# (same style as your Causal Forest Attrpred figure, for a direct visual
# comparison across chapters/models).

make_bcf_attrpred_figure <- function(
    summary_file = "bcf_explanation_summary.csv",
    figure_file = "bcf_attrpred_vs_signal.png"
) {
  if (!file.exists(summary_file)) stop("Summary file not found: ", summary_file)
  summary_df <- utils::read.csv(summary_file, stringsAsFactors = FALSE)
  summary_df$confound_type <- factor(
    summary_df$confound_type, levels = c("none", "prognostic", "predictive"),
    labels = c("None", "Prognostic", "Predictive")
  )
  
  plot_object <- ggplot(summary_df, aes(x = pred_scale, y = attrpred_mean,
                                        colour = confound_type, group = confound_type)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = pmax(0, attrpred_mean - attrpred_std),
                      ymax = pmin(1, attrpred_mean + attrpred_std)),
                  width = 0.025) +
    scale_x_continuous(breaks = c(0.05, 0.30, 1.00)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "BCF explanation fidelity (permutation importance)",
         subtitle = paste0("Mean Attrpred \u00b1 1 SD across ", length(SEEDS_EXPLAIN), " seeds"),
         x = "Predictive signal strength", y = "Attrpred",
         colour = "Confounding mechanism") +
    theme_minimal(base_size = 12)
  
  ggsave(filename = figure_file, plot = plot_object, width = 7.5, height = 5.2, dpi = 300)
  print(plot_object)
  invisible(plot_object)
}

# Figure 2: BCF vs Causal Forest built-in importance, side by side, for the
# Q4 comparison. Requires the Causal Forest summary CSV (cf_sammary.csv in
# your earlier files) alongside the BCF explanation summary produced above.

make_bcf_vs_cf_attrpred_figure <- function(
    bcf_summary_file = "bcf_explanation_summary.csv",
    cf_summary_file = "cf_sammary.csv",
    figure_file = "bcf_vs_cf_attrpred.png"
) {
  if (!file.exists(bcf_summary_file)) stop("BCF summary file not found: ", bcf_summary_file)
  if (!file.exists(cf_summary_file))  stop("Causal Forest summary file not found: ", cf_summary_file)
  
  bcf_df <- utils::read.csv(bcf_summary_file, stringsAsFactors = FALSE)
  cf_df  <- utils::read.csv(cf_summary_file, stringsAsFactors = FALSE)
  
  bcf_plot_df <- data.frame(
    pred_scale = bcf_df$pred_scale, confound_type = bcf_df$confound_type,
    model = "BCF (permutation)", attrpred_mean = bcf_df$attrpred_mean
  )
  cf_plot_df <- data.frame(
    pred_scale = cf_df$pred_scale, confound_type = cf_df$confound_type,
    model = "Causal Forest (built-in)", attrpred_mean = cf_df$cf_attrpred_mean
  )
  combined_df <- rbind(bcf_plot_df, cf_plot_df)
  combined_df$condition <- paste0(combined_df$pred_scale, "\n", combined_df$confound_type)
  
  plot_object <- ggplot(combined_df, aes(x = condition, y = attrpred_mean, fill = model)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Explanation fidelity: BCF vs Causal Forest",
         subtitle = "Mean Attrpred by experimental condition",
         x = "Predictive scale / confounding type", y = "Attrpred", fill = "Model") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 0))
  
  ggsave(filename = figure_file, plot = plot_object, width = 8.5, height = 5.2, dpi = 300)
  print(plot_object)
  invisible(plot_object)
}



#  Pilot: time ONE run at FULL settings 
# pilot_dataset <- generate_semi_synthetic_dataset(
#   covariates, pred_scale = 1.0, confound_type = "predictive",
#   confound_strength = 5.0, random_seed = 1
# )
# t0 <- Sys.time()
# pilot_result <- fit_bcf_and_explain(
#   pilot_dataset, random_seed = 1,
#   nburn = 1000, nsim = 1000, n_chains = 2, n_explain = 15
# )
# cat("Took", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "minutes\n")
# print(pilot_result)
#
# 
# run_bcf_explanation_grid_checkpointed(covariates = covariates)
#
#Summarise and plot 
# bcf_explanation_summary <- summarise_bcf_explanation_results()
# make_bcf_attrpred_figure()


