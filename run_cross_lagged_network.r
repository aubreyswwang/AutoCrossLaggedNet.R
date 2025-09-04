# run_cross-lagged_network.R
# Automated Cross-lagged Network Analysis Project
# Author: Aubrey WANG
# Purpose: Generate cross-lagged network analysis reports by modifying only the config_cross_lagged_network.csv file.
# =============================================================================
# Clear environment
rm(list = ls())

# Load required packages
cat("Loading R packages...\n")
pkgs <- c("dplyr", "readr", "readxl", "writexl", "ggplot2",
          "patchwork", "bootnet", "qgraph", "glmnet", "huge",
          "tibble")

# Install and load all packages
invisible(lapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}))

# Manually specify working directory
# setwd("/your/project/directory")  # Please replace with your project path
cat("Current working directory:\n", getwd(), "\n") # Please confirm your working directory is correct

# Read configuration file
cat("Reading configuration file...\n")
if (!file.exists("config_cross_lagged.csv")) {
  stop("Error: Configuration file not found. Please create configuration file first!")
}

# Load config
config_df <- read_csv("config_cross_lagged.csv") %>%
  filter(!is.na(key))

config <- config_df %>%
  select(key, value) %>%
  deframe() %>%
  as.list()

# Parse config
data_path       <- config$data_path
sheet_name      <- config$sheet_name
filter_enabled  <- as.logical(config$filter_enabled)
filters         <- config$filters
use_weight      <- as.logical(config$use_weight)
weight_col      <- config$weight_col
save_path       <- config$save_path

t0_start        <- as.integer(config$t0_start_col)
t0_end          <- as.integer(config$t0_end_col)
t1_start        <- as.integer(config$t1_start_col)
t1_end          <- as.integer(config$t1_end_col)

cov_cols        <- config$covariate_cols
if (is.na(cov_cols) || cov_cols == "" || cov_cols == "NA") {
  cov_cols <- character(0)
} else {
  cov_cols <- trimws(unlist(strsplit(cov_cols, ",")))
}

node_labels     <- trimws(unlist(strsplit(config$node_labels, ",")))
node_fullnames  <- setNames(
  trimws(unlist(strsplit(config$node_fullnames, ","))),
  node_labels
)
node_groups     <- trimws(unlist(strsplit(config$node_groups, ",")))

# Parse group_colors
color_list <- strsplit(config$group_colors, ";")[[1]]
group_colors <- setNames(
  sapply(strsplit(color_list, "="), `[`, 2),
  sapply(strsplit(color_list, "="), `[`, 1)
)

n_boots         <- as.integer(config$n_boots)
n_cores         <- as.integer(config$n_cores)
seed            <- as.integer(config$seed)

# Create output directory
dir.create(save_path, recursive = TRUE, showWarnings = FALSE)

# Load data
message("Loading data from: ", data_path)
data_raw <- read_excel(data_path, sheet = sheet_name)

# Function to apply filters
apply_filters <- function(df, filter_string) {
  if (is.na(filter_string) || is.null(filter_string) || trimws(filter_string) %in% c("", "NA")) {
    return(df)
  }

  conditions <- strsplit(trimws(filter_string), ";")[[1]]
  conditions <- conditions[conditions != ""]

  for (cond in conditions) {
    cond <- trimws(cond)
    op_match <- regexpr("==|!=|>=|<=|>|<|%in%", cond)
    if (op_match == -1) stop("Invalid operator in filter: ", cond)

    op <- regmatches(cond, op_match)
    parts <- regmatches(cond, regexec("(.+?)(==|!=|>=|<=|>|<|%in%)(.+)", cond))[[1]]

    col_name <- trimws(parts[2])
    operator <- parts[3]
    value_str <- trimws(parts[4])

    if (!col_name %in% colnames(df)) {
      stop("Filter failed: column not found - ", col_name)
    }

    lhs <- df[[col_name]]

    if (operator == "%in%") {
      rhs_values <- strsplit(value_str, ",")[[1]] %>% trimws()
      if (all(!is.na(suppressWarnings(as.numeric(rhs_values))))) {
        rhs_values <- as.numeric(rhs_values)
        lhs <- as.numeric(lhs)
      }
      mask <- lhs %in% rhs_values
    } else {
      rhs_val <- value_str
      if (!is.na(suppressWarnings(as.numeric(rhs_val)))) {
        rhs_val <- as.numeric(rhs_val)
        lhs <- as.numeric(lhs)
      } else {
        rhs_val <- as.character(rhs_val)
        lhs <- as.character(lhs)
      }

      switch(operator,
             "=="  = { mask <- lhs == rhs_val },
             "!="  = { mask <- lhs != rhs_val },
             ">"   = { mask <- lhs >  rhs_val },
             "<"   = { mask <- lhs <  rhs_val },
             ">="  = { mask <- lhs >= rhs_val },
             "<="  = { mask <- lhs <= rhs_val },
             stop("Unsupported operator: ", operator)
      )
    }
    df <- df[mask, , drop = FALSE]
  }
  return(df)
}

# Apply filtering
if (filter_enabled & !is.na(filters) && trimws(filters) != "" && filters != "NA") {
  message("Applying filters: ", filters)
  data_filtered <- apply_filters(data_raw, filters)
} else {
  message("Skipping data filtering.")
  data_filtered <- data_raw
}

if (nrow(data_filtered) == 0) stop("No data after filtering.")
message("Data loaded. N = ", nrow(data_filtered))

# Extract T0 and T1
if (max(t0_end, t1_end) > ncol(data_filtered)) {
  stop("Column index out of range.")
}

w0b <- data_filtered %>% select(t0_start:t0_end) %>% as.data.frame()
w1b <- data_filtered %>% select(t1_start:t1_end) %>% as.data.frame()

colnames(w0b) <- paste0(node_labels, "_t0")
colnames(w1b) <- paste0(node_labels, "_t1") 

w0b[] <- lapply(w0b, as.numeric)
w1b[] <- lapply(w1b, as.numeric)

# Prepare weight vector (optional)
weights_vec <- NULL
if (use_weight) {
  if (is.na(weight_col) || weight_col == "" || !(weight_col %in% colnames(data_filtered))) {
    stop("Weight column not found: ", weight_col)
  }
  weights_vec <- data_filtered[[weight_col]]
  if (any(is.na(weights_vec))) {
    warning("Missing values in weight column. Using 1 for NA.")
    weights_vec[is.na(weights_vec)] <- 1
  }
}

# Ensure covariates are numeric
data_filtered[cov_cols] <- lapply(data_filtered[cov_cols], function(x) {
  if (is.numeric(x)) return(x)
  if (is.factor(x)) x <- as.character(x)
  x_num <- suppressWarnings(as.numeric(x))
  if (all(is.na(x_num))) {
    warning("All values in column became NA after numeric conversion: ", deparse(substitute(x)))
  }
  return(x_num)
})

# Prepare data for cross-lagged model
w0_1 <- cbind(w0b, w1b)

if (length(cov_cols) > 0 && all(cov_cols %in% colnames(data_filtered))) {
  w0_1 <- cbind(w0_1, data_filtered[cov_cols])
} else if (length(cov_cols) > 0) {
  missing <- cov_cols[!cov_cols %in% colnames(data_filtered)]
  warning("Covariate columns not found, skipped: ", paste(missing, collapse = ", "))
}

if (use_weight) w0_1 <- cbind(w0_1, weights_vec)
W0_1 <- makeX(w0_1)

# Design matrix X: T0 + covariates
k <- ncol(w0b)
numCovar <- length(cov_cols)

X <- w0_1[, c(1:k, (k*2+1):(k*2+numCovar))]
Y <- w0_1[, (k+1):(k*2+numCovar)] # T1 variables

adjMat3 <- matrix(0, nrow = k + numCovar, ncol = k + numCovar)

set.seed(seed)
for (i in 1:k) {
  fit <- cv.glmnet(
    x = X,
    y = Y[, i],
    family = "gaussian",
    alpha = 1,
    standardize = TRUE,
    weights = weights_vec
  )
  lam <- fit$lambda.min
  coef_vec <- coef(fit, s = lam)
  adjMat3[, i] <- coef_vec[2:(k + numCovar + 1)]
}

adjMat3 <- getWmat(adjMat3, nNodes = k + numCovar, labels = c(node_labels, cov_cols), directed = TRUE)
adjMat3_final <- adjMat3[1:k, 1:k]

# Save Network Plots (No fixed size — user can modify later)
png(file.path(save_path, "network_with_self_loops.png"), width = 10, height = 8, units = "in", res = 300)
nwc1 <- qgraph(adjMat3_final, groups = node_groups, legend = TRUE, threshold = 0.05,
               labels = node_labels, colors = node_colors)
dev.off()
message("Saved: network_with_self_loops.png")

# No self-loops
diag(adjMat3_final) <- 0
png(file.path(save_path, "network_no_self_loops.png"), width = 10, height = 8, units = "in", res = 300)
qgraph(adjMat3_final, groups = node_groups, legend = TRUE, threshold = 0.05,
       labels = node_labels, colors = node_colors)
dev.off()
message("Saved: network_no_self_loops.png")

# Centrality plot
p_in <- centralityPlot(list("T0" = nwc1), weighted = TRUE, signed = TRUE,
                       scale = "z-scores", labels = node_labels,
                       include = "InExpectedInfluence", orderBy = "InExpectedInfluence")
p_out <- centralityPlot(list("T0" = nwc1), weighted = TRUE, signed = TRUE,
                        scale = "z-scores", labels = node_labels,
                        include = "OutExpectedInfluence", orderBy = "OutExpectedInfluence")
combined <- p_in + p_out + plot_layout(ncol = 2)
ggsave(file.path(save_path, "centrality_plot.png"), combined, width = 10, height = 6, dpi = 300)
message("Saved: centrality_plot.png")

# Bootstrap
boot_data <- data.frame(w0b, w1b)  # Only T0 for bootstrap design
if (length(cov_cols) > 0) boot_data <- cbind(boot_data, data_filtered[cov_cols])
if (use_weight) boot_data <- cbind(boot_data, ipw = weights_vec)
boot_data <- makeX(boot_data)

CLPN.OR <- function(data) {
  k <- ncol(w0b)
  X <- as.matrix(data[, 1:k])
  Y <- data[, (k + 1):(2 * k)]
  ipw <- if(use_weight) data[, "ipw"] else NULL
  adj <- matrix(0, nrow = k, ncol = k)
  for (i in 1:k) {
    fit <- cv.glmnet(
      x = X,
      y = Y[, i],
      family = "gaussian",
      alpha = 1,
      standardize = TRUE,
      weights = ipw
    )
    lam <- fit$lambda.min
    coef_vec <- coef(fit, s = lam)
    adj[, i] <- coef_vec[2:(k + 1)]
  }
  return(adj)
}

net.2 <- estimateNetwork(boot_data, fun = CLPN.OR, labels = node_labels, directed = TRUE)

# Save bootstrap network
png(file.path(save_path, "bootstrap_network.png"), width = 10, height = 8, units = "in", res = 300)
plot(net.2, groups = node_groups, legend = FALSE, colors = node_colors)
dev.off()
message("Saved: bootstrap_network.png")

# Run bootstrapping
set.seed(seed)
nonParBoot.w2 <- bootnet(net.2, type = "nonparametric", nBoots = n_boots, directed = TRUE,
                         statistics = c("edge", "outExpectedInfluence", "inExpectedInfluence"),
                         ncores = n_cores)

caseBoot.w2 <- bootnet(net.2, type = "case", nBoots = n_boots, directed = TRUE,
                       statistics = c("outExpectedInfluence", "inExpectedInfluence"))

# Extract results
cs_results <- corStability(caseBoot.w2)
boot_results <- nonParBoot.w2$sampleTable
boot_table <- nonParBoot.w2$bootTable

boot_stats <- boot_table %>%
  mutate(edge_name = paste(node1, "->", node2)) %>%
  group_by(edge_name) %>%
  summarise(
    se = sd(value, na.rm = TRUE),
    ci_lower = quantile(value, 0.025, na.rm = TRUE),
    ci_upper = quantile(value, 0.975, na.rm = TRUE),
    p_value = 2 * pmin(mean(value > 0, na.rm = TRUE), 1 - mean(value > 0, na.rm = TRUE)),
    .groups = 'drop'
  )

results_table <- boot_results %>%
  mutate(edge_name = paste(node1, "->", node2)) %>%
  left_join(boot_stats, by = "edge_name") %>%
  select(from = node1, to = node2, beta = value, type, name, se, ci_lower, ci_upper, p_value) %>%
  mutate(
    cs_inEI = cs_results["inExpectedInfluence"],
    cs_outEI = cs_results["outExpectedInfluence"]
  )

write_xlsx(results_table, path = file.path(save_path, "coefficients_with_stability.xlsx"))
message("Saved: coefficients_with_stability.xlsx")

# Save diagnostic plots (you can customize size later)
ggsave(file.path(save_path, "bootstrap_sample_network.png"), plot(nonParBoot.w2, order = 'sample'))
ggsave(file.path(save_path, "bootstrap_centrality.png"), plot(caseBoot.w2, statistics = c("outExpectedInfluence", "inExpectedInfluence")))
ggsave(file.path(save_path, "bootstrap_outEI_diff.png"), plot(nonParBoot.w2, 'outExpectedInfluence', plot = 'difference', order = 'sample'))
ggsave(file.path(save_path, "bootstrap_inEI_diff.png"), plot(nonParBoot.w2, 'inExpectedInfluence', plot = 'difference', order = 'sample'))
ggsave(file.path(save_path, "bootstrap_edge_diff.png"), plot(nonParBoot.w2, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample"))

# Top edges
top_idx <- order(adjMat3_final, decreasing = TRUE)[1:400]
pos <- arrayInd(top_idx, dim(adjMat3_final))
edges_df <- data.frame(
  from = node_labels[pos[,1]],
  to = node_labels[pos[,2]],
  value = adjMat3_final[top_idx],
  from_domain = node_groups[pos[,1]],
  to_domain = node_groups[pos[,2]]
)
write_xlsx(edges_df, path = file.path(save_path, "top_400_edges.xlsx"))
message("Saved: top_400_edges.xlsx")

message("All analysis completed successfully!")