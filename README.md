# Automated Cross-lagged Network Analysis
Just edit one CSV configuration file, run the provided R script, and the results will be automatically generated and saved into the result folder.

- Cross-lagged network estimation using LASSO regularization
- Bootstrap analysis for stability assessment
- Centrality analysis
- Network visualization
- Statistical testing with confidence intervals

**Note:** This script is designed specifically for continuous variables only (including covariates). Categorical variables are not supported.

## ✨ Features

- **Configuration-driven**: Modify only the CSV file to run different analyses
- **Automated workflow**: From data loading to result generation
- **Flexible filtering**: Advanced data filtering capabilities
- **Weight support**: Inverse probability weighting (IPW) integration

## 📦 Requirements

### R Version
- R ≥ 4.0.0

### Required R Packages
The script will automatically install missing packages:

```
"dplyr", "readr", "readxl", "writexl", "ggplot2",
"patchwork", "bootnet", "qgraph", "glmnet", "huge", "tibble"
```

## 📝 Usage

### Quick Start

1. **Prepare your data**: 
   - Excel format with T0 and T1 measurements
   - All variables must be continuous
   - Include any covariates as additional columns

2. **Configure the analysis**:
   - Edit `config_cross_lagged.csv` with your specific settings
   - Set data paths, variable ranges, and analysis parameters

3. **Run the script**:
```r
# Set working directory (uncomment and modify in the script)
# setwd("/path/to/your/project")

# Run the analysis
source("run_cross-lagged_network.R")
```

### Advanced Usage

The script supports complex filtering, multiple covariates, and weighted analysis. See the below for detailed options.

## ⚙️ Configuration Guide

All analysis parameters are controlled through `config_cross_lagged.csv`. Here's a detailed breakdown:

### Data Input Settings

| Parameter | Description | Example |
|-----------|-------------|---------|
| `data_path` | Path to Excel data file | `"data/your_data.xlsx"` |
| `sheet_name` | Excel sheet name | `"Sheet1"` |
| `filter_enabled` | Enable data filtering | `TRUE` or `FALSE` |
| `filters` | Filter conditions (see below) | `"Group==A;Age>=18"` |
| `use_weight` | Use inverse probability weights | `TRUE` or `FALSE` |
| `weight_col` | Weight column name | `"IPW_weight"` |
| `save_path` | Output directory | `"output/results"` |

### Variable Specification

| Parameter | Description | Example |
|-----------|-------------|---------|
| `t0_start_col` | First T0 variable column (1-based) | `5` |
| `t0_end_col` | Last T0 variable column | `14` |
| `t1_start_col` | First T1 variable column | `20` |
| `t1_end_col` | Last T1 variable column | `29` |
| `covariate_cols` | Covariate column names | `"Age,Gender,Education"` |

### Network Visualization

| Parameter | Description | Example |
|-----------|-------------|---------|
| `node_labels` | Short variable names | `"VAR1,VAR2,VAR3"` |
| `node_fullnames` | Descriptive names | `"Variable 1,Variable 2,Variable 3"` |
| `node_groups` | Domain categories | `"Domain A,Domain A,Domain B"` |
| `group_colors` | Color mapping | `"Domain A=#FF6B6B;Domain B=#4ECDC4"` |

### Bootstrap Settings

| Parameter | Description | Example |
|-----------|-------------|---------|
| `n_boots` | Bootstrap iterations | `1000` |
| `n_cores` | CPU cores for parallel processing | `4` |
| `seed` | Random seed | `123` |

### Filter Syntax

The filtering system supports multiple conditions:

```
# Single condition
"Group==A"

# Multiple conditions (AND logic)
"Group==A;Age>=18;Score<50"

# Operators supported
==    # Equal to
!=    # Not equal to
>     # Greater than
<     # Less than
>=    # Greater than or equal
<=    # Less than or equal
%in%  # Value in list

# Example with %in%
"Group %in% A,B,C"
```

## 📊 Output Files

The script generates outputs in the `save_path` directory:

### Network Visualizations
- `network_with_self_loops.png` - Network including autoregressions
- `network_no_self_loops.png` - Network without autoregressions  
- `bootstrap_network.png` - Bootstrap estimated network
- `centrality_plot.png` - In/Out expected influence centrality

### Bootstrap Diagnostics
- `bootstrap_sample_network.png` - Bootstrap sample variability
- `bootstrap_centrality.png` - Centrality stability
- `bootstrap_outEI_diff.png` - Outward expected influence differences
- `bootstrap_inEI_diff.png` - Inward expected influence differences
- `bootstrap_edge_diff.png` - Edge weight differences

### Statistical Results
- `coefficients_with_stability.xlsx` - Complete results table with:
  - Edge weights (beta coefficients)
  - Standard errors
  - 95% confidence intervals
  - P-values
  - Centrality stability coefficients
- `top_400_edges.xlsx` - Strongest connections ranked

## 📋 Data Requirements

### Data Structure
Your Excel file must contain:

1. **T0 variables**: Baseline measurements (continuous)
2. **T1 variables**: Follow-up measurements (continuous) 
3. **Covariates**: Additional variables (continuous)
4. **Weights**: If using IPW (continuous)

### Example Data Layout

| ID | Age | Gender | T0_Var1 | T0_Var2 | T0_Var3 | T1_Var1 | T1_Var2 | T1_Var3 | Weight | Group(Optional) |
|----|-----|--------|---------|---------|---------|---------|---------|---------|--------|-----------------|
| 1  | 25  | 1      | 2.3     | 1.8     | 3.1     | 2.5     | 2.0     | 3.0     | 1.2    | A               |
| 2  | 30  | 0      | 1.9     | 2.1     | 2.8     | 2.1     | 2.3     | 2.9     | 0.9    | B               |          

### Data Quality Requirements

- **No categorical variables**: All analysis variables must be continuous
- **Complete cases preferred**: Missing data in key variables may cause issues
- **Consistent scaling**: Consider standardizing variables if scales differ greatly
- **Sufficient sample size**: Minimum N=100 recommended for bootstrap stability

## 🔧 Troubleshooting

### Common Issues

**Error: "Column index out of range"**
- Check that `t0_end_col` and `t1_end_col` don't exceed your data columns
- Verify column indices are 1-based (not 0-based)

**Error: "Weight column not found"**
- Ensure `weight_col` exactly matches column name in data
- Check for spelling and case sensitivity

**Error: "No data after filtering"**
- Review filter syntax and conditions
- Check that filter column names exist in data
- Verify filter values match data values

**Error: "Configuration file not found"**
- Ensure `config_cross_lagged.csv` is in working directory
- Check file name spelling and extension

**Bootstrap fails or takes too long**
- Reduce `n_boots` for initial testing
- Decrease `n_cores` if memory issues occur
- Check for perfect correlations in data

### Performance Tips

- Start with smaller bootstrap samples (`n_boots=100`) for testing
- Use multiple cores (`n_cores=4` or higher) for large datasets
- Consider data filtering to reduce sample size if needed
