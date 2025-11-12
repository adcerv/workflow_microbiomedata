# RARE BACTERIA ANALYSIS AND CORRELATION WITH FOODS
# =================================================

# 1. BACTERIAL DATA PREPROCESSING
# --------------------------------

# Group taxonomically at Genus level and filter by prevalence
rare.bacteria.tax <- rare.bacteria.sex1 %>% 
  tax_glom("Genus") %>%
  phyloseq_filter_prevalence(prev.trh = 0.05)  # 5% prevalence threshold for large samples

# Convert to meco object and calculate abundances
dataset <- phyloseq2meco(rare.bacteria.tax)
sort(colnames(dataset$sample_table))
dataset$cal_abund()

# 2. DIFFERENTIAL ANALYSIS (LEfSe)
# --------------------------------

set.seed(185553)  # For reproducibility

# Perform LEFSE analysis to identify differential taxa
t2 <- trans_diff$new(
  dataset = dataset, 
  taxa_level = "Genus",
  transformation = "CPM",           # Counts per million
  method = "lefse",                 # LEFSE method
  p_adjust_method = "none",         # No p-value adjustment
  group = "hta_nueva",              # Grouping variable
  alpha = 0.05,                     # Significance level
  lefse_min_subsam = 40             # Minimum subsampling for LEFSE
)

# Visualize differential analysis results
t2$plot_diff_bar(threshold = 1.45,  # LDA threshold
                 xtext_size = 1.3) +
  theme(axis.text.y = element_text(size = 9, color = "black"))

# 3. SELECTION OF SIGNIFICANT TAXA
# --------------------------------

# Filter top 20 taxa for each group (highest LDA)
d1 <- t2$res_diff %>%
  filter(Group == "0") %>%          # Group 0
  arrange(LDA) %>%
  tail(20)                          # 20 taxa with highest LDA

d2 <- t2$res_diff %>%
  filter(Group == "1") %>%          # Group 1
  arrange(LDA) %>%
  tail(20)

# Combine results from both groups
d3 <- rbind(d1, d2)

# 4. CORRELATION ANALYSIS WITH FOOD DATA
# --------------------------------------

# NOTE:
# The metadata object must contain the variables of interest for correlation
# These can be selected using tidyverse 
# For reading, format by placing the sample ID in the row names

# Create object for environmental/correlation analysis
t1 <- trans_env$new(dataset = dataset, add_data = metadata) 

# Calculate correlations with selected taxa
t1$cal_cor(
  use_data = "other",               # Use additional data (foods)
  other_taxa = d3$Taxa,             # Selected taxa from LEFSE
  cor_method = "spearman"           # Spearman correlation method
  # p_adjust_method = "fdr"         # Optional: FDR adjustment
)

# 5. CORRELATION VISUALIZATION
# ----------------------------

# Correlation heatmap
t1$plot_cor() +
  theme(
    axis.text.y = element_text(size = 10, color = "black", face = "italic"),
    axis.text.x = element_text(size = 12, face = "bold", color = "black", angle = 90),
    legend.key.size = unit(0.9, 'cm'),
    legend.text = element_text(size = 11)
  ) +
  scale_fill_gradient2(
    low = "#0066FF",    # Blue for negative correlations
    high = "#FF0000",   # Red for positive correlations  
    mid = "white"       # White for zero correlation
  )

# 6. RESULTS EXTRACTION
# ---------------------

# Get correlation matrix
correlation_results <- t1$res_cor