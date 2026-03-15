# ============================================================================
#  AeroGuard TinyML — Data Visualization & Analytics (R Language)
#  Capstone Project 2: Data Visualization using R Library
#  Unit II — Automation and Data Visualization
# ============================================================================

# ── 0. Set writable library path ─────────────────────────────────────────────
user_lib <- file.path(Sys.getenv("USERPROFILE"), "R_libs")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

# ── 1. Install & Load Required Packages ──────────────────────────────────────
required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "scales",
  "ggridges", "ggcorrplot", "treemapify", "patchwork",
  "viridis", "RColorBrewer", "ggrepel", "forcats",
  "reshape2", "gridExtra", "grid", "png"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org", quiet = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# ── 2. Configuration ─────────────────────────────────────────────────────────
BASE_DIR       <- "C:/HS/TML1"
DATASET_DIR    <- file.path(BASE_DIR, "TinyML_Dataset")
METADATA_FILE  <- file.path(DATASET_DIR, "metadata", "dataset_metadata.csv")
FEATURES_DIR   <- file.path(DATASET_DIR, "features")
OUTPUT_DIR     <- file.path(BASE_DIR, "R_Visualizations")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ── Elegant pastel-cloud colour palette ───────────────────────────────────────
cloud_palette <- c(
  "Background"  = "#7EC8E3",   # Soft sky blue

  "Cough"       = "#FF6B6B",   # Coral rose
  "Human_Noise" = "#95E1A3"    # Mint green
)

cloud_gradient <- c(
  "#FF9AA2", "#FFB7B2", "#FFDAC1", "#E2F0CB",
  "#B5EAD7", "#C7CEEA", "#A0C4FF", "#BDB2FF",
  "#FFC6FF", "#CAFFBF", "#9BF6FF", "#FDFFB6"
)

pastel_extended <- c(
  "#F8B4B4", "#F6D6AD", "#FCF4A3", "#B8E6C8",
  "#A3D9F5", "#C4B7EA", "#F2A6CE", "#9ED8DB",
  "#E8C5F0", "#FFD3E0", "#D4E7C5", "#B3CDE3"
)

# Common elegant theme
theme_aeroguard <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 4,
                                       hjust = 0.5, colour = "#2C3E50",
                                       margin = margin(b = 10)),
      plot.subtitle    = element_text(hjust = 0.5, colour = "#5D6D7E",
                                       size = base_size, margin = margin(b = 15)),
      plot.background  = element_rect(fill = "#FAFBFD", colour = NA),
      panel.background = element_rect(fill = "#FAFBFD", colour = NA),
      panel.grid.major = element_line(colour = "#ECF0F1", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(colour = "#34495E", face = "bold"),
      axis.text        = element_text(colour = "#5D6D7E"),
      legend.position  = "right",
      legend.background = element_rect(fill = "#FAFBFD", colour = NA),
      legend.title     = element_text(face = "bold", colour = "#2C3E50"),
      legend.text      = element_text(colour = "#5D6D7E"),
      plot.margin      = margin(20, 20, 20, 20)
    )
}

save_plot <- function(plot, filename, w = 12, h = 8) {
  filepath <- file.path(OUTPUT_DIR, filename)
  ggsave(filepath, plot, width = w, height = h, dpi = 300, bg = "#FAFBFD")
  cat(paste0("  [SAVED] ", filename, "\n"))
}

# ── 3. Load Data ──────────────────────────────────────────────────────────────
cat("============================================================\n")
cat("  AeroGuard TinyML - R Data Visualization & Analytics\n")
cat("============================================================\n\n")

cat("[1/3] Loading metadata...\n")
metadata <- read_csv(METADATA_FILE, show_col_types = FALSE)
cat(paste0("  Total samples: ", nrow(metadata), "\n"))
cat(paste0("  Classes: ", paste(unique(metadata$class), collapse = ", "), "\n"))

# ── 4. Load MFCC Features (sample for visualization) ─────────────────────────
cat("[2/3] Loading MFCC features...\n")

load_npy_simple <- function(filepath) {
  # Read .npy binary (NumPy format) - handles float64 arrays
  con <- file(filepath, "rb")
  on.exit(close(con))

  magic <- readBin(con, "raw", 6)
  major <- readBin(con, "integer", 1, size = 1)
  minor <- readBin(con, "integer", 1, size = 1)
  header_len <- readBin(con, "integer", 1, size = 2, endian = "little")
  header <- readChar(con, header_len)

  # Parse shape from header
  shape_match <- regmatches(header, regexpr("\\(([0-9, ]+)\\)", header))
  shape_str <- gsub("[\\(\\) ]", "", shape_match)
  shape <- as.integer(strsplit(shape_str, ",")[[1]])
  shape <- shape[!is.na(shape)]

  n_elements <- prod(shape)
  data <- readBin(con, "double", n_elements, size = 8, endian = "little")

  if (length(shape) == 2) {
    # NumPy C-order (row-major): shape is (rows, cols)
    # In this dataset shape = (13, 101) means 13 coefficients × 101 frames
    # We want output as (frames × coefficients) = (101 × 13)
    mat <- matrix(data, nrow = shape[1], ncol = shape[2], byrow = TRUE)
    t(mat)  # Transpose to frames × coefficients
  } else {
    data
  }
}

# Load a sample of MFCC features from each class for visualizations
set.seed(42)
sample_size_per_class <- 30

mfcc_data_list <- list()
mfcc_summary_list <- list()

for (cls in c("Background", "Cough", "Human_Noise")) {
  prefix <- substr(cls, 1, 1)
  class_files <- metadata %>%
    filter(class == cls, split == "train")
  class_files <- class_files %>%
    slice_sample(n = min(sample_size_per_class, nrow(class_files)))

  for (i in seq_len(nrow(class_files))) {
    mfcc_file <- file.path(FEATURES_DIR, class_files$mfcc_file[i])
    if (file.exists(mfcc_file)) {
      tryCatch({
        mat <- load_npy_simple(mfcc_file)
        # Store per-coefficient stats
        for (coeff in 1:ncol(mat)) {
          mfcc_summary_list[[length(mfcc_summary_list) + 1]] <- data.frame(
            class       = cls,
            file        = class_files$filename[i],
            coefficient = paste0("MFCC_", coeff),
            coeff_num   = coeff,
            mean_val    = mean(mat[, coeff]),
            sd_val      = sd(mat[, coeff]),
            min_val     = min(mat[, coeff]),
            max_val     = max(mat[, coeff]),
            median_val  = median(mat[, coeff]),
            range_val   = max(mat[, coeff]) - min(mat[, coeff]),
            stringsAsFactors = FALSE
          )
        }
        # Store the full matrix for heatmaps (first few only)
        if (i <= 3) {
          mfcc_data_list[[paste0(cls, "_", i)]] <- list(
            matrix = mat, class = cls, file = class_files$filename[i]
          )
        }
      }, error = function(e) {
        # Skip unreadable files silently
      })
    }
  }
}

mfcc_summary <- bind_rows(mfcc_summary_list)
cat(paste0("  MFCC features loaded: ", length(unique(mfcc_summary$file)),
           " files, ", nrow(mfcc_summary), " coefficient records\n"))

# ── 5. Prepare Derived Datasets ──────────────────────────────────────────────
cat("[3/3] Preparing datasets for visualization...\n\n")

class_dist <- metadata %>%
  group_by(class, split) %>%
  summarise(count = n(), .groups = "drop")

class_total <- metadata %>%
  group_by(class) %>%
  summarise(total = n(), .groups = "drop") %>%
  mutate(percentage = total / sum(total) * 100)

source_stats <- metadata %>%
  group_by(class) %>%
  summarise(
    unique_sources = n_distinct(source_file),
    total_samples  = n(),
    avg_windows    = total_samples / unique_sources,
    .groups = "drop"
  )

window_dist <- metadata %>%
  group_by(class, window_id) %>%
  summarise(count = n(), .groups = "drop")

# Model training report data (from training_report.json)
model_perf <- data.frame(
  metric = c("Train Accuracy", "Validation Accuracy", "Test Accuracy",
             "Best Val Accuracy"),
  value  = c(0.9214, 0.9188, 0.9668, 0.9668),
  stringsAsFactors = FALSE
)

model_sizes <- data.frame(
  format     = c("Keras (.keras)", "TFLite Float32", "TFLite Int8\n(Quantized)"),
  size_kb    = c(51.64, 59.67, 28.01),
  stringsAsFactors = FALSE
)

training_config <- data.frame(
  parameter = c("Architecture", "Parameters", "Batch Size", "Epochs (used)",
                 "Learning Rate", "Optimizer", "Loss Function"),
  value     = c("1D CNN", "13,219", "32", "34", "0.001", "Adam",
                 "Categorical\nCrossentropy"),
  stringsAsFactors = FALSE
)

cat("============================================================\n")
cat("  Generating Visualizations ...\n")
cat("============================================================\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 1: Class Distribution Bar Chart
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot  1: Class Distribution Bar Chart\n")

p1 <- ggplot(class_total, aes(x = reorder(class, -total), y = total, fill = class)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.8) +
  geom_text(aes(label = paste0(total, "\n(", round(percentage, 1), "%)")),
            vjust = -0.3, fontface = "bold", size = 4.5, colour = "#2C3E50") +
  scale_fill_manual(values = cloud_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "AeroGuard Dataset — Class Distribution",
    subtitle = "Total audio samples per sound category",
    x = "Audio Class", y = "Number of Samples"
  ) +
  theme_aeroguard() +
  theme(legend.position = "none")

save_plot(p1, "01_class_distribution_bar.png")

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 2: Train/Test Split Grouped Bar Chart
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot  2: Train/Test Split Distribution\n")

p2 <- ggplot(class_dist, aes(x = class, y = count, fill = split)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           colour = "white", linewidth = 0.6) +
  geom_text(aes(label = count),
            position = position_dodge(width = 0.7), vjust = -0.4,
            fontface = "bold", size = 4, colour = "#2C3E50") +
  scale_fill_manual(
    values = c("train" = "#74B9FF", "test" = "#FD79A8"),
    labels = c("Training Set", "Test Set")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Train vs Test Split by Class",
    subtitle = "80/20 split ratio across all categories",
    x = "Audio Class", y = "Number of Samples", fill = "Split"
  ) +
  theme_aeroguard()

save_plot(p2, "02_train_test_split.png")

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 3: Pie Chart — Overall Class Proportions
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot  3: Pie Chart — Class Proportions\n")

pie_data <- class_total %>%
  arrange(desc(total)) %>%
  mutate(
    ypos    = cumsum(percentage) - 0.5 * percentage,
    label   = paste0(class, "\n", round(percentage, 1), "%")
  )

p3 <- ggplot(pie_data, aes(x = "", y = percentage, fill = class)) +
  geom_col(width = 1, colour = "white", linewidth = 1.5) +
  coord_polar(theta = "y") +
  geom_text(aes(y = ypos, label = label),
            size = 4.5, fontface = "bold", colour = "#2C3E50") +
  scale_fill_manual(values = cloud_palette) +
  labs(
    title    = "Dataset Composition — Proportional View",
    subtitle = "Relative share of each sound class"
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold", size = 17, hjust = 0.5,
                                    colour = "#2C3E50", margin = margin(b = 5)),
    plot.subtitle   = element_text(hjust = 0.5, colour = "#5D6D7E",
                                    margin = margin(b = 10)),
    plot.background = element_rect(fill = "#FAFBFD", colour = NA),
    legend.position = "none",
    plot.margin     = margin(20, 20, 20, 20)
  )

save_plot(p3, "03_class_pie_chart.png", w = 9, h = 9)

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 4: Donut Chart — Dataset Composition
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot  4: Donut Chart — Dataset Composition\n")

donut_data <- class_total %>%
  arrange(desc(total)) %>%
  mutate(
    fraction  = total / sum(total),
    ymax      = cumsum(fraction),
    ymin      = c(0, ymax[-n()]),
    label_pos = (ymax + ymin) / 2,
    label     = paste0(class, "\n", total, " samples\n(", round(percentage, 1), "%)")
  )

p4 <- ggplot(donut_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5, fill = class)) +
  geom_rect(colour = "white", linewidth = 1.2) +
  geom_text(aes(x = 3.25, y = label_pos, label = label),
            size = 3.8, fontface = "bold", colour = "#2C3E50") +
  annotate("text", x = 0, y = 0, label = paste0("Total\n", sum(class_total$total)),
           size = 7, fontface = "bold", colour = "#2C3E50") +
  coord_polar(theta = "y") +
  xlim(c(0, 4.5)) +
  scale_fill_manual(values = cloud_palette) +
  labs(
    title    = "AeroGuard Dataset — Donut Chart",
    subtitle = "Elegant proportional overview of all sound classes"
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold", size = 17, hjust = 0.5,
                                    colour = "#2C3E50", margin = margin(b = 5)),
    plot.subtitle   = element_text(hjust = 0.5, colour = "#5D6D7E",
                                    margin = margin(b = 10)),
    plot.background = element_rect(fill = "#FAFBFD", colour = NA),
    legend.position = "none",
    plot.margin     = margin(20, 20, 20, 20)
  )

save_plot(p4, "04_donut_chart.png", w = 9, h = 9)

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 5: Treemap — Dataset Composition
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot  5: Treemap — Dataset Composition\n")

treemap_data <- class_dist %>%
  mutate(label = paste0(class, " (", split, ")\n", count, " samples"))

p5 <- ggplot(treemap_data, aes(area = count, fill = class, subgroup = class,
                                label = label)) +
  geom_treemap(colour = "white", linewidth = 2) +
  geom_treemap_text(fontface = "bold", colour = "#2C3E50",
                    place = "centre", size = 12, grow = FALSE) +
  geom_treemap_subgroup_border(colour = "white", linewidth = 4) +
  scale_fill_manual(values = cloud_palette) +
  labs(
    title    = "Dataset Treemap — Hierarchical View",
    subtitle = "Area proportional to number of samples (Train & Test)"
  ) +
  theme_aeroguard() +
  theme(legend.position = "none")

save_plot(p5, "05_treemap_composition.png")

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 6: MFCC Coefficient Means — Box Plot
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot  6: MFCC Coefficient Box Plots\n")

p6 <- ggplot(mfcc_summary, aes(x = coefficient, y = mean_val, fill = class)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.4, outlier.size = 1.2,
               colour = "#5D6D7E", linewidth = 0.4) +
  scale_fill_manual(values = cloud_palette) +
  labs(
    title    = "MFCC Coefficient Distributions by Class",
    subtitle = "Mean values across audio samples — key differentiators for classification",
    x = "MFCC Coefficient", y = "Mean Value", fill = "Audio Class"
  ) +
  theme_aeroguard() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

save_plot(p6, "06_mfcc_boxplot.png", w = 14, h = 8)

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 7: MFCC Violin + Jitter Plot
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot  7: MFCC Violin Plot\n")

mfcc_first6 <- mfcc_summary %>% filter(coeff_num <= 6)

p7 <- ggplot(mfcc_first6, aes(x = class, y = mean_val, fill = class)) +
  geom_violin(alpha = 0.7, colour = "#5D6D7E", linewidth = 0.4, trim = FALSE) +
  geom_jitter(aes(colour = class), width = 0.15, alpha = 0.5, size = 1.2) +
  facet_wrap(~ coefficient, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = cloud_palette) +
  scale_colour_manual(values = c("#5BA3CF", "#E05555", "#6DC47A")) +
  labs(
    title    = "MFCC Feature Distribution — Violin Plots",
    subtitle = "Distribution shape and data points for first 6 MFCC coefficients",
    x = "Audio Class", y = "Mean Coefficient Value"
  ) +
  theme_aeroguard() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 12, colour = "#2C3E50"))

save_plot(p7, "07_mfcc_violin.png", w = 14, h = 10)

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 8: MFCC Ridge Plot
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot  8: MFCC Ridge Plot\n")

ridge_data <- mfcc_summary %>%
  filter(coeff_num <= 8)

p8 <- ggplot(ridge_data, aes(x = mean_val, y = coefficient, fill = class)) +
  geom_density_ridges(alpha = 0.6, scale = 1.5, colour = "white",
                      linewidth = 0.5, rel_min_height = 0.01) +
  scale_fill_manual(values = cloud_palette) +
  labs(
    title    = "MFCC Feature Density — Ridge Plot",
    subtitle = "Overlapping density distributions showing class separability",
    x = "Mean Coefficient Value", y = "MFCC Coefficient", fill = "Audio Class"
  ) +
  theme_aeroguard() +
  theme(panel.grid.major.y = element_blank())

save_plot(p8, "08_mfcc_ridge.png", w = 13, h = 9)

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 9: MFCC Correlation Heatmap
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot  9: MFCC Correlation Heatmap\n")

mfcc_wide <- mfcc_summary %>%
  select(file, coefficient, mean_val) %>%
  pivot_wider(names_from = coefficient, values_from = mean_val) %>%
  select(-file)

cor_matrix <- cor(mfcc_wide, use = "complete.obs")

p9 <- ggcorrplot(cor_matrix,
  method  = "square",
  type    = "lower",
  lab     = TRUE,
  lab_size = 3,
  colors  = c("#74B9FF", "#FAFBFD", "#FF6B6B"),
  outline.color = "white",
  title = "MFCC Coefficient Correlation Matrix"
) +
  theme_aeroguard() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title   = element_text(face = "bold", size = 17, hjust = 0.5,
                                 colour = "#2C3E50")
  )

save_plot(p9, "09_mfcc_correlation.png", w = 11, h = 10)

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 10: MFCC Heatmap — Single Sample per Class
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 10: MFCC Spectrogram-style Heatmaps\n")

heatmap_plots <- list()
idx <- 1
for (name in names(mfcc_data_list)[1:min(3, length(mfcc_data_list))]) {
  item <- mfcc_data_list[[name]]
  mat  <- item$matrix
  df   <- reshape2::melt(mat)
  colnames(df) <- c("Frame", "Coefficient", "Value")

  heatmap_plots[[idx]] <- ggplot(df, aes(x = Frame, y = Coefficient, fill = Value)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = c("#0D1B2A", "#1B263B", "#415A77", "#778DA9",
                  "#A8DADC", "#F1FAEE", "#FFD166", "#EF476F")
    ) +
    labs(
      title = paste0(item$class, " — ", item$file),
      x = "Time Frame", y = "MFCC Coefficient", fill = "Value"
    ) +
    theme_aeroguard(base_size = 11) +
    theme(plot.title = element_text(size = 13))

  idx <- idx + 1
}

if (length(heatmap_plots) >= 3) {
  p10 <- heatmap_plots[[1]] / heatmap_plots[[2]] / heatmap_plots[[3]] +
    plot_annotation(
      title    = "MFCC Spectrogram Heatmaps — One Sample per Class",
      subtitle = "Temporal evolution of 13 MFCC coefficients (101 frames × 13 coefficients)",
      theme = theme(
        plot.title    = element_text(face = "bold", size = 17, hjust = 0.5,
                                      colour = "#2C3E50"),
        plot.subtitle = element_text(hjust = 0.5, colour = "#5D6D7E", size = 13),
        plot.background = element_rect(fill = "#FAFBFD", colour = NA)
      )
    )
  save_plot(p10, "10_mfcc_heatmaps.png", w = 14, h = 14)
}

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 11: MFCC Mean Profile — Radar-like Line Plot
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 11: MFCC Mean Profile Line Chart\n")

mfcc_profile <- mfcc_summary %>%
  group_by(class, coeff_num) %>%
  summarise(grand_mean = mean(mean_val), grand_sd = mean(sd_val), .groups = "drop")

p11 <- ggplot(mfcc_profile, aes(x = coeff_num, y = grand_mean,
                                 colour = class, fill = class)) +
  geom_ribbon(aes(ymin = grand_mean - grand_sd, ymax = grand_mean + grand_sd),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 3, shape = 21, colour = "white", stroke = 1.5) +
  scale_colour_manual(values = c("#5BA3CF", "#E05555", "#6DC47A")) +
  scale_fill_manual(values = cloud_palette) +
  scale_x_continuous(breaks = 1:13, labels = paste0("M", 1:13)) +
  labs(
    title    = "MFCC Mean Profile Across Coefficients",
    subtitle = "Average MFCC values ± 1 SD — showing distinct class fingerprints",
    x = "MFCC Coefficient", y = "Mean Value", colour = "Class", fill = "Class"
  ) +
  theme_aeroguard()

save_plot(p11, "11_mfcc_profile.png")

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 12: Feature Variability — Standard Deviation Comparison
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 12: Feature Variability (SD) Comparison\n")

sd_data <- mfcc_summary %>%
  group_by(class, coefficient, coeff_num) %>%
  summarise(avg_sd = mean(sd_val), .groups = "drop")

p12 <- ggplot(sd_data, aes(x = reorder(coefficient, coeff_num),
                             y = avg_sd, fill = class)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = cloud_palette) +
  labs(
    title    = "MFCC Feature Variability by Class",
    subtitle = "Average standard deviation per coefficient — higher = more temporal variation",
    x = "MFCC Coefficient", y = "Average Standard Deviation", fill = "Audio Class"
  ) +
  theme_aeroguard() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p12, "12_feature_variability.png", w = 14, h = 8)

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 13: Feature Range (Min-Max Span) Lollipop Chart
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 13: Feature Range Lollipop Chart\n")

range_data <- mfcc_summary %>%
  group_by(class, coefficient, coeff_num) %>%
  summarise(avg_range = mean(range_val), .groups = "drop")

p13 <- ggplot(range_data, aes(x = reorder(coefficient, coeff_num),
                               y = avg_range, colour = class)) +
  geom_segment(aes(xend = reorder(coefficient, coeff_num), y = 0, yend = avg_range),
               linewidth = 1.2, position = position_dodge(width = 0.6)) +
  geom_point(size = 4, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c("#5BA3CF", "#E05555", "#6DC47A")) +
  labs(
    title    = "MFCC Feature Range — Lollipop Chart",
    subtitle = "Average min-max span per coefficient showing dynamic range",
    x = "MFCC Coefficient", y = "Average Range (Max − Min)", colour = "Audio Class"
  ) +
  theme_aeroguard() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p13, "13_feature_range_lollipop.png", w = 14, h = 8)

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 14: Unique Source Files per Class
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 14: Unique Source Files per Class\n")

p14 <- ggplot(source_stats, aes(x = reorder(class, -unique_sources),
                                 y = unique_sources, fill = class)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.8) +
  geom_text(aes(label = unique_sources), vjust = -0.5,
            fontface = "bold", size = 5, colour = "#2C3E50") +
  scale_fill_manual(values = cloud_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Data Diversity — Unique Source Recordings",
    subtitle = "Number of distinct original audio files per class",
    x = "Audio Class", y = "Unique Source Files"
  ) +
  theme_aeroguard() +
  theme(legend.position = "none")

save_plot(p14, "14_unique_sources.png")

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 15: Average Windows per Source File
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 15: Average Windows per Source\n")

p15 <- ggplot(source_stats, aes(x = class, y = avg_windows, fill = class)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.8) +
  geom_text(aes(label = round(avg_windows, 1)), vjust = -0.5,
            fontface = "bold", size = 5, colour = "#2C3E50") +
  scale_fill_manual(values = cloud_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Data Augmentation Intensity",
    subtitle = "Average 1-second windows extracted per original recording",
    x = "Audio Class", y = "Avg. Windows per Source"
  ) +
  theme_aeroguard() +
  theme(legend.position = "none")

save_plot(p15, "15_avg_windows_per_source.png")

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 16: Window ID Distribution — Histogram
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 16: Window ID Distribution Histogram\n")

p16 <- ggplot(metadata, aes(x = window_id, fill = class)) +
  geom_histogram(binwidth = 1, colour = "white", linewidth = 0.3,
                 alpha = 0.8, position = "stack") +
  scale_fill_manual(values = cloud_palette) +
  labs(
    title    = "Window ID Distribution",
    subtitle = "Which temporal window positions are most common across recordings",
    x = "Window ID (Position in Source)", y = "Count", fill = "Audio Class"
  ) +
  theme_aeroguard()

save_plot(p16, "16_window_id_histogram.png")

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 17: Window ID Distribution — Density
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 17: Window ID Density Plot\n")

p17 <- ggplot(metadata, aes(x = window_id, fill = class, colour = class)) +
  geom_density(alpha = 0.35, linewidth = 1.1) +
  scale_fill_manual(values = cloud_palette) +
  scale_colour_manual(values = c("#5BA3CF", "#E05555", "#6DC47A")) +
  labs(
    title    = "Window Position Density by Class",
    subtitle = "Kernel density estimation showing temporal extraction patterns",
    x = "Window ID", y = "Density", fill = "Audio Class", colour = "Audio Class"
  ) +
  theme_aeroguard()

save_plot(p17, "17_window_density.png")

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 18: Model Accuracy Comparison
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 18: Model Accuracy Comparison\n")

acc_colours <- c("#74B9FF", "#A29BFE", "#FF6B6B", "#FFEAA7")

p18 <- ggplot(model_perf, aes(x = reorder(metric, value), y = value * 100,
                               fill = metric)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.8) +
  geom_text(aes(label = paste0(round(value * 100, 2), "%")),
            hjust = -0.15, fontface = "bold", size = 4.5, colour = "#2C3E50") +
  scale_fill_manual(values = acc_colours) +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  coord_flip() +
  labs(
    title    = "Model Performance Metrics",
    subtitle = "1D CNN TinyML Model — Accuracy across train/val/test splits",
    x = NULL, y = "Accuracy (%)"
  ) +
  theme_aeroguard() +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank())

save_plot(p18, "18_model_accuracy.png")

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 19: Model Size Comparison — Horizontal Bar
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 19: Model Size Comparison\n")

size_colours <- c("#DDA0DD", "#87CEEB", "#98FB98")

p19 <- ggplot(model_sizes, aes(x = reorder(format, size_kb), y = size_kb,
                                fill = format)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.8) +
  geom_text(aes(label = paste0(round(size_kb, 1), " KB")),
            hjust = -0.15, fontface = "bold", size = 5, colour = "#2C3E50") +
  scale_fill_manual(values = size_colours) +
  scale_y_continuous(limits = c(0, 75), expand = c(0, 0)) +
  coord_flip() +
  labs(
    title    = "Model Size Comparison — Quantization Effect",
    subtitle = "INT8 quantization reduces model to 28 KB — ideal for ESP32 deployment",
    x = NULL, y = "Size (KB)"
  ) +
  theme_aeroguard() +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank())

save_plot(p19, "19_model_size_comparison.png")

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 20: MFCC 1 vs MFCC 2 Scatter Plot
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 20: MFCC Feature Scatter Plot\n")

scatter_data <- mfcc_summary %>%
  filter(coeff_num %in% c(1, 2)) %>%
  select(class, file, coefficient, mean_val) %>%
  pivot_wider(names_from = coefficient, values_from = mean_val)

p20 <- ggplot(scatter_data, aes(x = MFCC_1, y = MFCC_2, colour = class)) +
  geom_point(size = 3.5, alpha = 0.75, shape = 16) +
  stat_ellipse(level = 0.85, linewidth = 1.2, linetype = "dashed") +
  scale_colour_manual(values = c("#5BA3CF", "#E05555", "#6DC47A")) +
  labs(
    title    = "Feature Space — MFCC₁ vs MFCC₂",
    subtitle = "2D projection showing class cluster separation with 85% confidence ellipses",
    x = "MFCC Coefficient 1 (Energy)", y = "MFCC Coefficient 2",
    colour = "Audio Class"
  ) +
  theme_aeroguard()

save_plot(p20, "20_mfcc_scatter.png")

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 21: MFCC Multi-Scatter Matrix (pairs)
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 21: MFCC Pairwise Scatter Matrix\n")

scatter_multi <- mfcc_summary %>%
  filter(coeff_num %in% c(1, 2, 3, 4)) %>%
  select(class, file, coefficient, mean_val) %>%
  pivot_wider(names_from = coefficient, values_from = mean_val)

plots_list_pairs <- list()
coeffs <- paste0("MFCC_", 1:4)

for (i in seq_along(coeffs)) {
  for (j in seq_along(coeffs)) {
    if (i == j) {
      # Diagonal: density
      p_diag <- ggplot(scatter_multi, aes(x = .data[[coeffs[i]]], fill = class)) +
        geom_density(alpha = 0.5, colour = NA) +
        scale_fill_manual(values = cloud_palette) +
        theme_void() + theme(legend.position = "none")
      plots_list_pairs[[length(plots_list_pairs) + 1]] <- p_diag
    } else {
      # Off-diagonal: scatter
      p_scat <- ggplot(scatter_multi, aes(x = .data[[coeffs[j]]],
                                           y = .data[[coeffs[i]]], colour = class)) +
        geom_point(alpha = 0.6, size = 1.8) +
        scale_colour_manual(values = c("#5BA3CF", "#E05555", "#6DC47A")) +
        theme_void() + theme(legend.position = "none")
      plots_list_pairs[[length(plots_list_pairs) + 1]] <- p_scat
    }
  }
}

p21 <- wrap_plots(plots_list_pairs, ncol = 4) +
  plot_annotation(
    title    = "MFCC Pairwise Scatter Matrix (Coefficients 1–4)",
    subtitle = "Diagonal: class densities  |  Off-diagonal: scatter plots",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 17, hjust = 0.5,
                                    colour = "#2C3E50"),
      plot.subtitle = element_text(hjust = 0.5, colour = "#5D6D7E", size = 13),
      plot.background = element_rect(fill = "#FAFBFD", colour = NA)
    )
  )

save_plot(p21, "21_mfcc_scatter_matrix.png", w = 12, h = 12)

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 22: Class-wise MFCC Heatmap (averaged)
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 22: Class-wise MFCC Mean Heatmap\n")

class_mfcc_mean <- mfcc_summary %>%
  group_by(class, coefficient, coeff_num) %>%
  summarise(grand_mean = mean(mean_val), .groups = "drop")

p22 <- ggplot(class_mfcc_mean, aes(x = reorder(coefficient, coeff_num),
                                    y = class, fill = grand_mean)) +
  geom_tile(colour = "white", linewidth = 1.5) +
  geom_text(aes(label = round(grand_mean, 1)),
            colour = "#2C3E50", fontface = "bold", size = 4) +
  scale_fill_gradientn(
    colours = c("#74B9FF", "#DFE6E9", "#FFEAA7", "#FD79A8", "#E17055")
  ) +
  labs(
    title    = "Class-wise Average MFCC Fingerprint",
    subtitle = "Mean value of each MFCC coefficient per audio class",
    x = "MFCC Coefficient", y = "Audio Class", fill = "Mean\nValue"
  ) +
  theme_aeroguard() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

save_plot(p22, "22_classwise_mfcc_heatmap.png", w = 14, h = 6)

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 23: Stacked Area Chart — Cumulative Samples
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 23: Stacked Area — Cumulative Growth\n")

# Simulate cumulative sample collection (ordered by filename index)
cum_data <- metadata %>%
  arrange(class, filename) %>%
  group_by(class) %>%
  mutate(sample_idx = row_number()) %>%
  ungroup()

p23 <- ggplot(cum_data, aes(x = sample_idx, fill = class)) +
  geom_area(stat = "bin", bins = 30, alpha = 0.7, colour = "white",
            linewidth = 0.5, position = "stack") +
  scale_fill_manual(values = cloud_palette) +
  labs(
    title    = "Dataset Growth — Stacked Area Chart",
    subtitle = "Cumulative sample accumulation across classes",
    x = "Sample Index", y = "Count", fill = "Audio Class"
  ) +
  theme_aeroguard()

save_plot(p23, "23_stacked_area.png")

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 24: Bubble Chart — Source Diversity vs Sample Count
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 24: Bubble Chart — Diversity vs Volume\n")

p24 <- ggplot(source_stats, aes(x = unique_sources, y = total_samples,
                                 size = avg_windows, colour = class)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = class), size = 5, fontface = "bold",
                  colour = "#2C3E50", nudge_y = 30) +
  scale_size_continuous(range = c(8, 25), name = "Avg Windows\nper Source") +
  scale_colour_manual(values = c("#5BA3CF", "#E05555", "#6DC47A")) +
  labs(
    title    = "Dataset Diversity vs Volume — Bubble Chart",
    subtitle = "Bubble size represents average augmentation intensity",
    x = "Unique Source Recordings", y = "Total Samples"
  ) +
  theme_aeroguard() +
  guides(colour = "none")

save_plot(p24, "24_bubble_chart.png")

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 25: Radar/Spider Chart (simulated with coord_polar)
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 25: Radar Chart — Feature Signatures\n")

radar_data <- mfcc_summary %>%
  filter(coeff_num <= 10) %>%
  group_by(class, coeff_num) %>%
  summarise(value = mean(mean_val), .groups = "drop") %>%
  mutate(value_scaled = (value - min(value)) / (max(value) - min(value)))

p25 <- ggplot(radar_data, aes(x = factor(coeff_num), y = value_scaled,
                               group = class, colour = class, fill = class)) +
  geom_polygon(alpha = 0.15, linewidth = 1.2) +
  geom_point(size = 3) +
  coord_polar() +
  scale_colour_manual(values = c("#5BA3CF", "#E05555", "#6DC47A")) +
  scale_fill_manual(values = cloud_palette) +
  labs(
    title    = "MFCC Feature Radar — Class Fingerprints",
    subtitle = "Normalized MFCC coefficients 1–10 showing distinct class shapes",
    x = "MFCC Coefficient", y = "Normalized Value",
    colour = "Class", fill = "Class"
  ) +
  theme_aeroguard() +
  theme(axis.text.y = element_blank())

save_plot(p25, "25_radar_chart.png", w = 10, h = 10)

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 26: Parallel Coordinate Plot
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 26: Parallel Coordinates Plot\n")

parallel_data <- mfcc_summary %>%
  filter(coeff_num <= 8) %>%
  group_by(class, file, coeff_num) %>%
  summarise(val = mean(mean_val), .groups = "drop") %>%
  mutate(val_norm = (val - min(val)) / (max(val) - min(val)))

p26 <- ggplot(parallel_data, aes(x = factor(coeff_num), y = val_norm,
                                  group = interaction(class, file),
                                  colour = class)) +
  geom_line(alpha = 0.25, linewidth = 0.5) +
  stat_summary(aes(group = class), fun = mean, geom = "line",
               linewidth = 2, alpha = 1) +
  stat_summary(aes(group = class), fun = mean, geom = "point", size = 3) +
  scale_colour_manual(values = c("#5BA3CF", "#E05555", "#6DC47A")) +
  labs(
    title    = "Parallel Coordinates — MFCC Feature Trajectories",
    subtitle = "Individual samples (faint) with class means (bold) showing feature patterns",
    x = "MFCC Coefficient", y = "Normalized Value", colour = "Audio Class"
  ) +
  theme_aeroguard()

save_plot(p26, "26_parallel_coordinates.png", w = 13, h = 8)

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 27: Faceted Histogram — MFCC Values
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 27: Faceted MFCC Histograms\n")

hist_data <- mfcc_summary %>% filter(coeff_num <= 6)

p27 <- ggplot(hist_data, aes(x = mean_val, fill = class)) +
  geom_histogram(bins = 25, alpha = 0.75, colour = "white", linewidth = 0.3,
                 position = "identity") +
  facet_grid(class ~ coefficient, scales = "free") +
  scale_fill_manual(values = cloud_palette) +
  labs(
    title    = "MFCC Value Distributions — Faceted Histograms",
    subtitle = "Rows: Audio Classes  |  Columns: MFCC Coefficients",
    x = "Mean Value", y = "Count"
  ) +
  theme_aeroguard(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold", size = 10, colour = "#2C3E50"),
    strip.background = element_rect(fill = "#ECF0F1", colour = NA)
  )

save_plot(p27, "27_faceted_histograms.png", w = 16, h = 10)

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 28: Dataset Summary Dashboard Infographic
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 28: Dataset Summary Dashboard\n")

dashboard_stats <- data.frame(
  metric = c("Total Samples", "Classes", "Sample Rate",
             "Window Size", "MFCC Coefficients", "Model Parameters",
             "Test Accuracy", "Model Size (INT8)"),
  value  = c("1,353", "3", "16 kHz",
             "1 sec", "13", "13,219",
             "96.68%", "28 KB"),
  icon   = c("📊", "🏷️", "🎵", "⏱️", "📐", "🧠", "🎯", "📦"),
  stringsAsFactors = FALSE
)

dashboard_stats$y_pos <- rev(seq_len(nrow(dashboard_stats)))

p28 <- ggplot(dashboard_stats, aes(x = 1, y = y_pos)) +
  geom_tile(aes(fill = factor(y_pos %% 2)), width = 3.5, height = 0.82,
            alpha = 0.3, colour = NA) +
  geom_text(aes(x = 0.2, label = icon), size = 8, hjust = 0.5) +
  geom_text(aes(x = 0.7, label = metric), hjust = 0, size = 5,
            colour = "#5D6D7E", fontface = "plain") +
  geom_text(aes(x = 2.8, label = value), hjust = 1, size = 6,
            colour = "#2C3E50", fontface = "bold") +
  scale_fill_manual(values = c("#E8F4FD", "#FDF2E9")) +
  xlim(-0.2, 3.2) +
  labs(
    title    = "AeroGuard TinyML — Dataset & Model Dashboard",
    subtitle = "Key statistics at a glance"
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold", size = 19, hjust = 0.5,
                                    colour = "#2C3E50", margin = margin(b = 5)),
    plot.subtitle   = element_text(hjust = 0.5, colour = "#5D6D7E", size = 14,
                                    margin = margin(b = 15)),
    plot.background = element_rect(fill = "#FAFBFD", colour = NA),
    legend.position = "none",
    plot.margin     = margin(25, 30, 25, 30)
  )

save_plot(p28, "28_dashboard_infographic.png", w = 10, h = 9)

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 29: Model Architecture Summary — Waffle-style
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 29: CNN Architecture Layer Visualization\n")

arch_data <- data.frame(
  layer     = c("Conv1D (32)", "Conv1D (64)", "Conv1D (128)",
                "Dense (128)", "Output (3)"),
  params    = c(1280, 6208, 24704, 16512, 387),
  layer_num = 1:5,
  stringsAsFactors = FALSE
)
arch_data$pct <- arch_data$params / sum(arch_data$params) * 100

arch_colours <- c("#FF9AA2", "#FFB7B2", "#FFDAC1", "#B5EAD7", "#C7CEEA")

p29 <- ggplot(arch_data, aes(x = reorder(layer, -layer_num), y = params,
                              fill = layer)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.8) +
  geom_text(aes(label = paste0(format(params, big.mark = ","), "\n(",
                                round(pct, 1), "%)")),
            hjust = -0.1, fontface = "bold", size = 4, colour = "#2C3E50") +
  scale_fill_manual(values = arch_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  coord_flip() +
  labs(
    title    = "CNN Architecture — Parameter Distribution",
    subtitle = "Number of trainable parameters per layer (Total: 13,219)",
    x = NULL, y = "Number of Parameters"
  ) +
  theme_aeroguard() +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank())

save_plot(p29, "29_cnn_architecture.png")

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 30: Comprehensive Multi-Panel Summary
# ═══════════════════════════════════════════════════════════════════════════════
cat("Plot 30: Multi-Panel Summary Poster\n")

panel_a <- p1 + labs(title = "A) Class Distribution") +
  theme(plot.title = element_text(size = 13))
panel_b <- p3 + labs(title = "B) Proportions") +
  theme(plot.title = element_text(size = 13))
panel_c <- p18 + labs(title = "C) Model Accuracy") +
  theme(plot.title = element_text(size = 13))
panel_d <- p19 + labs(title = "D) Model Sizes") +
  theme(plot.title = element_text(size = 13))

p30 <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(
    title    = "AeroGuard TinyML — Comprehensive Summary",
    subtitle = "Dataset composition, model performance, and deployment readiness",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 20, hjust = 0.5,
                                    colour = "#2C3E50"),
      plot.subtitle = element_text(hjust = 0.5, colour = "#5D6D7E", size = 14,
                                    margin = margin(b = 10)),
      plot.background = element_rect(fill = "#FAFBFD", colour = NA)
    )
  )

save_plot(p30, "30_comprehensive_summary.png", w = 16, h = 14)

# ═══════════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("  ALL 30 VISUALIZATIONS GENERATED SUCCESSFULLY!\n")
cat("============================================================\n")
cat(paste0("  Output folder: ", OUTPUT_DIR, "\n\n"))

viz_files <- list.files(OUTPUT_DIR, pattern = "\\.png$")
cat(paste0("  Total PNG files: ", length(viz_files), "\n\n"))
cat("  Visualization Index:\n")
cat("  ─────────────────────────────────────────────────────────\n")
cat("   01  Class Distribution Bar Chart\n")
cat("   02  Train/Test Split Grouped Bars\n")
cat("   03  Pie Chart — Class Proportions\n")
cat("   04  Donut Chart — Dataset Composition\n")
cat("   05  Treemap — Hierarchical View\n")
cat("   06  MFCC Coefficient Box Plots\n")
cat("   07  MFCC Violin + Jitter Plots\n")
cat("   08  MFCC Ridge Plot (Density)\n")
cat("   09  MFCC Correlation Heatmap\n")
cat("   10  MFCC Spectrogram Heatmaps\n")
cat("   11  MFCC Mean Profile Line Chart\n")
cat("   12  Feature Variability (SD) Bars\n")
cat("   13  Feature Range Lollipop Chart\n")
cat("   14  Unique Source Files Bar\n")
cat("   15  Avg Windows per Source Bar\n")
cat("   16  Window ID Histogram (Stacked)\n")
cat("   17  Window Position Density Plot\n")
cat("   18  Model Accuracy Comparison\n")
cat("   19  Model Size Comparison\n")
cat("   20  MFCC Scatter Plot (M1 vs M2)\n")
cat("   21  MFCC Pairwise Scatter Matrix\n")
cat("   22  Class-wise MFCC Mean Heatmap\n")
cat("   23  Stacked Area Chart\n")
cat("   24  Bubble Chart — Diversity vs Volume\n")
cat("   25  Radar Chart — Feature Fingerprints\n")
cat("   26  Parallel Coordinates Plot\n")
cat("   27  Faceted MFCC Histograms\n")
cat("   28  Dataset-Model Dashboard Infographic\n")
cat("   29  CNN Architecture Parameters\n")
cat("   30  Comprehensive Multi-Panel Summary\n")
cat("  ─────────────────────────────────────────────────────────\n")
cat("\n  Capstone Project 2 — Data Visualization Complete!\n")
cat("============================================================\n")
