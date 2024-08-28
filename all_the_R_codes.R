# figure 4a
library(tidyverse)

methods <- c("XCMS", "Mzmine 3", "MS-DIAL", "MetCohort")
features <- c(16532, 40596, 22073, 15910)
hit <- c(425, 435, 383, 432)
labels <- c('425(97.5%)', '435(99.8%)', '383(87.8%)','432(99.1%)')

df <- data.frame(methods = factor(methods, levels = methods), features, hit, labels)

scaling_factor <- max(features) / max(hit) / 1.2

df_long <- data.frame(
  methods = factor(rep(methods, 2), levels = methods),
  measurement_type = factor(rep(c("Total features", "Hit number"), each = 4), 
                            levels = c("Total features", "Hit number")),
  value = c(features, hit * scaling_factor)
)

p <- ggplot(df_long, aes(x = methods, y = value, fill = measurement_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = c(features, labels)), 
            vjust = -0.9, 
            hjust = ifelse(df_long$measurement_type == "Total features", 1.15, 0.15),
            size = 2.2) +
  scale_y_continuous(name = "Total features", 
                     limits = c(0, max(features) * 1.1), 
                     sec.axis = sec_axis(~./scaling_factor, name = "Hit number")) +
  scale_fill_manual(values = c("Total features" = "#a6cee3", "Hit number" = "#1f78b4")) +
  labs(x = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.text.y = element_blank(),  # Remove y-axis labels
        # axis.ticks.y = element_blank(), # Remove y-axis ticks
        axis.text.x = element_text(colour = 'black', size = 8),
        legend.position = c(0.84, 0.99),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.background = element_rect(fill="white", linewidth=0.2),
        legend.margin = margin(0, 0.3, 0, 0.3),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black", size = 0))

print(p)
ggsave('D:/Metabolomics2022/文章/peakmat文章/peakmat article 20240513/pictures/figure4/mix39 detection rate_bar_2.tiff', plot = p, width = 4, height = 3, dpi = 300, compression = "lzw")

# figure 4b, 4c
library(extrafont)
library(VennDiagram)
# font_import()
# 1: MetCohort; 2: mzmine
venn.plot <- draw.pairwise.venn(area1 = 6965,
                                area2 = 7371,
                                cross.area = 4726,
                                category = c('MetCohort', 'MZmine 3'),
                                fill = c("#1f78b4", "#b2df8a"), 
                                col = c("#1f78b4", "#b2df8a"), 
                                fontfamily = "arial",
                                cat.fontfamily = "arial",
                                cat.pos = c(325, 35),
                                cat.dist = c(0.045,0.045)
)
tiff(
  filename = "D:/Metabolomics2022/文章/peakmat文章/peakmat article 20240513/pictures/figure4/venn_peak detection msv84790 MetCohort_mzmine-20240523.tiff",
  compression = "lzw",
  width = 600,
  height = 600,
  res = 200
)
grid.draw(venn.plot)
dev.off()

# 1: MetCohort; 2: xcms
venn.plot <- draw.pairwise.venn(area1 = 6965,
                                area2 = 4095,
                                cross.area = 3403,
                                category = c('MetCohort', 'XCMS'),
                                fill = c("#1f78b4", "#b2df8a"), 
                                col = c("#1f78b4", "#b2df8a"), 
                                fontfamily = "arial",
                                cat.fontfamily = "arial",
                                cat.pos = c(325, 35),
                                cat.dist = c(0.045,0.045)
)
tiff(
  filename = "D:/Metabolomics2022/文章/peakmat文章/peakmat article 20240513/pictures/figure4/venn_peak detection msv84790 MetCohort_xcms-20240523.tiff",
  compression = "lzw",
  width = 600,
  height = 600,
  res = 200
)
grid.draw(venn.plot)
dev.off()

# 1: MetCohort; 2: asari
venn.plot <- draw.pairwise.venn(area1 = 6965,
                                area2 = 1015,
                                cross.area = 706,
                                category = c('MetCohort', 'asari'),
                                fill = c("#1f78b4", "#b2df8a"), 
                                col = c("#1f78b4", "#b2df8a"), 
                                fontfamily = "arial",
                                cat.fontfamily = "arial",
                                cat.pos = c(325, 24),
                                cat.dist = c(0.045,0.045)
)
tiff(
  filename = "D:/Metabolomics2022/文章/peakmat文章/peakmat article 20240513/pictures/figure4/venn_peak detection msv84790 MetCohort_asari-20240523.tiff",
  compression = "lzw",
  width = 600,
  height = 600,
  res = 200
)
grid.draw(venn.plot)
dev.off()

# figure 4e, 4f (violin plot)
library(readxl)
df <- read_excel('D:\\Experiments\\peakmat_evaluation\\msv84790\\peakmat_vs_xcms_20240523.xlsx')

df$unique_file <- factor(df$unique_file, 
                         levels=c("{'xcms-4095', 'peakmat-v3-7066'}",
                                  "{'xcms-4095'}",
                                  "{'peakmat-v3-7066'}"))

p <- ggplot(df, aes(x=unique_file, y=`entropy index`, fill=unique_file)) + 
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1) +
  labs(y='Entropy Index', x=NULL) +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10, color="black"),  # Customize x-axis text size and color
        text = element_text(family="Arial")) +
  scale_x_discrete(labels=c("{'xcms-4095', 'peakmat-v3-7066'}"="MetCohort & XCMS",
                            "{'xcms-4095'}"="Only XCMS",
                            "{'peakmat-v3-7066'}"="Only MetCohort")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

print(p)
ggsave('D:/Metabolomics2022/文章/peakmat文章/peakmat article 20240513/pictures/figure4/MetCohort_xcms_violin_20240523.tiff',
       width = 4, height = 3, plot=p, dpi=300, bg='white')

df <- read_excel('D:\\Experiments\\peakmat_evaluation\\msv84790\\peakmat_vs_mzmine_20240523.xlsx')

df$unique_file <- factor(df$unique_file, 
                         levels=c("{'mzmine-7699', 'peakmat-v3-7066'}",
                                  "{'mzmine-7699'}",
                                  "{'peakmat-v3-7066'}"))
p <- ggplot(df, aes(x=unique_file, y=`entropy index`, fill=unique_file)) + 
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1) +
  labs(y='Entropy Index', x=NULL) +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10, color="black"),  # Customize x-axis text size and color
        text = element_text(family="Arial")) +
  scale_x_discrete(labels=c("{'mzmine-7699', 'peakmat-v3-7066'}"="MetCohort & MZmine",
                            "{'mzmine-7699'}"="Only MZmine",
                            "{'peakmat-v3-7066'}"="Only MetCohort")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

print(p)
ggsave('D:/Metabolomics2022/文章/peakmat文章/peakmat article 20240513/pictures/figure4/MetCohort_mzmine_violin_20240523.tiff',
       width = 4, height = 3, plot=p, dpi=300, bg='white')


df <- read_excel('D:\\Experiments\\peakmat_evaluation\\msv84790\\peakmat_vs_asari_20240523.xlsx')
df <- read_excel('D:\\Experiments\\peakmat_evaluation\\msv84790\\peakmat_vs_asari_20240523_remove_0_height.xlsx')

df$unique_file <- factor(df$unique_file, 
                         levels=c("{'asari_msv84790-1070', 'peakmat-v3-7066'}",
                                  "{'asari_msv84790-1070'}",
                                  "{'peakmat-v3-7066'}"))
p <- ggplot(df, aes(x=unique_file, y=`entropy index`, fill=unique_file)) + 
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1) +
  labs(y='Entropy Index', x=NULL) +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10, color="black"),  # Customize x-axis text size and color
        text = element_text(family="Arial")) +
  scale_x_discrete(labels=c("{'asari_msv84790-1070', 'peakmat-v3-7066'}"="MetCohort & asari",
                            "{'asari_msv84790-1070'}"="Only asari",
                            "{'peakmat-v3-7066'}"="Only MetCohort")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

print(p)
ggsave('D:/Metabolomics2022/文章/peakmat文章/peakmat article 20240513/pictures/figure4/MetCohort_asari_violin_20240523.tiff',
       width = 4, height = 3, plot=p, dpi=300, bg='white')
ggsave('D:/Metabolomics2022/文章/peakmat文章/peakmat article 20240513/pictures/figure4/MetCohort_asari_violin_20240523_remove_0_height.tiff',
       width = 4, height = 3, plot=p, dpi=300, bg='white')



## height distribution
df <- read_excel('D:\\Experiments\\peakmat_evaluation\\msv84790\\peakmat_vs_xcms_20240523.xlsx')

df$unique_file <- factor(df$unique_file, 
                         levels=c("{'xcms-4095', 'peakmat-v3-7066'}",
                                  "{'xcms-4095'}",
                                  "{'peakmat-v3-7066'}"))
df$log_height <- log(df$height)

p <- ggplot(df, aes(x=unique_file, y=log_height, fill=unique_file)) + 
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1) +
  labs(y='Log(Peak Intensity)', x=NULL) +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10, color="black"),  # Customize x-axis text size and color
        text = element_text(family="Arial")) +
  scale_x_discrete(labels=c("{'xcms-4095', 'peakmat-v3-7066'}"="MetCohort & XCMS",
                            "{'xcms-4095'}"="Only XCMS",
                            "{'peakmat-v3-7066'}"="Only MetCohort")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

print(p)
ggsave('D:/Metabolomics2022/文章/peakmat文章/peakmat article 20240513/pictures/figure4/MetCohort_xcms_int_violin_20240523.tiff',
       width = 4, height = 3, plot=p, dpi=300, bg='white')

df <- read_excel('D:\\Experiments\\peakmat_evaluation\\msv84790\\peakmat_vs_mzmine_20240523.xlsx')
df$log_height <- log(df$height)

df$unique_file <- factor(df$unique_file, 
                         levels=c("{'mzmine-7699', 'peakmat-v3-7066'}",
                                  "{'mzmine-7699'}",
                                  "{'peakmat-v3-7066'}"))
p <- ggplot(df, aes(x=unique_file, y=log_height, fill=unique_file)) + 
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1) +
  labs(y='Log(Peak intensity)', x=NULL) +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10, color="black"),  # Customize x-axis text size and color
        text = element_text(family="Arial")) +
  scale_x_discrete(labels=c("{'mzmine-7699', 'peakmat-v3-7066'}"="MetCohort & MZmine",
                            "{'mzmine-7699'}"="Only MZmine",
                            "{'peakmat-v3-7066'}"="Only MetCohort")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

print(p)
ggsave('D:/Metabolomics2022/文章/peakmat文章/peakmat article 20240513/pictures/figure4/MetCohort_mzmine_int_violin_20240523.tiff',
       width = 4, height = 3, plot=p, dpi=300, bg='white')


df <- read_excel('D:\\Experiments\\peakmat_evaluation\\msv84790\\peakmat_vs_asari_20240523.xlsx')
df <- read_excel('D:\\Experiments\\peakmat_evaluation\\msv84790\\peakmat_vs_asari_20240523_remove_0_height.xlsx')

df$unique_file <- factor(df$unique_file, 
                         levels=c("{'asari_msv84790-1070', 'peakmat-v3-7066'}",
                                  "{'asari_msv84790-1070'}",
                                  "{'peakmat-v3-7066'}"))
df$log_height <- log(df$height)
p <- ggplot(df, aes(x=unique_file, y=log_height, fill=unique_file)) + 
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1) +
  labs(y='Log(Peak intensity)', x=NULL) +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10, color="black"),  # Customize x-axis text size and color
        text = element_text(family="Arial")) +
  scale_x_discrete(labels=c("{'asari_msv84790-1070', 'peakmat-v3-7066'}"="MetCohort & asari",
                            "{'asari_msv84790-1070'}"="Only asari",
                            "{'peakmat-v3-7066'}"="Only MetCohort")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

print(p)
ggsave('D:/Metabolomics2022/文章/peakmat文章/peakmat article 20240513/pictures/figure4/MetCohort_asari_int_violin_20240523.tiff',
       width = 4, height = 3, plot=p, dpi=300, bg='white')
ggsave('D:/Metabolomics2022/文章/peakmat文章/peakmat article 20240513/pictures/figure4/MetCohort_asari_int_violin_20240523_remove_0_height.tiff',
       width = 4, height = 3, plot=p, dpi=300, bg='white')


# supplementary figure 15
# SP61 data, venn diagram between Peakmat and XCMS
venn.plot <- draw.pairwise.venn(area1 = 5172,
                                area2 = 6807,
                                cross.area = 3584,
                                category = c('MetCohort', 'XCMS'),
                                fill = c("#1f78b4", "#b2df8a"), 
                                col = c("#1f78b4", "#b2df8a"), 
                                fontfamily = "arial",
                                cat.fontfamily = "arial",
                                cat.pos = c(325, 35),
                                cat.dist = c(0.045,0.045)
)
tiff(
  filename = "D:/Metabolomics2022/文章/peakmat文章/peakmat article 20240513/pictures/figure5/venn_peak detection sp61 MetCohort_xcms-20240527.tiff",
  compression = "lzw",
  width = 600,
  height = 600,
  res = 200
)
grid.draw(venn.plot)
dev.off()


# intensity statistics
library(tidyverse)
library(readxl)
df <- read_excel('D:\\Experiments\\peakmat_evaluation\\msv84790\\peakmat_vs_xcms_20240523.xlsx')

df$unique_file <- factor(df$unique_file, 
                         levels=c("{'xcms-4095', 'peakmat-v3-7066'}",
                                  "{'xcms-4095'}",
                                  "{'peakmat-v3-7066'}"))
median_heights <- df %>%
  group_by(unique_file) %>%
  summarize(median_height = median(height, na.rm = TRUE))

df <- read_excel('D:\\Experiments\\peakmat_evaluation\\msv84790\\peakmat_vs_mzmine_20240523.xlsx')

df$unique_file <- factor(df$unique_file, 
                         levels=c("{'mzmine-7699', 'peakmat-v3-7066'}",
                                  "{'mzmine-7699'}",
                                  "{'peakmat-v3-7066'}"))
median_heights <- df %>%
  group_by(unique_file) %>%
  summarize(median_height = median(height, na.rm = TRUE))
