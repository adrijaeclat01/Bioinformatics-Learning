

setwd("C:/Users/Adrija/OneDrive/Desktop/Research Works/PSM Documents/converted_file.xlsx")
getwd

df <- read.delim(file.choose(), header = TRUE)
head(df)
colnames(df)
dim(df)
summary(df)

summary(df$X..Change)
head(df$X..Change)


install.packages(c("dplyr","ggplot2","pheatmap","readr","tidyr","matrixStats"), dependencies=TRUE)

library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(tidyr)
library(matrixStats)


df$X..Change <- as.numeric(as.character(df$X..Change))  

summary(df$X..Change)
hist(df$X..Change, breaks=80, main="Distribution of X..Change", xlab="Change")
boxplot(df$X..Change, horizontal = TRUE)

up_cutoff <- 30
down_cutoff <- -30
upregulated <- df %>% filter(X..Change >= up_cutoff)
downregulated <- df %>% filter(X..Change <= down_cutoff)

nrow(upregulated)
nrow(downregulated)

top_up <- df %>% arrange(desc(X..Change)) %>% slice_head(n=20)
top_up

top_down <- df %>% arrange(X..Change) %>% slice_head(n=20)
top_down

boxplot(df$X..Change, horizontal = TRUE,
        main = "Boxplot of X..Change with thresholds", xlab = "X..Change")
abline(v = up_cutoff, col = "red", lwd = 2, lty = 2)
abline(v = down_cutoff, col = "blue", lwd = 2, lty = 2)

library(ggplot2)
ggplot(df, aes(x = X..Change)) +
  geom_histogram(bins = 120) +
  geom_vline(xintercept = c(down_cutoff, up_cutoff), color = "red", linetype = "dashed") +
  labs(title = "Histogram of change with thresholds", x = "X..Change")

df_rank <- df %>% mutate(absChange = abs(X..Change)) %>% arrange(desc(absChange)) %>% mutate(rank = row_number())
ggplot(df_rank %>% slice_head(n = 500), aes(x = rank, y = X..Change)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = c(down_cutoff, up_cutoff), color = "red", linetype = 2) +
  labs(title = "Top 500 proteins by absolute change", x = "Rank", y = "X..Change")

if(!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)




if(!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
if(!requireNamespace("ggrepel", quietly=TRUE)) install.packages("ggrepel")
library(ggplot2); library(ggrepel); library(dplyr)


change_col <- "X..Change"
p_col      <- "Pvalue"
q_col      <- "Qvalue"   # we'll prefer q_col for significance
up_cutoff  <- 30
down_cutoff <- -30
q_thresh   <- 0.05       # FDR threshold


df[[change_col]] <- as.numeric(as.character(df[[change_col]]))
df[[p_col]]      <- as.numeric(as.character(df[[p_col]]))
df[[q_col]]      <- as.numeric(as.character(df[[q_col]]))


volc <- df %>%
  filter(!is.na(.data[[change_col]])) %>%
  mutate(neglog10p = -log10(.data[[p_col]]),
         sig = case_when(
           (.data[[change_col]] >= up_cutoff) & (.data[[q_col]] < q_thresh) ~ "Up_sig",
           (.data[[change_col]] <= down_cutoff) & (.data[[q_col]] < q_thresh) ~ "Down_sig",
           TRUE ~ "NotSig"
         ))


cat("Points used:", nrow(volc), "\n")
cat("Significant Up:", sum(volc$sig=="Up_sig"), "  Significant Down:", sum(volc$sig=="Down_sig"), "\n")


p <- ggplot(volc, aes_string(x = change_col, y = "neglog10p")) +
  geom_point(aes(color = sig), alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("Up_sig"="red","Down_sig"="blue","NotSig"="grey70")) +
  geom_vline(xintercept = c(down_cutoff, up_cutoff), linetype="dashed") +
  geom_hline(yintercept = -log10(q_thresh), linetype="dashed") +
  labs(title = "Volcano plot (Qvalue < 0.05 & |Change| >= 30)",
       x = paste0(change_col, " (difference)"),
       y = "-log10(Pvalue)") +
  theme_minimal()


top10 <- volc %>% arrange(desc(neglog10p * abs(.data[[change_col]]))) %>% slice_head(n = 10)
label_col <- if("ProteinGroups" %in% colnames(df)) "ProteinGroups" else if("Genes" %in% colnames(df)) "Genes" else NA

if(!is.na(label_col)){
  p <- p + geom_text_repel(data = top10, aes_string(label = label_col), size = 3)
} else {
  p <- p + geom_text_repel(data = top10, aes(label = rownames(top10)), size = 3)
}

print(p)
ggsave("volcano_labeled_top10.png", p, width = 7, height = 6, dpi = 300)
cat("Saved volcano_labeled_top10.png to", getwd(), "\n")

install.packages("gprofiler2")
library(gprofiler2)

up_genes <- df %>%
  filter(X..Change >= 30 & Qvalue < 0.05) %>%
  pull(Genes) %>%
  unique() %>%
  na.omit()

down_genes <- df %>%
  filter(X..Change <= -30 & Qvalue < 0.05) %>%
  pull(Genes) %>%
  unique() %>%
  na.omit()

go_up <- gost(query = up_genes,
              organism = "hsapiens",
              sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"))

go_up

head(go_up$result, 10)

setwd("C:/Users/Adrija/OneDrive/Desktop/Research Works/PSM Documents")
install.packages("readxl")
library(readxl)

df <- read_excel("converted_file.xlsx")

head(df)
colnames(df)
dim(df)
summary(df)

library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(tidyr)
library(matrixStats)
library(ggrepel)

df$Ratio   <- as.numeric(as.character(df$Ratio))
df$Pvalue  <- as.numeric(as.character(df$Pvalue))
df$Qvalue  <- as.numeric(as.character(df$Qvalue))

summary(df$Ratio)
head(df$Ratio)

hist(df$Ratio, breaks = 80,
     main = "Distribution of Ratio", xlab = "Ratio")

boxplot(df$Ratio, horizontal = TRUE,
        main = "Boxplot of Ratio", xlab = "Ratio")

up_cutoff   <- 30
down_cutoff <- -30

upregulated   <- df %>% filter(Ratio >= up_cutoff)
downregulated <- df %>% filter(Ratio <= down_cutoff)

nrow(upregulated)
nrow(downregulated)

top_up <- df %>% arrange(desc(Ratio)) %>% slice_head(n = 20)
top_down <- df %>% arrange(Ratio) %>% slice_head(n = 20)


df_rank <- df %>%
  mutate(absRatio = abs(Ratio)) %>%
  arrange(desc(absRatio)) %>%
  mutate(rank = row_number())

ggplot(df_rank %>% slice_head(n = 500),
       aes(x = rank, y = Ratio)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = c(down_cutoff, up_cutoff),
             linetype = 2, color = "red") +
  labs(title = "Top 500 Proteins by Absolute Ratio Change",
       x = "Rank", y = "Ratio")

library(dplyr)
library(ggplot2)
library(ggrepel)

df$`AVG Log2 Ratio` <- as.numeric(df$`AVG Log2 Ratio`)
df$Pvalue <- as.numeric(df$Pvalue)
df$Qvalue <- as.numeric(df$Qvalue)

logFC_col <- "AVG Log2 Ratio"

volc <- df %>%
  filter(!is.na(.data[[logFC_col]])) %>%
  mutate(
    neglog10p = -log10(Pvalue),
    sig = case_when(
      (.data[[logFC_col]] >= 1 & Qvalue < 0.05) ~ "Up_sig",
      (.data[[logFC_col]] <= -1 & Qvalue < 0.05) ~ "Down_sig",
      TRUE ~ "NotSig"
    )
  )

p <- ggplot(volc, aes(x = `AVG Log2 Ratio`, y = neglog10p)) +
  geom_point(aes(color = sig), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Up_sig"="red", "Down_sig"="blue", "NotSig"="grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype="dashed") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  labs(title = "Volcano Plot", x = "AVG Log2 Ratio (log2 fold change)", y = "-log10(P-value)") +
  theme_minimal()

print(p)


