
## Load libraries
library(ggplot2)
library(ggrepel)

## Load data
dea_res <- readRDS(file = "/path/to/data/output/dea_res.rds")
vst_counts_long <- readRDS(file = "/path/to/data/output/vst_transformed_counts_long.rds")

# ----------------------
# Volcano Plot
# ----------------------
p1 <- ggplot() +
  geom_point(data = dea_res, 
             aes(x = Bcl2.gsvsBcl2.const_logFC, 
                 y = -log10(Bcl2.gsvsBcl2.const_LIMMA.Adj.P.value),
                 color = ifelse(Gene == "Bcl2", "highlight", "background")),
             size = 0.5) +
  scale_color_manual(breaks = c("highlight", "background"),
                     values = c("red", "gray")) +
  geom_text_repel(data = subset(dea_res, Gene == "Bcl2"),
                  aes(label = Gene),
                  size = 4, color = "steelblue") +
  labs(x = expression(Log[2]*"(Inactive Gene Sensor / Constitutive)"),
       y = expression(-Log[10]*"(Adj. P-Value)")) +
  scale_x_continuous(limits = c(-7.5, 7.5)) +
  scale_y_continuous(limits = c(0, 30),
                     expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none")
ggsave(plot = p1, 
       filename = "/path/to/data/output/volcano.png", 
       width = 5, height = 3, dpi = 300)

# ----------------------
# Jitter Plot
# ----------------------
p2 <- ggplot() +
  geom_jitter(data = subset(vst_counts_long, GeneSymbol == "Bcl2"), 
              aes(x = Group, y = Expression, color = Group), 
              width = 0.1) +
  scale_x_discrete(labels = c("WT" = "Wild Type",
                              "Bcl2.const" = "Constitutive",
                              "Bcl2.gs" = "Inactive Gene Sensor")) +
  labs(x = "", y = "Normalized Expression (VST)") +
  ggtitle("Bcl2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave(plot = p2, 
       filename = "/path/to/data/output/jitter.png", 
       width = 5, height = 3, dpi = 300)
