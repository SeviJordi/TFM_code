# Figuras del TFM

# Librerias 
library(tidyverse)
library(scales)
library(plotly)
library(DT)
library(marginaleffects)
library(UpSetR)
library(scales)
library(ggpubr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# aesthetics
tools_colors = c(
  clipkit = "#DC3977",
  fixcore = "#7C1D6F",
  hmmcleaner = "#F0746E",
  snippy = "#FFA500",
  subset = "gray60",
  treeshrink = "#107D55",
  trimal = "#045275")
  
# Figure 4 ####
## Data
data <- data.frame(
  ST = c(
    "ST258",
    "ST258",
    "ST258",
    "ST258",
    "ST307",
    "ST307",
    "ST307",
    "ST307"
    ),
  tool = c("Panaroo default",
           "Panaroo",
           "Panacota",
           "Roary default",
           "Panaroo default",
           "Panaroo",
           "Panacota",
           "Roary default"
           ),
  genes = c(4806, 4859, 4673, 4573, 4772, 4792, 4571, 4654)
)

core_sizes <- ggplot(data) + 
  aes(x = tool, y = genes, fill = tool) + 
  geom_col() + 
  scale_fill_viridis_d(
    option = "H",
    labels = c("PanACoTA", "Panaroo-PanACoTA", "Panaroo", "Roary"),
    guide = "none"
    ) +
  geom_text(aes(label = genes, y = genes + 150)) +
  facet_grid(~ST,
             labeller = labeller(ST = c(ST258 = "ST258", ST307 = "ST307"))) + 
  theme_bw() + 
  scale_x_discrete(
    labels = c("PanACoTA", "Panaroo-PanACoTA", "Panaroo", "Roary")
    ) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(
    fill = "Pangenome pipeline",
    x = "Pangenome pipeline",
    y = "Core genes"
    )

ggsave(plot = core_sizes,
       filename = "../../../projects/tfm_figures/Figure4.svg",
       device = "svg", 
       width = 15.5, units = "cm", height = 11
       )

# Figure 5 ####
data_sub <- 
  read_tsv("JSF_TFM_alignment_stats.tsv", show_col_types = F)

data_sub <- data_sub %>%
  rowwise() %>%
  mutate(
    Alignment_name = gsub(".aln","", Alignment_name),
    ST = strsplit(Alignment_name,"_")[[1]][1],
    pantool = strsplit(Alignment_name,"_")[[1]][2],
    curation = strsplit(Alignment_name,"_")[[1]][3]
  ) %>%
  ungroup() %>%
  mutate(
    curation = case_when(curation == "raw" ~ "subset", T  ~ curation)
    )


figure5 <- data_sub %>%
  filter(pantool != "parsnp", pantool != "roary") %>%
  group_by(ST, pantool) %>%
  mutate(snippy_var = Proportion_variable_sites[4],
         snippy_mis = Missing_percent[4]) %>%
  ungroup() %>%
  mutate(
    Proportion_variable_sites = Proportion_variable_sites -snippy_var,
    Missing_percent = Missing_percent -snippy_mis
  ) %>% 
  filter(curation != "snippy") %>%
  ggplot() +
  aes(y = Proportion_variable_sites, x = Missing_percent/100, color = curation) + 
  geom_point(size = 3) + 
  facet_grid(ST~pantool, 
             labeller = labeller(
               ST = c(ST258 = "ST258"), 
               pantool = c(panacota = "PanACoTA", 
                           panaroo = "Panaroo-PanACoTA", 
                           panaroodefault = "Panaroo", 
                           roarydefault = "Roary"))) + 
  scale_shape(guide = "none") +
  scale_color_manual(
    values = tools_colors,
    labels = c("ClipKit", "FixCore", "HmmCleaner", "Raw", "TreeShrink", "TrimAI")
    ) + 
  scale_y_continuous(labels = percent) + 
  scale_x_continuous(labels = percent) +
  theme_bw() + 
  labs(y = "Δ Perc. variable sites",
       x = " Δ Perc. Missing",
       color = "Filtering tool"
       ) + 
  geom_hline(yintercept = 0, colour = "gray30") + 
  geom_vline(xintercept = 0, colour = "gray30") +
  theme(legend.position = "bottom")

ggsave(plot = figure5,
       filename = "../../../projects/tfm_figures/Figure5.svg",
       device = "svg", 
       width = 15.5, units = "cm", height = 10
)

# Figure 6 ####
summary <- read_csv("JSF_TFM_alignment_diversities.csv")

tool_order <- 
  c("snippy","subset" ,"fixcore", "hmmcleaner", "treeshrink", "trimal", "clipkit")


st258 <- summary %>%
  filter(ST=="ST258") %>%
  group_by(ST, pan) %>%
  mutate(pi = pi - pi[4]) %>%
  ungroup() %>%
  filter(tool != "snippy") %>%
  ggplot() +
  aes(x = factor(tool, tool_order), y = pi, color = tool) + 
  geom_point(size = 3) + 
  facet_grid(~pan, scales = "free_y", labeller = labeller(
    pan = c(panacota = "PanACoTA", 
            panaroo = "Panaroo-PanACoTA", 
            panaroodefault = "Panaroo", 
            roarydefault = "Roary"))) +
  #  geom_errorbar(aes(ymin = lim_inf, ymax = lim_sup)) +
  scale_shape(guide = "none") +
  scale_color_manual(values = tools_colors, guide = "none") + 
  theme_bw() + 
  scale_x_discrete(labels = c("Raw", "FixCore", "HmmCleaner", "TreeShrink", "TrimAI", "ClipKit")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(x = "Filtering tool", y = "Δπ")

st307 <- summary %>%
  filter(ST=="ST307") %>%
  group_by(ST, pan) %>%
  mutate(pi = pi - pi[4]) %>%
  ungroup() %>%
  filter(tool != "snippy") %>%
  ggplot() +
  aes(x = factor(tool, tool_order), y = pi, color = tool) + 
  geom_point(size = 3) + 
  facet_grid(~pan, scales = "free_y", labeller = labeller(
    pan = c(panacota = "PanACoTA", 
            panaroo = "Panaroo-PanACoTA", 
            panaroodefault = "Panaroo", 
            roarydefault = "Roary"))) +
  #  geom_errorbar(aes(ymin = lim_inf, ymax = lim_sup)) +
  scale_color_manual(values = tools_colors, guide = "none") + 
  scale_shape(guide = "none") +
  theme_bw() + 
  scale_x_discrete(labels = c("Raw", "FixCore", "HmmCleaner", "TreeShrink", "TrimAI", "ClipKit")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(x = "Filtering tool", y = "Δπ")

figure6 <- ggarrange(st258, st307, align = "h", nrow = 2, labels = c("A)", "B)"))
ggsave(plot = figure6,
       filename = "../../../projects/tfm_figures/Figure6.svg",
       device = "svg", 
       width = 15.5, units = "cm", height = 11
)

# Figure 7 ####
## ST307 ####

data <- read_tsv("JSF_TFM_SNPs_ST307.tsv")

data_l <- data %>%
  filter(snippy != "*") %>%
  pivot_longer(cols = -c("gene", "seq", "len", "POS", "ALT", "reference", "snippy"), names_to = "pantool", values_to = "nt") %>%
  mutate(
    nt = case_when(nt != "*" ~ nt)
  ) %>%
  filter(!is.na(nt), snippy != nt) %>%
  mutate(
    pan = case_when(
      startsWith(pantool, "pc")  ~ "panacota",
      startsWith(pantool, "prf")  ~ "panaroodefault",
      startsWith(pantool, "pr")  ~ "panaroo",
      startsWith(pantool, "rf")  ~ "roarydefault"
    ),
    tool = case_when(
      endsWith(pantool, "_c") ~ "clipkit",
      endsWith(pantool, "_t") ~ "trimal",
      endsWith(pantool, "_ts") ~ "treeshrink",
      endsWith(pantool, "_h") ~ "hmmcleaner",
      endsWith(pantool, "_f") ~ "fixcore",
      T ~ "subset"
    )
  ) %>%
  group_by(gene,seq, POS, pan, nt) %>%
  filter("subset" %in% tool) %>%
  ungroup()

snps_307 <- data_l %>% 
  group_by(pan, tool) %>%
  summarise(n = n()) %>% 
  ggplot() + 
  aes(x = factor(tool, tool_order), y = n, fill = tool) + 
  geom_col() + 
  theme_bw() + 
  facet_wrap(~pan, nrow = 1, labeller = labeller(
    pan = c(panacota = "PanACoTA", 
            panaroo = "Panaroo-PanACoTA", 
            panaroodefault = "Panaroo", 
            roarydefault = "Roary"), ST = c(ST258 = "ST258"))) + 
  scale_x_discrete( labels = c("Raw","FixCore", "HmmCleaner","TreeShrink", "TrimAI", "ClipKit")) + 
  scale_fill_manual(values = tools_colors, guide = "none") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(y = "Nº SNPs", x = "Filtering tool")


## ST258 ####
data_258 <- read_tsv("JSF_TFM_SNPs_ST258.tsv")

data_258_l <- data_258 %>%
  filter(snippy != "*") %>%
  pivot_longer(
    cols = -c("gene", "seq", "len", "POS", "ALT", "reference", "snippy"),
    names_to = "pantool", values_to = "nt") %>%
  mutate(
    nt = case_when(nt != "*" ~ nt)
  ) %>%
  filter(!is.na(nt), snippy != nt) %>%
  mutate(
    pan = case_when(
      startsWith(pantool, "pc")  ~ "panacota",
      startsWith(pantool, "prf")  ~ "panaroodefault",
      startsWith(pantool, "pr")  ~ "panaroo",
      startsWith(pantool, "rf")  ~ "roarydefault"
    ),
    tool = case_when(
      endsWith(pantool, "_c") ~ "clipkit",
      endsWith(pantool, "_t") ~ "trimal",
      endsWith(pantool, "_ts") ~ "treeshrink",
      endsWith(pantool, "_h") ~ "hmmcleaner",
      endsWith(pantool, "_f") ~ "fixcore",
      T ~ "subset"
    )
  ) %>%
  group_by(gene,seq, POS, pan, nt) %>%
  filter("subset" %in% tool) %>%
  ungroup()

snps_258 <- data_258_l %>% 
  group_by(pan, tool) %>%
  summarise(n = n()) %>% 
ggplot() + 
  aes(x = factor(tool, tool_order), y = n, fill = tool) + 
  geom_col() + 
  theme_bw() + 
  facet_wrap(~pan, nrow = 1, labeller = labeller(
    pan = c(panacota = "PanACoTA", 
            panaroo = "Panaroo-PanACoTA", 
            panaroodefault = "Panaroo", 
            roarydefault = "Roary"), ST = c(ST258 = "ST258"))) + 
  scale_x_discrete( labels = c("Raw","FixCore", "HmmCleaner","TreeShrink", "TrimAI", "ClipKit")) +
  scale_y_continuous(labels =c("0", "100000", "200000"), breaks = c(0,100000,200000)) +
  scale_fill_manual(values = tools_colors, guide = "none") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(y = "Nº SNPs", x = "Filtering tool")


figure7 <- ggarrange(snps_258, snps_307, align = "h", nrow = 2, labels = c("A)", "B)"))
ggsave(plot = figure7,
       filename = "../../../projects/tfm_figures/Figure7.svg",
       device = "svg", 
       width = 15.5, units = "cm", height = 13
)

# Figura 8 ####

gene_div <- 
  read_tsv("JSF_TFM_gene_diversities.tsv",
           col_names = c("ST", "pan", "tool", "gene", "pi"))

gene_div <- gene_div %>%
  mutate(
    pan = case_when(
      pan == "panaroodef" ~ "panaroodefault",
      pan == "roarydef" ~ "roarydefault", T ~ pan)
  ) %>%
  group_by(ST,pan) %>%
  pivot_wider(names_from = "tool", values_from = "pi") %>%
  ungroup() %>%
  filter(!is.na(raw), !is.na(snippy)) %>%
  mutate(diff = raw - snippy)

## ST258 ####
names_258 <- read_tsv("JSF_TFM_ST258_gene_families.tsv") %>%
  pivot_longer(-c(gene_name), names_to = "pan", values_to = "gene" ) %>%
  mutate(ST = "ST258",
         gene = gsub("ST258", "ST258-", gene))

divs_258 <- gene_div %>%
  filter(ST == "ST258") %>%
  left_join(names_258) %>%
  select(pan, gene_name, snippy) %>%
  rename(gene = gene_name,
         diversity = snippy) %>%
  group_by(pan) %>%
  mutate(variable = case_when(
    diversity > median(diversity) ~ "diverse", T ~ "No_diverse"
  )) %>%
  ungroup()

var_snps_258 <- data_258_l %>%
  filter(gene != "gene_216", gene != "gene_1528") %>%
  left_join(divs_258) %>%
  group_by(pan,variable) %>%
  mutate(n_var = n_distinct(gene)) %>%
  ungroup() %>%
  group_by(pan, tool, variable, n_var) %>%
  summarise(n = n()) %>%
  ggplot() +
  aes(x = factor(tool, tool_order), y = n/n_var, fill = tool, alpha = variable) +
  geom_col(position = "dodge") +
  theme_bw() +
  facet_grid(cols = vars(pan), labeller = labeller(
    pan = c(panacota = "PanACoTA",
            panaroo = "Panaroo-PanACoTA",
            panaroodefault = "Panaroo",
            roarydefault = "Roary"), ST = c(ST258 = "ST258"))) +
  scale_alpha_manual("Variable \ngene family",values = c(1,0.7), labels = c("Yes", "No")) +
  scale_x_discrete( labels = c("FixCore", "HmmCleaner","Raw", "TreeShrink", "TrimAI", "ClipKit")) +
  scale_fill_manual(values = tools_colors, guide = "none") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(y = "Nº SNPs/gene family", x = "Filtering tool")



pct_var_258 <- data_258_l %>%
  filter(gene != "gene_216", gene != "gene_1528") %>%
  left_join(divs_258) %>%
  group_by(pan,variable) %>%
  mutate(n_var = n_distinct(gene)) %>%
  ungroup() %>%
  group_by(pan, tool, variable, n_var) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(pan, variable) %>%
  summarise(
    tool = tool,
    n_pf = n,
    fc = 1 - (n_pf/max(n_pf))
  ) %>%
  ungroup() %>%
  select(-n_pf) %>%
  pivot_wider(names_from = variable, values_from = fc) %>% 
  ggplot() +
  aes(x = diverse, y = No_diverse, color = tool, shape = pan) +
  geom_point(size = 3) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "darkred" ) +
  scale_color_manual(values = tools_colors, labels = c("ClipKit","FixCore", "HmmCleaner","Raw", "TreeShrink", "TrimAI")) +
  scale_shape(labels = c(panacota = "PanACoTA",
                         panaroo = "Panaroo-PanACoTA",
                         panaroodefault = "Panaroo",
                         roarydefault = "Roary")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlim(c(0,1)) + ylim(c(0,1)) +
  scale_x_continuous(labels = percent, limits = c(0,1)) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  labs(
    y = "%SNPs reduction no variable genes",
    x  = "%SNPs reduction variable genes",
    color = "FIltering tool",
    shape = "Pangenome approach")

## ST307 ####
names_307 <- read_tsv("JSF_TFM_ST307_gene_families.tsv") %>%
  pivot_longer(-c(gene_name), names_to = "pan", values_to = "gene" ) %>%
  mutate(ST = "ST307",
         gene = gsub("ST307", "ST307-", gene))

divs_307 <- gene_div %>%
  filter(ST == "ST307") %>%
  left_join(names_307) %>%
  select(pan, gene_name, snippy) %>%
  rename(gene = gene_name,
         diversity = snippy) %>%
  group_by(pan) %>%
  mutate(variable = case_when(
    diversity > median(diversity) ~ "diverse", T ~ "No_diverse"
  )) %>%
  ungroup()

var_snps_307 <- data_l %>%
  filter(gene != "gene_289") %>%
  left_join(divs_307) %>%
  group_by(pan,variable) %>%
  mutate(n_var = n_distinct(gene)) %>%
  ungroup() %>%
  group_by(pan, tool, variable, n_var) %>%
  summarise(n = n()) %>%
  ggplot() +
  aes(x = factor(tool, tool_order), y = n/n_var, fill = tool, alpha = variable) +
  geom_col(position = "dodge") +
  theme_bw() +
  facet_grid(cols = vars(pan), labeller = labeller(
    pan = c(panacota = "PanACoTA",
            panaroo = "Panaroo-PanACoTA",
            panaroodefault = "Panaroo",
            roarydefault = "Roary"), ST = c(ST258 = "ST258"))) +
  scale_alpha_manual("Variable \ngene family",values = c(1,0.7), labels = c("Yes", "No")) +
  scale_x_discrete( labels = c("FixCore", "HmmCleaner","Raw", "TreeShrink", "TrimAI", "ClipKit")) +
  scale_fill_manual(values = tools_colors, guide = "none") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(y = "Nº SNPs/gene family", x = "Filtering tool")


pct_var_307 <- data_l %>%
  filter(gene != "gene_289") %>%
  left_join(divs_307) %>%
  group_by(pan,variable) %>%
  mutate(n_var = n_distinct(gene)) %>%
  ungroup() %>%
  group_by(pan, tool, variable, n_var) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(pan, variable) %>%
  summarise(
    tool = tool,
    n_pf = n,
    fc = 1 - (n_pf/max(n_pf))
  ) %>%
  ungroup() %>%
  select(-n_pf) %>%
  pivot_wider(names_from = variable, values_from = fc) %>%
  ggplot() +
  aes(x = diverse, y = No_diverse, color = tool, shape = pan) +
  geom_point(size = 3) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "darkred" ) +
  scale_color_manual(values = tools_colors, labels = c("ClipKit", "FixCore", "HmmCleaner","Raw", "TreeShrink", "TrimAI")) +
  scale_shape(labels = c(panacota = "PanACoTA",
                         panaroo = "Panaroo (PanACoTA)",
                         panaroodefault = "Panaroo",
                         roarydefault = "Roary")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlim(c(0,1)) + ylim(c(0,1)) +
  scale_x_continuous(labels = percent, limits = c(0,1)) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  labs(
    y = "%SNPs reduction no variable genes",
    x  = "%SNPs reduction variable genes",
    color = "FIltering tool",
    shape = "Pangenome approach")


pa <- ggarrange(pct_var_258, pct_var_307, align = "h", common.legend = T, legend = "bottom", labels = c("C)", "D)"))
pr <- ggarrange(var_snps_258, var_snps_307, align = "h", common.legend = T, legend = "bottom", labels = c("A)", "B)"), ncol = 1)

figure8 <- ggarrange(pr, pa, ncol = 1, heights = c(1,0.7))

ggsave(plot = figure8,
       filename = "../../../projects/tfm_figures/Figure8.svg",
       device = "svg", 
       width = 15.5, units = "cm", height = 20
)


# Figure 9 ####
data_cuttree <- read_csv("JSF_TFM_clustering.csv", show_col_types = F)

figure9 <- data_cuttree %>% filter(pan != "roary", pan != "parsnp") %>%
  ggplot() + 
  aes(x = threshold, y = n_clusters, color = tool, shape = pan) + 
  geom_point() + 
  geom_line() + 
  facet_grid(ST~pan, scales = "free_y", labeller = labeller(
    pan = c(panacota = "PanACoTA", 
            panaroo = "Panaroo-PanACoTA", 
            panaroodefault = "Panaroo", 
            roarydefault = "Roary"), ST = c(ST258 = "ST258"))) + 
  scale_color_manual(values = tools_colors,  labels = c("ClipKit", "FixCore", "HmmCleaner", "Mapping", "Raw", "TreeShrink", "TrimAI")) + 
  theme_bw() + 
  labs(x = "SNPs/Mb", y = "Nº clusters", color = "Filtering Tool", shape = "Pangenome tool") + 
  scale_shape(guide = "none") + 
  theme(legend.position = "bottom")

ggsave(plot = figure9,
       filename = "../../../projects/tfm_figures/Figure9.svg",
       device = "svg", 
       width = 15.5, units = "cm", height = 10
)

# Figure 10 ####
dists <- read.csv("JSF_TFM_ST258_wRF.csv") %>%
  mutate(
    pantool = case_when(pantool == "trimmal" ~ "trimal", T ~ pantool)
  )


st258_rf <- ggplot(filter(dists, pantool != "roary", tool1 != "snippy")) +
  aes(x = factor(tool1, tool_order), y = rfdist, color = tool1) +
  geom_point(size = 3) +
  scale_color_manual(values = tools_colors, guide = "none") + 
  scale_shape(guide = "none") + 
  theme_bw() + 
  facet_wrap(~pantool, nrow = 1, labeller = labeller(
    pantool = c(panacota = "PanACoTA", 
                panaroo = "Panaroo-PanACoTA", 
                panaroo_def = "Panaroo", 
                roary_def = "Roary"))) + 
  scale_x_discrete(labels = c("Raw", "FixCore", "HmmCleaner", "TreeShrink", "TrimAI", "ClipKit")) +
  labs(x = "Filtering tool", y = "wRF distance", color = "Tool", shape = "Pangenome tool" ) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))



dists_307 <- read.csv("JSF_TFM_ST307_wRF.csv") %>%
  mutate(
    tool1 = case_when(tool1 == "trimmal" ~ "trimal", T ~ tool1)
  )

st307_rf <- ggplot(filter(dists_307, pantool != "roary", tool1 != "snippy")) +
  aes(x = factor(tool1, tool_order), y = rfdist, color = tool1) +
  geom_point(size = 3) +
  scale_color_manual(values = tools_colors, guide = "none") + 
  scale_shape(guide = "none") + 
  theme_bw() + 
  facet_wrap(~pantool, nrow = 1, labeller = labeller(
    pantool = c(panacota = "PanACoTA", 
                panaroo = "Panaroo-PanACoTA", 
                panaroo_def= "Panaroo", 
                roaryd_def = "Roary"))) + 
  scale_x_discrete(labels = c("Raw", "FixCore", "HmmCleaner", "TreeShrink", "TrimAI", "ClipKit")) +
  labs(x = "Filtering tool", y = "wRF distance", color = "Tool", shape = "Pangenome tool" ) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

figure10 <- ggarrange(st258_rf,
                      st307_rf,
                      ncol = 1,
                      align = "h",
                      labels = c("A)", "B)"))

ggsave(plot = figure10,
       filename = "../../../projects/tfm_figures/Figure10.svg",
       device = "svg", 
       width = 15.5, units = "cm", height = 11
)

# Figure 11 ####
gene_div <- 
  read_tsv("JSF_TFM_gene_diversities.tsv",
           col_names = c("ST", "pan", "tool", "gene", "pi"))

gene_div <- gene_div %>%
  mutate(
    pan = case_when(
      pan == "panaroodef" ~ "panaroodefault",
      pan == "roarydef" ~ "roarydefault", T ~ pan)
    ) %>%
  group_by(ST,pan) %>%
  pivot_wider(names_from = "tool", values_from = "pi") %>%
  ungroup() %>%
  filter(!is.na(raw), !is.na(snippy)) %>%
  mutate(diff = raw - snippy)


genes_st307 <- gene_div %>%
  filter(pan != "roary", ST == "ST307", pan != "parsnp") %>%
  mutate(outlier = case_when(diff > mean(diff) + 3*sd(diff) ~ "yes", T ~ "no")) %>%
  ggplot() + 
  aes(y = raw, x = snippy,color = outlier) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c(yes = "#1A73E8", no = "black"), guide = "none") +
  theme_bw() + 
  facet_grid(cols = vars(pan), labeller = labeller(
    pan = c(panacota = "PanACoTA", 
            panaroo = "Panaroo-PanACoTA", 
            panaroodefault = "Panaroo", 
            roarydefault = "Roary"), ST = c(ST258 = "CG258"))) + 
  theme(aspect.ratio = 1) + 
  geom_abline(aes(slope = 1, intercept = 0), color = "red", linetype = 2) +
  labs(x = "Mapping diversity", y = "Pangenome diversity") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  scale_shape(guide = "none")

genes_st258 <- gene_div %>%
  filter(pan != "roary", ST == "ST258", pan != "parsnp") %>%
  mutate(outlier = case_when(diff > mean(diff) + 3*sd(diff) ~ "yes", T ~ "no")) %>%
  ggplot() + 
  aes(y = raw, x = snippy,color = outlier) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c(yes = "#1A73E8", no = "black"), guide = "none") +
  theme_bw() + 
  facet_grid(cols = vars(pan), labeller = labeller(
    pan = c(panacota = "PanACoTA", 
            panaroo = "Panaroo-PanACoTA", 
            panaroodefault = "Panaroo", 
            roarydefault = "Roary"), ST = c(ST258 = "CG258"))) + 
  geom_abline(aes(slope = 1, intercept = 0), color = "red", linetype = 2) +
  theme(aspect.ratio = 1) + 
  labs(x = "Mapping diversity", y = "Pangenome diversity") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  scale_shape(guide = "none")

figure11 <- ggarrange(genes_st258,
                      genes_st307,
                      ncol = 1,
                      align = "h",
                      labels = c("A)", "B)"))

ggsave(plot = figure11,
       filename = "../../../projects/tfm_figures/Figure11.svg",
       device = "svg", 
       width = 15.5, units = "cm", height = 12
)




# Figura 12 ####
pos_snps_258 <- data_258_l %>%
  filter(gene != "gene_216", gene != "gene_1528") %>%
  mutate(rel_len = POS/len,
         rel_len = cut(rel_len,breaks = seq(0,1,0.1), include.lowest = T)) %>%
  ggplot() +
  aes(x = rel_len, fill = factor(tool, tool_order)) +
  geom_bar() +
  theme_bw() +
  facet_grid(cols = vars(pan), labeller = labeller(
    pan = c(panacota = "PanACoTA",
            panaroo = "Panaroo-PanACoTA",
            panaroodefault = "Panaroo",
            roarydefault = "Roary"), ST = c(ST258 = "ST258"))) +
  scale_fill_manual(values = tools_colors, labels = c("Raw","FixCore", "HmmCleaner","TreeShrink", "TrimAI", "ClipKit")) +
  scale_x_discrete(labels = paste(seq(0,100,10), "%", sep = "")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(y = "Nº SNPs", x = "Relative gene position", fill = "Filtering tool")


pos_snps_307 <- data_l %>%
  filter(gene != "gene_289") %>%
  mutate(rel_len = POS/len,
         rel_len = cut(rel_len,breaks = seq(0,1,0.1), include.lowest = T)) %>%
  ggplot() +
  aes(x = rel_len, fill = factor(tool,tool_order)) +
  geom_bar() +
  theme_bw() +
  facet_grid(cols = vars(pan), labeller = labeller(
    pan = c(panacota = "PanACoTA",
            panaroo = "Panaroo-PanACoTA",
            panaroodefault = "Panaroo",
            roarydefault = "Roary"), ST = c(ST258 = "ST258"))) +
  scale_fill_manual(values = tools_colors, labels = c("Raw","FixCore", "HmmCleaner","TreeShrink", "TrimAI", "ClipKit")) +
  scale_x_discrete(labels = paste(seq(0,100,10), "%", sep = "")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(y = "Nº SNPs", x = "Relative gene position", fill = "Filtering tool")

figure12 <- ggarrange(pos_snps_258,
                pos_snps_307,
                align = "h",
                common.legend = T,
                legend = "bottom",
                labels = c("A)", "B)"),
                ncol = 1)

ggsave(plot = figure12,
       filename = "../../../projects/tfm_figures/Figure12.svg",
       device = "svg", 
       width = 15.5, units = "cm", height = 12
)
