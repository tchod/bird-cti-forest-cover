# Bird communities warm at the same rate
# but by different pathways in open and forested landscapes

# Data -------------------------------------------------------------------------

library(tidyverse)
library(segmented)
library(glmmTMB)
library(emmeans)
library(patchwork)

site_year <- read_csv("site_year.csv", show_col_types = FALSE) %>%
  mutate(id      = factor(id),
         habitat = factor(habitat, levels = c("open", "mixed", "closed")),
         yr      = year - 2000,
         fyear   = factor(year),
         period  = factor(if_else(year <= 2012, "early", "late"),
                          levels = c("early", "late")))

nrow(site_year)                                    # plot-years
n_distinct(site_year$id)                           # plots
site_year %>% count(year) %>% summarise(min = min(n), max = max(n), mean = mean(n))
site_year %>% distinct(id, habitat) %>% count(habitat)

# Plot mean CTI along the forest-cover gradient --------------------------------

cols_hab <- c(open = "#7570b3", mixed = "#d95f02", closed = "#1b9e77")
hab_labs <- c(open = "Open", mixed = "Mixed", closed = "Closed")

plot_cti <- site_year %>%
  group_by(id, habitat) %>%
  summarise(forest = mean(forest), CTI = mean(CTI), n_years = n(), .groups = "drop")

mean(plot_cti$CTI)
plot_cti %>%
  group_by(habitat) %>%
  summarise(n_plots = n(), mean_CTI = mean(CTI), sd_CTI = sd(CTI))

ggplot(plot_cti, aes(forest, CTI, colour = habitat)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_vline(xintercept = c(25, 75), linewidth = 0.6, colour = "grey60") +
  geom_hline(yintercept = mean(plot_cti$CTI), linewidth = 0.6, colour = "grey60") +
  scale_colour_manual(values = cols_hab, labels = hab_labs, name = NULL) +
  labs(x = "Forest cover in transect buffer (%)", y = "Plot mean CTI (°C)") +
  theme_bw()

# Breakpoint in the annual mean CTI --------------------------------------------

cti_year <- site_year %>%
  group_by(year) %>%
  summarise(n_plots = n(), se = sd(CTI) / sqrt(n()), CTI = mean(CTI), .groups = "drop")

m_lin <- lm(CTI ~ year, data = cti_year)
m_seg <- segmented(m_lin, seg.Z = ~ year, psi = 2012)
summary(m_seg)
confint(m_seg)
slope(m_seg)

# Models -----------------------------------------------------------------------
# 9 responses x 2 resolutions (all plots, forest-cover classes) = 18 models
# CTI, CTI_cold, CTI_warm: Gaussian, trends in degC per year
# abundance: negative binomial, richness: Poisson, trends on the log scale

ctrl <- glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))

# distribution of abundance and richness chosen by AIC
m_ABU_pois <- glmmTMB(Ntot ~ yr * period + (1 | id), data = site_year, family = poisson())
m_ABU_nb   <- glmmTMB(Ntot ~ yr * period + (1 | id), data = site_year, family = nbinom2())
AIC(m_ABU_pois, m_ABU_nb)

m_SPE_pois <- glmmTMB(S ~ yr * period + (1 | id), data = site_year, family = poisson())
m_SPE_nb   <- glmmTMB(S ~ yr * period + (1 | id), data = site_year, family = nbinom2())
AIC(m_SPE_pois, m_SPE_nb)

# plot-years without any individual of a group drop out of that group's CTI model
d_warm <- site_year %>% filter(!is.na(CTI_warm))
d_cold <- site_year %>% filter(!is.na(CTI_cold))

# all species
m_CTI     <- glmmTMB(CTI  ~ yr * period           + (1 | id), data = site_year, control = ctrl)
m_CTI_hab <- glmmTMB(CTI  ~ yr * period * habitat + (1 | id), data = site_year, control = ctrl)
m_ABU     <- glmmTMB(Ntot ~ yr * period           + (1 | id), data = site_year, family = nbinom2(), control = ctrl)
m_ABU_hab <- glmmTMB(Ntot ~ yr * period * habitat + (1 | id), data = site_year, family = nbinom2(), control = ctrl)
m_SPE     <- glmmTMB(S    ~ yr * period           + (1 | id), data = site_year, family = poisson(), control = ctrl)
m_SPE_hab <- glmmTMB(S    ~ yr * period * habitat + (1 | id), data = site_year, family = poisson(), control = ctrl)

# cold-associated species
m_CTI_cold     <- glmmTMB(CTI_cold ~ yr * period           + (1 | id), data = d_cold, control = ctrl)
m_CTI_cold_hab <- glmmTMB(CTI_cold ~ yr * period * habitat + (1 | id), data = d_cold, control = ctrl)
m_ABU_cold     <- glmmTMB(N_cold   ~ yr * period           + (1 | id), data = site_year, family = nbinom2(), control = ctrl)
m_ABU_cold_hab <- glmmTMB(N_cold   ~ yr * period * habitat + (1 | id), data = site_year, family = nbinom2(), control = ctrl)
m_SPE_cold     <- glmmTMB(S_cold   ~ yr * period           + (1 | id), data = site_year, family = poisson(), control = ctrl)
m_SPE_cold_hab <- glmmTMB(S_cold   ~ yr * period * habitat + (1 | id), data = site_year, family = poisson(), control = ctrl)

# warm-associated species
m_CTI_warm     <- glmmTMB(CTI_warm ~ yr * period           + (1 | id), data = d_warm, control = ctrl)
m_CTI_warm_hab <- glmmTMB(CTI_warm ~ yr * period * habitat + (1 | id), data = d_warm, control = ctrl)
m_ABU_warm     <- glmmTMB(N_warm   ~ yr * period           + (1 | id), data = site_year, family = nbinom2(), control = ctrl)
m_ABU_warm_hab <- glmmTMB(N_warm   ~ yr * period * habitat + (1 | id), data = site_year, family = nbinom2(), control = ctrl)
m_SPE_warm     <- glmmTMB(S_warm   ~ yr * period           + (1 | id), data = site_year, family = poisson(), control = ctrl)
m_SPE_warm_hab <- glmmTMB(S_warm   ~ yr * period * habitat + (1 | id), data = site_year, family = poisson(), control = ctrl)

# Trends per period ------------------------------------------------------------
# all plots: ~ period; forest-cover classes: ~ period | habitat
# Gaussian models (CTI) may return lower.CL/upper.CL/t.ratio instead of
# asymp.LCL/asymp.UCL/z.ratio, depending on the emmeans version; rename them

new_names <- c(asymp.LCL = "lower.CL", asymp.UCL = "upper.CL", z.ratio = "t.ratio")

slopes_CTI     <- summary(emtrends(m_CTI, ~ period, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  rename(any_of(new_names)) %>%
  mutate(indicator = "CTI", group = "all species", habitat = "all"); slopes_CTI
slopes_CTI_hab <- summary(emtrends(m_CTI_hab, ~ period | habitat, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  rename(any_of(new_names)) %>%
  mutate(indicator = "CTI", group = "all species", habitat = as.character(habitat)); slopes_CTI_hab
slopes_ABU     <- summary(emtrends(m_ABU, ~ period, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Abundance", group = "all species", habitat = "all")
slopes_ABU_hab <- summary(emtrends(m_ABU_hab, ~ period | habitat, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Abundance", group = "all species", habitat = as.character(habitat))
slopes_SPE     <- summary(emtrends(m_SPE, ~ period, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Species richness", group = "all species", habitat = "all")
slopes_SPE_hab <- summary(emtrends(m_SPE_hab, ~ period | habitat, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Species richness", group = "all species", habitat = as.character(habitat))

slopes_CTI_cold     <- summary(emtrends(m_CTI_cold, ~ period, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  rename(any_of(new_names)) %>%
  mutate(indicator = "CTI", group = "cold", habitat = "all")
slopes_CTI_cold_hab <- summary(emtrends(m_CTI_cold_hab, ~ period | habitat, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  rename(any_of(new_names)) %>%
  mutate(indicator = "CTI", group = "cold", habitat = as.character(habitat))
slopes_ABU_cold     <- summary(emtrends(m_ABU_cold, ~ period, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Abundance", group = "cold", habitat = "all")
slopes_ABU_cold_hab <- summary(emtrends(m_ABU_cold_hab, ~ period | habitat, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Abundance", group = "cold", habitat = as.character(habitat))
slopes_SPE_cold     <- summary(emtrends(m_SPE_cold, ~ period, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Species richness", group = "cold", habitat = "all")
slopes_SPE_cold_hab <- summary(emtrends(m_SPE_cold_hab, ~ period | habitat, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Species richness", group = "cold", habitat = as.character(habitat))

slopes_CTI_warm     <- summary(emtrends(m_CTI_warm, ~ period, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  rename(any_of(new_names)) %>%
  mutate(indicator = "CTI", group = "warm", habitat = "all")
slopes_CTI_warm_hab <- summary(emtrends(m_CTI_warm_hab, ~ period | habitat, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  rename(any_of(new_names)) %>%
  mutate(indicator = "CTI", group = "warm", habitat = as.character(habitat))
slopes_ABU_warm     <- summary(emtrends(m_ABU_warm, ~ period, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Abundance", group = "warm", habitat = "all")
slopes_ABU_warm_hab <- summary(emtrends(m_ABU_warm_hab, ~ period | habitat, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Abundance", group = "warm", habitat = as.character(habitat))
slopes_SPE_warm     <- summary(emtrends(m_SPE_warm, ~ period, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Species richness", group = "warm", habitat = "all")
slopes_SPE_warm_hab <- summary(emtrends(m_SPE_warm_hab, ~ period | habitat, var = "yr"), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Species richness", group = "warm", habitat = as.character(habitat))

slopes <- bind_rows(slopes_CTI, slopes_CTI_hab, slopes_ABU, slopes_ABU_hab, slopes_SPE, slopes_SPE_hab,
                    slopes_CTI_cold, slopes_CTI_cold_hab, slopes_ABU_cold, slopes_ABU_cold_hab,
                    slopes_SPE_cold, slopes_SPE_cold_hab,
                    slopes_CTI_warm, slopes_CTI_warm_hab, slopes_ABU_warm, slopes_ABU_warm_hab,
                    slopes_SPE_warm, slopes_SPE_warm_hab)

# Change in trend between periods (late - early) --------------------------------

contrast_CTI     <- summary(pairs(emtrends(m_CTI, ~ period, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  rename(any_of(new_names)) %>%
  mutate(indicator = "CTI", group = "all species", habitat = "all")
contrast_CTI_hab <- summary(pairs(emtrends(m_CTI_hab, ~ period | habitat, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  rename(any_of(new_names)) %>%
  mutate(indicator = "CTI", group = "all species", habitat = as.character(habitat))
contrast_ABU     <- summary(pairs(emtrends(m_ABU, ~ period, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Abundance", group = "all species", habitat = "all")
contrast_ABU_hab <- summary(pairs(emtrends(m_ABU_hab, ~ period | habitat, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Abundance", group = "all species", habitat = as.character(habitat))
contrast_SPE     <- summary(pairs(emtrends(m_SPE, ~ period, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Species richness", group = "all species", habitat = "all")
contrast_SPE_hab <- summary(pairs(emtrends(m_SPE_hab, ~ period | habitat, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Species richness", group = "all species", habitat = as.character(habitat))

contrast_CTI_cold     <- summary(pairs(emtrends(m_CTI_cold, ~ period, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  rename(any_of(new_names)) %>%
  mutate(indicator = "CTI", group = "cold", habitat = "all")
contrast_CTI_cold_hab <- summary(pairs(emtrends(m_CTI_cold_hab, ~ period | habitat, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  rename(any_of(new_names)) %>%
  mutate(indicator = "CTI", group = "cold", habitat = as.character(habitat))
contrast_ABU_cold     <- summary(pairs(emtrends(m_ABU_cold, ~ period, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Abundance", group = "cold", habitat = "all")
contrast_ABU_cold_hab <- summary(pairs(emtrends(m_ABU_cold_hab, ~ period | habitat, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Abundance", group = "cold", habitat = as.character(habitat))
contrast_SPE_cold     <- summary(pairs(emtrends(m_SPE_cold, ~ period, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Species richness", group = "cold", habitat = "all")
contrast_SPE_cold_hab <- summary(pairs(emtrends(m_SPE_cold_hab, ~ period | habitat, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Species richness", group = "cold", habitat = as.character(habitat))

contrast_CTI_warm     <- summary(pairs(emtrends(m_CTI_warm, ~ period, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  rename(any_of(new_names)) %>%
  mutate(indicator = "CTI", group = "warm", habitat = "all")
contrast_CTI_warm_hab <- summary(pairs(emtrends(m_CTI_warm_hab, ~ period | habitat, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  rename(any_of(new_names)) %>%
  mutate(indicator = "CTI", group = "warm", habitat = as.character(habitat))
contrast_ABU_warm     <- summary(pairs(emtrends(m_ABU_warm, ~ period, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Abundance", group = "warm", habitat = "all")
contrast_ABU_warm_hab <- summary(pairs(emtrends(m_ABU_warm_hab, ~ period | habitat, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Abundance", group = "warm", habitat = as.character(habitat))
contrast_SPE_warm     <- summary(pairs(emtrends(m_SPE_warm, ~ period, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Species richness", group = "warm", habitat = "all")
contrast_SPE_warm_hab <- summary(pairs(emtrends(m_SPE_warm_hab, ~ period | habitat, var = "yr"), reverse = TRUE), infer = TRUE) %>% as.data.frame() %>%
  mutate(indicator = "Species richness", group = "warm", habitat = as.character(habitat))

contrasts <- bind_rows(contrast_CTI, contrast_CTI_hab, contrast_ABU, contrast_ABU_hab, contrast_SPE, contrast_SPE_hab,
                       contrast_CTI_cold, contrast_CTI_cold_hab, contrast_ABU_cold, contrast_ABU_cold_hab,
                       contrast_SPE_cold, contrast_SPE_cold_hab,
                       contrast_CTI_warm, contrast_CTI_warm_hab, contrast_ABU_warm, contrast_ABU_warm_hab,
                       contrast_SPE_warm, contrast_SPE_warm_hab)

# CTI trends, all plots and classes (Table 1) -----------------------------------

hab_levels <- c("all", "open", "mixed", "closed")

table1 <- bind_rows(slopes_CTI, slopes_CTI_hab) %>%
  dplyr::select(habitat, period, yr.trend, asymp.LCL, asymp.UCL) %>%
  pivot_wider(names_from = period, values_from = c(yr.trend, asymp.LCL, asymp.UCL)) %>%
  left_join(bind_rows(contrast_CTI, contrast_CTI_hab) %>%
              dplyr::select(habitat, d_beta = estimate, z.ratio, p.value),
            by = "habitat") %>%
  arrange(factor(habitat, levels = hab_levels))
table1

# Period boundary shifted by one year (Table S3) --------------------------------

site_year <- site_year %>%
  mutate(period11 = factor(if_else(year <= 2011, "early", "late"), levels = c("early", "late")),
         period13 = factor(if_else(year <= 2013, "early", "late"), levels = c("early", "late")))

m_CTI_2011 <- glmmTMB(CTI ~ yr * period11 + (1 | id), data = site_year)
summary(emtrends(m_CTI_2011, ~ period11, var = "yr"), infer = TRUE)
summary(pairs(emtrends(m_CTI_2011, ~ period11, var = "yr"), reverse = TRUE), infer = TRUE)

summary(emtrends(m_CTI, ~ period, var = "yr"), infer = TRUE)
summary(pairs(emtrends(m_CTI, ~ period, var = "yr"), reverse = TRUE), infer = TRUE)

m_CTI_2013 <- glmmTMB(CTI ~ yr * period13 + (1 | id), data = site_year)
summary(emtrends(m_CTI_2013, ~ period13, var = "yr"), infer = TRUE)
summary(pairs(emtrends(m_CTI_2013, ~ period13, var = "yr"), reverse = TRUE), infer = TRUE)

# All trends and their change between periods (Table S4) ------------------------

slopes %>%
  dplyr::select(indicator, group, habitat, period, yr.trend, asymp.LCL, asymp.UCL) %>%
  pivot_wider(names_from = period, values_from = c(yr.trend, asymp.LCL, asymp.UCL)) %>%
  left_join(contrasts %>%
              dplyr::select(indicator, group, habitat, d_beta = estimate,
                            d_lwr = asymp.LCL, d_upr = asymp.UCL, z.ratio, p.value),
            by = c("indicator", "group", "habitat")) %>%
  arrange(group, indicator, factor(habitat, levels = hab_levels)) %>%
  print(n = Inf, width = Inf)

# Classes compared within each period, Tukey (Table S5) -------------------------
# CTI, CTI_cold and abundance of cold-associated species

summary(pairs(emtrends(m_CTI_hab, ~ habitat | period, var = "yr")), infer = TRUE)
summary(pairs(emtrends(m_CTI_cold_hab, ~ habitat | period, var = "yr")), infer = TRUE)
summary(pairs(emtrends(m_ABU_cold_hab, ~ habitat | period, var = "yr")), infer = TRUE)

# CTI across all plots (Figure 2B) ----------------------------------------------
# points: annual means from a model with year as a factor; lines: period trends

m_CTI_year <- glmmTMB(CTI ~ fyear + (1 | id), data = site_year, control = ctrl)

i_CTI <- emmeans(m_CTI_year, ~ fyear) %>%
  as.data.frame() %>%
  rename(any_of(new_names)) %>%
  mutate(year = as.integer(as.character(fyear)))

pred_CTI <- tibble(yr = seq(0, 24, by = 0.25)) %>%
  mutate(year   = yr + 2000,
         period = factor(if_else(year <= 2012, "early", "late"), levels = c("early", "late")),
         id     = site_year$id[1])

pr <- predict(m_CTI, newdata = pred_CTI, re.form = NA, se.fit = TRUE)
pred_CTI <- pred_CTI %>%
  mutate(fit = pr$fit,
         lwr = pr$fit - 1.96 * pr$se.fit,
         upr = pr$fit + 1.96 * pr$se.fit)

ggplot() +
  geom_vline(xintercept = 2012.5, linetype = 2, colour = "grey40") +
  geom_ribbon(data = pred_CTI, aes(year, ymin = lwr, ymax = upr, group = period), alpha = 0.15) +
  geom_linerange(data = i_CTI, aes(year, ymin = asymp.LCL, ymax = asymp.UCL), colour = "grey60", linewidth = 0.4) +
  geom_point(data = i_CTI, aes(year, emmean), size = 1.6, shape = 21,
             fill = "white", colour = "black", stroke = 0.7) +
  geom_line(data = pred_CTI, aes(year, fit, group = period), linewidth = 1) +
  labs(x = "Year", y = "Community Temperature Index (°C)") +
  scale_x_continuous(limits = c(2000, 2024), breaks = seq(2000, 2024, 8)) +
  theme_bw()

# CTI trajectories by forest-cover class (Figure 3) -----------------------------

pred_CTI_hab <- expand_grid(yr = seq(0, 24, by = 0.25), habitat = c("open", "mixed", "closed")) %>%
  mutate(year    = yr + 2000,
         habitat = factor(habitat, levels = c("open", "mixed", "closed")),
         period  = factor(if_else(year <= 2012, "early", "late"), levels = c("early", "late")),
         id      = site_year$id[1])

pr <- predict(m_CTI_hab, newdata = pred_CTI_hab, re.form = NA, se.fit = TRUE)
pred_CTI_hab <- pred_CTI_hab %>%
  mutate(fit = pr$fit,
         lwr = pr$fit - 1.96 * pr$se.fit,
         upr = pr$fit + 1.96 * pr$se.fit)

# annual means per class (descriptive)
raw_CTI_hab <- site_year %>%
  group_by(habitat, year) %>%
  summarise(CTI = mean(CTI), .groups = "drop")

ggplot() +
  geom_vline(xintercept = 2012.5, linetype = 2, colour = "grey40") +
  geom_ribbon(data = pred_CTI_hab, aes(year, ymin = lwr, ymax = upr, fill = habitat,
                                       group = interaction(period, habitat)),
              alpha = 0.12, colour = NA) +
  geom_point(data = raw_CTI_hab, aes(year, CTI, colour = habitat), size = 1.1, alpha = 0.45) +
  geom_line(data = pred_CTI_hab, aes(year, fit, colour = habitat,
                                     group = interaction(period, habitat)), linewidth = 1) +
  scale_colour_manual(values = cols_hab, labels = hab_labs, name = NULL) +
  scale_fill_manual(values = cols_hab, guide = "none") +
  labs(x = "Year", y = "Community Temperature Index (°C)") +
  scale_x_continuous(limits = c(2000, 2024), breaks = seq(2000, 2024, 8)) +
  theme_bw()

# Trends of all indices, 3 x 3 panels (Figure 4) --------------------------------
# rows: all species (A-C), cold-associated (D-F), warm-associated (G-I)
# columns: CTI (CTI, CTI_cold, CTI_warm), abundance, species richness

period_labs <- c(early = "2000-2012", late = "2013-2024")
cols_fig    <- c(all = "grey35", cols_hab)
shapes_fig  <- c(all = 18, open = 16, mixed = 16, closed = 16)
labs_fig    <- c(all = "All plots", hab_labs)

fig_data <- slopes %>%
  mutate(row     = factor(group, levels = c("all species", "cold", "warm"),
                          labels = c("All species", "Cold-associated", "Warm-associated")),
         col     = factor(indicator, levels = c("CTI", "Abundance", "Species richness")),
         habitat = factor(habitat, levels = hab_levels),
         period  = factor(period, levels = c("early", "late")))

tags <- expand_grid(row = levels(fig_data$row), col = levels(fig_data$col)) %>%
  mutate(tag = LETTERS[1:9],
         row = factor(row, levels = levels(fig_data$row)),
         col = factor(col, levels = levels(fig_data$col)))

p_col_CTI <- fig_data %>%
  filter(col == "CTI") %>%
  ggplot(aes(period, yr.trend, colour = habitat, shape = habitat)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey55") +
  geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge(width = 0.6),
                  size = 0.45, linewidth = 0.5) +
  geom_text(data = filter(tags, col == "CTI"), aes(x = -Inf, y = Inf, label = tag),
            inherit.aes = FALSE, hjust = -0.4, vjust = 1.4, fontface = 2, size = 3.8) +
  facet_wrap(~ row, ncol = 1) +
  scale_colour_manual(values = cols_fig, labels = labs_fig, name = NULL) +
  scale_shape_manual(values = shapes_fig, labels = labs_fig, name = NULL) +
  scale_x_discrete(labels = period_labs) +
  labs(title = "CTI", x = NULL, y = "Trend (°C/year)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8.5),
        panel.grid.major.x = element_blank())

# row names only in the first column
p_col_ABU <- fig_data %>%
  filter(col == "Abundance") %>%
  ggplot(aes(period, yr.trend, colour = habitat, shape = habitat)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey55") +
  geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge(width = 0.6),
                  size = 0.45, linewidth = 0.5) +
  geom_text(data = filter(tags, col == "Abundance"), aes(x = -Inf, y = Inf, label = tag),
            inherit.aes = FALSE, hjust = -0.4, vjust = 1.4, fontface = 2, size = 3.8) +
  facet_wrap(~ row, ncol = 1) +
  scale_colour_manual(values = cols_fig, labels = labs_fig, name = NULL) +
  scale_shape_manual(values = shapes_fig, labels = labs_fig, name = NULL) +
  scale_x_discrete(labels = period_labs) +
  labs(title = "Abundance", x = NULL, y = "Trend (β, log scale)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major.x = element_blank())

p_col_SPE <- fig_data %>%
  filter(col == "Species richness") %>%
  ggplot(aes(period, yr.trend, colour = habitat, shape = habitat)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey55") +
  geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge(width = 0.6),
                  size = 0.45, linewidth = 0.5) +
  geom_text(data = filter(tags, col == "Species richness"), aes(x = -Inf, y = Inf, label = tag),
            inherit.aes = FALSE, hjust = -0.4, vjust = 1.4, fontface = 2, size = 3.8) +
  facet_wrap(~ row, ncol = 1) +
  scale_colour_manual(values = cols_fig, labels = labs_fig, name = NULL) +
  scale_shape_manual(values = shapes_fig, labels = labs_fig, name = NULL) +
  scale_x_discrete(labels = period_labs) +
  labs(title = "Species richness", x = NULL, y = "Trend (β, log scale)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major.x = element_blank())

(p_col_CTI | p_col_ABU | p_col_SPE) +
  plot_layout(guides = "collect") &
  theme(aspect.ratio = 1,
        text = element_text(size = 9),
        axis.text = element_text(size = 7.5),
        panel.spacing.x = unit(4, "mm"),
        plot.margin = margin(2, 4, 0, 4),
        legend.position = "bottom",
        legend.box.margin = margin(-6, 0, 0, 0),
        legend.margin = margin(0, 0, 0, 0))

# Decomposition of the CTI trend (Figure 5) -------------------------------------
# CTI = p_w * CTI_warm + p_c * CTI_cold, so the CTI trend splits into
#   within cold group = mean p_c * trend in CTI_cold
#   within warm group = mean p_w * trend in CTI_warm
#   between groups    = the rest of the CTI trend
# mean shares are averaged over plot-years, across all plots and per class;
# the between-group term also absorbs the small joint effect of simultaneous
# changes in shares and group indices

cols_thermal <- c(warm = "#b2182b", cold = "#2166ac")

shares <- bind_rows(
  site_year %>% summarise(p_w = mean(p_w)) %>% mutate(habitat = "all"),
  site_year %>% group_by(habitat) %>% summarise(p_w = mean(p_w), .groups = "drop") %>%
    mutate(habitat = as.character(habitat))
) %>%
  mutate(p_c = 1 - p_w)
shares

trend_CTI <- bind_rows(slopes_CTI, slopes_CTI_hab) %>%
  transmute(habitat, period = as.character(period), CTI_trend = yr.trend)
trend_cold <- bind_rows(slopes_CTI_cold, slopes_CTI_cold_hab) %>%
  transmute(habitat, period = as.character(period), cold_trend = yr.trend)
trend_warm <- bind_rows(slopes_CTI_warm, slopes_CTI_warm_hab) %>%
  transmute(habitat, period = as.character(period), warm_trend = yr.trend)

decomp <- trend_CTI %>%
  left_join(trend_cold, by = c("habitat", "period")) %>%
  left_join(trend_warm, by = c("habitat", "period")) %>%
  left_join(shares, by = "habitat") %>%
  mutate(within_cold = p_c * cold_trend,
         within_warm = p_w * warm_trend,
         between     = CTI_trend - within_cold - within_warm)

decomp %>%
  dplyr::select(habitat, period, CTI_trend, within_cold, within_warm, between) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

decomp_long <- decomp %>%
  pivot_longer(c(within_cold, within_warm, between), names_to = "component", values_to = "value") %>%
  mutate(habitat   = factor(habitat, levels = hab_levels, labels = labs_fig[hab_levels]),
         period    = factor(period, levels = c("early", "late"), labels = period_labs),
         component = factor(component, levels = c("between", "within_warm", "within_cold"),
                            labels = c("Between groups", "Within warm group", "Within cold group")))

decomp_total <- decomp %>%
  mutate(habitat = factor(habitat, levels = hab_levels, labels = labs_fig[hab_levels]),
         period  = factor(period, levels = c("early", "late"), labels = period_labs))

ggplot(decomp_long, aes(habitat, value, fill = component)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_col(width = 0.65) +
  geom_errorbar(data = decomp_total, aes(habitat, ymin = CTI_trend, ymax = CTI_trend),
                inherit.aes = FALSE, width = 0.8, linewidth = 0.9) +
  facet_wrap(~ period) +
  scale_fill_manual(values = c("Between groups"    = "grey65",
                               "Within warm group" = cols_thermal[["warm"]],
                               "Within cold group" = cols_thermal[["cold"]]),
                    name = NULL) +
  labs(x = NULL, y = "Contribution to CTI trend (°C/year)") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank())
