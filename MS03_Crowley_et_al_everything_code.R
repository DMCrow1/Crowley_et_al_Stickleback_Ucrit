# ============================================================
# MASTER ANALYSIS SCRIPT
# ============================================================
# ==================================================================================================================
#----------------- PART 1  -----------------#
#----------------- Full burst analysis -----------------#
# ==================================================================================================================
library(car)
library(DHARMa)
library(emmeans)
require(gamlss)
library(ggplot2)
library(glmmTMB)
library(here)
library(lme4)
library(lmerTest)
library(readxl)
library(sjPlot)
library(splines)
library(tidyverse)
library(mgcViz)

setwd("/path/to/your/data")
# Load the data set
fish.data <- read_excel("Crowley_et_all_Lab_data.xlsx", sheet = "Burst_data")

# Rename 't' to avoid conflict with base R's transpose
fish.data <- fish.data %>% rename(time_sec = t)%>%
  mutate(across(c(body_length, time_sec, x, water_speed),
                ~ as.numeric(gsub(",", ".", .)))) %>%
  mutate(b_g = tolower(b_g))

## ----------------------------------------------------------------------------------
## Remove burst_IDs that cross multiple water speeds and with only 2 b's & no g"s
## ----------------------------------------------------------------------------------

# Step 1: Bursts that cross multiple water speeds
bad_bursts_speeds <- fish.data %>%
  group_by(burst_ID) %>%
  summarise(
    n_speeds = n_distinct(water_speed),
    .groups = "drop"
  ) %>%
  filter(n_speeds > 1) %>%
  pull(burst_ID)

# Step 2: Bursts with 2 or fewer 'b' and 0 'g'
bad_bursts_b_only <- fish.data %>%
  group_by(burst_ID) %>%
  summarise(
    n_b = sum(b_g == "b"),
    n_g = sum(b_g == "g"),
    .groups = "drop"
  ) %>%
  filter(n_b <= 2 & n_g == 0) %>%
  pull(burst_ID)

# Step 3: Combine all bad bursts
all_bad_bursts <- union(bad_bursts_speeds, bad_bursts_b_only)

# Step 4: Remove bad bursts
df_clean <- fish.data %>%
  filter(!(burst_ID %in% all_bad_bursts)) %>%
  group_by(fish_ID, burst_ID) %>%
  arrange(time_sec, .by_group = TRUE) %>%
  mutate(
    dx    = x - lag(x),
    dt    = time_sec - lag(time_sec),
    dx_BL = dx / body_length
  ) %>%
  ungroup()

# Optional: check removed bursts
removed_bursts <- fish.data %>%
  filter(burst_ID %in% all_bad_bursts) %>%
  distinct(fish_ID, burst_ID)

print(paste("Number of removed bursts:", nrow(removed_bursts)))

# ----------------------------------------------------------
# Burst-level summary table
# ----------------------------------------------------------
DF <- df_clean %>%
  group_by(fish_ID, exp_condition, housing, water_speed, burst_ID, tunnel) %>%
  summarise(
    burst_distance = sum(abs(dx[b_g == "b"]), na.rm = TRUE),
    max_vx = max(abs(vx[b_g == "b"]), na.rm = TRUE),
    max_ax = max(abs(ax[b_g == "b"]), na.rm = TRUE),
    avg_vx = mean(abs(vx[b_g == "b"]), na.rm = TRUE),
    avg_ax = mean(abs(ax[b_g == "b"]), na.rm = TRUE),
    avg_burst_distance = mean(abs(dx[b_g == "b"]), na.rm = TRUE),
    .groups = "drop"
  )

# ==================================================================================================================
#----------------- Models  -----------------#
#------------------------------------------------------------------------------#

DF$housing       <- factor(DF$housing)
DF$exp_condition <- factor(DF$exp_condition)
DF$fish_ID       <- factor(DF$fish_ID)
DF$tunnel        <- factor(DF$tunnel)
DF$BurstNr       <- as.numeric(
  unlist(lapply(strsplit(DF$burst_ID, "_"), function(x){x[[3]]}))
)

# ----------------------------------------------------------
# GAM 1: avg_burst_distance
# ----------------------------------------------------------
GAM <- gam(
  log2(avg_burst_distance) ~ exp_condition *housing  +
    tunnel +
    s(water_speed, k = 3) +
    s(BurstNr, k = 3) +
    s(fish_ID, bs = "re"),
  data = DF,
  family = scat(),
  method = "REML"
)
plot(GAM, residuals = TRUE)
summary(GAM)
gam.check(GAM)
qq.gam(GAM)

# More important
simres <- simulateResiduals(GAM, n = 1e3)
plot(simres, quantreg = TRUE)

# Formal comparison
EMGAM1 <- emmeans(GAM,pairwise ~ exp_condition | housing,
                  type = "response")
EMGAM1$contrasts
EM2GAM1 <- emmeans(GAM,pairwise ~ housing | exp_condition,
                   type = "response")
EM2GAM1$contrasts

leveneTest(log2(avg_burst_distance) ~ housing,data = DF)

# Kolmogorov-Smirnov test by condition
ks.test(DF$avg_ax[DF$housing == "constant"],
        DF$avg_ax[DF$housing == "standard"])

# This specifically tests if the effect of condition changes by house
contrast(emmeans(GAM, ~ exp_condition * housing), interaction = "pairwise",type = "response")

# ----------------------------------------------------------
# GAM 2: avg_ax
# ----------------------------------------------------------
GAM2 <- gam(
  log2(avg_ax) ~ exp_condition * housing +
    tunnel +
    s(water_speed, k = 3)+
    s(fish_ID, bs = "re"),
  data = DF,
  family = scat()
)
plot(GAM2, residuals = TRUE)
summary(GAM2)
gam.check(GAM2)
qq.gam(GAM2)

# More important
simres <- simulateResiduals(GAM2, n = 1e3)
plot(simres, quantreg = TRUE)

# Formal comparison
EMGAM2 <- emmeans(GAM2,pairwise ~ exp_condition | housing,
                  type = "response")
EMGAM2$contrasts
EM2GAM2 <- emmeans(GAM2,pairwise ~ housing | exp_condition,
                   type = "response")
EM2GAM2$contrasts

leveneTest(log2(avg_ax) ~ housing,data = DF)

# Kolmogorov-Smirnov test by condition
ks.test(DF$avg_ax[DF$housing == "constant"],
        DF$avg_ax[DF$housing == "standard"])

# This specifically tests if the effect of condition changes by house
contrast(emmeans(GAM2, ~ exp_condition * housing), interaction = "pairwise",type = "response")

# ----------------------------------------------------------
# GAM 3: avg_vx
# ----------------------------------------------------------
GAM3 <- gam(
  log(avg_vx) ~ exp_condition * housing +
    tunnel +
    s(fish_ID, bs = "re"),
  data = DF,
  family = scat()
)
plot(GAM3, residuals = TRUE)
summary(GAM3)
gam.check(GAM3)
qq.gam(GAM3)

# More important
simres <- simulateResiduals(GAM3, n = 1e3)
plot(simres, quantreg = TRUE)

# Formal comparison
EMGAM3 <- emmeans(GAM3,pairwise ~ exp_condition | housing,
                  type = "response")
EMGAM3$contrasts
EM2GAM3 <- emmeans(GAM3,pairwise ~ housing | exp_condition,
                   type = "response")
EM2GAM3$contrasts

leveneTest(log2(avg_vx) ~ housing,data = DF)

# Kolmogorov-Smirnov test by condition
ks.test(DF$avg_vx[DF$housing == "constant"],
        DF$avg_vx[DF$housing == "standard"])

# This specifically tests if the effect of condition changes by house
contrast(emmeans(GAM3, ~ exp_condition * housing), interaction = "pairwise",type = "response")

emmeans(GAM3, ~ exp_condition * housing,
        type = "response")

# Check means within housing 
DF %>%
  group_by(exp_condition, housing) %>%
  summarise(
    n = n(),
    mean_log = mean(log(avg_vx), na.rm = TRUE),
    sd_log = sd(log(avg_vx), na.rm = TRUE)
  )

ggplot(DF, aes(x = exp_condition, y = avg_vx, colour = housing)) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA)

# ----------------------------------------------------------
# Total number of bursts per fish per condition
# ----------------------------------------------------------
total_bursts <- df_clean %>%
  distinct(fish_ID, housing, tunnel, exp_condition, burst_ID) %>%
  count(fish_ID, housing, exp_condition)

total_bursts <- df_clean %>%
  group_by(fish_ID, housing,tunnel, exp_condition) %>%
  summarize(mean_n = mean(burst_count, na.rm = TRUE), .groups = "drop")

total_bursts <- total_bursts %>%
  arrange(fish_ID, tunnel, exp_condition) %>%
  group_by(fish_ID, housing) %>%
  mutate(
    direction = ifelse(mean_n[2] > mean_n[1], "up", "down")
  ) %>%
  ungroup()

total_bursts$exp_condition <- as.factor(total_bursts$exp_condition)
LMM1 <- lmer (mean_n ~ exp_condition * housing + tunnel+ (1 | fish_ID), data = total_bursts)
# Residuals vs fitted
plot(LMM)
simres <- simulateResiduals(LMM)
plot(simres, quantreg = FALSE)
summary(LMM1)
levels(total_bursts$exp_condition)

# ----------------------------------------------------------
# Frequency of bursts per minute
# ----------------------------------------------------------

# Step 1: Calculate duration per burst and glide
seq_summary <- df_clean %>%
  group_by(fish_ID, housing, exp_condition, burst_count,
           burst_ID, b_g, water_speed,tunnel) %>%
  summarise(
    duration = max(time_sec) - min(time_sec),
    .groups = "drop"
  )

# Step 2: Summarise number of bursts per fish × condition × speed
burst_counts <- seq_summary %>%
  filter(b_g == "b") %>%  # only bursts
  group_by(fish_ID, housing, exp_condition, water_speed,tunnel) %>%
  summarise(
    n_bursts = n(),
    duration_min = 5,  # each observation window is 5 minutes
    rate_bursts_per_min = n_bursts / duration_min,
    .groups = "drop"
  )

# Quick check of burst count distribution
# Summary statistics
summary(burst_counts$n_bursts)
cat("Mean:", mean(burst_counts$n_bursts), "Variance:", var(burst_counts$n_bursts), "\n")

# Histogram
hist(burst_counts$n_bursts, breaks = 20, 
     main = "Distribution of burst counts", 
     xlab = "Number of bursts", col = "lightblue")

# Step 3:  Group means & variances using more granular groups
burst_counts$group <- interaction(burst_counts$exp_condition,
                                  burst_counts$housing,
                                  burst_counts$water_speed,
                                  burst_counts$tunnel,
                                  drop = TRUE)

grpMeans <- tapply(burst_counts$n_bursts, burst_counts$group, mean)
grpVars  <- tapply(burst_counts$n_bursts, burst_counts$group, var)

# Remove groups with only 1 observation (variance = NA)
keep <- !(is.na(grpVars) | is.na(grpMeans))
grpMeans <- grpMeans[keep]
grpVars  <- grpVars[keep]

# Fit mean–variance models
lm_qp <- lm(grpVars ~ grpMeans - 1)
phi_fit <- coef(lm_qp)

lm_nb <- lm(grpVars ~ I(grpMeans^2) + offset(grpMeans) - 1)
k_fit <- 1 / coef(lm_nb)

# Loess with automatic span (safe even for small samples)
loess_fit <- loess(grpVars ~ grpMeans,
                   span = 0.75,    # prevent small-span errors
                   degree = 1)

# Plot variance per water speed 
plot(grpVars ~ grpMeans,
     xlab = "Group means",
     ylab = "Group variances",
     main = "Mean–variance diagnostics for burst counts")

# Poisson line
abline(a = 0, b = 1, lty = 2)
# Quasi-Poisson
curve(phi_fit * x, col = 2, add = TRUE)
# Negative binomial
curve(x * (1 + x/k_fit), col = 4, add = TRUE)
# Loess smoother (safe)
mvec <- seq(min(grpMeans), max(grpMeans), length.out = 200)
lines(mvec, predict(loess_fit, data.frame(grpMeans = mvec)),
      col = 5, lwd = 2)
ggplot(burst_counts, aes(x = water_speed, y = n_bursts)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  facet_wrap(~ exp_condition + housing)

# Step 4: Fit negative binomial GLMM on burst counts due to variance outcome
model_tunnel <- glmmTMB(
  n_bursts ~ exp_condition * housing +
    poly(water_speed, 2) +
    tunnel +
    offset(log(duration_min)) +
    (1 | fish_ID),
  family = nbinom2,
  data = burst_counts
)
summary(model_tunnel)

# Step 5: Post-hoc tests for condition or home differences
emmeans(model_tunnel, pairwise ~ exp_condition * housing, type = "response",adjust = "tukey")

# ----------------------------------------------------------
# Compute relative x-location (relative to average location)
# ----------------------------------------------------------
df_clean$exp_condition <- factor(df_clean$exp_condition)
ref_x <- tapply(df_clean$x, list(df_clean$fish_ID, df_clean$exp_condition, df_clean$burst_count), mean)

df_clean$rel_x <- numeric(nrow(df_clean))
for(i in 1:nrow(ref_x)){
  cat("fish", i, "\n")
  for(j in 1:2){
    current <-  df_clean[df_clean$fish_ID == rownames(ref_x)[i] & 
                           df_clean$exp_condition == levels(df_clean$exp_condition)[j], ]
    n_burst <- length(unique(current$burst_count))
    current$burst_count <- factor(current$burst_count)
    for(k in 1:n_burst){
      wh <- df_clean$fish_ID == rownames(ref_x)[i] & 
        df_clean$exp_condition == levels(df_clean$exp_condition)[j] &
        df_clean$burst_count == levels(current$burst_count)[k]
      df_clean$rel_x[wh] <- df_clean$x[wh] - ref_x[i, j, k]
    }
  }
}

for(i in 1:nrow(df_clean)){
  if(i / 1000 == round(i / 1000)){
    cat(round(i/nrow(DF) * 100, 1), "%\n")
  }
  xi <- df_clean$fish_ID[i]
  xj <- df_clean$exp_condition[i]
  xk <- df_clean$burst_count[i]
  df_clean$rel_x[i] <- df_clean$x[i] - ref_x[xi, xj, xk]
}

# Check that it worked
summary(tapply(df_clean$rel_x, list(df_clean$fish_ID, df_clean$exp_condition, df_clean$burst_count), mean))

# Spread in x-location per fish
A <- tapply(df_clean$rel_x, list(df_clean$fish_ID, df_clean$exp_condition, df_clean$burst_count), sd)
DF_2 <- as.data.frame.table(A)
DF_2 <- DF_2[!is.na(DF_2$Freq), ]
names(DF_2) <- c("ID", "Exp", "Burst", "SDx")
DF_2$House <- factor(substr(DF_2$ID, 2, 2))

DF_2$Tunnel <- character(nrow(DF_2))
DF_2$WaterSpeed <- numeric(nrow(DF_2))
for(i in 1:nrow(DF_2)){
  if(i / 100 == round(i / 100)){
    cat(round(i/nrow(DF_2) * 100, 1), "%\n")
  }
  wh <- which(df_clean$fish_ID == DF_2$ID[i] &
                df_clean$exp_condition == DF_2$Exp[i] &
                df_clean$burst_count == DF_2$Burst[i])[1]
  DF_2$Tunnel[i]     <- df_clean$tunnel[wh]
  DF_2$WaterSpeed[i] <- df_clean$water_speed[wh]
}
DF_2$Tunnel <- factor(DF_2$Tunnel)

LMM <- lmer(log2(SDx) ~ Exp * House + ns(WaterSpeed, 2) + Tunnel + (1 | ID), data = DF_2)

# Residuals vs fitted
plot(LMM)
simres <- simulateResiduals(LMM)
plot(simres, quantreg = FALSE)

plot_model(LMM, type = "pred", terms = c("WaterSpeed", "Exp", "House"))

summary(LMM)

# Formal comparison
EMM1 <- emmeans(LMM, pairwise ~ Exp | House,
                type = "response")
EMM1$contrasts
EMM2 <- emmeans(LMM, pairwise ~ House | Exp,
                type = "response")
EMM2$contrasts
# ==================================================================================================================
#----------------- Graphs  -----------------#
#------------------------------------------------------------------------------#
# ----------------------------------------------------------
# Avg_burst_distance, Avg_Ax & Avg_vx
# ----------------------------------------------------------
p <- plot_model(GAM3, type = "pred", terms = c("housing", "exp_condition"))

# Create background rectangles as a separate ggplot layer
background <- list(
  annotate(
    "rect",
    xmin = 0.5, xmax = 1.5,
    ymin = -Inf, ymax = Inf,
    fill = "#56B4E9",
    alpha = 0.3
  ),
  annotate(
    "rect",
    xmin = 1.5, xmax = 2.5,
    ymin = -Inf, ymax = Inf,
    fill = "#E69F00",
    alpha = 0.3
  )
)
# Put background layers before the existing plot layers
p_final <- p
p_final$layers <- c(background, p$layers)
p_final <- p_final +
  # Colour palette
  scale_color_manual(values = c("darkgrey", "lightyellow")) +
  scale_fill_manual(values = c("darkgrey", "lightyellow")) +
  # Theme customization
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    # Reduce spacing between facets
    panel.spacing = unit(0, "cm"),
    # Move x-axis title closer to axis
    axis.title.x = element_text(
      margin = margin(t = 0.2, unit = "cm")
    ),
    # Reduce bottom plot margin
    plot.margin = unit(c(1, 1, 0.5, 1), "cm")
  ) +
  # Expand y-scale slightly at bottom
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    x = "Housing condition",
    y = "Average burst distance (cm)",
    colour = "Experimental light condition",
    fill = "Experimental light condition"
  )
p_final
ggsave(
  filename = "avg_burst.png", plot = p_final, width = 5, height = 5, dpi = 300)

# ----------------------------------------------------------
# Total number of bursts per fish per condition
# ----------------------------------------------------------
p3 <- ggplot(total_bursts, aes(x = exp_condition, 
                               y = mean_n, 
                               fill = housing)) +
  # Background rectangles
  geom_rect(data = NULL, aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.1) +
  geom_rect(data = NULL, aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf),
            fill = "lightyellow", alpha = 0.1) +
  
  # Boxplots
  geom_boxplot(alpha = 0.8) +
  
  # Lines colored by direction
  geom_line(aes(group = fish_ID, color = direction), linewidth = 0.6) +
  
  # Points
  geom_point(aes(color = ), position = position_dodge(width = .5), size = 2) +
  
  # Facets
  facet_wrap(~ housing) +
  
  # Manual color scale for lines/points
  scale_color_manual(values = c("down" = "black", "up" = "red")) +
  
  # Manual fill scale for boxplots
  scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
  
  # Labels
  labs(
    x = "Experimental light condition",
    y = "Mean number of bursts",
  ) +
  
  # Theme
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
p3
ggsave(
  filename = "avg_number_burst.png", plot = p3, width = 8, height = 5, dpi = 300)

# ----------------------------------------------------------
# Burst frequency
# ----------------------------------------------------------
# Summarise mean rate per fish × condition × housing
burst_summary <- burst_counts %>%
  group_by(fish_ID, housing, exp_condition) %>%
  summarise(
    mean_rate = mean(rate_bursts_per_min, na.rm = TRUE),
    .groups = "drop"
  )
burst_summary <- burst_summary %>%
  mutate(mean_rate = as.numeric(mean_rate))

ggplot(burst_summary, aes(x = exp_condition, 
                          y = mean_rate, 
                          fill = housing)) +
  geom_rect(data=NULL,aes(xmin=0.5,xmax=1.5,ymin=-Inf,ymax=Inf),
            fill="lightgrey", alpha=0.1) +
  geom_rect(data=NULL,aes(xmin=1.5,xmax=2.5,ymin=-Inf,ymax=Inf),
            fill="lightyellow", alpha=0.1) +
  geom_boxplot(alpha=.8) +
  
  # # Lines connecting each fish (paired data)
  # # Lines connecting each fish (all black)
  # geom_line(aes(group = fish_ID), color = "black", linewidth = 0.8)+
  # Points for each fish
  geom_point(size = 2) +
  facet_wrap(~ housing) +
  scale_fill_manual(values = c("#56B4E9","#E69F00")) +
  
  # Labels and theme
  labs(
    x = "Experimental light condition",
    y = "Mean burst rate (bursts/min)",
    title = "Paired mean burst rate across light conditions by housing"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
ggsave("burstfreq_plot.png", width = 8, height = 6, units = "in", dpi = 300)

# ==================================================================================================================
#----------------- PART 2  -----------------#
#----------------- Ucrit and Ugt analysis and visulisation  -----------------#
# ==================================================================================================================
require("car")
require("DHARMa")
require("effects")
require("lmerTest")
require("readxl")
require("sjPlot")
require(emmeans)
require(ggplot2)

setwd("/path/to/your/data")
DF2 <- read_excel("Crowley_et_all_Lab_data.xlsx", sheet = "Capacity_data")
# ==================================================================================================================
#----------------- Models Ucrit -----------------#
#------------------------------------------------------------------------------#
## ----------------------------------------------------------------------------------
# Start condition: did we fix the start condition from the field?
## ----------------------------------------------------------------------------------
LMM <- lmer(
  Ugt ~ Condition * Home_light + Start_Condition + (1 | Fish_ID),
  data = DF2
)

# Some diagnostic plots
plot(LMM)
qqPlot(resid(LMM))
simres <- simulateResiduals(LMM)
plot(simres)

summary(LMM)

# YES! true for both Ucrit and Ugt

## ----------------------------------------------------------------------------------
# Including variables of interest 
# Swim tunnel: Does swim tunnel matter?
## ----------------------------------------------------------------------------------
as.numeric(DF2$Swim_tunnel)
T1 <- subset(DF2, Swim_tunnel == "1")
T7 <- subset(DF2, Swim_tunnel == "7")
t.test(T1$Ucrit,T7$Ucrit)

# YES = must be in all models Ucrit + Ugt

LMM2 <- lmer(
  Ucrit ~ Condition * Home_light + Swim_tunnel +(1 | Fish_ID),
  data = DF2
)
summary(LMM2)

EM1 <- emmeans(LMM2, pairwise ~ Condition | Home_light,type = "response")
EM1$contrasts

EM2 <- emmeans(LMM2,pairwise ~  Home_light | Condition,
               type = "response")
EM2$contrasts
sjPlot::plot_model(LMM2, type = "int")

#----------------- Model Ugt -----------------#
#------------------------------------------------------------------------------#
DF2$Condition <- as.factor(DF2$Condition)
LMM3 <- lmer(
  Ugt ~ Condition * Home_light + Swim_tunnel + (1 | Fish_ID),
  data = DF2
)
summary(LMM3)
levels(DF2$Condition)
# Some diagnostic plots
plot(LMM3)
qqPlot(resid(LMM3))
simres <- simulateResiduals(LMM4)
plot(simres)

# ==================================================================================================================
#----------------- Graphs -----------------#
#------------------------------------------------------------------------------#
## ----------------------------------------------------------------------------------
# Ucrit
## ----------------------------------------------------------------------------------
## Paired data graph with box plot ##
Ucrit <- ggplot(DF2, aes(x = Condition, y = Ucrit, fill = Home_light)) +
  ylab("Ucrit (cm s-1)") +
  scale_fill_manual(values = c("#56B4E9","#E69F00"))+
  geom_rect(data=NULL,aes(xmin=0.5,xmax=1.5,ymin=-Inf,ymax=Inf),
            fill="lightgrey", alpha=0.1) +
  geom_rect(data=NULL,aes(xmin=1.5,xmax=2.5,ymin=-Inf,ymax=Inf),
            fill="lightyellow", alpha=0.1) +
  geom_boxplot(alpha=.8) +
  geom_line(aes(group = Fish_ID,color = Ucrit_trend)) +
  scale_color_manual(values=c("#FF3030", "#333333","blue"))+
  geom_point(size = 2) +
  facet_wrap(~ Home_light) +
  scale_x_discrete("Experimental light condtion") +
  theme_minimal() +
  theme(legend.position = "", 
        panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

Ucrit
ggsave("Ucrit_new.png", plot = Ucrit, device = "png", dpi = 300, width = 10, height = 8, units = "in")

## ----------------------------------------------------------------------------------
# Ugt
## ----------------------------------------------------------------------------------
## Paired data graph with box plot ##
Ugt <- ggplot(DF2, aes(x = Condition, y = Ugt, fill = Home_light)) +
  ylab("Ugt (cm s-1)") +
  scale_fill_manual(values = c("#56B4E9","#E69F00"))+
  geom_rect(data=NULL,aes(xmin=0.5,xmax=1.5,ymin=-Inf,ymax=Inf),
            fill="lightgrey", alpha=0.1) +
  geom_rect(data=NULL,aes(xmin=1.5,xmax=2.5,ymin=-Inf,ymax=Inf),
            fill="lightyellow", alpha=0.1) +
  geom_boxplot(alpha=.8) +
  geom_line(aes(group = Fish_ID,color = Ugt_trend)) +
  scale_color_manual(values=c("#FF3030", "#333333","blue"))+
  geom_point(size = 2) +
  facet_wrap(~ Home_light) +
  scale_x_discrete("Experimental light condtion") +
  theme_minimal() +
  theme(legend.position = "", 
        panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
Ugt
ggsave("Ugt.png", plot = Ugt, device = "png", dpi = 300, width = 10, height = 8, units = "in")
