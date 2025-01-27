#############################################################################
# SCRIPT: Ontogenetic mechanisms of differentiation in two Salvia species 
# with different pollinators. 
# Davies, Agustín; Benitez-Vieyra, Santiago
############################################################################

library(car)
library(ggeffects)
library(ggplot2)
library(lme4)
library(lmerTest)
library(patchwork)
library(performance)

###############################
### ALLOMETRIC TRAJECTORIES ###
###############################

# data input and preparation
dat1 <- read.csv("trayectorias_ont.csv", header = T, dec = ".", sep = ",")
dat1$ind <- paste(dat1$species, dat1$ID, sep ="_")
dat1$ind <- as.factor(dat1$ind)
dat1$adu <- factor(as.factor(dat1$bud), labels = c("p", "p", "p", "p",
                                                    "p", "p", "p", "p", 
                                                    "p", "a"))
dat1$species <- as.factor(dat1$species)

# variables are centered to assure convergence
dat1$slogULL <- scale(log(dat1$ULL), scale = F)
dat1$slogCTL <- scale(log(dat1$CTL), scale = F)

# data exploration
ggplot(dat1) + geom_point(aes(x = log(ULL), y = log(CTL), color = ind))

# model selection - random component
pit1 <- lmer(slogCTL ~ slogULL*species + (1|ind), data = dat1)
pit2 <- lmer(slogCTL ~ slogULL*species + (1|ind) + (0 + slogULL|ind), 
             data = dat1)

anova(pit1, pit2)    
compare_performance(pit1, pit2)         

# model selection - fixed component
pit3 <- lmer(slogCTL ~ slogULL + species + (1|ind) + (0 + slogULL|ind), 
             data = dat1)
pit4 <- lmer(slogCTL ~ slogULL + (1|ind) + (0 + slogULL|ind), 
             data = dat1)
anova(pit2, pit3) 
anova(pit3, pit4) 
compare_performance(pit2, pit3, pit4)

summary(pit4)

# model diagnosis
check_model(pit2)

## refit model 3 and 4 with uncentered variables
dat1$logULL <- log(dat1$ULL)
dat1$logCTL <- log(dat1$CTL)
fm1 <- lmer(logCTL ~ logULL+species + (1|ind) + (0 + logULL|ind), 
            data = dat1)
fm2 <- lmer(logCTL ~ logULL + (1|ind) + (0 + logULL|ind), 
                     data = dat1)
anova(fm1, fm2)
compare_performance(fm1, fm2)

# Get predictions
# (change type to "random" to include variation due to random effects,
# in such case, they will be 95% prediction intervals)
pred <- ggpredict(fm1, terms = c("logULL", "species"), type = "fixed")

# Plot using ggplot2
t1 <- ggplot(pred, aes(x = x, y = predicted, color = group)) +
  geom_line(alpha = 1, linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.2, linewidth = 0.3) +
  geom_point(data = dat1, aes(x = logULL, y = logCTL, color = species, 
                              shape = adu), alpha = 0.5) +
  scale_color_manual(values = c("#1C1089", "#EC990C")) +
  scale_shape_manual(values = c(19, 1)) +
  xlab("ln(corolla upper lip length)") + 
  ylab("ln(corolla tube length)") + 
  theme_bw() +
  theme(legend.position = "none", axis.title=element_text(size=8))
t1

# further details were edited in inkscape
svg(file = "FIG3.svg", width = 3.15, height = 3)
t1
dev.off()


######################################
### DURATION OF FLOWER DEVELOPMENT ###
######################################

# check relation between calix and corolla size
cc <- read.csv("calyx_corolla.csv", header = TRUE)
cc$log_calyx <- log(cc$calyx_length)
cc$log_tube <- log(cc$tube_length)

ggplot(cc, aes(x = log_calyx, y = log_tube, color = sp)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("#1C1089", "#EC990C")) +
  theme_minimal() +
  labs(x = "log(calix length)",
       y = "log(corolla tube length)")

# model to predict tube length from calyx length
fit <- lm(log_tube ~ sp* log_calyx, data = cc)
summary(fit) ## R2  = 0.9212
anova(fit)

######################################
# relation between length and time to anthesis
dat <- read.csv("growth.csv", header = TRUE, stringsAsFactors = TRUE)
dat$days <- dat$hours_pre/24
dat$IDS <- paste(dat$sp, dat$ID, sep = "_")
summary(dat)

# data exploration
ggplot(dat, aes(x = days, y = long_calyx, 
color = IDS, shape = sp)) +
  geom_point() + geom_line() +
  theme_minimal() +
  labs(x = "days pre antesis", y = "calyx mm") 

# predict tube length from calix length
# (only in closed buds)
newdat <- data.frame(
  sp = dat$sp, 
  log_calyx = log(dat$long_calyx)
)

dat$tube <- ifelse(dat$cat == "f",
  dat$long_tube,
  exp(predict(fit, newdata = newdat))
)

# data exploration
ggplot(dat, aes(x = days, y = tube, 
color = IDS, shape = sp)) +
  geom_point() + geom_line() +
  theme_minimal() +
  labs(x = "days pre anthesis", y = "tube mm") 

# exclude data from open flowers
buds <- subset(dat, dat$cat != "f")

ggplot(buds, aes(x = days, y = tube, 
color = IDS, shape = sp)) +
  geom_point() + geom_line() +
  theme_minimal() +
  labs(x = "days pre anthesis", y = "tube mm") 

# model
m1 <- lmer(days ~ tube + (1|IDS), data = buds)
m2 <- lmer(days ~ sp + tube + (1|IDS), data = buds)
m3 <- lmer(days ~ sp * tube + (1|IDS), data = buds)
compare_performance(m1, m2, m3)
anova(m2, m3)

summary(m3)

new_data2 <- expand.grid(
  tube = seq(0,
    max(buds$tube, na.rm = TRUE),
    length.out = 100
  ),
  sp = levels(buds$sp)
)

# Function to generate predictions
predict_fun <- function(model) {
  predict(model, newdata = new_data2, re.form = NA) # Fixed effects only
}

# Bootstrap with bootMer
set.seed(123) # For reproducibility
boot_res <- bootMer(m3, predict_fun, nsim = 1000, type = "parametric")

# Extract 95% CI
new_data2$pred <- predict(m3, new_data2, re.form = NA) # Predicted values
new_data2$lower <- apply(boot_res$t, 2, quantile, probs = 0.025) # Lower 95% CI
new_data2$upper <- apply(boot_res$t, 2, quantile, probs = 0.975) # Upper 95% CI

# inverse plot (tube in y)
tube <- ggplot(buds, aes(x = days, y = tube, col = sp)) +
  geom_line(data = new_data2, aes(x = pred, y = tube), linewidth = 0.5) +
  geom_line(data = new_data2, aes(x = upper, y = tube), size = 0.5, linetype = "solid") +
  geom_line(data = new_data2, aes(x = lower, y = tube), size = 0.5, linetype = "solid") +
  geom_point(aes(group = IDS), size = 0.5, alpha = 0.5) + # Second plot's points
  geom_line(aes(group = IDS), alpha = 0.5, linetype = "solid", linewidth = 0.25) + # Second plot's lines
  scale_color_manual(values = c("#1C1089", "#EC990C")) +
  xlim(-16, 0) +
  labs(x = "days pre anthesis", y = "tube length (mm)") +
  guides(color = "none", shape = "none") +
  theme_bw()
tube

svg(file = "FIG4.svg", width = 3.15, height = 3)
tube # further details were edited in inkscape
dev.off()

## predict days at tube = 0
### valores de la predicción
stach_0 <- predict(m3, newdata = data.frame(tube = 0, sp = "S.stachydifolia"), 
re.form = NA, se.fit = TRUE)
stach_0$fit # -12.4685 days
stach_0$se # 0.3401 se

guar_0 <- predict(m3, newdata = data.frame(tube = 0, sp = "S.guaranitica"), 
re.form = NA, se.fit = TRUE)
guar_0$fit # -14.5707 days
guar_0$se # 0.3613


####################################################
## transforming measures of tube length in stages ##
####################################################

inv_m3 <- lmer(tube ~ sp * days + (1|IDS), data = buds)
summary(inv_m3)

calc_intervalos <- function(dia){
diax <- bootMer(inv_m3, FUN = function(x) {
  predict(x,
    newdata = expand.grid(
      days = dia,
      sp = c("S.guaranitica", "S.stachydifolia")
    ),
    re.form = ~0
  )
}, nsim = 1000)

ci_lower <- apply(diax$t, 2, quantile, probs = 0.025)
ci_upper <- apply(diax$t, 2, quantile, probs = 0.975)
return(c(ci_lower[1], ci_upper[1], ci_lower[2], ci_upper[2]))
}

set.seed(2468) # for reproducibility
intervalos <- t(sapply(-1:-10, calc_intervalos))
colnames(intervalos) <- c("guar_lower", "guar_upper", "stach_lower", "stach_upper")
int <- intervalos*1000

##########################
### Cell proliferation ###
##########################

# data input and preparation
dat2 <- read.csv("celulas_num.csv", header = TRUE)
dat2 <- na.omit(dat2)
dat2$logTL <- log(dat2$TL)
dat2$ID <- as.factor(paste(dat2$species, dat2$ID, sep = "_"))
dat2$species <- as.factor(dat2$species)
dat2$adul <- factor(as.factor(dat2$bud), labels = c("p", "p", "p", "p", "a"))

aggregate(dat2$CN, by = list(dat2$species, dat2$adul), FUN = mean)
aggregate(dat2$CN, by = list(dat2$species, dat2$adul), FUN = sd)

### create categories based on confidence intervals
dat2$categoria <- ifelse(dat2$species == "S.guaranitica" &
  dat2$TL <= int[1,2] & dat2$TL >= int[4,2], "d1-3",
ifelse(dat2$species == "S.guaranitica" &
  dat2$TL <= int[4,2] & dat2$TL >= int[7,2], "d4-6",
ifelse(dat2$species == "S.guaranitica" &
  dat2$TL <= int[7,2] & dat2$TL >= int[10,2], "d7-9",
ifelse(dat2$species == "S.guaranitica" &
  dat2$TL <= int[10,2], "d10plus",
ifelse(dat2$species == "S.guaranitica" &
  dat2$adul == "p" &
  dat2$TL >= int[1,2], "d0-1",
ifelse(dat2$adul == "a", "flower",
  ifelse(dat2$species == "S.stachydifolia" &
    dat2$TL <= int[1,4] & dat2$TL >= int[4,4], "d1-3",
  ifelse(dat2$species == "S.stachydifolia" &
    dat2$TL <= int[4,4] & dat2$TL >= int[7,4], "d4-6",
  ifelse(dat2$species== "S.stachydifolia" &
    dat2$TL <= int[7,4] & dat2$TL >= int[10,4], "d7-9",
  ifelse(dat2$species == "S.stachydifolia" &
    dat2$TL >= int[1,4], "d0-1",
  ifelse(dat2$species == "S.stachydifolia" &
    dat2$TL <= int[10,4], "d10plus", NA)
  ))))))))))

dat2$categoria <- factor(dat2$categoria, 
levels = c("d10plus","d7-9","d4-6","d1-3", "d0-1", "flower"), 
labels = c("10+", "7-9", "4-6", "1-3", "< 1", "flower"))
table(dat2$categoria, dat2$species)

num_cel <- ggplot(data = dat2, aes(x = categoria, y = CN, color = species)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.75), width = 0.2) +
  stat_summary(fun = mean, geom = "line", aes(group = species), position = position_dodge(width = 0.75)) +
  scale_color_manual(values = c("#1C1089", "#EC990C")) +
  labs(x = "stages", y = "cell number") +
  guides(color = "none") +
  theme_bw()
num_cel

svg(file = "FIG5.svg", width = 3.15, height = 3)
num_cel
dev.off()

fit_nc1 <- glmer(CN ~ categoria*species + (1|ID), data = dat2, family = poisson())
fit_nc2 <- glmer(CN ~ categoria+species + (1|ID), data = dat2, family = poisson())

overdisp_fun <- function(m) { # from Ben Bolker
  rdf <- df.residual(m)
  rp <- residuals(m,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(fit_nc1)

anova(fit_nc1, fit_nc2) # no interaction
compare_performance(fit_nc1, fit_nc2)

############################
### Cell growth patterns ###
############################

# data input and preparation
dat3 <- read.csv("celulas_tam.csv", header = TRUE, dec = ".", 
                 sep = ",", stringsAsFactors = TRUE)

dat3$region <- factor(dat3$region, levels = c("basal", "medium", "distal"))
stac <- subset(dat3, dat3$species == "S.stachydifolia")
guar <- subset(dat3, dat3$species == "S.guaranitica")

dat3$logTL <- log(dat3$TL)
dat3$logW <- log(dat3$W)
dat3$logL <- log(dat3$L)
dat3$propLW <- dat3$L/dat3$W

dat3$ind <- as.factor(paste(dat3$species, dat3$ID, sep = "_"))
dat3$reg_sp <- as.factor(paste(dat3$species, dat3$region, sep = "_"))
dat3$stage <- as.factor(ifelse(dat3$bud == 5, "flower", "bud"))

dat3 <- na.omit(dat3)

summary(dat3)

# create categories based on confidence intervals
dat3$categoria <- ifelse(dat3$species == "S.guaranitica" &
  dat3$TL <= int[1,2] & dat3$TL >= int[4,2], "d1-3",
ifelse(dat3$species == "S.guaranitica" &
  dat3$TL <= int[4,2] & dat3$TL >= int[7,2], "d4-6",
ifelse(dat3$species == "S.guaranitica" &
  dat3$TL <= int[7,2] & dat3$TL >= int[10,2], "d7-9",
ifelse(dat3$species == "S.guaranitica" &
  dat3$TL <= int[10,2], "d10plus",
ifelse(dat3$species == "S.guaranitica" &
  dat3$stage == "bud" &
  dat3$TL >= int[1,2], "d0-1",
ifelse(dat3$stage == "flower", "flower",
  ifelse(dat3$species == "S.stachydifolia" &
    dat3$TL <= int[1,4] & dat3$TL >= int[4,4], "d1-3",
  ifelse(dat3$species == "S.stachydifolia" &
    dat3$TL <= int[4,4] & dat3$TL >= int[7,4], "d4-6",
  ifelse(dat3$species== "S.stachydifolia" &
    dat3$TL <= int[7,4] & dat3$TL >= int[10,4], "d7-9",
  ifelse(dat3$species == "S.stachydifolia" &
    dat3$TL >= int[1,4], "d0-1",
  ifelse(dat3$species == "S.stachydifolia" &
    dat3$TL <= int[10,4], "d10plus", NA)
  ))))))))))

dat3$categoria <- factor(dat3$categoria, 
levels = c("d10plus","d7-9","d4-6","d1-3", "d0-1", "flower"), 
labels = c("10+", "7-9", "4-6", "1-3", "< 1", "flower"))
table(dat3$categoria, dat3$species, dat3$region)

## models and plots
# cell length

largo_cel <- ggplot(data = dat3, aes(x = categoria, y = logL, color = species)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5),
    alpha = 0.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
    position = position_dodge(width = 0.5), width = 0.5, linewidth = 0.8) +
  stat_summary(fun = mean, geom = "line", aes(group = species),
    position = position_dodge(width = 0.2), linewidth = 0.8) +
  scale_color_manual(values = c("#1C1089", "#EC990C")) +
  labs(x = "stages", y = "log(cell length)") +
  facet_grid(. ~ region) +
  guides(color = "none") +
  theme_bw()
largo_cel

fit_lb1 <- lmer(logL ~ categoria * species + (1 | ID),
  data = subset(dat3, dat3$region == "basal"))
fit_lb2 <- lmer(logL ~ categoria + species + (1 | ID),
  data = subset(dat3, dat3$region == "basal"))
fit_lb3 <- lmer(logL ~ categoria + (1 | ID),
  data = subset(dat3, dat3$region == "basal"))

compare_performance(fit_lb1, fit_lb2, fit_lb3)
anova(fit_lb1, fit_lb2, test = "Chi") ## best model with interaction

fit_lm1 <- lmer(logL ~ categoria * species + (1 | ID),
  data = subset(dat3, dat3$region == "medium"))
fit_lm2 <- lmer(logL ~ categoria + species + (1 | ID),
  data = subset(dat3, dat3$region == "medium"))
fit_lm3 <- lmer(logL ~ categoria + (1 | ID),
  data = subset(dat3, dat3$region == "medium"))

compare_performance(fit_lm1, fit_lm2, fit_lm3)
anova(fit_lm1, fit_lm2, test = "Chi") ## best model with interaction

fit_ld1 <- lmer(logL ~ categoria * species + (1 | ID),
  data = subset(dat3, dat3$region == "distal"))
fit_ld2 <- lmer(logL ~ categoria + species + (1 | ID),
  data = subset(dat3, dat3$region == "distal"))
fit_ld3 <- lmer(logL ~ categoria + (1 | ID),
  data = subset(dat3, dat3$region == "distal"))

compare_performance(fit_ld1, fit_ld2, fit_ld3)
anova(fit_ld1, fit_ld2, test = "Chi")

## cell width
ancho_cel <- ggplot(data = dat3, aes(x = categoria, y = logW, color = species)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5),
    alpha = 0.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
    position = position_dodge(width = 0.5), width = 0.5, linewidth = 0.8) +
  stat_summary(fun = mean, geom = "line", aes(group = species),
    position = position_dodge(width = 0.2), linewidth = 0.8) +
  scale_color_manual(values = c("#1C1089", "#EC990C")) +
  facet_grid(. ~ region) +
  labs(x = "stages", y = "log(cell width)") +
  guides(color = "none") +
  theme_bw()
ancho_cel

fit_wb1 <- lmer(logW ~ categoria * species + (1 | ID),
  data = subset(dat3, dat3$region == "basal"))
fit_wb2 <- lmer(logW ~ categoria + species + (1 | ID),
  data = subset(dat3, dat3$region == "basal"))
fit_wb3 <- lmer(logW ~ categoria + (1 | ID),
  data = subset(dat3, dat3$region == "basal"))

compare_performance(fit_wb1, fit_wb2, fit_wb3)
anova(fit_wb1, fit_wb2, test = "Chi") ## best model with interaction

fit_wm1 <- lmer(logW ~ categoria * species + (1 | ID),
  data = subset(dat3, dat3$region == "medium"))
fit_wm2 <- lmer(logW ~ categoria + species + (1 | ID),
  data = subset(dat3, dat3$region == "medium"))
fit_wm3 <- lmer(logW ~ categoria + (1 | ID),
  data = subset(dat3, dat3$region == "medium"))

compare_performance(fit_wm1, fit_wm2, fit_wm3)
anova(fit_wm1, fit_wm2, test = "Chi") ## best model with interaction

fit_wd1 <- lmer(logW ~ categoria * species + (1 | ID),
  data = subset(dat3, dat3$region == "distal"))
fit_wd2 <- lmer(logW ~ categoria + species + (1 | ID),
  data = subset(dat3, dat3$region == "distal"))
fit_wd3 <- lmer(logW ~ categoria + (1 | ID),
  data = subset(dat3, dat3$region == "distal"))

compare_performance(fit_wd1, fit_wd2, fit_wd3)
anova(fit_wd1, fit_wd2, test = "Chi") ## best model with interaction

## prop L/W

prop_cel <- ggplot(data = dat3, aes(x = categoria, y = propLW, color = species)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5),
    alpha = 0.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
    position = position_dodge(width = 0.5), width = 0.5, linewidth = 0.8) +
  stat_summary(fun = mean, geom = "line", aes(group = species),
    position = position_dodge(width = 0.2), linewidth = 0.8) +
  scale_color_manual(values = c("#1C1089", "#EC990C")) +
  labs(x = "stages", y = "cell proportion length/width") +
  facet_grid(. ~ region) +
  guides(color = "none") +
  theme_bw()
prop_cel

fit_pb1 <- lmer(propLW ~ categoria*species + (1|ID), 
  data = subset(dat3, dat3$region == "basal"))
fit_pb2 <- lmer(propLW ~ categoria+species + (1|ID), 
  data = subset(dat3, dat3$region == "basal"))
fit_pb3 <- lmer(propLW ~ categoria + (1|ID), 
  data = subset(dat3, dat3$region == "basal"))

compare_performance(fit_pb1, fit_pb2, fit_pb3)
anova(fit_pb1, fit_pb2, test = "Chi") #best model with interaction

fit_pm1 <- lmer(propLW ~ categoria*species + (1|ID), 
data = subset(dat3, dat3$region == "medium"))
fit_pm2 <- lmer(propLW ~ categoria+species + (1|ID), 
data = subset(dat3, dat3$region == "medium"))
fit_pm3 <- lmer(propLW ~ categoria + (1|ID), 
data = subset(dat3, dat3$region == "medium"))

compare_performance(fit_pm1, fit_pm2, fit_pm3)
anova(fit_pm1, fit_pm2, test = "Chi") # best model with interaction

fit_pd1 <- lmer(propLW ~ categoria*species + (1|ID), 
data = subset(dat3, dat3$region == "distal"))
fit_pd2 <- lmer(propLW ~ categoria+species + (1|ID), 
data = subset(dat3, dat3$region == "distal"))
fit_pd3 <- lmer(propLW ~ categoria + (1|ID), 
data = subset(dat3, dat3$region == "distal"))

compare_performance(fit_pd1, fit_pd2, fit_pd3) 
anova(fit_pd1, fit_pd2, test = "Chi") # best model with interaction

svg(file = "FIG6.svg", width = 7, height = 7)
largo_cel / ancho_cel / prop_cel
dev.off()


###########
### END ###
###########

