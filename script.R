#############################################################################
# SCRIPT: Ontogenetic mechanisms of differentiation in two Salvia species 
# with different pollinators. 
# Davies, Agustín; Benitez-Vieyra, Santiago
############################################################################

library(cAIC4)
library(ggeffects)
library(ggplot2)
library(gamm4) 
library(lme4)
library(lmerTest)
library(patchwork)
library(tidymv) 

################################
### Ontogenetic trajectories ###
################################

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
cAIC(pit1)
cAIC(pit2)

# model selection - fixed component
pit3 <- lmer(slogCTL ~ slogULL + species + (1|ind) + (0 + slogULL|ind), 
             data = dat1)
pit4 <- lmer(slogCTL ~ slogULL + (1|ind) + (0 + slogULL|ind), 
             data = dat1)

anova(pit2, pit3, pit4) 
anova(pit3, pit4) 

# diagnostics
plot(predict(pit2), residuals(pit2)) 
qqnorm(residuals(pit2))

## refit model 3 and 4 with uncentered variables
dat1$logULL <- log(dat1$ULL)
dat1$logCTL <- log(dat1$CTL)
fm1 <- lmer(logCTL ~ logULL+species + (1|ind) + (0 + logULL|ind), 
            data = dat1)
fm2 <- lmer(logCTL ~ logULL + (1|ind) + (0 + logULL|ind), 
                     data = dat1)
anova(fm1, fm2)

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
svg(file = "FIG2.svg", width = 3.15, height = 3)
t1
dev.off()



##########################
### Cell proliferation ###
##########################

# data input and preparation

dat2 <- read.csv("celulas_num.csv", header = T, dec = ".", sep = ",")
dat2 <- na.omit(dat2)
stac2 <- subset(dat2, dat2$species == "sta")
guar2 <- subset(dat2, dat2$species == "gua")

dat2$logTL <- log(dat2$TL)
dat2$ID <- as.factor(paste(dat2$species, dat2$ID, sep = "_"))
dat2$species <- as.factor(dat2$species)
dat2$adul <- factor(as.factor(dat2$bud), labels = c("p", "p", "p", "p", "a"))

# Exploratory graphics
n1 <- ggplot(data = stac2, aes(x = log(TL), y = CN)) + geom_point() +
  stat_smooth() + theme_bw()
n2 <- ggplot(data = guar2, aes(x = log(TL), y = CN)) + geom_point() +
  stat_smooth() + theme_bw()
n1 | n2


# model selection
mit1 <- gamm4(CN ~ species + s(logTL, by= species), random =~ (1|ID), 
              data = dat2, family = poisson(), REML = T)
mit2 <- gamm4(CN ~ species + s(logTL), random =~ (1|ID), 
              data = dat2, family = poisson(), REML = T)

cAIC(mit1) #383.32
cAIC(mit2) #417.66

summary(mit1$gam)

# Plot
f4 <- plot_smooths(mit1$gam, series = logTL, comparison = species, transform = exp) + 
  geom_point(data = dat2, aes(x = logTL, y = CN, color = species, shape = adul), 
             size = 1.5, alpha = 0.5) +
  scale_shape_manual(values = c(19, 1)) + 
  scale_fill_manual(values = c("#1C1089", "#EC990C")) + 
  scale_color_manual(values = c("#1C1089", "#EC990C")) + theme_bw() +
  xlab("log(LSBI)") + ylab("cell number") +
  theme(legend.position = "none", axis.title=element_text(size=8))
f4

# further details were edited in inkscape
svg(file = "FIG3.svg", width = 3.15, height = 3)
f4
dev.off()

############################
### Cell growth patterns ###
############################

# data input and preparation
dat3 <- read.csv("celulas_tam.csv", header = T, dec = ".", 
                 sep = ",", stringsAsFactors = T)

dat3$region <- factor(dat3$region, levels = c("basal", "medium", "distal"))
stac <- subset(dat3, dat3$species == "S. stachydifolia")
guar <- subset(dat3, dat3$species == "S. guaranitica")

dat3$logTL <- log(dat3$TL)
dat3$logW <- log(dat3$W)
dat3$logL <- log(dat3$L)
dat3$propLW <- dat3$L/dat3$W

dat3$ind <- as.factor(paste(dat3$species, dat3$ID, sep = "_"))
dat3$reg_sp <- as.factor(paste(dat3$species, dat3$region, sep = "_"))
dat3$stage <- ifelse(dat3$bud == 5, "flower", "bud")

dat3 <- na.omit(dat3)

# Exploratory plots

# S stachydifolia
g1 <- ggplot(data = stac, aes(x = log(TL), y = log(W), color = region)) + 
  geom_point() + stat_smooth() + scale_color_viridis_d() + 
  theme_bw() + guides(color = "none")
g2 <- ggplot(data = stac, aes(x = log(TL), y = log(L), color = region)) + 
  geom_point() + stat_smooth() + scale_color_viridis_d() + 
  theme_bw() + guides(color = "none")
g3 <- ggplot(data = stac, aes(x = log(TL), y = W/L, color = region)) + 
  geom_point() + stat_smooth() + scale_color_viridis_d() + 
  theme_bw() + guides(color = "none")

(g1 | g2 | g3)

# S guaranitica
g4 <- ggplot(data = guar, aes(x = log(TL), y = log(W), color = region)) + 
  geom_point() + stat_smooth() + scale_color_viridis_d() + 
  theme_bw() + guides(color = "none")
g5 <- ggplot(data = guar, aes(x = log(TL), y = log(L), color = region)) + 
  geom_point() + stat_smooth() + scale_color_viridis_d() + 
  theme_bw() + guides(color = "none")
g6 <- ggplot(data = guar, aes(x = log(TL), y = W/L, color = region)) + 
  geom_point() + stat_smooth() + scale_color_viridis_d() + 
  theme_bw() + guides(color = "none")

(g4 | g5 | g6)


# 1. model selection - cell length
lit1 <- gamm4(logL ~ region * species + s(logTL, by= interaction(region,species)),
              random =~ (1|ind), data = dat3, REML = T)
lit2 <- gamm4(logL ~ region * species + s(logTL, by= region),
              random =~ (1|ind), data = dat3, REML = T)
lit3 <- gamm4(logL ~ region * species + s(logTL, by= species),
              random =~ (1|ind), data = dat3, REML = T)
lit4 <- gamm4(logL ~ region * species + s(logTL),
              random =~ (1|ind), data = dat3, REML = T)

cAIC(lit1) # -9.16
cAIC(lit2) #
cAIC(lit3) # 100.52
cAIC(lit4) # 111.36

summary(lit1$gam)

# Plot 
f1 <- plot_smooths(lit1$gam, series = logTL, comparison = species, facet_terms = region) + 
  geom_point(data = dat3, aes(x = logTL, y = logL, color = species, shape = stage), 
             alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(19, 1)) + 
  scale_fill_manual(values = c("#1C1089", "#EC990C")) + 
  scale_color_manual(values = c("#1C1089", "#EC990C")) + theme_bw() +
  xlab("log(LSBI)") + ylab("log(cell length)") +
  theme(legend.position = "none", axis.title=element_text(size=8))
f1 



# 2. model selection - cell width
fit1 <- gamm4(logW ~ region * species + s(logTL, by = interaction(region, species)),
              random =~ (1|ind), data = dat3, REML = T)
fit1b <- gamm4(logW ~ region * species + s(logTL, by = reg_sp), 
              random =~ (1|ind), data = dat3, REML = T)
fit2 <- gamm4(logW ~ region * species + s(logTL, by= region),
              random =~ (1|ind), data = dat3, REML = T)
fit3 <- gamm4(logW ~ region * species + s(logTL, by= species),
              random =~ (1|ind), data = dat3, REML = T)
fit4 <- gamm4(logW ~ region * species + s(logTL),
              random =~ (1|ind), data = dat3, REML = T)

cAIC(fit1)  # -539.98
cAIC(fit1b) # -538.76
cAIC(fit2)  # -467.83
cAIC(fit3)  # -331.90
cAIC(fit4)  # -289.14

summary(fit1$gam)

# Plot 
f2 <- plot_smooths(fit1$gam, series = logTL, facet_terms = region, comparison = species) + 
  geom_point(data = dat3, aes(x = logTL, y = logW, color = species, shape = stage), 
             alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(19, 1)) + 
  scale_fill_manual(values = c("#1C1089", "#EC990C")) + 
  scale_color_manual(values = c("#1C1089", "#EC990C")) + theme_bw() +
  xlab("log(LSBI)") + ylab("log(cell width)") +
  theme(legend.position = "none", axis.title=element_text(size=8))
f2 


## 3. model selection - cell proportion
kit1 <- gamm4(propLW ~ region * species + s(logTL, by= interaction(region,species)),
              random =~ (1|ind), data = dat3, REML = T)
kit2 <- gamm4(propLW ~ region * species + s(logTL, by= region),
              random =~ (1|ind), data = dat3, REML = T)
kit3 <- gamm4(propLW ~ region * species + s(logTL, by= species),
              random =~ (1|ind), data = dat3, REML = T)
kit4 <- gamm4(propLW ~ region * species + s(logTL),
              random =~ (1|ind), data = dat3, REML = T)

cAIC(kit1) # 1327.09
cAIC(kit2) # 1477.68
cAIC(kit3) # 1510.47
cAIC(kit4) # 1571.78

summary(kit1$gam)

# Plot 
f3 <- plot_smooths(kit1$gam, series = logTL, comparison = species, facet_terms = region) + 
  geom_point(data = dat3, aes(x = logTL, y = propLW, color = species, shape = stage), 
             alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(19, 1)) + 
  scale_fill_manual(values = c("#1C1089", "#EC990C")) + 
  scale_color_manual(values = c("#1C1089", "#EC990C")) + theme_bw() +
  ylim(c(-2.5, 10)) +
  xlab("log(LSBI)") + ylab("cell length/cell width") +
  theme(legend.position = "none", axis.title=element_text(size=8))
f3 

# further details were edited in inkscape
svg(file = "FIG4.svg", width = 7, height = 7)
f1/f2/f3
dev.off()

###########
### END ###
###########

