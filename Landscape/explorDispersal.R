library(dplyr)
library(ggplot2)
library(ggfortify)
library(qpcR)

## Dispersal 
# d10.20.10 es igual a: noBosque = 10,  Plantaciones = 20,  Nativo 1 + Slope = 10

f <- factor(c(1,10, 20), levels = c(1,10, 20))
exp <- expand.grid(f, f, f)
di <- rbind(exp[exp[, 1] == 20,], exp[exp[, 2] == 20,], exp[exp[, 3] == 20,]) %>% distinct()
attach(di)
di[order(Var3, Var2),]
apply(di, 1, paste, collapse = ".") %>% paste("d", ., sep ="")

pc <- read.csv("C:/Trabajos/Manuscritos/Carpintero/3_Dispersal footprint/pcSpat.csv", header = T)
t <- read.csv("C:/Trabajos/Manuscritos/Carpintero/3_Dispersal footprint/GEE/GEEdispersal.csv", header = T)
t <- apply(t, 2, function(x) x / max(x))
# t <- apply(t, 2, function(x) -log(x))
pc <- cbind(pc, t)
pc$picNu <- as.factor(pc$picNu)
for(i in 1:nrow(pc)) if(is.na(pc$dr[i]) == T) pc$dr[i] <- 0
for(i in 1:nrow(pc)) if(is.na(pc$det[i]) == T){ pc$dr[i] <- NA }else{ if(pc$det[i] == 1) pc$dr[i] <- 1}

## glm logit ----

# Detections (USE "dr"!! - see above)
sum <- data.frame(mod = colnames(pc)[9:ncol(pc)], aic = NA)
for(i in 9:ncol(pc)){
  fit <- glm(pc$dr ~ pc[, i], family = "binomial")
  # resid(fit) %>% plot %>% abline(0, 0)
  sum$aic[i-8] <- (fit %>% summary)$aic
}
sum <- sum[order(sum$aic),] 
min <- sum[1, 2]
sum %>% mutate(aic = as.numeric(aic), delta = aic - min)

# New pecking
NewPecking <- list()
sum <- data.frame(mod = colnames(pc)[9:ncol(pc)], aic = NA)
for(i in 9:ncol(pc)){
  fit <- glm(pc$picNu ~ pc[, i], family = "binomial")
  NewPecking[i-8] <- fit %>% summary
  # resid(fit) %>% plot %>% abline(0, 0)
  sum$aic[i-8] <- (fit %>% summary)$aic
}
sum <- sum[order(sum$aic),] 
min <- sum[1, 2]
np <- sum %>% mutate(aic = as.numeric(aic), delta = aic - min)
npPar <- sub("d", "", np$mod)
npPar <- scan(text = npPar, sep = ".", quiet = TRUE) %>% matrix(nrow = length(npPar), ncol = 3, byrow = T) %>% data.frame
colnames(npPar) <- c("nontree", "plant", "slope")
np <- np %>% bind_cols(npPar)
npSel <- np[np$delta < 3, ]
npSel <- npSel %>% mutate(weight = akaike.weights(npSel$aic)$weights)
npSel <- npSel %>% mutate(nontreeW = nontree * weight, plantW = plant * weight, slopeW = slope * weight)
apply(npSel %>% select(nontreeW, plantW, slopeW), 2, sum)

glm(pc$picNu ~ pc[, which(pc %>% colnames == np[1,1])], family = "binomial") %>% summary


# Old pecking 
OldPecking <- list()
sum <- data.frame(mod = colnames(pc)[9:ncol(pc)], aic = NA)
for(i in 9:ncol(pc)){
  fit <- glm(pc$picVi ~ pc[, i], family = "binomial")
  OldPecking[i-8] <- fit %>% summary
  # resid(fit) %>% plot %>% abline(0, 0)
  sum$aic[i-8] <- (fit %>% summary)$aic
}
sum <- sum[order(sum$aic),] 
min <- sum[1, 2]
op <- sum %>% mutate(aic = as.numeric(aic), delta = aic - min)
opPar <- sub("d", "", op$mod)
opPar <- scan(text = opPar, sep = ".", quiet = TRUE) %>% matrix(nrow = length(opPar), ncol = 3, byrow = T) %>% data.frame
colnames(opPar) <- c("nontree", "plant", "slope")
op <- op %>% bind_cols(opPar)
opSel <- op[op$delta < 3, ]
opSel <- opSel %>% mutate(weight = akaike.weights(opSel$aic)$weights)
opSel <- opSel %>% mutate(nontreeW = nontree * weight, plantW = plant * weight, slopeW = slope * weight)
apply(opSel %>% select(nontreeW, plantW, slopeW), 2, sum)

glm(pc$picVi ~ pc[, which(pc %>% colnames == op[1,1])], family = "binomial") %>% summary

colnames(pc[9:ncol(pc)])
seq(1, 20, by = 0.2)

# Recycle
pc.a <- pc %>% select(d1.1.1, d1.1.20, d1.20.1, d10.1.1) 
pc.a <- apply(pc.a, 2, as.numeric)
autoplot(prcomp(pc.a), data = pc.a, loadings = T, loadings.label = T)
corrgram(pc, lower.panel = panel.shade, upper.panel = panel.pts)

# Plots exploratorios 
ggplot(pc, aes(x = -log(d10.20.1), y = -log(d1.1.10))) + 
  geom_point(aes(color = picNu))

ggplot(pc, aes(x = picNu, y =  1/log(d1.20.1))) + 
  geom_boxplot(aes(color = picNu))

ggplot(pc, aes(x = picNu, y = d1.1.1)) + 
  geom_violin(aes(color = picNu))

ggplot(pc, aes(x = picNu, y = -log(d1.1.1))) + 
  geom_boxplot(aes(color = picNu))

ggplot(pc, aes(x = picNu, y = -log(d1.1.1))) + 
  geom_violin(aes(color = picNu)) + 
  geom_boxplot(width=0.1)

pc$d1.10.1 %>% log %>% hist 
