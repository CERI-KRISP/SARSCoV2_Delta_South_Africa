# Analyzing the spread of SARS-CoV-2 variants of concern (VOC) in South Africa using multinomial logistic regression
# Christian L. Althaus, 25 September 2021

# Load libraries
library(tidyverse)
library(RColorBrewer)
library(lubridate)
library(nnet)
library(splines)
library(emmeans)

# Define colors
cols <- brewer.pal(9, "Set1")
cols <- cols[c(3, 2, 1, 4)]
t.cols <- cols
for(i in 1:length(cols)) {
  x <- col2rgb(cols[i])
  t.cols[i] <- rgb(x[1, ], x[2, ], x[3, ], alpha = 125, maxColorValue = 255)
}

# Read data (download variant surveillance data from GISAID)
country <- "South Africa"
gisaid_full <- read_tsv("variant_surveillance.tsv")
gisaid_full <- as.data.frame(gisaid_full)
gisaid_country <- subset(gisaid_full, grepl(country, Location))
gisaid <- gisaid_country
table(gisaid$`Pango lineage`)

# Set variants and start date
variants <- c("B.1.1.7", "B.1.351", "B.1.351.2", "B.1.351.3", "P.1", "P.1.1", "P.1.2", "B.1.617.2", "AY.1", "AY.2", "AY.3", "AY.4", "AY.5", "AY.6", "AY.7", "AY.8", "AY.9", "AY.10", "AY.11", "AY.12")
variant_names <- c("Alpha", "Beta", "Gamma", "Delta")
begin <- ymd(20200801)

# Clean data (only consider variants when >= 25 samples)
gisaid$date <- ymd(gisaid$`Collection date`)
gisaid <- subset(gisaid, !is.na(date) & date >= begin)
gisaid$variant <- gisaid$`Pango lineage`
w <- which(table(gisaid$variant) >= 25)
variants <- intersect(variants, names(w))
w <- is.element(gisaid$variant, variants)
sort(table(gisaid$variant))
gisaid$variant[is.element(gisaid$variant, "B.1.1.7")] <- "Alpha"
gisaid$variant[is.element(gisaid$variant, c("B.1.351", "B.1.351.2", "B.1.351.3"))] <- "Beta"
gisaid$variant[is.element(gisaid$variant, c("P.1", "P.1.1", "P.1.2"))] <- "Gamma"
gisaid$variant[is.element(gisaid$variant, c("B.1.617.2", "AY.1", "AY.2", "AY.3", "AY.4", "AY.5", "AY.6", "AY.7", "AY.8", "AY.9", "AY.10", "AY.11", "AY.12"))] <- "Delta"
gisaid$variant[!w] <- "Other"
gisaid <- subset(gisaid, select = c("date", "variant"))
gisaid$variant <- as.factor(gisaid$variant)
gisaid$variant2 <- relevel(gisaid$variant, ref = "Beta")
table(gisaid$variant)
table(gisaid$variant2)

# Multinomial logistic regression
gisaid_fit <- gisaid
gisaid_fit$date <- as.numeric(gisaid_fit$date) - as.numeric(min(gisaid_fit$date))
fit <- multinom(variant2 ~ ns(date, 2), data = gisaid_fit)
summary(fit)

# Convert to weekly data
gisaid <- cbind(gisaid, week = round_date(gisaid$date, unit = "week", week_start = 4))
week1 <- min(gisaid$week)
week2 <- max(gisaid$week)
week_y <- seq(week1, week2, 7)
var_prop <- as.data.frame(matrix(NA, nrow = length(week_y), ncol = 3 + length(variant_names)))
names(var_prop) <- c("date", "samples", variant_names, "Other")
var_prop$date <- week_y

for(i in week_y) {
  x <- subset(gisaid, week == i)
  var_prop[var_prop$date == i, 2] <- length(x$variant2)
  var_prop[var_prop$date == i, 3] <- length(x$variant2[x$variant2 == "Alpha"])/length(x$variant2)
  var_prop[var_prop$date == i, 4] <- length(x$variant2[x$variant2 == "Beta"])/length(x$variant2)
  var_prop[var_prop$date == i, 5] <- length(x$variant2[x$variant2 == "Gamma"])/length(x$variant2)
  var_prop[var_prop$date == i, 6] <- length(x$variant2[x$variant2 == "Delta"])/length(x$variant2)
  var_prop[var_prop$date == i, 7] <- length(x$variant2[x$variant2 == "Other"])/length(x$variant2)
}
var_prop

# Estimating the proportion of variants over time
t1 <- min(gisaid_fit$date)
t2 <- max(gisaid_fit$date)
prediction_date <- as_date(min(gisaid$date):max(gisaid$date))
prediction <- data.frame(emmeans(fit, ~ variant2, mode = "prob", by = "date", at = list(date = t1:t2)))

# Estimating the growth difference between Delta and Beta
t3 <- max(prediction$date[prediction$variant2 == "Delta" & prediction$prob < 0.5])
fit_emtrends <- emtrends(fit, trt.vs.ctrl ~ variant2, var = "date", mode = "latent", at = list(date = t3))
fit_emtrends
confint(fit_emtrends)$contrasts

# Plot results
pdf(paste0("variants_multinomial_", country, ".pdf"), width = 12, height = 4)
par(mfrow = c(1, 2))
plot(var_prop$date, var_prop$sample, ty = "n",
     xlim = c(min(gisaid$date), min(gisaid$date) + t2), ylim = c(0, 1),
     xlab = NA, ylab = "Proportion variant", axes = FALSE, frame = FALSE)
abline(h = seq(0, 1, 0.2), col = "lightgray", lty = 3)
abline(v = ymd(20200901) + 0:6*months(2), col = "lightgray", lty = 3)
axis(1, ymd(20200901) + 0:6*months(2), c("Sep\n2020", "Nov\n2020", "Jan\n2021", "Mar\n2021", "May\n2021", "Jul\n2021", "Sep\n2021"), padj = 0.5)
axis(2)

for(i in 1:length(levels(prediction$variant2))) {
  x <- subset(prediction, variant2 == levels(prediction$variant2)[i])
  x$date <- min(gisaid$date) + x$date
  polygon(c(x$date, rev(x$date)),
          c(x$lower.CL, rev(x$upper.CL)),
          col = t.cols[i], border = NA)
  lines(x$date, x$prob, col = cols[i])
}

for(i in 1:length(levels(prediction$variant2))) {
  points(var_prop$date, var_prop[, levels(prediction$variant2)[i]], pch = 19, cex = 1.5*var_prop$samples/max(var_prop$samples), col = cols[i])
}

legend("left", inset = 0, levels(prediction$variant2), col = cols, pch = 19, bty = "n")

plot(var_prop$date, var_prop$sample, ty = "n",
     xlim = c(min(gisaid$date), min(gisaid$date) + t2), ylim = c(1e-2, 1), log = "y",
     xlab = NA, ylab = "Proportion variant (log10)", axes = FALSE, frame = FALSE)
abline(h = 10^(-3:0), col = "lightgray", lty = 3)
abline(v = ymd(20200901) + 0:6*months(2), col = "lightgray", lty = 3)
axis(1, ymd(20200901) + 0:6*months(2), c("Sep\n2020", "Nov\n2020", "Jan\n2021", "Mar\n2021", "May\n2021", "Jul\n2021", "Sep\n2021"), padj = 0.5)
axis(2)

for(i in 1:length(levels(prediction$variant2))) {
  x <- subset(prediction, variant2 == levels(prediction$variant2)[i])
  x$date <- min(gisaid$date) + x$date
  polygon(c(x$date, rev(x$date)),
          c(x$lower.CL, rev(x$upper.CL)),
          col = t.cols[i], border = NA)
  lines(x$date, x$prob, col = cols[i])
}

for(i in 1:length(levels(prediction$variant2))) {
  points(var_prop$date, var_prop[, levels(prediction$variant2)[i]], pch = 19, cex = 1.5*var_prop$samples/max(var_prop$samples), col = cols[i])
}
dev.off()
