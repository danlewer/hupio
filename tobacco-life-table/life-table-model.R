options(scipen = 999)

# ===============================
# libraries and general functions
# -------------------------------

library(data.table)
library(splines)
library(RColorBrewer)
library(devEMF) # for enhanced metafile (Windows vector graphic)

# functions to assist plotting

roundup <- function (x, dig = 2) {
  l <- floor(log10(x)) + 1
  r <- ceiling(x / (10 ^ (l-dig)))
  r * 10 ^ (l-dig)
}

yax <- function (x, tickabove = F, ntick = 5) { # make y-axis
  l <- c(c(1, 2, 4, 5, 25) %o% 10^(0:8))
  d <- l[which.min(abs(x/ntick - l))]
  d <- 0:(ntick+1) * d
  i <- findInterval(x, d)
  if (tickabove) {i <- i + 1}
  d[seq_len(i)]
}

add.alpha <- function (cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)

# all-cause life table

life.table <- function (mx, cohort = 100000, EX = T) { # if EX is true, just return life expectancy
  n <- length(mx) + 1
  qx <- 2 * mx / (2 + mx)
  qx <- c(qx, 1) # forced method - mortality rate max age + 1 is 100%
  lx <- c(1, cumprod(1 - qx)) * cohort
  dx <- -c(diff(lx), lx[n] * qx[n])
  t <- (lx + c(lx[-1], 0)) / 2
  Tx <- rev(cumsum(rev(t)))
  ex <- Tx / lx
  lt <- data.frame(lx = lx, dx = dx, t = t, Tx = Tx, ex = ex)[1:n,]
  if (EX) ex[1] else lt
}

# cause-specific life table function
# this works by creating a life table based on the all-cause mortality rates
# it then applies the proportions of deaths at each age
# the 'mat' argument is a matrix of cause-specific mortality rates

cs.life.table <- function (mat, cohort = 100000) {
  # all-cause life table
  mx <- rowSums(mat)
  n <- length(mx) + 1
  qx <- 2 * mx / (2 + mx)
  qx <- c(qx, 1) # forced method - mortality rate max age + 1 is 100%
  lx <- c(1, cumprod(1 - qx)) * cohort
  dx <- -c(diff(lx), lx[n] * qx[n])
  # make life table data frame
  lt <- data.frame(lx = lx, dx = dx)[1:n,]
  # matrix of proportions of deaths by cause at each age
  cm <- rbind(mat, 1) / rowSums(rbind(mat, 1))
  cs.dx <- lt$dx * cm
  colnames(cs.dx) <- paste0('dx_', colnames(mat))
  # outputs
  cbind(lt, cs.dx)[-nrow(lt),]
}

# ==============================
# read data and aggregate by sex
# ------------------------------

# data is from the 'HUPIO' study of people who use illicit opioids
# Codelist: https://wellcomeopenresearch.org/articles/5-282
# Causes of death: https://www.thelancet.com/journals/lanpub/article/PIIS2468-2667(21)00254-1/fulltext

m <- fread("https://raw.githubusercontent.com/danlewer/hupio/main/mortality/single-year-of-age-rates/hupio_rates_drugs_vs_smoking_25aug2023.csv", drop = 'gender', col.names = c('opioids', 'age', 'smoking', 'cvd', 'other_cancers', 'drug_related', 'hep', 'other', 'follow_up', 'all_cause'))
m <- m[opioids == T & age <= 70, lapply(.SD, sum), by = 'age'][, -'opioids']

# =====================================
# estimate age-specific mortality rates
# -------------------------------------

causes <- c('all_cause', 'smoking', 'drug_related', 'hep', 'cvd', 'other_cancers', 'other')
titles <- c('All cause', 'Smoking\n-specific\n(respiratory cancers\nand COPD)', 'Drug\npoisoning', 'Viral\nhepatitis', 'Cardiovascular\ndiseases', 'Other\ncancers\n(non-respiratory)', 'Other\nunderlying causes\nof death')

# estimate age-specific mortality rates using poisson model and cubic splines

pf <- function (outcome, nd = data.table(age = 18:70, follow_up = 1), ci = T, data = m) {
  f <- as.formula(paste0(outcome, '~ ns(age, 3) + offset(log(follow_up))'))
  model <- glm(f, data = data, family = 'poisson')
  inv <- model$family$linkinv
  p <- predict(model, newdata = nd, se.fit = T)
  if (ci == T) {
    cbind(nd, 
          pred = inv(p$fit), 
          lower = inv(p$fit - qnorm(0.975) * p$se.fit),
          upper = inv(p$fit + qnorm(0.975) * p$se.fit))
  } else {
    inv(p$fit)
  }
}
pm <- lapply(causes, pf)
names(pm) <- causes

# does the sum of modelled cause specific rates approximate the all-cause rate?

sum_cause_specific <- rowSums(sapply(pm[causes[-1]], function (x) x$pred))
plot(18:70, pm$all_cause$pred * 1e5, type = 'l', col = 'red', xlab = 'Age', ylab = 'Modelled mortality rate / 100,000', main = 'All-cause (red) vs.\nsum of cause-specific (blue)')
lines(18:70, sum_cause_specific * 1e5, col = 'blue')

# visualise mortality rates by age with confidence intervals

cols <- brewer.pal(3, 'Set1')[3:1]

emf('FigA.emf', height = 6.5, width = 8, units = 'in')

par(mfrow = c(2, 3), mar = c(2, 3, 1, 1), oma = c(3, 3, 0, 0))
lapply(1:6, function (i) {
  y <- m[, get(causes[i])] / m$follow_up * 1e5
  ymax <- roundup(max(y), 2)
  plot(1, type = 'n', xlab = NA, ylab = NA, xlim = c(18, 70), ylim = c(0, ymax), axes = F)
  rect(18, 0, 70, ymax, col = 'grey98')
  points(18:70, y)
  axis(1, c(18, 3:7 * 10), pos = 0)
  axis(2, yax(ymax), pos = 18, las = 2)
  text(20, ymax * 0.98, titles[i], adj = c(0, 1))
  with(pm[causes[i]][[1]], {
    polygon(c(age, rev(age)), c(lower, rev(pmin(upper, ymax/1e5))) * 1e5, border = NA, col = add.alpha(cols[1], 0.3))
    lines(18:70, pred * 1e5, col = cols[1])
  })
})
mtext('Age', side = 1, line = 1, outer = T, cex = 0.8)
mtext('Deaths per 100,000 person-years', side = 2, line = 1, outer = T, cex = 0.8)

dev.off()

# visualise mortality rates by age in single plot

causes2 <- causes[-c(1, 7)]
titles <- c('Respiratory cancers\nand COPD', 'Drug\npoisoning', 'Viral\nhepatitis', 'Cardiovascular\ndiseases', 'Cancers other\nthan repiratory')
cols <- brewer.pal(5, 'Set1')
pchs <- 0:4
final_vals <- sapply(pm[causes2], function (x) x$pred)[53,]

png('Fig1.png', height = 6, width = 9, units = 'in', res = 300)

par(mar = c(5, 5, 0, 13), xpd = NA)
plot(1, type = 'n', xlim = c(18, 70), ylim = c(0, 1300), axes = F, xlab = NA, ylab = 'Mortality rate per 100,000 pys')
for (i in seq_along(causes2)) {
  y <- m[, get(causes2[i])] / m$follow_up * 100000
  y[y > 1300] <- NA
  points(m$age, y, col = cols[i], cex = 0.7, pch = pchs[i])
  y <- pm[[causes2[i]]]$pred * 100000
  lines(m$age, pmin(y, 1300), col = cols[i])
}
axis(1, 2:7 * 10, pos = 0)
axis(1, c(18, 70), labels = F, pos = 0)
axis(2, 0:13 * 100, pos = 18, las = 2)
rect(18, 0, 70, 1300)
title(xlab = 'Age', line = 2)
ys <- seq(400, 1200, length.out = 5)
segments(72, ys, 78, col = cols[order(final_vals)])
points(rep(75, 5), ys, pch = pchs[order(final_vals)], col = cols[order(final_vals)])
text(79, ys, titles[order(final_vals)], adj = 0)

dev.off()

# ===========================
# Attributable risk fractions
# ---------------------------

fractions <- c(smoking = 1, cvd = 0.5, other_cancers = 0.3, drug_related = 1, hep = 1)

# ==============================================
# estimate deaths and YLLs in simulated datasets
# ----------------------------------------------

# simulate deaths & YLLs in given scenario

deaths_ylls <- function (sims = 100, 
                         smoking = 1, 
                         drugs = 1,
                         point = F,
                         DATA = m,
                         summary_only = T) {
  lapply (1 : if (point) 1 else sims, function(x) {
    if (x %% 100 == 0) print (x)
    m_sim <- sapply(DATA[, causes, with = F], function (x) rpois(length(x), x))
    m_sim <- cbind(DATA[, c('age', 'follow_up')], m_sim)
    mr <- cbind(age = 18:70, as.data.frame(`names<-`(lapply(causes, pf, ci = F, data = if (point) DATA else m_sim), causes))) # modelled rates
    mat <- cbind(smoking_smoking = mr$smoking * fractions['smoking'] * smoking,
                 smoking_notSmoking = mr$smoking * (1-fractions['smoking']),
                 cvd_smoking = mr$cvd * fractions['cvd'] * smoking,
                 cvd_notSmoking = mr$cvd * (1-fractions['cvd']),
                 otherCancers_smoking = mr$other_cancers * fractions['other_cancers'] * smoking,
                 otherCancers_notSmoking = mr$other_cancers * (1-fractions['other_cancers']),
                 hep_drugs = mr$hep * fractions['hep'] * drugs,
                 hep_notDrugs = mr$hep * (1-fractions['hep']),
                 drugRelated_drugs = mr$drug_related * fractions['drug_related'] * drugs,
                 drugRelated_notDrugs = mr$drug_related * (1-fractions['drug_related']),
                 other = mr$other)
    cs.lt <- cs.life.table(mat)
    r <- list(life_table = cs.lt,
              summary = rbind(deaths = colSums(cs.lt[-53,])[-(1:2)],
                              ylls = colSums(cs.lt[-53,] * (69 - 18:69 + 0.5))[-(1:2)]))
    if (summary_only) r$summary else r
  })
}

# generate results in different scenarios of smoking and drugs

set.seed(562)
B <- 1000 # number of sims

Ap <- deaths_ylls(point = T, smoking = 1, drugs = 1)
Am <- deaths_ylls(sims = B, smoking = 1, drugs = 1)
Bp <- deaths_ylls(point = T, smoking = 1, drugs = 0)
Bm <- deaths_ylls(sims = B, smoking = 1, drugs = 0)
Cp <- deaths_ylls(point = T, smoking = 0, drugs = 1)
Cm <- deaths_ylls(sims = B, smoking = 0, drugs = 1)
Dp <- deaths_ylls(point = T, smoking = 0, drugs = 0)
Dm <- deaths_ylls(sims = B, smoking = 0, drugs = 0)

# function to summarise results of simulations 

smoking_vars <- c('dx_smoking_smoking', 'dx_cvd_smoking', 'dx_otherCancers_smoking')
drug_vars <- c('dx_hep_drugs', 'dx_drugRelated_drugs')
s <- function (d, i = 1) { # 1 = deaths, 2 = YLLS
  if (length(d) == 1) {
    x <- d[[1]][i,]
    x <- c(all = sum(x), smoking = sum(x[smoking_vars]), drugs = sum(x[drug_vars]))
    x <- c(x, other = x[1] - (x[2] + x[3]))
    return(c(x, pc_smoking = x[2] / x[1] * 100, pc_drugs = x[3] / x[1] * 100, pc_other = x[4] / x[1] * 100))
  }
  x <- t(sapply(d, function (x) x[i,]))
  x <- cbind(all = rowSums(x), smoking = rowSums(x[, smoking_vars]), drugs = rowSums(x[, drug_vars]))
  x <- cbind(x, other = x[,1] - (x[,2] + x[,3]))
  cbind(x, pc_smoking = x[,2] / x[,1] * 100, pc_drugs = x[,3] / x[,1] * 100, pc_other = x[,4] / x[,1] * 100)
}

# create deaths results

f <- function (x, digs = 1, units = c(1000, 1000, 1000, 1000, 1, 1, 1)) {
  x <- t(t(x) / units)
  x <- formatC(round(x, digs), big.mark = ',', format = 'f', digits = digs)
  `names<-`(paste0(x[1,], ' (', x[3,], ', ', x[4,], ')'), colnames(x))
}
deaths_scenarios <- rbind(cbind(observed = f(rbind(point = s(Ap), apply(s(Am), 2, quantile, probs = c(0.5, 0.025, 0.975)))),
                                eliminate_drugs = f(rbind(point = s(Bp), apply(s(Bm), 2, quantile, probs = c(0.5, 0.025, 0.975)))),
                                eliminate_smoking = f(rbind(point = s(Cp), apply(s(Cm), 2, quantile, probs = c(0.5, 0.025, 0.975)))),
                                eliminate_both = f(rbind(point = s(Dp), apply(s(Dm), 2, quantile, probs = c(0.5, 0.025, 0.975))))),
                          cbind(observed = NA,
                                eliminate_drugs = f(rbind(point = s(Bp) - s(Ap), apply(s(Bm) - s(Am), 2, quantile, probs = c(0.5, 0.025, 0.975)))),
                                eliminate_smoking = f(rbind(point = s(Cp) - s(Ap), apply(s(Cm) - s(Am), 2, quantile, probs = c(0.5, 0.025, 0.975)))),
                                eliminate_both = f(rbind(point = s(Dp) - s(Ap), apply(s(Dm) - s(Am), 2, quantile, probs = c(0.5, 0.025, 0.975))))))

# create YLL results

f2 <- function (x) f(x, digs = 2, units = c(100000, 100000, 100000, 100000, 1, 1, 1))
yll_scenarios <- rbind(cbind(observed = f2(rbind(point = s(Ap, 2), apply(s(Am, 2), 2, quantile, probs = c(0.5, 0.025, 0.975)))),
                             eliminate_drugs = f2(rbind(point = s(Bp, 2), apply(s(Bm, 2), 2, quantile, probs = c(0.5, 0.025, 0.975)))),
                             eliminate_smoking = f2(rbind(point = s(Cp, 2), apply(s(Cm, 2), 2, quantile, probs = c(0.5, 0.025, 0.975)))),
                             eliminate_both = f2(rbind(point = s(Dp, 2), apply(s(Dm, 2), 2, quantile, probs = c(0.5, 0.025, 0.975))))),
                       cbind(observed = NA,
                             eliminate_drugs = f2(rbind(point = s(Bp, 2) - s(Ap, 2), apply(s(Bm, 2) - s(Am, 2), 2, quantile, probs = c(0.5, 0.025, 0.975)))),
                             eliminate_smoking = f2(rbind(point = s(Cp, 2) - s(Ap, 2), apply(s(Cm, 2) - s(Am, 2), 2, quantile, probs = c(0.5, 0.025, 0.975)))),
                             eliminate_both = f2(rbind(point = s(Dp, 2) - s(Ap, 2), apply(s(Dm, 2) - s(Am, 2), 2, quantile, probs = c(0.5, 0.025, 0.975))))))

# save table

scenarios <- rbind(deaths_scenarios, yll_scenarios)
write.csv(scenarios, 'scenarios_table.csv')

# =========================
# plot of cumulative deaths
# -------------------------

cuml_deaths <- deaths_ylls(point = T, summary_only = F)[[1]]$life_table
cuml_deaths <- sapply(cuml_deaths, cumsum)
cuml_deaths <- subset(cuml_deaths, select = -c(dx, lx))
cuml_deaths <- cuml_deaths[, colSums(cuml_deaths) > 0]
cuml_deaths <- cuml_deaths[, c('dx_smoking_smoking', 'dx_cvd_smoking', "dx_otherCancers_smoking", "dx_hep_drugs", "dx_drugRelated_drugs", 'dx_cvd_notSmoking', "dx_otherCancers_notSmoking", "dx_other")]
ynums <- apply(cuml_deaths, 1, cumsum)
ynums <- rbind(0, ynums)
cols <- c(brewer.pal(3, 'Oranges'),
          brewer.pal(3, 'Purples')[c(2, 3)],
          brewer.pal(3, 'Greens'))
sdo70 <- sapply(split(cuml_deaths[53,], c(1, 1, 1, 2, 2, 3, 3, 3)), sum)
sdo69 <- sapply(split(cuml_deaths[52,], c(1, 1, 1, 2, 2, 3, 3, 3)), sum)
labs3 <- c('Attributable to tobacco', 'Attributable to illicit drugs', 'Other deaths')
labs3 <- paste0(labs3,
                "\n", paste0(round(sdo69 / 1000, 1), '%', ' cumulative risk'),
                "\n", paste0(round(sdo69 / sum(sdo69) * 100, 1), '%', ' of premature deaths'))

png('Fig2.png', height = 7, width = 10, res = 300, units = 'in')
par(mar = c(5, 6, 0, 15), xpd = NA)
plot(1, type = 'n', xlim = c(18, 70), ylim = c(0, 100000), xlab = NA, ylab = NA, axes = F)
for (i in 1:8) {
  polygon(c(18:70, 70:18), y = c(ynums[i,], rev(ynums[i+1,])), col = cols[i], border = 'black')
}
rect(18, 0, 70, 100000)
axis(1, pos = 0)
axis(2, 0:5 * 20000, paste0(0:5 * 20, '%'), las = 2, pos = 18)
title(xlab = 'Age', line = 2)
title(ylab = 'Cumulative risk of death in a cohort of\npeople who use illicit opioids', line = 4)
ys <- seq(20000, 97500, length.out = 9)
rect(20, ys[-length(ys)], 23, ys[-1], col = cols)
text(25, ys[-length(ys)] + diff(ys)/2, c('Respiratory cancers\nand COPD', 'Cardiovascular diseases\nattributable to tobacco', 'Cancers other than respiratory\nattributable to tobacco', 'Viral\nhepatitis', 'Drug\npoisoning', 'Cardiovascular diseases\nnot attributable to tobacco', 'Cancers not attributable\nto tobacco', 'All other\ndeaths'), adj = 0)
arrows(73, c(0, cumsum(sdo70)[-3]) + 500, y1 = cumsum(sdo70) - 500, code = 3, length = 0.1, angle = 45, col = cols[c(3, 5, 8)])
text(75, c(0, cumsum(sdo70))[-4] + diff(c(0, cumsum(sdo70)))/2, labs3, adj = 0)
dev.off()

# =============================================
# bar plot of attributable and prevented deaths
# ---------------------------------------------

x <- rbind(s(Ap), s(Bp), s(Cp), s(Dp)) / 1e5
#x <- x[, c(1, 4:2)]
y1 <- t(apply(cbind(0, x[,-1]), 1, cumsum))
xl <- 0:3 * 4
delta <- -t(x[1,] - t(x))
waterfall_targets <-t(apply(delta[-1,-1], 1, function (y) x[1,1] + cumsum(y)))

# waterfall chart

waterfall <- function (start, targets, xleft, width = 0.5, cols = 1:3) {
  yvals <- cbind(c(0, start, targets), c(start, targets, 0))[2:(length(targets)+1),]
  mapply(rect,
         xleft = xleft,
         ybottom = apply(yvals, 1, min),
         xright = xleft + width,
         ytop = apply(yvals, 1, max),
         col = cols)
  segments(x0 = xleft, y0 = targets, x1 = xleft + width * 2, lty = 3)
}

cols <- brewer.pal(3, 'Set2')
ys <- seq(0.4, 0.7, length.out = 4)

png('Fig3.png', height = 5, width = 8, units = 'in', res = 300)
par(mar = c(6, 5, 0, 10), xpd = NA)
plot(1, type = 'n', xlim = c(0, 12.75), ylim = c(0, 0.7), axes = F, xlab = NA, ylab = NA)
rect(-0.25, 0:6/10, 12.75, 1:7/10, col = rep(c('white', 'grey93'), 4), border = NA)
rect(-0.25, 0, 12.75, 0.7)
segments(-0.25, x[1,1], x1 = 12.75, lty = 3)
mapply(rect,
       xleft = xl,
       ybottom = split(y1[,-4], f = 1:4),
       xright = xl + 0.5,
       ytop = split(y1[,-1], f = 1:4),
       col = rep(list(cols), each = 4))
waterfall(start = x[1,1], targets = waterfall_targets[1,], xleft = 1:3, cols = cols)
waterfall(start = x[1,1], targets = waterfall_targets[2,], xleft = 5:7, cols = cols)
waterfall(start = x[1,1], targets = waterfall_targets[3,], xleft = 9:11, cols = cols)
segments(c(-0.25, 0.75, 4.75, 8.75, 12.75), y0 = 0, y1 = 0.7)
axis(2, 0:7 / 10, paste0(0:7 * 10, '%'), las = 2, pos = -0.25)
title(ylab = 'Risk of death before age 70')
text(0.25, -0.01, 'Observed\nmortality rates', srt = 90, adj = 1)
text(c(2.75, 6.75, 10.75), -0.01, c('Eliminating\nillicit\ndrugs', 'Eliminating\ntobacco\nsmoking', 'Eliminating illicit\ndrugs and\ntobacco\nsmoking'), adj = c(0.5,1))
rect(13, ys[-length(ys)], 13.5, ys[-1], col = cols)
text(13.75, ys[-length(ys)] + diff(ys)/2, c('Attributable to\nillicit drugs', 'Attributable to\ntobacco smoking', 'Attributable to\nother causes'), adj = 0)
dev.off()

# ==================================
# comparison with general population
# ----------------------------------

# gen pop data

mgp <- fread("https://raw.githubusercontent.com/danlewer/hupio/main/mortality/single-year-of-age-rates/hupio_rates_drugs_vs_smoking_25aug2023.csv", drop = 'gender', col.names = c('opioids', 'age', 'smoking', 'cvd', 'other_cancers', 'drug_related', 'hep', 'other', 'follow_up', 'all_cause'))
mgp <- mgp[opioids == F & age <= 70, lapply(.SD, sum), by = 'age'][, -'opioids']
s(deaths_ylls(point = T, DATA = mgp))
s(deaths_ylls(point = T, DATA = mgp), i = 2) / 100000
