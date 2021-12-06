options(scipen = 999)

library(data.table)
library(RColorBrewer)
library(data.table)
library(devEMF)
library(stringi)

vpt <- function(x, t, digs = 2, ...) {
  f <- function(xl, tl) {
    y <- poisson.test(xl, tl, ...)
    c(xl/tl, y$conf.int[1:2])
  }
  mapply(f, xl = x, tl = t)
}
vf <- function(res, digs = 2, ci_only = F) {
  res <- format(round(res, digs), digits = digs, nsmall = digs)
  res <- if (ci_only) {
    apply(res, 2, function(x) paste0(x[2], '-', x[3]))
  } else {
    apply(res, 2, function(x) paste0(x[1], '(', x[2], '-', x[3], ')'))
  }
  res <- gsub(' ', '', res)
  gsub('\\(', ' (', res)
}
add.alpha <- function (cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)
rnd <- function(x, dig = 1) format(round(x, digits = dig), nsmall = dig, digits = dig)
yax <- function(x, tickabove = F, ntick = 5) {
  l <- c(c(1, 2, 4, 5, 25) %o% 10^(0:8))
  d <- l[which.min(abs(x/ntick - l))]
  d <- 0:(ntick+1) * d
  i <- findInterval(x, d)
  if (tickabove) {i <- i + 1}
  d[seq_len(i)]
}

cs_smr <- read.csv('https://raw.githubusercontent.com/danlewer/hupio/main/mortality/lph_fig1.csv')
periods <- read.csv('https://raw.githubusercontent.com/danlewer/hupio/main/mortality/lph_fig2.csv', col.names = c('exposure', 'cy4', 'cause', 's', 'l', 'u'))
ages <- read.csv('https://raw.githubusercontent.com/danlewer/hupio/main/mortality/lph_fig3.csv', col.names = c('exposure', 'ag', 'cause', 's', 'l', 'u'))
ages_uam <- read.csv('https://raw.githubusercontent.com/danlewer/hupio/main/mortality/lph_fig4.csv', col.names = c('syear', 'cause', 's', 'l', 'u', 'top', 'bottom', 'col'))
setDT(cs_smr); setDT(periods); setDT(ages); setDT(ages_uam)


year_lims4 <- 2001 + 0:5 * 3
cause_titles <- c('All cause', 'Drug-related\ndeaths', 'Cardiovascular\ndiseases', 'COPD', 'Other respiratory\ndiseases', 'Respiratory\ncancers', 'Other cancers', 'Liver\ndisease', 'External\ncauses')
causes <- c('drd', 'cvd', 'copd', 'other_resp', 'respiratory_cancer', 'other_cancers', 'liver', 'external')
causes2 <- c(list(causes), causes)
causes_with_other <- c('drd', 'external', 'cvd', 'copd', 'other_resp', 'respiratory_cancer', 'other_cancers', 'liver', 'other')
cause_titles_w_other <- c('Drug-related\ndeaths', 'External\ncauses', 'Cardiovascular\ndiseases', 'COPD', 'Other respiratory\ndiseases', 'Respiratory\ncancers', 'Other cancers', 'Liver\ndisease')

#  ::::::::
#  figure 1
#  ........

cs_smr_all <- cs_smr[gender == 'both']

cairo_pdf('lph_fig1.pdf', height = 5.5, width = 7, family = 'Corbel')

par(mar = c(3, 0, 0, 0), oma = c(0, 0, 3, 0), xpd = NA)
layout(matrix(1:2, ncol = 2), widths = c(2, 1))

plot(1, type = 'n', xlim = c(0, 5.2), ylim = c(0, 25), axes = F, xlab = NA, ylab = NA)

xs <- c(1.5, 2:5)
xs <- xs - 0.5
xs <- c(1.1, 1.6, 2.6, 3.6, 4.4)
xs <- xs + 0.2
gap <- 0.05
ys <- nrow(cs_smr_all):1

text(xs[1], ys, cs_smr_all$label, adj = 1, cex = 0.7)

text(xs[2], ys, cs_smr_all$d_f, adj = 1, cex = 0.7)
text(xs[2] + gap, ys, cs_smr_all$pc_f, adj = 0, cex = 0.7)

text(xs[3], ys, cs_smr_all$e_f, adj = 1, cex = 0.7)
text(xs[3] + gap, ys, cs_smr_all$e_pc_f, adj = 0, cex = 0.7)

text(xs[4], ys, cs_smr_all$exc_f, adj = 1, cex = 0.7)
text(xs[4] + gap, ys, cs_smr_all$exc_pc_f, adj = 0, cex = 0.7)

text(xs[5], ys, cs_smr_all$smr_f, adj = 1, cex = 0.7)
text(xs[5] + gap, ys, cs_smr_all$cs_bs, adj = 0, cex = 0.7)

mapply(text, x = `[<-`(xs, i = 2:5, xs[2:5] + gap),
       y = max(ys) + 1,
       labels = c('Underlying cause\nof death', 'Observed\ndeaths (%)', 'Expected\ndeaths (%)', 'Excess\ndeaths (%)', 'SMR\n(95% CI)'),
       adj = list(c(1, 0), c(0.5, 0), c(0.5, 0), c(0.5, 0), c(0.5, 0)),
       cex = 0.7)

plot(1, type = 'n', xlim = c(0, 6), ylim = c(0, 25), axes = F, ylab = NA, xlab = 'SMR')
points(log(cs_smr_all$smr), ys, pch = 18, cex = 0.7)
arrows(log(cs_smr_all$bs_lower), ys, log(cs_smr_all$bs_upper), ys, code = 3, angle = 90, length = 0.03)
xx <- c(1, 2, 5, 10, 25, 50, 100, 200, 400)
axis(1, log(xx), xx, pos = 0.5, labels = F)
text(log(xx), -1, xx, srt = 45, adj = 1, cex = 0.7)
segments(0, 0, 0, 24.5)
text(log(max(xx)) / 2, -3, 'SMR (95% CI)', cex = 0.7)
segments(-13, c(0.5, 3.5, 5.5, 8.5, 11.5, 12.5, 18.5, 21.5, 24.5), 6, lty = c(1, rep(3, 7), 1))

dev.off()

#  ::::::::
#  figure 2
#  ........

cols <- brewer.pal(3, 'Set1')[1:2]
cols2 <- add.alpha(cols, 0.2)

nx <- length(year_lims4)

ymaxes <- c(1750, 750, rep(200, 7))

cairo_pdf('lph_fig2.pdf', height = 7, width = 8, family = 'Corbel')

par(mar = rep(0, 4), oma = rep(4, 4), xpd = NA)
layout(matrix(c(1, 4, 7, 2, 5, 8, 3, 6, 9, 10, 10, 10), ncol = 4), widths = c(2, 2, 2, 1))
for (i in seq_along(causes2)) {
  pd <- periods[cause == c('all_cause', causes)[i]]
  plot(1, type = 'n', xlim = c(0.25, nx), ylim = c(0, ymaxes[i]), axes = F, xlab = NA, ylab = NA)
  rect(1, 0, nx, ymaxes[i], border = NA, col = 'grey95')
  box()
  with(pd[exposure == T], polygon(c(1:nx, nx:1), c(l, rev(u)), col = cols2[1], border = NA))
  with(pd[exposure == F], polygon(c(1:nx, nx:1), c(l, rev(u)), col = cols2[2], border = NA))
  lines(1:nx, pd[exposure == T, s], col = cols[1])
  lines(1:nx, pd[exposure == F, s], col = cols[2])
  text(nx/2, ymaxes[i] * 0.93, cause_titles[i], font = 2)
  axis(2, tcl = 0.5, mgp = c(0, -2, 0), las = 2)
  if (i %in% 7:9) axis(1, 1:nx, labels = F)
  if (i %in% 7:9) text(1:nx, -ymaxes[i] * 0.1, paste0(year_lims4, '-', stri_sub(c(year_lims4[-1], 2019) -1, -2)), srt = 45, adj = 1)
}
plot(1, type = 'n', xlab = NA, ylab = NA, axes = F, xlim = c(0, 10), ylim = c(0, 10))
segments(x0 = c(1, 1), y0 = c(4, 6), x1 = c(3, 3), col = rev(cols))
polygon(c(1, 3, 3, 1), c(3.5, 3.5, 4.5, 4.5), col = cols2[2], border = NA)
polygon(c(1, 3, 3, 1), c(5.5, 5.5, 6.5, 6.5), col = cols2[1], border = NA)
text(3.5, c(4, 6), c('Comparison\ngroup', 'History of using\nillicit opioids'), adj = 0)
title(ylab = 'Standardised mortality rate per 100,000 person-years', outer = T, line = 2)

dev.off()

#  :::::::::::::::::::::::::::::::::
#  figure 3 version in supplementary
#  .................................

agesw <- dcast(ages, exposure + ag ~ cause, value.var = 's')
agesw <- agesw[, c('exposure', 'ag', causes_with_other), with = F]

pv <- agesw[, causes_with_other, with = F]
pv <- cbind(bl = 0, pv)
pv <- as.matrix(pv)
pvy <- t(apply(pv, 1, cumsum))

pvy_e <- pvy[agesw$exposure == T,]
pvc_e <- pvy_e / pvy_e[, ncol(pvy_e)]

pvy_u <- pvy[agesw$exposure == F,]
pvc_u <- pvy_u / pvy_u[, ncol(pvy_u)]

x1 <- 18:65

cols <- brewer.pal(9, 'Set3')
xx <- c(18, seq(25, 65, 10))

emf('cause_by_age_comparison_11june2021.emf', height = 6, width = 8, family = 'Corbel')

par(mar = rep(0, 4), oma = c(2, 5, 5, 5), xpd = NA)
layout(matrix(c(1, 3, 6, 2, 4, 6, 5, 5, 7), ncol = 3), widths = c(2, 2, 1), heights = c(6, 6, 1.5))

plot(1, type = 'n', xlim = c(18, 65), ylim = c(0, 6000), axes = F, xlab = NA, ylab = NA)
rect(18, 0, 65, 6000)
for (i in 1:(ncol(pvy)-1)) polygon(c(x1, rev(x1)), c(pvy_e[,i], rev(pvy_e[, i+1])), col = cols[i], lwd = 0.5)
axis(2, 0:6 * 1000, prettyNum(0:6 * 1000, big.mark = ','), pos = 18, las = 2)
text(18, 6000 * 1.1, 'Standardised mortality rate\nper 100,000 person-years', adj = 0)
text(20, 6000 * 0.9, 'History of using\nillicit opioids', adj = 0)

plot(1, type = 'n', xlim = c(18, 65), ylim = c(0, 1), axes = F, xlab = NA, ylab = NA)
rect(18, 0, 65, 1)
for (i in 1:(ncol(pvy)-1)) polygon(c(x1, rev(x1)), c(pvc_e[,i], rev(pvc_e[, i+1])), col = cols[i], lwd = 0.5)
axis(4, 0:5/5, paste0(0:5 * 20, '%'), pos = 65, las = 2)
text(65, 1 * 1.1, 'Proportion of\ndeaths', adj = 1)

plot(1, type = 'n', xlim = c(18, 65), ylim = c(0, 1200), axes = F, xlab = NA, ylab = NA)
rect(18, 0, 65, 1200)
for (i in 1:(ncol(pvy)-1)) polygon(c(x1, rev(x1)), c(pvy_u[,i], rev(pvy_u[, i+1])), col = cols[i], lwd = 0.5)
axis(2, 0:6 * 200, prettyNum(0:6 * 200, big.mark = ','), pos = 18, las = 2)
axis(1, xx, pos = 0, labels = F)
text(xx, 1200*-0.1, xx, srt = 45, adj = 1)
text(20, 1200 * 0.9, 'Matched comparison\ngroup', adj = 0)

plot(1, type = 'n', xlim = c(18, 65), ylim = c(0, 1), axes = F, xlab = NA, ylab = NA)
rect(18, 0, 65, 1)
for (i in 1:(ncol(pvy)-1)) polygon(c(x1, rev(x1)), c(pvc_u[,i], rev(pvc_u[, i+1])), col = cols[i], lwd = 0.5)
axis(4, 0:5/5, paste0(0:5 * 20, '%'), pos = 65, las = 2)
axis(1, xx, pos = 0, labels = F)
text(xx, -0.1, xx, srt = 45, adj = 1)

par(mar = c(0, 3, 0, 0))
plot(1, type = 'n', xlab = NA, ylab = NA, axes = F, xlim = c(0, 10), ylim = c(0, 10))
ys <- seq(1, 9, length.out = ncol(pvy))
rect(3, ys[-1], 5, ys[-length(ys)], col = cols)
text(6, ys[-length(ys)] + diff(ys)/2, c(cause_titles_w_other, 'Other'), adj = 0)
text(3, max(ys) * 1.05, 'Underlying cause\nof death', adj = 0)

par(mar = rep(0, 4))
plot(1, type = 'n', xlab = NA, ylab = NA, axes = F, xlim = c(0, 10), ylim = c(0, 10))
text(5, 0, 'Age')

dev.off()

#  ::::::::
#  figure 3
#  ........

cairo_pdf('lph_fig3.pdf', height = 4, width = 8, family = 'Corbel')

par(mar = rep(0, 4), oma = c(2, 5, 5, 5), xpd = NA)
layout(matrix(c(1, 4, 2, 4, 3, 5), ncol = 3), widths = c(2, 2, 1), heights = c(9, 1))

plot(1, type = 'n', xlim = c(18, 65), ylim = c(0, 6000), axes = F, xlab = NA, ylab = NA)
for (i in 1:(ncol(pvy)-1)) polygon(c(x1, rev(x1)), c(pvy_e[,i], rev(pvy_e[, i+1])), col = cols[i], lwd = 0.5)
rect(18, 0, 65, 6000)
axis(2, 0:6 * 1000, prettyNum(0:6 * 1000, big.mark = ','), pos = 18, las = 2)
text(41.5, 6000 * 1.15, 'Standardised mortality rate\nper 100,000 person-years', adj = c(0.5, 1))
axis(1, xx, pos = 0, labels = F)
text(xx, 6000*-0.05, xx, srt = 45, adj = 1)

plot(1, type = 'n', xlim = c(18, 65), ylim = c(0, 1), axes = F, xlab = NA, ylab = NA)
for (i in 1:(ncol(pvy)-1)) polygon(c(x1, rev(x1)), c(pvc_e[,i], rev(pvc_e[, i+1])), col = cols[i], lwd = 0.5)
rect(18, 0, 65, 1)
axis(4, 0:5/5, paste0(0:5 * 20, '%'), pos = 65, las = 2)
text(41.5, 1 * 1.15, 'Proportion of\ndeaths', adj = c(0.5, 1))
axis(1, xx, pos = 0, labels = F)
text(xx, -0.05, xx, srt = 45, adj = 1)

par(mar = c(0, 3, 0, 0))
plot(1, type = 'n', xlab = NA, ylab = NA, axes = F, xlim = c(0, 10), ylim = c(0, 10))
ys <- seq(0, 10, length.out = ncol(pvy))
rect(3, ys[-1], 5, ys[-length(ys)], col = cols)
text(6, ys[-length(ys)] + diff(ys)/2, c(cause_titles_w_other, 'Other'), adj = 0)
text(3, max(ys) * 1.15, 'Underlying cause\nof death', adj = c(0, 1))

par(mar = rep(0, 4))
plot(1, type = 'n', xlab = NA, ylab = NA, axes = F, xlim = c(0, 10), ylim = c(0, 10))
text(5, 0, 'Age')

dev.off()

#  ::::::::
#  figure 4
#  ........

cols <- brewer.pal(9, 'Set3')

cairo_pdf('lph_fig4.pdf', height = 7, width = 9, family = 'Corbel')

par(mar = c(5, 5, 0, 10), xpd = NA)

plot(1, type = 'n', xlim = c(2001, 2019), ylim = c(0, 1300), axes = F, xlab = NA, ylab = NA)
rect(18, 0, 65, 6000)
with(ages_uam[cause != 'all_cause'], rect(syear, bottom, syear + 1, top, col = col))
axis(2, 0:5 * 250, prettyNum(0:5 * 250, big.mark = ','), pos = 2001, las = 2)
title(ylab = 'Standardised mortality rate\nper 100,000 person-years', line = 3)
title(xlab = 'Year', line = 2)

axis(1, 2001:2019, pos = 0, labels = F)
xx <- seq(2001, 2018, 3)
text(xx + 0.5, -35, xx, srt = 45, adj = 1)

ys <- seq(0, 1200, length.out = 10)
rect(2020, ys[-1], 2021, ys[-length(ys)], col = cols)
text(2021.5, ys[-length(ys)] + diff(ys)/2, c(cause_titles_w_other, 'Other'), adj = 0)
text(2020, max(ys) * 1.07, 'Underlying cause\nof death', adj = 0)

dev.off()

#  :::::::
#  table 2
#  .......

table2 <- rbind(ages[exposure == T & cause == 'all_cause' & ag %in% seq(20, 55, 5), .(cause = 'all_cause', s = sum(s), l = sum(l), u = sum(u)), ag],
                ages[exposure == T & cause %in% c('liver', 'cvd', 'copd', 'other_resp', 'respiratory_cancer', 'other_cancers') & ag %in% seq(20, 55, 5), .(cause = 'ncd', s = sum(s), l = sum(l), u = sum(u)), ag],
                ages[exposure == T & cause == 'drd' & ag %in% seq(20, 55, 5), .(cause = 'drd', s = sum(s), l = sum(l), u = sum(u)), ag])
table2 <- table2[cause == 'all_cause', .(ag = ag, ac = s)][table2, on = 'ag']
table2[, pc := paste0(format(round(s / ac * 100, 1), digits = 1, nsmall = 1), ' %')]
table2[, c('s', 'l', 'u') := lapply(.SD, round, digits = 0), .SDcols = c('s', 'l', 'u')]
table2[, c('s', 'l', 'u') := lapply(.SD, function(x) pmax(0, x)), .SDcols = c('s', 'l', 'u')]
table2[, c('s', 'l', 'u') := lapply(.SD, formatC, big.mark = ','), .SDcols = c('s', 'l', 'u')]
table2[, rate := paste0(s, '(', l, '-', u, ')')]
table2[, rate := gsub('\\(', ' (', gsub(' ', '', rate))]
table2 <- dcast(table2, ag ~ cause, value.var = 'rate')[dcast(table2, ag ~ cause, value.var = 'pc'), on = 'ag']
table2 <- table2[, c('ag', 'all_cause', 'drd', 'i.drd', 'ncd', 'i.ncd')]
setnames(table2, c('i.drd', 'i.ncd'), c('pc_drd', 'pc_ncd'))
fwrite(table2, 'table2_4dec2021.csv')

#  ::::::::::::::::::::::::
#  tables for supplementary
#  ........................

# SMRs stratified by sex

fut <- structure(list(gender = c("male", "female", "both"), fu = c(668834.056125941, 291877.305954825, 960711.362080766)), row.names = c(NA, -3L), class = c("data.table", "data.frame"))
supp_t6 <- copy(cs_smr)
supp_t6 <- fut[supp_t6, on = 'gender']
supp_t6[, crude_rate := vf(vpt(d, fu / 100000, digs = 1), digs = 1)]
supp_t6[, sex := factor(gender, c('both', 'female', 'male'), c('All', 'Female', 'Male'))]
supp_t6[, cause := factor(cause, fread('label_lookup.csv')$cause)]
supp_t6 <- supp_t6[order(cause, sex)]
supp_t6[, smr_f2 := paste0(smr_f, ' ', cs_bs)]
supp_t6 <- supp_t6[, c('label', 'sex', 'd_f', 'crude_rate', 'e_f', 'smr_f2')]
supp_t6[, label := fifelse(label == shift(label), '', as.character(label))]
supp_t6$label[1] <- 'All cause'
fwrite(supp_t6, 'supplementary_table_smr_by_sex_4dec2021.csv')

# fig2

fig2_table <- copy(periods)
fig2_table[, ncd := cause %in% c('liver', 'other_cancers', 'respiratory_cancer', 'other_resp', 'copd', 'cvd')]
fig2_table <- rbind(fig2_table, fig2_table[ncd == T, .(s = sum(s), l = sum(l), u = sum(u), cause = 'ncd'), c('cy4', 'exposure')], fill = T)
fig2_table[, cause := factor(cause, c('all_cause', causes_with_other, 'ncd'))]
fig2_table[, exposure := factor(exposure, c('TRUE', 'FALSE'))]
fig2_table[, c('s', 'l', 'u') := lapply(.SD, round, digits = 0), .SDcols = c('s', 'l', 'u')]
fig2_table[, c('s', 'l', 'u') := lapply(.SD, function(x) pmax(0, x)), .SDcols = c('s', 'l', 'u')]
fig2_table[, c('s', 'l', 'u') := lapply(.SD, formatC, big.mark = ','), .SDcols = c('s', 'l', 'u')]
fig2_table[, rate := paste0(s, '(', l, '-', u, ')')]
fig2_table[, rate := gsub('\\(', ' (', gsub(' ', '', rate))]
fig2_table <- dcast(fig2_table, exposure + cy4 ~ cause, value.var = 'rate')
fwrite(fig2_table, 'fig2_table_4dec2021.csv')

# fig3

fig3_table <- copy(ages)
fig3_table[, ncd := cause %in% c('liver', 'other_cancers', 'respiratory_cancer', 'other_resp', 'copd', 'cvd')]
fig3_table <- rbind(fig3_table, fig3_table[ncd == T, .(s = sum(s), l = sum(l), u = sum(u), cause = 'ncd'), c('exposure', 'ag')], fill = T)
fig3_table[, cause := factor(cause, c('all_cause', causes_with_other, 'ncd'))]
fig3_table[, exposure := factor(exposure, c('TRUE', 'FALSE'))]
fig3_table[, c('s', 'l', 'u') := lapply(.SD, round, digits = 0), .SDcols = c('s', 'l', 'u')]
fig3_table[, c('s', 'l', 'u') := lapply(.SD, function(x) pmax(0, x)), .SDcols = c('s', 'l', 'u')]
fig3_table[, c('s', 'l', 'u') := lapply(.SD, formatC, big.mark = ','), .SDcols = c('s', 'l', 'u')]
fig3_table[, rate := paste0(s, '(', l, '-', u, ')')]
fig3_table[, rate := gsub('\\(', ' (', gsub(' ', '', rate))]
fig3_table <- dcast(fig3_table, exposure + ag ~ cause, value.var = 'rate')
fwrite(fig3_table, 'fig3_table_4dec2021.csv')

# fig4

fig4_table <- copy(ages_uam)
fig4_table[, ncd := cause %in% c('liver', 'other_cancers', 'respiratory_cancer', 'other_resp', 'copd', 'cvd')]
fig4_table <- rbind(fig4_table, fig4_table[ncd == T, .(s = sum(s), l = sum(l), u = sum(u), cause = 'ncd'), 'syear'], fill = T)
fig4_table[, cause := factor(cause, c('all_cause', causes_with_other, 'ncd'))]
fig4_table[, c('s', 'l', 'u') := lapply(.SD, round, digits = 0), .SDcols = c('s', 'l', 'u')]
fig4_table[, c('s', 'l', 'u') := lapply(.SD, function(x) pmax(0, x)), .SDcols = c('s', 'l', 'u')]
fig4_table[, c('s', 'l', 'u') := lapply(.SD, formatC, big.mark = ','), .SDcols = c('s', 'l', 'u')]
fig4_table[, rate := paste0(s, '(', l, '-', u, ')')]
fig4_table[, rate := gsub('\\(', ' (', gsub(' ', '', rate))]
fig4_table <- dcast(fig4_table, syear ~ cause, value.var = 'rate')
fwrite(fig4_table, 'fig4_table_4dec2021.csv')
