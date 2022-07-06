library(data.table)
library(devEMF)
library(stringi)
library(RColorBrewer)

d <- read.csv(url('https://raw.githubusercontent.com/danlewer/hupio/main/copd/copd_treatment_adverse_outcomes_30sept2021.csv'))
setDT(d)
d[, label := factor(outcome, c('dod_resp', 'dod_all', 'resp_hosp', 'ae_hosp', 'sc', 'meds', 'pulm_rehab', 'pnu', 'flu'),
                    c('Death with underlying respiratory cause',
                      'All-cause death',
                      'Emergency admission for respiratory disease',
                      'Hospitalised exacerbation of COPD',
                      'Support with smoking cessation',
                      'Bronchodilators or corticosteroids',
                      'Pulmonary rehabilitation',
                      'Pneumococcal vaccine',
                      'Seasonal influenza vaccine'))]
d <- d[order(label, sens)]

extract_num <- function (x) {
  pe <- stri_sub(x, to = stri_locate_first_fixed(x, " ")[,1] - 1)
  lower <- stri_sub(x, from = stri_locate_first_fixed(x, "(")[,1]+1, to = stri_locate_first_fixed(x, "-")[,1]-1)
  upper <- stri_sub(x, from = stri_locate_first_fixed(x, "-")[,1]+1, to = stri_locate_first_fixed(x, ")")[,1]-1)
  r <- cbind(pe, lower, upper)
  `class<-`(r, 'numeric')
}

unadj <- extract_num(d$unadj)
adj <- extract_num(d$adj)

unadj <- cbind(d[, -'adj'], unadj)
adj <- cbind(d[, -'unadj'], adj)
names(unadj)[4] <- 'est'
names(adj)[4] <- 'est'

ygap <- 0.5
ys <- c(1:5, 6:12 + ygap)
ys2 <- ys[c(1:4, 6:10)]
ydisp <- 0.1

xs <- c(-5.5, -2.25, -1.25)

pf <- function (pd1, pd2, heads = c('Unadjusted HR\n(95% CI)', 'Adjusted HR\n(95% CI)')) {
  
  cols <- brewer.pal(3, 'Set1')[1:2]
  
  rr <- c(0.5, 0.75, 1, 1.5, 2, 3, 5)
  
  par(xpd = NA)
  
  plot(1, type = 'n', xlim = c(min(xs), log(max(rr))), ylim = c(0, max(ys)), axes = F, xlab = NA, ylab = NA)
  
  axis(1, log(rr), labels = rr, pos = 0)
  segments(0, 0, y1 = 12)
  
  with(pd1, {
    points(log(pe), ys2 + ydisp, pch = 19, col = cols[1])
    arrows(log(lower), ys2 + ydisp, x1 = log(upper), code = 3, angle = 90, length = 0.03, col = cols[1])
  })
  with(pd2, {
    points(log(pe), ys2 - ydisp, pch = 19, col = cols[2])
    arrows(log(lower), ys2 - ydisp, x1 = log(upper), code = 3, angle = 90, length = 0.03, col = cols[2])
  })
  
  text(xs[1], ys2, pd1$label, adj = 0)
  text(xs[1], ys[c(5, 11)], c('Adverse events', 'Evidence-based treatment'), adj = 0, font = 2)
  text(xs[2], ys2, pd1$est, col = cols[1])
  text(xs[3], ys2, pd2$est, col = cols[2])
  text(xs[2:3], max(ys), heads, col = cols)
  
  text(0, -2, 'Hazard ratio')
  
}


cairo_pdf('Fig2_v2.pdf', height = 6, width = 11, family = 'Franklin Gothic Book')

par(mar = c(4, 0, 4, 0))

pf(pd1 = unadj[sens == 'no'], pd2 = adj[sens == 'no'])
#text(-5.5, 14.5, 'Figure 2: Hazard ratios of evidence-based treatment and adverse events after diagnosis of COPD, comparing participants with a\nhistory of illicit opioid use to those without', adj = 0)
#abline(h = 13.5)

dev.off()

emf('sensitivity_current_smokers.emf', height = 6, width = 11, family = 'Franklin Gothic Book')

pf(pd1 = adj[sens == 'no'], pd2 = adj[sens == 'yes'], heads = c('Main\nanalysis', 'Restricted to\ncurrent smokers'))

dev.off()
