library(data.table)
library(RColorBrewer)
library(devEMF)

# read deaths data
# source: https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/deathsrelatedtodrugpoisoningenglandandwalesreferencetable

h <- `names<-`(read.csv(url("https://raw.githubusercontent.com/danlewer/hupio/main/general_pop_rates/ons_heroin_morphine_deaths.csv")), c('year', 'Under 20', '20-29', '30-39', '40-49', '50-69', '70 and over'))
setDT(h)
h <- melt(h, id.vars = 'year', variable.name = 'age', value.name = 'd')
year_lims <- c(2001 + 0:5 * 3)
h[, cy4 := findInterval(year, year_lims)]
h <- h[year %in% 2001:2019, .(deaths = sum(d)), c('cy4', 'age')]

# read population data
# source: nomisweb

p <- read.csv(url("https://raw.githubusercontent.com/danlewer/hupio/main/general_pop_rates/nomis_population_estimates.csv"))
names(p) <- stri_replace_all_fixed(stri_replace_all_fixed(names(p), '...', '-'), '.', ' ')
setDT(p)
p <- melt(p, id.vars = c('year', 'sex'), variable.name = 'age', value.name = 'pop')
p[, age := gsub('Aged|year|and| |s', '', age)]
p <- p[!(age %in% c('under1', '1-4', '5-9', '10-14', '80-84', '85over')) & year %in% 2001:2018]
p[, cy4 := findInterval(year, year_lims)]
p <- p[, .(pop = sum(pop)), c('sex', 'age', 'cy4')]

p[, age2 := factor(age, 
                   c('15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79'),
                   c('Under 20', '20-29', '20-29', '30-39', '30-39', '40-49', '40-49', '50-69', '50-69', '50-69', '50-69', '70 and over', '70 and over'))]
p2 <- p[, .(pop = sum(pop)), c('cy4', 'age2')]
setnames(p2, 'age2', 'age')

# combine and standardise

h <- p2[h, on = c('cy4', 'age')]
std <- p2[, .(std = round(sum(pop)/10, 0)), age]

h <- std[h, on = 'age']
h[, r := deaths / pop]
h[, e := r * std]

# standardised rate per million

srpm <- h[, .(so = sum(e) / sum(std) * 1000000), cy4][order(cy4)]

# add results from hupio study

hupio <- data.table(cy4 = 1:6,
                    hupio_drd = c(479.7, 532.6, 484.1, 391.1, 550.2, 599.7))

srpm <- srpm[hupio, on = 'cy4']
ymax <- 40
ymax2 <- 800
srpm[, hupio_rescale := hupio_drd * ymax/ymax2]

cols <- brewer.pal(3, 'Set1')[1:2]

emf('drd_hupio_vs_gen_pop_6june2021.emf', width = 6, height = 4, family = 'Corbel')

par(mar = c(5, 6, 0, 6), xpd = NA)
plot(1, type = 'n', xlim = c(0.5, 6.5), ylim = c(0, 40), axes = F, xlab = NA, ylab = NA)
with(srpm, lines(cy4, so, col = cols[1]))
with(srpm, lines(cy4, hupio_rescale, col = cols[2]))
with(srpm, points(cy4, so, col = cols[1], pch = 19))
with(srpm, points(cy4, hupio_rescale, col = cols[2], pch = 19))
axis(2, pos = 0.8, las = 2, col = cols[1], col.axis = cols[1])
ys <- 0:4 * 200
axis(4, ys * ymax / ymax2, ys, pos = 6.2, las = 2, col = cols[2], col.axis = cols[2])
xlabs <- paste0(year_lims, '-', stri_sub(year_lims+3, -2))
axis(1, 1:6, labels = F, pos = -2)
text(1:6, -5, paste0(year_lims, '-', stri_sub(year_lims+3, -2)), srt = 45, adj = 1)
text(-1, 20, 'Standardised rate of heroin or\nmorphine-related death per million\n(general population)', srt = 90, col = cols[1])
text(8, 20, 'Standardised rate of drug-\nrelated death per 100,000\nin study', srt = 270, col = cols[2])

dev.off()
