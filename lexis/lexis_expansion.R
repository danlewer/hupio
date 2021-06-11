library(data.table)
library(lubridate)

#  lexis expansion function
#  ........................

lexis <- function(data, id = 'id', start_date = 'entry', end_date = 'exit', baseline = 'exposure_baseline', opioids = 'becomes_exposed', dob = 'dob', event = 'death_date', age_groups = c(0, 18, seq(25, 100, 5)), period_starts = c(10227, 11323, 12418, 13514, 14610, 15706, 16801), ng = 50) {
  nr <- seq_len(nrow(data))
  ind <- findInterval(nr, quantile(nr, 0:(ng-1)/ng))
  x <- split(data, ind)
  l <- lapply(seq_along(x), function(y) {
    days <- x[[y]][, .(day = seq.int(from = get(start_date), to = get(end_date))), c(id, dob, event, opioids, baseline)]
    names(days) <- c('id', 'dob', 'death_date', 'opioids', 'exposure', 'day')
    days[, age := (day - dob) / 365.25]
    days[, age_group := findInterval(age, age_groups)]
    days[, opioids := exposure == 1 | (!is.na(opioids) & day >= opioids)]
    days[, event := death_date == day & !is.na(death_date)]
    days[, period := findInterval(day, period_starts)]
    print(paste0(y, '/', ng))
    days[, .(follow_up = .N, event = any(event)), c('id', 'age_group', 'period', 'opioids')]
  })
  rbindlist(l)
}

#  generate example data
#  .....................

# inputs

study_start <- as.integer(as.Date('1998-01-01'))
study_end <- as.integer(as.Date('2018-10-30'))
earliest_dob <- as.integer(as.Date('1960-01-01'))
latest_dob <- as.integer(as.Date('1995-12-31'))
youngest_entrant <- as.integer(19 * 365.25)
study_duration <- study_end - study_start
opioid_mortality_ratio <- 10
proportion_comparison_become_exposed <- 0.03
number_opioids <- 1e4
ratio <- 3

# make data

set.seed(31)
n <- number_opioids * (ratio + 1)
dat <- data.table(id = seq_len(n),
                  exposure_baseline = rep(1:0, c(number_opioids, number_opioids * ratio)),
                  entry = sample(study_start:study_end, n, replace = T),
                  follow_up = sample(seq_len(study_duration), n, replace = T),
                  dob = sample(earliest_dob:latest_dob, n, replace = T),
                  becomes_exposed = rbinom(n, 1, proportion_comparison_become_exposed))
dat[, entry := pmax(entry, dob + youngest_entrant)]
dat[, exit := pmin(entry + follow_up, study_end)]
dat[, follow_up := exit - entry]
dat[, age_entry := (entry - dob)/365.25]
dat[, death_risk := (age_entry^3)/5e6]
dat[, death_risk := death_risk * c(1, opioid_mortality_ratio)[exposure_baseline + 1]]
dat[, death := rbinom(n, 1, death_risk)]
dat[, death_date := (entry + as.integer(runif(n) * follow_up)) * death]
dat$death_date[dat$death_date == 0] <- NA_integer_
dat[, exit := pmin(exit, death_date, na.rm = T)]
dat[, follow_up := exit - entry]
dat[, becomes_exposed := (entry + as.integer(runif(n) * follow_up)) * becomes_exposed]
dat$becomes_exposed[dat$becomes_exposed == 0 | dat$exposure_baseline == 1] <- NA_integer_

#  create lexis-expanded data
#  ..........................

age_lims <- c(0, 18, seq(25, 100, 5))
period_starts <- seq(1998, 2018, 3)
period_starts <- make_date(year = period_starts, 1, 1)

lexis_data <- lexis(dat, id = 'id', start_date = 'entry', end_date = 'exit', baseline = 'exposure_baseline', opioids = 'becomes_exposed', dob = 'dob', event = 'death_date', age_groups = age_lims, period_starts = as.integer(period_starts))
lexis_data[, age_group := factor(age_group, seq_along(age_lims), age_lims)]
lexis_data[, period := factor(period, seq_along(period_starts), year(period_starts))]
lexis_data <- droplevels(lexis_data)

# check total follow-up duration (should be TRUE)

dat[, sum(exit - entry + 1)] == lexis_data[, sum(follow_up)]

#  smr calculation
#  ...............

smr <- lexis_data[, .(person_years = sum(follow_up)/365.25, deaths = sum(event)), c('age_group', 'opioids')]
ref <- smr[opioids == F]
ref[, ref_rate := deaths / person_years]
smr <- ref[, c('age_group', 'ref_rate')][smr[opioids == T, -'opioids'], on = 'age_group']
smr[, expected := ref_rate * person_years]
smr[, .(observed = sum(deaths), expected = sum(expected), smr = sum(deaths) / sum(expected))]

#  poisson model
#  .............

model <- glm(event ~ opioids + age_group + period + offset(log(follow_up)), data = lexis_data, family = 'poisson')
exp(cbind(coef(model), confint(model)))[-1,]

