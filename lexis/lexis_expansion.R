library(data.table)

#  lexis expansion function
#  ........................

lexis <- function(data, id, start_date, end_date, time_varying, dob, event, age_groups = 0:20 * 5, ng = 50) {
  nr <- seq_len(nrow(data))
  ind <- findInterval(nr, quantile(nr, 0:(ng-1)/ng))
  x <- split(data, ind)
  l <- lapply(seq_along(x), function(y) {
    t0 <- proc.time()
    days <- x[[y]][, .(day = seq.int(from = get(start_date), to = get(end_date))), c(id, dob, event, time_varying)]
    names(days) <- c('id', 'dob', 'event_date', 'time_varying', 'day')
    days[, age := (day - dob) / 365.25]
    days[, age_group := findInterval(age, age_groups)]
    days[, tv := !is.na(time_varying) & day >= time_varying]
    days[, event := event_date == day & !is.na(event_date)]
    print(paste0(y, '/', ng))
    days[, .(follow_up = .N, event = any(event)), c('id', 'age_group', 'tv')]
  })
  rbindlist(l)
}

#  generate example data
#  .....................

study_start <- as.integer(as.Date('1980-01-01'))
study_end <- as.integer(as.Date('2020-12-31'))
earliest_dob <- as.integer(as.Date('1950-01-01'))
latest_dob <- as.integer(as.Date('1970-12-31'))
study_duration <- study_end - study_start
event_probability <- 0.1
time_varying_probability <- 0.3
n <- 1e5

set.seed(33)
dat <- data.table(id = seq_len(n), 
                  entry = sample(study_start:study_end, n, replace = T),
                  follow_up = sample(seq_len(study_duration), n, replace = T),
                  event = rbinom(n, 1, event_probability),
                  dob = sample(earliest_dob:latest_dob, n, replace = T),
                  time_varying = rbinom(n, 1, time_varying_probability))
dat[, exit := entry + follow_up]
dat[, exit := pmin(exit, study_end)]
dat[, follow_up := exit - entry]
dat[, event_date := (entry + floor(runif(n) * follow_up)) * event]
dat$event_date[dat$event_date == 0] <- NA_integer_
dat[, exit := pmin(exit, event_date, na.rm = T)]
dat[, follow_up := exit - entry]
dat[, time_varying := (entry + floor(runif(n) * follow_up)) * time_varying]
dat$time_varying[dat$time_varying == 0] <- NA_integer_

#  create lexis-expanded data
#  ..........................

lexis_data <- lexis(dat, id = 'id', start_date = 'entry', end_date = 'exit', event = 'event_date', time_varying = 'time_varying', dob = 'dob')
