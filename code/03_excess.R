# EDSD 2025
# Analysis of mortality disturbances course
# Instructor: Enrique Acosta (CED)
# Lab 3: Very basic introduction to excess mortality estimation

rm(list=ls())
source("code/00_setup.R")

# 3 baselines: 3 colors
cols <- brewer.pal(3, "Dark2")

# loading mortality data ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
dt <- read_rds("data_input/canada_weekly_deaths_pop_2015_2023.rds")

dt %>% 
  ggplot()+
  geom_line(aes(date, dts), linewidth = 1)+
  theme_bw()

# ######################################
# estimating the baseline mortality ====
# ######################################

# ~~~~~~~~~~~~~~~~~
# average-week ====
# ~~~~~~~~~~~~~~~~~
# the simplest approach (maybe the worst)

# estimating the average weekly deaths for the whole training period
# as the baseline
w_av <- 
  dt %>% 
  filter(year <= 2019) %>%
  summarise(bsn = mean(dts, na.rm = TRUE)) %>% 
  pull(bsn)

dt_w_av <- 
  dt %>% 
  mutate(bsn = w_av)

dt_w_av %>% 
  ggplot()+
  geom_line(aes(date, dts), linewidth = 1)+
  geom_line(aes(date, bsn), linewidth = 1, col = cols[1])+
  geom_vline(xintercept = ymd("2020-03-15"), linetype = "dashed")+
  theme_bw()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Week-specific average ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
ws_av <- 
  dt %>% 
  filter(year <= 2019) %>%
  group_by(week) %>% 
  summarise(bsn = mean(dts, na.rm = TRUE)) %>% 
  ungroup()

dt_ws_av <- 
  dt %>% 
  left_join(ws_av)

dt_ws_av %>% 
  ggplot()+
  geom_line(aes(date, dts), linewidth = 1)+
  geom_line(aes(date, bsn), linewidth = 1, col = cols[2])+
  geom_vline(xintercept = ymd("2020-03-15"), linetype = "dashed")+
  theme_bw()


# it seems not bad at all!!

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Excess mortality using a Poisson model ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# let's adjust a bit the data to take advantage of information on the 
# population exposed to risk 
dt2 <- 
  dt %>% # merging with weekly population
  mutate(exposure = pop / 52, # exposure in person-weeks
         t = 1:n(), # a variable for the secular trend 
         w = ifelse(date <= "2020-03-15", 1, 0)) # a variable for the weights 

# Fitting a GAM model to estimate expected mortality ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# GAM models allow us to include parametric and semi-parametric terms together 
# the "mgcv" package is very convenient for working with GAM models
# https://cran.r-project.org/web/packages/mgcv/mgcv.pdf

# fitting the model
gam_model <- 
  gam(dts ~
        # linear term (exponential outside Poisson) for the secular trend
        t +
        # cyclic spline term for the seasonal trend
        s(week, bs = 'cp') + 
        # controlling for population changes over time
        offset(log(exposure)), 
      # to avoid including the pandemic period in the baseline estimation
      weights = w, 
      data = dt2, 
      # using a quasipoisson distribution to account for overdispersion
      family = "quasipoisson") 

# how does it look?
summary(gam_model)

# weekly slope (t) = 0.00012142
(exp(0.00012657)-1)*100
# 0.013% increase every week

# example for predicting estimates
bsn_poi <- predict(gam_model, newdata = dt2, type = "response", se.fit = TRUE)

# obtaining estimates for the three models
bsn <- 
  dt2 %>% 
  mutate(bsn = bsn_poi$fit,
         se = bsn_poi$se.fit,
         ll = bsn - 1.96*se,
         ul = bsn + 1.96*se)

bsn %>% 
  ggplot()+
  geom_ribbon(aes(date, ymin = ll, ymax = ul), alpha = 0.5, fill = cols[3])+
  geom_line(aes(date, dts), linewidth = 1)+
  geom_line(aes(date, bsn), linewidth = 1, col = cols[3])+
  theme_bw()

# excess estimation
exc <- 
  bsn %>% 
  filter(date >= "2020-03-15",
         date <= "2023-12-31") %>% 
  mutate(exc = dts - bsn) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# comparing excess estimation approaches ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# binding together the three estimates

bsns <- 
  bind_rows(dt_w_av %>% 
              select(year, week, date, dts, bsn) %>% 
              mutate(type = "weekly_average"),
            dt_ws_av %>% 
              select(year, week, date, dts, bsn) %>% 
              mutate(type = "weekly_spc_average"),
            bsn %>% 
              select(year, week, date, dts, bsn) %>% 
              mutate(type = "poisson_model"))
unique(bsns$type)

# plotting the three baselines
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p_bsns <- 
  bsns %>% 
  ggplot()+
  geom_line(aes(date, dts), linewidth = 1)+
  geom_line(aes(date, bsn, col = type), linewidth = 1)+
  scale_color_manual(values = cols)+
  theme_bw()
p_bsns

# estimating excess = observed - counterfactual scenario
excs <- 
  bsns %>% 
  mutate(exc = dts - bsn,
         psc = dts / bsn) %>% 
  filter(date >= "2020-03-15",
         date <= "2023-12-31")

# visualizing weekly excess deaths
excs %>% 
  ggplot()+
  geom_line(aes(date, exc, col = type), linewidth = 1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_color_manual(values = cols)+
  theme_bw()

# visualizing cumulative excess
excs %>% 
  drop_na(exc) %>% 
  group_by(type) %>% 
  mutate(exc_cum = cumsum(exc)) %>% 
  ggplot()+
  # geom_line(aes(date, exc))+
  geom_line(aes(date, exc_cum, col = type), linewidth = 1)+
  theme_bw()


# obtaining annual excess
yr_exc <- 
  excs %>% 
  filter(date >= "2020-03-15" & date <= "2023-12-31") %>% 
  group_by(year, type) %>% 
  summarise(exc = sum(exc, na.rm = TRUE)) %>% 
  ungroup()

yr_exc %>% 
  spread(year, exc)

tots <- 
  yr_exc %>% 
  spread(type, exc) %>% 
  summarise(poisson_model = sum(poisson_model),
            weekly_average = sum(weekly_average),
            weekly_spc_average = sum(weekly_spc_average),
            year = "total")

tots

t1 <- round(tots$poisson_model)
t2 <- round(tots$weekly_average)
t3 <- round(tots$weekly_spc_average)

d2 <- paste0(round((t2 / t1 - 1)*100, 1), "%")
d3 <- paste0(round((t3 / t1 - 1)*100, 1), "%")

# plotting excess estimates by year
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
yr_exc %>% 
  ggplot(aes(type, exc, fill = factor(year))) + 
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(x = "poisson_model", y = 0), label = paste0(t1, "\nreference"), vjust = -1)+
  geom_text(aes(x = "weekly_average", y = 0), label = paste0(t2, "\n", d2), vjust = -1)+
  geom_text(aes(x = "weekly_spc_average", y = 0), label = paste0(t3, "\n", d3), vjust = -1)+
  labs(fill = "year")+
  coord_cartesian(expand = 0)+
  theme_bw()


# Looking at the four years 2020-2023, estimates using averages overestimate in 
# 84% those obtained from the Poisson model!!!

# ~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~

# a function for analyzing different countries:

# let's compare Canada with the estimates we already obtained 
excess <-
  obtain_excess(cd = "CAN", sx = "b", ag = "TOT", ymin = 2015)

excess[[1]]
excess[[2]]
excess[[3]] %>% spread(type, exc)

# Now, let's look at other countries 
excess <-
  obtain_excess(cd = "NZL", sx = "b", ag = "TOT", ymin = 2015)

excess[[1]]
excess[[2]]
excess[[3]] %>% spread(type, exc)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Now, let's look at disturbances beyond mortality

# monthly births in Brazil, Colombia, and Mexico
bts <- 
  read_csv("data_input/births_bra_col_mx_2015_2022.csv") 

bts %>% 
  ggplot(aes(date, bts))+
  geom_point()+
  geom_line()+
  facet_grid(country~.,
             scales = "free_y")+
  theme_bw()

bra <- 
  bts %>% 
  filter(country == "BRA")

# plotting monthly births
bra %>% 
  ggplot(aes(date, bts))+
  geom_point()+
  geom_line()+
  theme_bw()

# excluding covid period (> March 2020)
covid <- seq(ymd('2020-03-15'),ymd('2021-12-15'), by = '1 month')

bra2 <- 
  bra %>% 
  mutate(per = case_when(date %in% covid ~ "covid",
                         TRUE ~ "typical"),
         # weights with value 0 during COVID
         w = ifelse(date %in% covid, 0, 1))

# fitting a GAM model with loglinear and cyclical terms
md <- 
  gam(bts ~ t + s(mth, bs = 'cp'), 
      weights = w,
      data = bra2,
      family = "quasipoisson")

# predicting the model
p <- predict(md, newdata = bra2, type = "response", se.fit = TRUE)

# obtaining the baseline and confidence intervals
bra3 <- 
  bra2 %>% 
  mutate(bsn = p$fit,
         ul = p$fit + (2 * p$se.fit),
         ll = p$fit - (2 * p$se.fit))

# plotting
bra3 %>% 
  ggplot()+
  geom_ribbon(aes(date, ymin = ll, ymax = ul), fill = "red", alpha = 0.3)+
  geom_point(aes(date, bts))+
  geom_line(aes(date, bts))+
  geom_line(aes(date, bsn), col = "red")+
  theme_bw()


# excluding periods affected by Zika
# National emergency Nov 2015 - Nov 2016
# seems like 7-month lag for the effect
zica <- seq(ymd('2016-05-15'),ymd('2017-05-15'), by = '1 month')

bra2 <- 
  bra %>% 
  mutate(per = case_when(date %in% zica ~ "zica",
                         date %in% covid ~ "covid",
                         TRUE ~ "typical"),
         # weights with value 0 during Zika and COVID
         w = ifelse(date %in% c(zica, covid), 0, 1))

md <- 
  gam(bts ~ t + 
        s(mth, bs = 'cp'), 
      weights = w,
      data = bra2,
      family = "quasipoisson")

p <- predict(md, newdata = bra2, type = "response", se.fit = TRUE)

bra3 <- 
  bra2 %>% 
  mutate(bsn = p$fit,
         ul = p$fit + (2 * p$se.fit),
         ll = p$fit - (2 * p$se.fit))

bra3 %>% 
  ggplot()+
  geom_ribbon(aes(date, ymin = ll, ymax = ul), fill = "red", alpha = 0.3)+
  geom_point(aes(date, bts))+
  geom_line(aes(date, bts))+
  geom_line(aes(date, bsn), col = "red")+
  theme_bw()

# now, let's summarize the disturbances during each crisis
bra3 %>% 
  summarise(bts = sum(bts),
            bsn = round(sum(bsn)),
            exc = bts - bsn,
            psc = 100*exc/bsn,
            mts = n(),
            .by = c(per))


# Assignment: what is the impact of the C19 pandemic on births in 
# Colombia and Mexico? 
