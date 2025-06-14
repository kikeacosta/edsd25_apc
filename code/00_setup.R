# EDSD 2025
# Analysis of mortality disturbances course
# Instructor: Enrique Acosta (CED)
# Preparing environment

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# installing missing packages ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
libs <- c("tidyverse", 
          "lubridate",
          "readr",
          "viridisLite",
          "viridis", 
          "mgcv",
          "ISOweek",
          "countrycode",
          "patchwork",
          "isoband",
          "RColorBrewer")

# installing from CRAN
for (i in libs){
  if (!(i %in% rownames(installed.packages()))) install.packages(i)
}

# installing MortalitySmooth (Camarda, 2012) from a mirror copy in my GitHub repository
options(timeout = 600)
if (!("MortalitySmooth" %in% rownames(installed.packages()))) remotes::install_github("kikeacosta/MortalitySmooth")

# Loading required packages 
lapply(libs, require, character.only = T)
library("MortalitySmooth")

# avoiding scientific notation
options(scipen=999)
# let's keep the same seed for reproducibility of results
set.seed(2019) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# redefining the Lexis shape specs ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
redo_lexis_shape <- 
  function(pmin, pmax, amin, amax){
    lexis_shape <<-  
      list(
        geom_vline(xintercept = seq(pmin, pmax, 10), 
                   linewidth = 0.2, linetype = "dashed", 
                   alpha = 0.8, color = "grey30"),
        geom_hline(yintercept = seq(amin, amax, 10), 
                   linewidth = 0.2, linetype = "dashed", 
                   alpha = 0.8, color = "grey30"),
        # adding cohorts
        geom_abline(intercept = seq(-pmax, -(pmin-amax), 10), slope = 1, 
                    linetype = "dashed", color = "grey30", 
                    linewidth = .2, alpha = 0.8),
        coord_equal(expand = 0),
        # adding proper labels to both axis
        scale_x_continuous(breaks = seq(pmin, pmax, 10)),
        scale_y_continuous(breaks = seq(amin, amax, 10)),     
        # adding axis titles
        labs(y = "Age", x = "Period"),
        theme_bw()
      )
  }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot a Lexis surface of mortality change over time ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_change <- function(c, s, amin, amax, ymin, ymax){
  # filtering data, and computing rates and log_rates 
  db2 <- 
    hmd %>% 
    mutate(deaths = dts + 1,
           Mx = 100000 * dts / pop, 
           log_m = log(Mx)) %>% 
    filter(code == c,
           sex == s,
           age %in% amin:amax, 
           year %in% ymin:ymax)
  
  amin2 <- db2 %>% pull(age) %>% min()
  amax2 <- db2 %>% pull(age) %>% max()
  ymin2 <- db2 %>% pull(year) %>% min()
  ymax2 <- db2 %>% pull(year) %>% max()
  
  ylist <- unique(db2$year) %>% sort() # all periods in data, ordered
  alist <- unique(db2$age) %>% sort() # all ages in data, ordered
  
  # mortality
  deaths <- 
    matrix(db2$dts, nrow = length(alist), ncol = length(ylist), byrow = F)
  colnames(deaths) <- ylist
  rownames(deaths) <- alist
  
  # population at risk (population/exposure)
  exposure <- 
    matrix(db2$pop, nrow = length(alist), ncol = length(ylist), byrow=F)
  colnames(exposure) <- ylist
  rownames(exposure) <- alist
  
  # smoothing mortality with best AIC (that is method=2, also possible best BIC with method=1)
  fit <- 
    Mort2Dsmooth(x = alist, y = ylist, Z = deaths, offset = log(exposure),
                 overdispersion = TRUE, method = 2)
  
  # transforming them to smoothed deaths
  mx_smooth <- (exp(fit$logmortality) * 100000)
  
  # from matrix to tidy form (adjusting ages to start in 0)
  smt <- 
    mx_smooth %>% 
    as_tibble(rownames = "age") %>%
    gather(-age, key = year, value = Mx) %>% 
    mutate(type = "m_smoothed",
           age = as.integer(age),
           year = as.integer(year))
  
  # replacing missing values and estimating log_rates (/100k)
  smt2 <- 
    smt %>% 
    replace_na(list(Mx = 0)) %>% 
    mutate(log_m = log(Mx))
  
  db_per <- 
    smt2 %>% 
    group_by(age) %>% 
    mutate(ch = ((Mx / lag(Mx)) - 1) * 100) %>% 
    ungroup() %>% 
    drop_na()
  
  qt <- 25
  # one for decrease of mortality (from green to blue), 
  # another for mortality increases (from yellow to red)
  col_scale <- c(colorRampPalette(c("royalblue", "springgreen"), space = "Lab")(qt),
                 colorRampPalette(c("yellow", "red"), space = "Lab")(qt))
  
  # Definition of brackets for the scale of change
  val <- unique(db_per$ch)
  # separate negative and positive values
  pval <- val[val>0] # all positive values of change (mortality deterioration)
  nval <- val[val<0] # all negative values of change (mortality improvement)
  # identification of brackets for positive values (minimum + 23 quantiles + maximum)
  pcop <-c(min(pval), quantile(pval, prob=1/(qt-1)*(1:(qt-2))), max(pval)*1.01)
  # the same as above but for negative values (minimum + 23 quantiles + maximum)
  ncop <-c(min(nval)*1.01, quantile(nval, prob=1/(qt-1)*(1:(qt-2))), max(nval)*1.01) 
  # chain of brackets 25 ranges for negative changes, central value of no change (0), 
  # and 25 ranges for positive changes
  breaks_mc <- c(ncop, 0, pcop) 
  
  # adding to each value of change (continuous) the corresponding bracket (a discrete interval)
  db_per2 <- 
    db_per %>% 
    mutate(ch_cut = cut(db_per$ch, breaks = breaks_mc)) 
  
  # all categories
  cuts <- db_per2 %>% pull(ch_cut) %>% unique() %>% sort()
  
  # assigning color to each category
  col_values <- setNames(col_scale, cuts)
  
  # Plot of mortality change over periods
  p_change <- 
    db_per2 %>% 
    ggplot(aes(year, age, z = ch_cut)) +
    geom_tile(aes(fill = ch_cut))+
    # adding the color palette constructed above
    scale_fill_manual(values = col_values, 
                      breaks = cuts[c(1, 10, 20, 25, 30, 40, 50)], 
                      name = "Mortality\nchange %")+ 
    #adding contour lines when the slope is 0
    geom_contour(aes(z = ch), breaks = 0, col="black", 
                 alpha=.7, linewidth = .3)+ 
    scale_x_continuous(expand = c(0,0), breaks = seq(ymin2, ymax2, 10)) +
    scale_y_continuous(expand = c(0,0), breaks = seq(amin2, amax2, 10))+
    # adding the grid of the Lexis diagram
    geom_vline(xintercept = seq(ymin2, ymax2, 10), linetype = "dashed", 
               color = "grey30", linewidth = .20, alpha = 0.8) +
    geom_hline(yintercept = seq(amin2, amax2, 10), linetype = "dashed", 
               color = "grey30", linewidth = .20, alpha = 0.8) +
    geom_abline(intercept = seq(-2020, -(ymin2-amax2), 10), slope = 1, 
                linetype = "dashed", color = "grey30", 
                linewidth = .2, alpha = 0.8)+
    labs(x="Period", y="Age")+
    # aesthetic details
    coord_equal() +
    labs(title = paste0(c, "_", s),
         x="Period", y="Age")+
    theme_minimal()+
    # deleting default discrete legend
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      plot.title = element_text(size = 12),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(colour = NA),
      panel.grid.minor = element_line(colour = NA),
      plot.background = element_rect(fill = "white", colour = "transparent")
    )
  
  return(p_change)
}


# extracting coefficients from an APC model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extract_coeffs <- function(mod){
  tp1 <- 
    coef(summary(mod)) %>% 
    as_tibble(rownames = "coeff") %>% 
    mutate(tdim = case_when(str_detect(coeff, "\\(A\\)") ~ "Age",
                            str_detect(coeff, "\\(P\\)") ~ "Period",
                            str_detect(coeff, "\\(C\\)") ~ "Cohort",
                            str_detect(coeff, "I\\(C") ~ "Drift",
                            str_detect(coeff, "I\\(P") ~ "Drift"),
           effect = exp(Estimate),
           ll = exp(Estimate - 1.96*`Std. Error`),
           ul = exp(Estimate + 1.96*`Std. Error`)) %>% 
    separate(coeff, c("trash", "value"), sep = "\\)") %>% 
    mutate(value = value %>% as.double()) %>% 
    select(tdim, value, effect, ll, ul)
  
  ps_model <- tp1 %>% filter(tdim == "Period") %>% pull(value)
  ps_data <- unique(dt2$P)
  ps_miss <- ps_data[!(ps_data %in% ps_model)]
  
  cs_model <- tp1 %>% filter(tdim == "Cohort") %>% pull(value)
  cs_data <- unique(dt2$C)
  cs_miss <- cs_data[!(cs_data %in% cs_model)]
  
  tp2 <- 
    tp1 %>% 
    bind_rows(tibble(tdim = "Period", value = ps_miss, effect = 1, ll = 1, ul = 1),
              tibble(tdim = "Cohort", value = cs_miss, effect = 1, ll = 1, ul = 1)) %>% 
    mutate(tdim = factor(tdim, levels = c("Age", "Period", "Cohort"))) %>% 
    arrange(tdim, value) %>% 
    filter(tdim != "Drift")
  
  drift <- 
    tp1 %>% 
    filter(tdim == "Drift") %>% pull(effect)
  
  out <- 
    list("coeffs" = tp2, "drift" = drift)
  return(out)
}

# extracting the drift from an APC model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_drift <- 
  function(mdl){
    
    pcoeff <- 
      coef(summary(mdl))[,1] %>%
      as_tibble(rownames = "coeff") %>%
      filter(coeff %in% c("P", "C")) %>%
      pull(value)
    
    # this is in log scale, let's transform it in percentage value
    # change of mortality between periods (%)
    drift <- (exp(pcoeff) - 1)*100
    drift
    
    # interpretation
    cat(paste0(round(drift, 3),
               "% change in mortality each period/cohort category"))
  }


# plotting Carstensen APC model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_carst <- 
  function(mod){
    
    age_eff <- 
      mod$Age %>% 
      as_tibble() %>% 
      mutate(dim = "Age")
    
    per_eff <- 
      mod$Per %>% 
      as_tibble() %>% 
      rename(Year = 1,
             RR = 2) %>% 
      mutate(dim = "Period")
    
    coh_eff <- 
      mod$Coh %>% 
      as_tibble() %>% 
      rename(Year = 1,
             RR = 2) %>% 
      mutate(dim = "Cohort")
    
    bind_rows(coh_eff, per_eff)
    
    p_a_eff <- 
      age_eff %>% 
      ggplot()+
      geom_line(aes(Age, Rate))+
      scale_y_log10()+
      facet_wrap(~dim, scales = "free_x")+
      theme_bw()
    
    p_pc_eff <- 
      bind_rows(coh_eff, per_eff) %>% 
      mutate(dim = factor(dim, levels = c("Period", "Cohort"))) %>% 
      ggplot()+
      geom_ribbon(aes(Year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.3)+
      geom_line(aes(Year, RR))+
      scale_y_log10(breaks = c(0.2, 0.5, 0.7, 0.8, 1, 1.2, 1.5, 2, 4, 5))+
      scale_x_continuous(breaks = seq(1800, 2020, 10))+
      theme_bw()+
      facet_grid(~dim, scales = "free_x", space = "free_x")+
      geom_hline(yintercept = 1, linetype = "dashed")
    
    p_a_eff+p_pc_eff+ 
      plot_layout(widths = c(1, 2.5))
  }


# cd <- "CAN"
# sx <- "b"
# ag <- "TOT"
# ymin <- 2015

obtain_excess <- 
  function(cd, sx, ag, ymin){
    
    zip_files <- unzip("data_input/STMFinput.zip", list = TRUE)
    
    file_name <- zip_files %>% filter(str_detect(Name, cd)) %>% pull(Name)
    
    dt <- 
      read_csv(unz("data_input/STMFinput.zip", file_name)) %>% 
      mutate(Week = as.double(Week))
    
    # adding date to each ISO week, using the package ISOweek
    # we need the week in this format: 2000-W01-7
    dt2 <- 
      dt %>% 
      # renaming all variablees to lower case
      rename_with(str_to_lower) %>% 
      # adding the date
      mutate(isoweek = paste0(year, "-W", sprintf("%02d", week), "-7"),
             date = ISOweek2date(isoweek)) %>% 
      filter(age == ag,
             sex == sx,
             year >= ymin) %>% 
      select(code = popcode, year, week, dts = deaths, date)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # estimating the baseline mortality ====
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # 1. average-week
    # ~~~~~~~~~~~~~~~
    # the simplest approach (maybe the worst)
    
    # estimating the average weekly deaths for the whole training period
    # as the baseline
    w_av <- 
      dt2 %>% 
      filter(year <= 2019) %>%
      summarise(bsn = mean(dts, na.rm = TRUE)) %>% 
      pull(bsn)
    
    dt_w_av <- 
      dt2 %>% 
      mutate(bsn = w_av)
    
    
    # 2. Week-specific average
    # ~~~~~~~~~~~~~~~~~~~~~~~~
    ws_av <- 
      dt2 %>% 
      filter(year <= 2019) %>%
      group_by(week) %>% 
      summarise(bsn = mean(dts, na.rm = TRUE)) %>% 
      ungroup()
    
    dt_ws_av <- 
      dt2 %>% 
      left_join(ws_av)
    
    
    
    # 3. Excess mortality using a Poisson model
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # exposures
    # ~~~~~~~~~
    # issue: population counts only available annually
    # solution: interpolation
    
    # loading total population counts from WPP 
    pop <- read_rds("data_input/wpp2022_pop.rds")
    
    # selecting the Spanish population
    pop2 <- 
      pop %>% 
      mutate(code = countrycode(name, origin = "country.name",
                                destination = "iso3c")) %>% 
      filter(code == cd)
    
    # from annual to weekly exposures
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # wee need weekly exposures between 2010 and 2023
    # then, we can interpolate between 2009 and 2024
    
    # first, create a grid with all weeks that we need to populate with 
    # weekly exposures
    pop_empty <- 
      # 52 weeks per year
      expand_grid(year = 2009:2024, week = 1:52) %>% 
      # adding extra weeks for leap years 2009, 2015, and 2020
      bind_rows(tibble(year = c(2009, 2015, 2020), week = 53)) %>% 
      # arrange it chronologically
      arrange(year, week)
    
    # preparing data for interpolation
    pop_inter <- 
      pop_empty %>%  
      # adding annual population in midyear to week 26 (the midyear week!)
      left_join(pop2 %>%
                  select(year, pop) %>% 
                  mutate(week = 26)) %>% 
      mutate(t = 1:n(), # creating a continuous variable for week sequence 
             # transforming year-week to date, using ISOweek2date() function
             isoweek = paste0(year, "-W", sprintf("%02d", week), "-7"),
             date = ISOweek2date(isoweek))
    
    # extracting values for interpolation
    xs <- pop_inter %>% drop_na() %>% pull(t) # weeks with data
    ys <- pop_inter %>% drop_na() %>% pull(pop) # available pop data
    ts <- pop_inter %>% pull(t) # all weeks in which we need estimates
    
    # smoothing using cubic splines
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # the "spline()" function allows to interpolate by constraining the curve 
    # to match the original data. In other words, it is strictly interpolation 
    # rather than smoothing  
    
    # extracting predictions from the model
    interpolated_pop <- spline(xs, ys, xout = ts)
    
    pop_inter2 <- 
      pop_inter %>% 
      mutate(pop2 = interpolated_pop$y)
    
    # visualizing annual values and weekly estimates
    pop_inter2 %>% 
      ggplot()+
      geom_line(aes(date, pop2))+
      geom_point(aes(date, pop), col = "red")+
      theme_bw()
    
    # weekly exposures to use
    pop3 <- 
      pop_inter2 %>% 
      select(-pop, -t, -isoweek) %>% 
      rename(pop = pop2)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # merging deaths and exposures
    dt3 <- 
      dt2 %>% 
      left_join(pop3) %>% # merging with weekly population
      mutate(exposure = pop / 52, # exposure in person-weeks
             t = 1:n(), # a variable for the secular trend 
             w = ifelse(date <= "2020-03-15", 1, 0)) # a variable for the weights 
    
    
    # GAM models allow us to include parametric and semiparametric terms together 
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
          data = dt3, 
          # using a quasipoisson distribution to account for overdispersion
          family = "quasipoisson") 
    
    
    # example for predicting estimates
    bsn_poi <- predict(gam_model, newdata = dt3, type = "response", se.fit = TRUE)
    
    # obtaining estimates for the three models
    bsn <- 
      dt3 %>% 
      mutate(bsn = bsn_poi$fit,
             se = bsn_poi$se.fit,
             ll = bsn - 1.96*se,
             ul = bsn + 1.96*se)
    
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
    
    # plotting the three baselines
    p_bsns <- 
      bsns %>% 
      ggplot()+
      geom_line(aes(date, dts), linewidth = 1)+
      geom_line(aes(date, bsn, col = type), linewidth = 1)+
      scale_color_manual(values = cols)+
      labs(title = cd)+
      theme_bw()
    p_bsns
    
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
    
    # obtaining annual excess
    yr_exc <- 
      excs %>% 
      filter(date >= "2020-03-15",
             date <= "2023-12-31") %>% 
      group_by(year, type) %>% 
      summarise(exc = sum(exc, na.rm = TRUE)) %>% 
      ungroup()
    
    tots <- 
      yr_exc %>% 
      spread(type, exc) %>% 
      summarise(poisson_model = sum(poisson_model),
                weekly_average = sum(weekly_average),
                weekly_spc_average = sum(weekly_spc_average),
                year = "total")
    
    t1 <- round(tots$poisson_model)
    t2 <- round(tots$weekly_average)
    t3 <- round(tots$weekly_spc_average)
    
    d2 <- paste0(round((t2 / t1 - 1)*100, 1), "%")
    d3 <- paste0(round((t3 / t1 - 1)*100, 1), "%")
    
    # plotting excess estimates by year
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    p_tots <- 
      yr_exc %>% 
      ggplot(aes(type, exc, fill = factor(year))) + 
      geom_bar(position="stack", stat="identity")+
      geom_text(aes(x = "poisson_model", y = 0), label = paste0(t1, "\nreference"), vjust = -1)+
      geom_text(aes(x = "weekly_average", y = 0), label = paste0(t2, "\n", d2), vjust = -1)+
      geom_text(aes(x = "weekly_spc_average", y = 0), label = paste0(t3, "\n", d3), vjust = -1)+
      geom_hline(yintercept = 0, lty = "dashed")+
      labs(fill = "year", title = cd)+
      coord_cartesian(expand = 0)+
      theme_bw()
    
    
    out_f <- list(plot_bsns <- p_bsns,
                  plot_tots <- p_tots,
                  yr_exc <- yr_exc)
    
    return(out_f)
  }
