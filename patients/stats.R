library(tidyverse)

dice_dt <- read_tsv('dice_stats.txt', col_names = c('patient', 'week', 'dice')) %>%
    mutate(patient = gsub("_.*", "", patient)) %>%
    mutate(week.no = as.numeric(gsub("cbct_", "", week)))

ggplot(data = dice_dt, aes(x = week.no, y = dice, color = patient)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    facet_wrap(~ patient, scales = 'free')

ggplot(data = dice_dt, aes(x = patient, y = dice, fill = factor(patient))) +
    geom_boxplot() +
    geom_jitter(alpha = 0.5, width = 0.15)

dice_stats <- dice_dt %>%
    group_by(patient) %>%
    summarize(mean_dice = mean(dice),
              sd_dice = sd(dice),
              min_dice = min(dice),
              max_dixe = max(dice))

#######################
vol_dt <- read_tsv('vol_stats.txt', col_names = c('patient', 'week', 'nvox', 'dx', 'dy', 'dz')) %>%
    mutate(patient = gsub("_.*", "", patient)) %>%
    mutate(week.no = as.numeric(gsub("cbct_", "", week))) %>%
    mutate(vol = nvox*dx*dy*dz/1000)

base_vols <- vol_dt %>%
    filter(week.no == 0) %>% mutate(pat = patient, v = vol) %>% select(pat, v)
vol_dt <- vol_dt %>%
    filter(week.no > 0) %>%
    group_by(patient) %>%
    mutate(vol_ratio = vol/filter(base_vols, pat == patient)$v)

ggplot(data = vol_dt, aes(x = week.no, y = vol_ratio, color = patient)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    facet_wrap(~ patient, scales = 'free')

ggplot(data = vol_dt, aes(x = patient, y = vol_ratio, fill = factor(patient))) +
    geom_boxplot() +
    geom_jitter(alpha = 0.5, width = 0.15)



    
vol_stats <- vol_dt %>%
    group_by(patient) %>%
    summarize(mean_vol = mean(vol_ratio),
              sd_vol = sd(vol_ratio),
              min_vol = min(vol_ratio),
              max_vol = max(vol_ratio))
