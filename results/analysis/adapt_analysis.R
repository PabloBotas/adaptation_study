library(tidyverse)
library(data.table)
library(stringr)
library(lettercase)
library(gridExtra)
library(grid)
library(lattice)
library(zoo)

cbPalette <- c('#1f78b4', '#a6cee3', '#b2df8a', '#33a02c', '#fb9a99', '#000000', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#b15928',
               '#1f78b4', '#a6cee3', '#b2df8a', '#33a02c', '#fb9a99', '#000000', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#b15928')
if (Sys.info()["nodename"] == "gmoc-Precision-T7600") {
    setwd('~/Desktop/pablo/adaptive_project/results/analysis')
}

patients <- c('P01', 'P02', 'P03', 'P04', 'P05', 'P07', 'P10', 'P14', 'P15', 'P16')
# patients <- c('P01', 'P02', 'P03', 'P05', 'P07', 'P14', 'P15')
# patients <- c('P05')
data_location <- '../data/dvh'
# plots_dir <- 'plots_adapted'
plots_dir <- 'adapted'
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
target_dose <- 60
dose_fraction <- 60
plan_fractions <- 1

source('adapt_parse_data.R')
dvhs <- readData(patients, data_location, plan_fractions, target_dose, FALSE)
### Filter current analysis ---------------------------------
# dvhs <- filter(dvhs, week.no < 7)
# dvhs <- filter(dvhs, week.no != 2)
# dvhs <- filter(dvhs, method != 'Geometric')
# dvhs <- filter(dvhs, constraint != 'Range shifter')
# dvhs <- filter(dvhs, constraint != 'Isocenter - Range shifter') %>%
    # mutate(constraint = factor(constraint))
# dvhs$method <- mapvalues(dvhs$method, from = c("Weights"), to = c("Geometric"))

## Flip Submandible glands for patient 2
dvhs <- dvhs %>%
    mutate(struct = ifelse(patient == "Patient 2",
                           plyr::mapvalues(as.character(struct), c('L. Submandible gland', 'R. Submandible gland'),
                                         c('R. Submandible gland', 'L. Submandible gland')),
                           as.character(struct)))
# summary(dvhs)

### Get stats and normalize Plan dose ------------------
source('adapt_plans_utils.R')
#print_summaries(dvhs)
dvhs_back <- dvhs
dvhs <- dvhs_back
norm_expression <- 'V98 = 98'
# temp <- calculateStats(dvhs, dose_fraction, norm_expression)
temp <- calculateStats(dvhs, dose_fraction)
dt.stats <- temp$stats
dvhs <- temp$mod_dvhs
rm(temp)
summary(dt.stats)

adapt_list <- select(filter(dt.stats, constraint == 'None', struct == 'CTV', V95 < 95), pat.week)
dt.stats <- dt.stats %>% mutate(need_adapt = (pat.week %in% adapt_list$pat.week) | (week.no == 0))
dvhs <- dvhs %>% mutate(need_adapt = (pat.week %in% adapt_list$pat.week) | (week.no == 0))

### Patient DVHs --------------------------------
for (p in levels(dvhs$patient_orig)) {
    # p <- "P16"
    print(paste('Creating patient', p, 'DVH'))
    plt <- ggplot(filter(dvhs, patient_orig == p, !(method %in% c('Plan', 'None', 'Robust'))),
           aes(x = dose_pct, y = vol, color = struct)) +
        geom_path(data = select(filter(dvhs, patient_orig == p, struct == 'CTV', method == 'Plan'), -week.name, -method),
                  color = 'green', alpha = 1, size = 1) +
        geom_path(data = select(filter(dvhs, patient_orig == p, struct == 'CTV', method == 'None'), -method),
                  color = 'red', alpha = 1, size = 1) +
        geom_path(data = select(filter(dvhs, patient_orig == p, struct == 'CTV', method == 'Robust'), -method),
                  color = 'orange', alpha = 1, size = 1) +
        geom_path(aes(linetype = constraint), alpha = 1) +
        facet_grid(method ~ week.name) +
        scale_x_continuous(limits = c(0, 125), breaks = seq(0, 100, by = 10)) +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
        scale_linetype_manual('Method', values = c('solid', 'dashed', 'dotted', 'longdash', '12345678', 'solid', 'twodash', 'F4')) +
        scale_colour_manual('Contour', values = cbPalette) +
        theme_bw() +
        theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
        labs(x = "Dose (%)", y = "Contour volume (%)", title = paste(p, "DVH evolution"))
    print(plt)
    ggsave(paste0(plots_dir, '/DVHs_', p, '.pdf'), width = 50, height = 14, units = "cm")
}

## Single DVH week
plt <- ggplot(filter(dvhs, patient_orig == p, !(method %in% c('Plan', 'None', 'Robust'))),
              aes(x = dose_pct, y = vol, color = struct)) +
    geom_path(data = select(filter(dvhs, patient_orig == p, struct == 'CTV', method == 'Plan'), -week.name, -method),
              color = 'green', alpha = 1, size = 1) +
    geom_path(data = select(filter(dvhs, patient_orig == p, struct == 'CTV', method == 'None'), -method),
              color = 'red', alpha = 1, size = 1) +
    geom_path(aes(linetype = constraint), alpha = 1) +
    facet_grid(method ~ week.name) +
    scale_x_continuous(limits = c(0, 125), breaks = seq(0, 100, by = 10)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
    scale_linetype_manual('Method', values = c('solid', 'dashed', 'dotted', 'longdash', '12345678', 'solid', 'twodash', 'F4')) +
    scale_colour_manual('Contour', values = cbPalette) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = "Dose (%)", y = "Contour volume (%)", title = paste(p, "DVH evolution"))
print(plt)
ggsave(paste0(plots_dir, '/DVHs_V98_98_', p, '.pdf'), width = 50, height = 14, units = "cm")

### Useful functions -----------------------
plot_week_parameter <- function(df, s, par, lab, plots_dir, additional_dir = '', to_file=FALSE) {
    # Deduce file name and y axis label
    str_struct <- ifelse(str_detect(s, 'TV'), s, str_lowercase(s))
    str_par <- paste0(toupper(substr(par, 1, 1)), substr(par, 2, nchar(par)))
    capt <- bquote(bold("Fig.:") ~ .(str_par) ~ "of plan and adaptations across fractions in" ~ .(str_struct) ~ ".")
    
    n_patients <- nrow(unique(select(filter(dt.stats, struct == s), patient)))
    text_angle <- ifelse(n_patients %in% c(3,5,6), 90, 0)

    p <- ggplot(data = filter(df, struct == s, method %in% c("Geometric", "Weights")),
                     aes(x = week.name, y = get(par))) +
        geom_point(data = select(filter(df, struct == s, method == "Plan"), -method), size = 2, color = "red") +
        geom_hline(data = select(filter(df, struct == s, method == "Plan"), -method),
                   aes(yintercept = get(par)), color = 'forestgreen', alpha = 0.65, linetype = 'dashed') +
        geom_line(data = select(filter(df, struct == s, method == "None"), -method),
                  group = 1, color = 'black', alpha = 0.6, size = 2) +
        geom_line(aes(group = case, color = constraint, linetype = method)) +
        geom_point(aes(group = case, color = constraint, shape = constraint), size = 2) +
        facet_wrap( ~ patient, scales = "free_y") +
        scale_fill_manual(values = c("black"), labels = c("")) +
        scale_x_discrete(drop = FALSE) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(angle = text_angle, size = 8),
              legend.position = 'top') +
        guides(linetype = guide_legend(title = "Method", title.position = "top", title.hjust = 0.5),
               color = guide_legend(title = "Shift constraint", title.position = "top", title.hjust = 0.5),
               shape = guide_legend(title = "Shift constraint", title.position = "top", title.hjust = 0.5),
               fill = guide_legend("Plan at week"), title.position = "top", title.hjust = 0.5) +
        labs(y = lab, caption = capt) +
        theme(plot.caption = element_text(hjust = 0.5))
    
    print(p)
    
    if (to_file) {
        directory <- ifelse(filter(dt.stats, struct == s)$oar[1], "OARs", "target")
        directory <- paste0(plots_dir, '/', directory, additional_dir)
        dir.create(directory, showWarnings = FALSE, recursive = TRUE)
        path <- paste0(directory, '/adapted.')
        outfile <- paste0(path, gsub( " .*$", "", lab), ".", s, ".week")
        ggsave(paste0(outfile, '.pdf'), width = 20, height = 12.5, units = "cm")
        ggsave(paste0(outfile, '.jpg'), width = 20, height = 12.5, units = "cm")
    }
}

# dt.stats <- dt.stats %>%
#     mutate(patient = as.character(patient)) %>%
#     filter(patient != "Patient 4") %>%
#     mutate(patient = ifelse(patient == "Patient 5", "Patient 4", patient)) %>%
#     mutate(patient = ifelse(patient == "Patient 6", "Patient 5", patient)) %>%
#     mutate(patient = factor(patient))
# 
# ## FOR HARALD's PLOTS
# dt.stats <- dt.stats %>%
#     filter(constraint %in% c('Free', 'Isocenter', 'None', 'Plan')) %>%
#     filter(week.no != 7) %>%
#     mutate_if(is.factor, factor)
# 
# ggplot(data = filter(dt.stats, struct == 'CTV',
#                      method %in% c('Weights', 'None', 'Plan'),
#                      constraint %in% c('Free', 'Isocenter', 'None', 'Plan')),
#        aes(x = case, y = shape, fill = case)) +
#     geom_boxplot() +
#     facet_wrap(~ patient)


### CTV Target Plots ------------------------------------------------
# Mean dose per week CTV
plot_week_parameter(dt.stats, 'CTV', 'mean', "Mean Dose (Gy(RBE))", plots_dir, '_V98_98', FALSE)
plot_week_parameter(dt.stats, 'CTV', 'max', "Max Dose (Gy(RBE))", plots_dir, '_V98_98', FALSE)
plot_week_parameter(dt.stats, 'CTV', 'min', "Min Dose (Gy(RBE))", plots_dir, '_V98_98', FALSE)
plot_week_parameter(dt.stats, 'CTV', 'D98', "D98 (Gy(RBE))", plots_dir, '_V98_98', FALSE)
plot_week_parameter(dt.stats, 'CTV', 'D95', "D95 (Gy(RBE))", plots_dir, '_V98_98', FALSE)
plot_week_parameter(dt.stats, 'CTV', 'V100', "V100 (%)", plots_dir, '_V98_98', FALSE)
plot_week_parameter(dt.stats, 'CTV', 'V98', "V98 (%)", plots_dir, '_V98_98', FALSE)
plot_week_parameter(dt.stats, 'CTV', 'V95', "V95 (%)", plots_dir, '_V98_98', FALSE)
plot_week_parameter(dt.stats, 'CTV', 'V107', "V107 (%)", plots_dir, '_V98_98', FALSE)
plot_week_parameter(dt.stats, 'CTV', 'V110', "V110 (%)", plots_dir, '_V98_98', FALSE)
plot_week_parameter(dt.stats, 'CTV', 'D5_D95', "D5-D95 (Gy(RBE))", plots_dir, '_V98_98', FALSE)
plot_week_parameter(dt.stats, 'CTV', 'D2_D98', "D2-D98 (Gy(RBE))", plots_dir, '_V98_98', FALSE)

ctv.averages <- dt.stats %>% filter(oar == FALSE) %>% group_by(patient, method, constraint) %>%
    # summarise_at(.vars = c("V98", "D2_D98"), .funs = c(mean = "mean", sd = "sd")) %>%
    summarise_at(.vars = c("V98", "D2_D98", "mean"), .funs = c(mean = "mean")) %>%
    group_by(patient) %>%
    mutate(mean_mean_ratio_plan = 100*(mean_mean - mean_mean[method == "Plan"])/mean_mean[method == "Plan"]) %>%
    mutate(V98_mean_ratio_plan = 100*(V98_mean - V98_mean[method == "Plan"])/V98_mean[method == "Plan"]) %>%
    mutate(D2_D98_mean_ratio_plan = 100*(D2_D98_mean - D2_D98_mean[method == "Plan"])/D2_D98_mean[method == "Plan"]) %>%
    mutate(mean_mean_ratio_frac = 100*(mean_mean - mean_mean[method == "None"])/mean_mean[method == "None"]) %>%
    mutate(V98_mean_ratio_frac = 100*(V98_mean - V98_mean[method == "None"])/V98_mean[method == "None"]) %>%
    mutate(D2_D98_mean_ratio_frac = 100*(D2_D98_mean - D2_D98_mean[method == "None"])/D2_D98_mean[method == "None"])
write.table(ctv.averages, 'CTV.averages.AAPM_submission.txt', sep = "\t")

free_weights <- filter(select(ctv.averages, patient, method, constraint, D2_D98_mean, D2_D98_mean_ratio_plan, D2_D98_mean_ratio_frac),
                       method == 'Weights' & constraint == 'Free')
a <- mean(free_weights$D2_D98_mean)
sa <- sd(free_weights$D2_D98_mean)

none <- filter(select(ctv.averages, patient, method, constraint, D2_D98_mean, D2_D98_mean_ratio_plan, D2_D98_mean_ratio_frac),
               method == 'None' & constraint == 'None')
b <- mean(none$D2_D98_mean)
sb <- sd(none$D2_D98_mean)

s <- sqrt(100*sa*sa/(b*b) + (-100*a/(b*b)*sb)^2)
s

### OARs Plots ------------------------------------------------
plot_week_parameter(dt.stats, 'L. Submandible gland', 'mean', "Mean Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)
plot_week_parameter(dt.stats, 'R. Submandible gland', 'mean', "Mean Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)
plot_week_parameter(dt.stats, 'Oral cavity', 'mean', "Mean Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)
plot_week_parameter(dt.stats, 'Spinal Cord', 'mean', "Mean Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)
plot_week_parameter(dt.stats, 'Larynx', 'mean', "Mean Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)
plot_week_parameter(dt.stats, 'L. Parotid', 'mean', "Mean Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)
plot_week_parameter(dt.stats, 'R. Parotid', 'mean', "Mean Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)
plot_week_parameter(dt.stats, 'Mandible', 'mean', "Mean Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)

plot_week_parameter(dt.stats, 'L. Submandible gland', 'max', "Max Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)
plot_week_parameter(dt.stats, 'R. Submandible gland', 'max', "Max Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)
plot_week_parameter(dt.stats, 'Oral cavity', 'max', "Max Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)
plot_week_parameter(dt.stats, 'Spinal Cord', 'max', "Max Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)
plot_week_parameter(dt.stats, 'Larynx', 'max', "Max Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)
plot_week_parameter(dt.stats, 'L. Parotid', 'max', "Max Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)
plot_week_parameter(dt.stats, 'R. Parotid', 'max', "Max Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)
plot_week_parameter(dt.stats, 'Mandible', 'max', "Max Dose (Gy(RBE))", plots_dir, '_V98_98', TRUE)

lsubmandible.averages <- dt.stats %>% filter(struct == 'L. Submandible gland') %>% group_by(patient, method, constraint) %>%
    summarise_at(.vars = c("mean", "max"), .funs = c(mean = "mean")) %>%
    group_by(patient) %>%
    mutate(mean_mean_ratio_plan = 100*(mean_mean - mean_mean[method == "Plan"])/mean_mean[method == "Plan"]) %>%
    mutate(max_mean_ratio_plan = 100*(max_mean - max_mean[method == "Plan"])/max_mean[method == "Plan"]) %>%
    mutate(mean_mean_ratio_frac = 100*(mean_mean - mean_mean[method == "None"])/mean_mean[method == "None"]) %>%
    mutate(max_mean_ratio_frac = 100*(max_mean - max_mean[method == "None"])/max_mean[method == "None"])
write.table(lsubmandible.averages, 'lsubmandible.averages.AAPM_submission.txt', sep = "\t")

oral.averages <- dt.stats %>% filter(struct == 'Oral cavity') %>% group_by(patient, method, constraint) %>%
    summarise_at(.vars = c("mean", "max"), .funs = c(mean = "mean")) %>%
    group_by(patient) %>%
    mutate(mean_mean_ratio_plan = 100*(mean_mean - mean_mean[method == "Plan"])/mean_mean[method == "Plan"]) %>%
    mutate(max_mean_ratio_plan = 100*(max_mean - max_mean[method == "Plan"])/max_mean[method == "Plan"]) %>%
    mutate(mean_mean_ratio_frac = 100*(mean_mean - mean_mean[method == "None"])/mean_mean[method == "None"]) %>%
    mutate(max_mean_ratio_frac = 100*(max_mean - max_mean[method == "None"])/max_mean[method == "None"])
write.table(oral.averages, 'larinx.averages.AAPM_submission.txt', sep = "\t")

mandible.averages <- dt.stats %>% filter(struct == 'Mandible') %>% group_by(patient, method, constraint) %>%
    summarise_at(.vars = c("mean", "max"), .funs = c(mean = "mean")) %>%
    group_by(patient) %>%
    mutate(mean_mean_ratio_plan = 100*(mean_mean - mean_mean[method == "Plan"])/mean_mean[method == "Plan"]) %>%
    mutate(max_mean_ratio_plan = 100*(max_mean - max_mean[method == "Plan"])/max_mean[method == "Plan"]) %>%
    mutate(mean_mean_ratio_frac = 100*(mean_mean - mean_mean[method == "None"])/mean_mean[method == "None"]) %>%
    mutate(max_mean_ratio_frac = 100*(max_mean - max_mean[method == "None"])/max_mean[method == "None"])
write.table(mandible.averages, 'mandible.averages.AAPM_submission.txt', sep = "\t")

mean(filter(select(mandible.averages, patient, method, constraint, mean_mean, mean_mean_ratio_plan, mean_mean_ratio_frac),
            method == 'Weights' & constraint == 'Isocenter')$mean_mean)
sd(filter(select(mandible.averages, patient, method, constraint, mean_mean, mean_mean_ratio_plan, mean_mean_ratio_frac),
            method == 'Weights' & constraint == 'Free')$mean_mean)


mean(filter(select(mandible.averages, patient, method, constraint, max_mean, max_mean_ratio_plan, max_mean_ratio_frac),
            method == 'Weights' & constraint == 'Free')$max_mean)
sd(filter(select(lsubmandible.averages, patient, method, constraint, max_mean, max_mean_ratio_plan, max_mean_ratio_frac),
            method == 'None' & constraint == 'None')$max_mean)

get_ribbon_data <- function(input) {
    require(matrixStats)
    
    output <- input %>%
        mutate(facet_name = "All weeks") %>%
        select(patient, struct, vol, dose_pct, method, facet_name, week.name) %>%
        spread(week.name, vol)

    other_nms <- c("patient", "struct", "dose_pct", "method", "facet_name")
    nweeks <- ncol(output) - length(other_nms)
    week.names <- NULL
    for (i in 1:nweeks) {
        week.names <- c(week.names, paste('week', i, sep = '.'))
    }
    nms <- c(other_nms, week.names)
    names(output) <- nms
    output <- output %>%
        group_by(patient, struct) %>%
        mutate(week.1 = if (length(na.trim(week.1, sides = "right")) < n())
            c(na.trim(week.1, sides = "right"), rep(0, n() - length(na.trim(week.1, sides = "right"))))
            else week.1) %>%
        mutate(week.2 = if (length(na.trim(week.2, sides = "right")) < n())
            c(na.trim(week.2, sides = "right"), rep(0, n() - length(na.trim(week.2, sides = "right"))))
            else week.2) %>%
        mutate(week.3 = if (length(na.trim(week.3, sides = "right")) < n())
            c(na.trim(week.3, sides = "right"), rep(0, n() - length(na.trim(week.3, sides = "right"))))
            else week.3) %>%
        mutate(week.4 = if (length(na.trim(week.4, sides = "right")) < n())
            c(na.trim(week.4, sides = "right"), rep(0, n() - length(na.trim(week.4, sides = "right"))))
            else week.4) %>%
        mutate(week.5 = if (length(na.trim(week.5, sides = "right")) < n())
            c(na.trim(week.5, sides = "right"), rep(0, n() - length(na.trim(week.5, sides = "right"))))
            else week.5)
    
    if ('week.6' %in% names(output)) {
        output <- output %>%
            mutate(week.6 =
                if (all(is.na(week.6)))
                    NA
                else if (length(na.trim(week.6, sides = "right")) < n())
                    c(na.trim(week.6, sides = "right"), rep(0, n() - length(na.trim(week.6, sides = "right"))))
                else week.6)
    }
    if ('week.7' %in% names(output)) {
        output <- output %>%
            mutate(week.7 =
                if (all(is.na(week.7)))
                    NA
                else if (length(na.trim(week.7, sides = "right")) < n())
                    c(na.trim(week.7, sides = "right"), rep(0, n() - length(na.trim(week.7, sides = "right"))))
                else week.7)
    }

        
    output <- output %>%
        mutate(week.1 = na.approx(week.1, na.rm = FALSE)) %>%
        mutate(week.2 = na.approx(week.2, na.rm = FALSE)) %>%
        mutate(week.3 = na.approx(week.3, na.rm = FALSE)) %>%
        mutate(week.4 = na.approx(week.4, na.rm = FALSE)) %>%
        mutate(week.5 = na.approx(week.5, na.rm = FALSE))
    if ('week.6' %in% names(output)) {
        output <- mutate(output, week.6 = if (!all(is.na(week.6))) na.approx(week.6, na.rm = FALSE) else NA)
    }
    if ('week.7' %in% names(output)) {
        output <- mutate(output, week.7 = if (!all(is.na(week.7))) na.approx(week.7, na.rm = FALSE) else NA)
    }

    output <- output %>%
        ungroup() %>%
        mutate(mean = rowMeans(.[(length(other_nms) + 1):(length(other_nms) + nweeks)], na.rm = TRUE)) %>%
        mutate(sd = rowSds(as.matrix(.[(length(other_nms) + 1):(length(other_nms) + nweeks)]), na.rm = TRUE)) %>%
        rowwise() %>%
        mutate(ribbon_min = mean - sd) %>%
        mutate(ribbon_min = max(0, ribbon_min)) %>%
        mutate(ribbon_max = mean + sd) %>%
        mutate(ribbon_max = min(100, ribbon_max))

    return(output)
}

### Adaptations vs plan all patients ----------------------------
### 
library(foreach)
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)

cases <- unique(data.table(method = as.character(dt.stats$method), constraint = as.character(dt.stats$constraint)))
cases <- filter(cases, method == 'Weights')
# cases <- filter(cases, method == 'Weights' | method == 'Geometric')
strt <- Sys.time()
foreach(i = 1:nrow(cases), .packages = c('dplyr', 'tidyr', 'ggplot2', 'zoo', 'gridExtra')) %dopar% {
# for (i in 1:nrow(cases)) {
    meth <- cases[i, 'method']
    const <- cases[i, 'constraint']
    
    main_data <- filter(dvhs, method == meth, constraint == const)
    plan_data <- select(filter(dvhs, method == 'Plan'), -week.name, -method)
    none_data <- filter(dvhs, method == 'None')

    p <- ggplot(main_data, aes(x = dose_pct, y = vol, color = struct)) +
        geom_path(aes(linetype = 'solid')) +
        geom_path(data = plan_data, aes(linetype = 'dashed')) +
        geom_path(data = none_data, aes(linetype = 'dotted')) +
        facet_grid(week.name ~ patient) +
        scale_x_continuous(limits = c(0, 125), breaks = seq(0, 125, by = 20)) +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
        scale_colour_manual('Contour', values = cbPalette, guide = FALSE) +
        scale_linetype_manual('Case', labels = c("Plan", "Non Adapt", "Adapt"),
                              values = c('solid', 'dashed', 'dotted')) +
        theme_bw() +
        theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(face = c(rep('plain',5), 'bold', rep('plain', 1))),
              legend.position = "top",
              strip.text.x = element_text(size = 8, face = "bold"),
              strip.text.y = element_text(size = 8, face = "bold")) +
        labs(x = "", y = "Contour volume (%)",
             title = paste("Patient evolution:", meth, "method,", const, "strategy"))

    ribbon_none_data <- get_ribbon_data(none_data)
    ribbon_none_data <- mutate(ribbon_none_data, facet_name = "Non Adapt")
    q <- ggplot() +
        geom_ribbon(data = ribbon_none_data, aes(x = dose_pct, ymin = ribbon_min,
                                                 ymax = ribbon_max, fill = struct),
                    colour = 'black', size = 0.2, alpha = 0.6) +
        geom_path(data = plan_data, aes(x = dose_pct, y = vol, color = struct)) +
        facet_grid(facet_name ~ patient) +
        scale_x_continuous(limits = c(0, 125), breaks = seq(0, 125, by = 20)) +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
        scale_colour_manual('Contour', values = cbPalette) +
        scale_fill_manual('Contour', values = cbPalette) +
        theme_bw() +
        theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(face = c(rep('plain',5), 'bold', rep('plain', 1))),
              legend.position = "None",
              strip.text.y = element_text(size = 8, face = "bold"),
              strip.text.x = element_blank(),
              axis.title.x = element_text(color = "white"),
              axis.title.y = element_text(color = "white")) +
        guides(color = guide_legend(ncol = 6))
    
    ribbon_data <- get_ribbon_data(main_data)
    ribbon_data <- mutate(ribbon_data, facet_name = "Adapt")
    r <- ggplot() +
        geom_ribbon(data = ribbon_data, aes(x = dose_pct, ymin = ribbon_min,
                                            ymax = ribbon_max, fill = struct),
                    colour = 'black', size = 0.2, alpha = 0.6) +
        geom_path(data = plan_data, aes(x = dose_pct, y = vol, color = struct)) +
        facet_grid(facet_name ~ patient) +
        scale_x_continuous(limits = c(0, 125), breaks = seq(0, 125, by = 20)) +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
        scale_colour_manual('Contour', values = cbPalette) +
        scale_fill_manual('Contour', values = cbPalette) +
        theme_bw() +
        labs(x = "Dose (%)", y = "Contour volume (%)") +
        theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(face = c(rep('plain',5), 'bold', rep('plain', 1))),
              legend.position = "bottom",
              strip.text.y = element_text(size = 8, face = "bold"),
              strip.text.x = element_blank(),
              axis.title.y = element_text(color = "white")) +
        guides(color = guide_legend(ncol = 6))
    
    pq <- grid.arrange(p, q, r, ncol = 1, heights = c(2.25, 0.5, 1))
    ggsave(paste0(plots_dir, '/DVHs_', meth, '_', const, '.pdf'), pq, width = 25.0, height = 25.0, units = "cm")
    
    rm(p, ribbon_data, ribbon_none_data, main_data, plan_data, none_data, pq, meth, const)
}
stopCluster(cl)
print(Sys.time() - strt)

### Adaptations vs plan single case ----------------------------
### 
cases <- distinct(mutate_all(select(dt.stats, method, constraint), as.character))
cases <- filter(cases, method == 'Weights', constraint %in% c('Range shifter', 'Isocenter - Range shifter'))

for (pat.id in levels(dt.stats$patient_orig)) {
# for (pat.id in 'P15') {
    for (i in 1:nrow(cases)) {
        # for (i in 1:nrow(cases)) {
        meth <- cases[i, 'method']
        const <- cases[i, 'constraint']
        
        main_data <- filter(dvhs, method == meth, constraint == const, patient_orig == pat.id, week.no != 6)
        plan_data <- select(filter(dvhs, method == 'Plan', patient_orig == pat.id, week.no != 6), -week.name, -method)
        none_data <- filter(dvhs, method == 'None', patient_orig == pat.id, week.no != 6)
        
        p <- ggplot(main_data, aes(x = dose_pct, y = vol, color = struct)) +
            geom_path(aes(linetype = 'solid')) +
            geom_path(data = plan_data, aes(linetype = 'dashed')) +
            geom_path(data = none_data, aes(linetype = 'dotted')) +
            facet_wrap(~ week.name) +
            scale_x_continuous(limits = c(0, 125), breaks = seq(0, 125, by = 20)) +
            scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
            scale_colour_manual('Contour', values = cbPalette) +
            scale_linetype_manual('Case', labels = c("Plan", "Non Adapt", "Adapt"),
                                  values = c('solid', 'dashed', 'dotted')) +
            theme_bw() +
            theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
                  axis.text.x = element_text(face = c(rep('plain',5), 'bold', rep('plain', 1))),
                  legend.position = "top",
                  strip.text.x = element_text(size = 8, face = "bold"),
                  strip.text.y = element_text(size = 8, face = "bold")) +
            labs(x = "", y = "Contour volume (%)",
                 title = paste("Patient evolution:", meth, "method,", const, "strategy"))
        
        ribbon_none_data <- get_ribbon_data(none_data)
        ribbon_none_data <- mutate(ribbon_none_data, facet_name = "Non Adapted")
        q <- ggplot() +
            geom_ribbon(data = ribbon_none_data, aes(x = dose_pct, ymin = ribbon_min,
                                                     ymax = ribbon_max, fill = struct),
                        colour = 'black', size = 0.2, alpha = 0.6) +
            geom_path(data = plan_data, aes(x = dose_pct, y = vol, color = struct)) +
            facet_grid(patient ~ facet_name) +
            scale_x_continuous(limits = c(0, 125), breaks = seq(0, 125, by = 20)) +
            scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
            scale_colour_manual('Contour', values = cbPalette) +
            scale_fill_manual('Contour', values = cbPalette) +
            theme_bw() +
            labs(x = "Dose (%)", y = "Contour volume (%)") +
            theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
                  axis.text.x = element_text(face = c(rep('plain',5), 'bold', rep('plain', 1))),
                  legend.position = "None",
                  strip.text.x = element_text(size = 18, face = "bold"),
                  strip.text.y = element_blank()) +
            guides(color = guide_legend(ncol = 6))
        
        ribbon_data <- get_ribbon_data(main_data)
        ribbon_data <- mutate(ribbon_data, facet_name = "Adapted")
        r <- ggplot() +
            geom_ribbon(data = ribbon_data, aes(x = dose_pct, ymin = ribbon_min,
                                                ymax = ribbon_max, fill = struct),
                        colour = 'black', size = 0.2, alpha = 0.6) +
            geom_path(data = plan_data, aes(x = dose_pct, y = vol, color = struct)) +
            facet_grid(~ facet_name) +
            scale_x_continuous(limits = c(0, 125), breaks = seq(0, 125, by = 20)) +
            scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
            scale_colour_manual('Contour', values = cbPalette) +
            scale_fill_manual('Contour', values = cbPalette) +
            theme_bw() +
            labs(x = "Dose (%)", y = "Contour volume (%)") +
            theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
                  axis.text.x = element_text(face = c(rep('plain',5), 'bold', rep('plain', 1))),
                  legend.position = "None",
                  strip.text.x = element_text(size = 18, face = "bold"),
                  strip.text.y = element_blank()) +
            guides(color = guide_legend(ncol = 6))
        
        lay <- rbind(c(1,1),
                     c(2,3))
        pqr <- grid.arrange(p, q, r, heights = c(2, 1), layout_matrix = lay)
        ggsave(paste0(plots_dir, '/DVHs_normalized_', meth, '_', const, '_', pat.id, '.pdf'),
               pqr, width = 24.0, height = 24.0, units = "cm")
        # ggsave('~/Desktop/DVHs_pat14.pdf', pqr, width = 24.0, height = 24.0, units = "cm")
        # ggsave('~/Desktop/DVHs_pat1_nonadapt.pdf', q, width = 15.0, height = 12.0, units = "cm")
        # ggsave('~/Desktop/DVHs_pat1_adapt.pdf', r, width = 15.0, height = 10.0, units = "cm")
        
        rm(p, q, r, ribbon_data, ribbon_none_data, main_data, plan_data, none_data, pqr, meth, const)
    }
}


### Adaptations vs robust plan ----------------------------
cases <- unique(data.table(method = as.character(dt.stats$method), constraint = as.character(dt.stats$constraint)))
cases <- filter(cases, method == 'Weights' | method == 'Geometric')
dvhs <- filter(dvhs, patient != 'Patient 4')
dvhs <- filter(dvhs, patient == 'Patient 1')
dvhs$patient <- factor(dvhs$patient)
for (i in 1:nrow(cases)) {
    meth <- cases[i, 'method']
    const <- cases[i, 'constraint']
    
    adapt_data <- filter(dvhs, method == meth, constraint == const)
    robust_data <- filter(dvhs, method == 'Robust')
    main_data <- filter(dvhs, as.character(method) %in% c(meth, 'Robust'), as.character(constraint) %in% c(const, 'None'))
    none_data <- filter(dvhs, method == 'None')

    p <- ggplot(main_data, aes(x = dose_pct, y = vol, color = struct, linetype = method)) +
        geom_path() +
        facet_wrap(patient ~ week.name) +
        scale_x_continuous(limits = c(0, 125), breaks = seq(0, 125, by = 20)) +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
        scale_colour_manual('Contour', values = cbPalette) +
        scale_linetype_manual('Method', values = c('solid', 'dotted')) +
        theme_bw() +
        theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(face = c(rep('plain',10), 'bold', rep('plain', 2))),
              legend.position = "top",
              strip.text.x = element_text(size = 12, face = "bold"),
              strip.text.y = element_text(size = 12, face = "bold")) +
        labs(x = "Dose (%)", y = "Contour volume (%)",
             title = paste("Patient evolution versus ROBUST OPTIMIZATION:", meth, "method,", const, "strategy"))

    ribbon_adapt_data <- get_ribbon_data(adapt_data)
    ribbon_plan_data <- get_ribbon_data(robust_data)
    ribbon_data <- bind_rows(ribbon_adapt_data, ribbon_plan_data)
    ribbon_data$method <- factor(ribbon_data$method)
    ribbon_data$facet_name <- factor(ribbon_data$facet_name)
    rm(adapt_data, robust_data, ribbon_adapt_data, ribbon_plan_data)
    
    q <- ggplot(data = ribbon_data) +
        geom_ribbon(aes(x = dose_pct, ymin = ribbon_min,
                        ymax = ribbon_max, fill = struct,
                        linetype = method), alpha = 0.6, size = 0.3, colour = 'black') +
        facet_grid(facet_name ~ patient) +
        scale_x_continuous(limits = c(0, 125), breaks = seq(0, 125, by = 20)) +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
        scale_fill_manual('Contour', values = cbPalette) +
        scale_linetype_manual('Method', values = c('solid', 'dotted')) +
        theme_bw() +
        labs(x = "Dose (%)", y = "Contour volume (%)") +
        theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(face = c(rep('plain',10), 'bold', rep('plain', 2))),
              legend.position = "right",
              strip.text.y = element_text(size = 12, face = "bold"),
              strip.text.x = element_blank())
    
    pq <- grid.arrange(p, q, ncol = 1, heights = c(2.25, 1))
    ggsave(paste0(plots_dir, '/DVHs_V98_98_robust_', meth, '_', const, '.pdf'), pq, width = 21.5, height = 27.9, units = "cm")
    
    # rm(adapt_data, ribbon_adapt_data, ribbon_plan_data, robust_data, none_data, lay, pq, meth, const)
}


### Patient evolution ------------------------------
weekly_evolution_box0 <- function(df, x, y, fill, d, label) {
    ## color
    gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
    }
    cols <- gg_color_hue(3)
    if (str_detect(y, '_plan')) {
        cols <- cols[c(1,3,1)]
    } else {
        cols <- cols[c(2,1,3)]
    }
    ## filename
    filestr <- paste0(x, '_', gsub('_plan', '_diff', y))
    ## plot
    p <- ggplot(data = df, aes(x = get(x), y = get(y), fill = get(fill))) +
        geom_boxplot() +
        # geom_point() +
        geom_text(aes(label = patient.no), position = position_jitter(width = 0.2)) +
        geom_jitter(aes(color = patient.no), position = position_jitter(width = 0.2),
                    alpha = 1, stroke = 0.1) +
        scale_fill_manual('Stage', values = cols) +
        theme(legend.position = 'none') +
        labs(y = label, x = element_blank())
    if (str_detect(filestr, 'diff')) {
        p <- p + geom_hline(yintercept = 0)
    }
    ## save
    print(p)
    ggsave(plot = p, paste0(d, '/target/plan_evolution/plan_evolution_', filestr, '.pdf'),
           width = 13, height = 8, units = "cm")
}

weekly_evolution_comp <- function(df, y, outdir, label, label_thres = NA, label_margin = -100000) {
    df <- filter(df, struct == 'CTV', method %in% 'None')
    
    ## filename
    filestr <- paste0('week.name_', gsub('_plan', '_diff', y))
    ## plot
    p <- ggplot(data = df, aes(x = week.name, y = get(y))) +
        geom_hline(yintercept = 0, alpha = 0.6) +
        geom_boxplot() +
        geom_vline(xintercept = 1.5, alpha = 0.6)
    ## If point highlight
    if (!is.na(label_thres)) {
        p <- p + ggrepel::geom_label_repel(data = filter(df, get(y) < label_thres),
                                           aes(x = week.name, y = get(y), label = patient.no),
                                           nudge_y = label_margin, box.padding = 0.15,
                                           segment.size  = 0.3, force = 2,
                                           segment.color = "grey50", size = 2,
                                           arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"))
    }
    p <- p +
        geom_line(data = filter(df, stage == 'Weekly'),
                  aes(x = week.name, y = get(y), group = patient.no, color = patient.no, linetype = patient.no)) +
        geom_point(aes(fill = patient.no), pch = 21, stroke = 0.3) +
        guides(fill = guide_legend(title = "Pat. no.", keywidth = 2),
               linetype = guide_legend(title = "Pat. no.", keywidth = 2),
               color = guide_legend(title = "Pat. no.", keywidth = 2)) +
        theme(legend.position = 'right') +
        labs(y = label, x = element_blank())
    
    
        
    ## save
    print(p)
    for (ext in c('.pdf', '.png')) {
        f <- paste0(outdir, '/target/plan_evolution/plan_evolution_', filestr, ext)
        ggsave(plot = p, f, width = 15, height = 8, units = "cm", dpi = 600)
    }
}
weekly_evolution_comp(dt.stats, 'min_plan', plots_dir, 'Min dose difference [Gy(RBE)]', -15) 

patient_evolution_box <- function(df, y, color, d, label) {
    df <- filter(df, struct == 'CTV', method %in% 'None', stage %in% c('Cum.', 'Weekly'))
    ## color
    gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
    }
    cols <- gg_color_hue(3)
    cols <- cols[c(1,3,1)]
    print(cols)
    cols2 <- c('black', "#619CFF")

    ## filename
    filestr <- paste0('patient.no_', gsub('_plan', '_diff', y))
    ## plot
    p <- ggplot(data = df, aes(x = patient.no, y = get(y), color = get(color), fill = get(color))) +
        geom_boxplot() +
        scale_color_manual('Type', values = cols2) +
        scale_fill_manual('Type', values = cols) +
        labs(y = label, x = 'Patient number')
    if (str_detect(filestr, 'diff')) {
        p <- p + geom_hline(yintercept = 0)
    }
    ## save
    print(p)
    ggsave(plot = p, paste0(d, '/target/plan_evolution/plan_evolution_', filestr, '.pdf'),
           width = 13, height = 8, units = "cm", dpi = 600)
}

# == mean =================================================== ===
weekly_evolution_box(filter(dt.stats, struct == 'CTV', method %in% c('Plan', 'None')),
                  'week.name', 'mean', 'stage', plots_dir, 'Mean dose [Gy(RBE)]')
weekly_evolution_comp(dt.stats, 'mean_plan', plots_dir, 'Mean dose difference [Gy(RBE)]', -1.75)

patient_evolution_box(dt.stats, 'V95_plan', plots_dir, 'Mean dose difference (Gy(RBE))')

ggplot(data = filter(dt.stats, struct == 'CTV', method %in% c('Plan', 'None'), stage %in% c('Cum.')),
        aes(x = patient.no, y = V98_plan, fill = stage)) +
    geom_col(position = 'dodge')

ggplot(data = filter(dt.stats, struct == 'CTV', method %in% c('Plan', 'None'), stage %in% c('Cum.', 'Weekly')),
       aes(x = patient.no, y = D98_plan, color = stage)) +
    geom_boxplot()

# == max =================================================== ===
weekly_evolution_box(filter(dt.stats, struct == 'CTV', method %in% c('Plan', 'None')),
                  'week.name', 'max', 'stage', plots_dir, 'Max dose [Gy(RBE)]')
weekly_evolution_comp(dt.stats, 'max_plan', plots_dir, 'Max dose difference [Gy(RBE)]')

# == min =================================================== ===
weekly_evolution_box(filter(dt.stats, struct == 'CTV', method %in% c('Plan', 'None')),
                  'week.name', 'min', 'stage', plots_dir, 'Min dose [Gy(RBE)]')
weekly_evolution_comp(dt.stats, 'min_plan', plots_dir, 'Min dose difference [Gy(RBE)]', -15) 

# == D2_D98 =================================================== ===
weekly_evolution_box(filter(dt.stats, struct == 'CTV', method %in% c('Plan', 'None')),
                  'week.name', 'D2_D98', 'stage', plots_dir, 'D2-D98 [Gy(RBE)]')
weekly_evolution_comp(dt.stats, 'D2_D98_plan', plots_dir, 'D2-D98 difference [Gy(RBE)]')

# == V95 =================================================== ===
weekly_evolution_box(filter(dt.stats, struct == 'CTV', method %in% c('Plan', 'None')),
                  'week.name', 'V95', 'stage', plots_dir, 'V95 (%)')
weekly_evolution_comp(dt.stats, 'V95_plan', plots_dir, expression('V95 difference ('*Delta*'%)'), -12.5)

# == V98 =================================================== ===
weekly_evolution_box(filter(dt.stats, struct == 'CTV', method %in% c('Plan', 'None')),
                  'week.name', 'V98', 'stage', plots_dir, 'V98 (%)')
weekly_evolution_comp(dt.stats, 'V98_plan', plots_dir, expression('V98 difference ('*Delta*'%)'), -25)

# == V100 =================================================== ===
weekly_evolution_box(filter(dt.stats, struct == 'CTV', method %in% c('Plan', 'None')),
                  'week.name', 'V100', 'stage', plots_dir, 'V100 (%)')
weekly_evolution_comp(dt.stats, 'V100_plan', plots_dir, expression('V100 difference ('*Delta*'%)'), -35)

# == V102 =================================================== ===
weekly_evolution_box(filter(dt.stats, struct == 'CTV', method %in% c('Plan', 'None')),
                  'week.name', 'V102', 'stage', plots_dir, 'V102 (%)')
weekly_evolution_comp(dt.stats, 'V102_plan', plots_dir, expression('V102 difference ('*Delta*'%)'))

# == V105 =================================================== ===
weekly_evolution_box(filter(dt.stats, struct == 'CTV', method %in% c('Plan', 'None')),
                  'week.name', 'V105', 'stage', plots_dir, 'V105 (%)')
weekly_evolution_comp(dt.stats, 'V105_plan', plots_dir, expression('V105 difference ('*Delta*'%)'))

# == V107 =================================================== ===
weekly_evolution_box(filter(dt.stats, struct == 'CTV', method %in% c('Plan', 'None')),
                  'week.name', 'V107', 'stage', plots_dir, 'V107 (%)')
weekly_evolution_comp(dt.stats, 'V107_plan', plots_dir, expression('V107 difference ('*Delta*'%)'))


### Method comparison: target geometric -----------------------------

# ggplot(data = filter(dt.stats, struct == 'CTV',
#                      method %in% c("None", "Plan", "Geometric")),
#        aes(x = week.name, y = V98, group = interaction(case, week.name),
#            fill = constraint)) +
#     geom_boxplot() +
#     theme(legend.position = 'bottom') + 
#     labs(y = 'V98 (%)', x = element_blank())
# ggsave(paste0(plots_dir, '/target/','target_V98_geometric_all_weekly', '.pdf'),
#        width = 14, height = 8, units = "cm")
# 
# ggplot(data = filter(dt.stats, struct == 'CTV',
#                      method %in% c("None", "Plan", "Geometric"),
#                      need_adapt == TRUE),
#        aes(x = week.name, y = V98, group = interaction(case, week.name),
#            fill = constraint)) +
#     geom_boxplot() +
#     theme(legend.position = 'bottom') + 
#     labs(y = 'V98 (%)', x = element_blank())
# ggsave(paste0(plots_dir, '/target/','target_V98_geometric_need_adapt_weekly', '.pdf'),
#        width = 14, height = 8, units = "cm")


ggplot(data = filter(dt.stats, struct == 'CTV',
                     method %in% c("None", "Plan", "Geometric")),
       aes(x = short_const, y = V98, group = interaction(case, need_adapt),
           fill = need_adapt)) +
    geom_boxplot() +
    theme(legend.position = 'bottom') + 
    labs(y = 'V98 (%)', x = element_blank())
ggsave(paste0(plots_dir, '/target/','target_V98_geometric_need_adapt', '.pdf'),
       width = 16, height = 8, units = "cm")


ggplot(data = filter(dt.stats, struct == 'CTV',
                     method %in% c("None", "Plan", "Weights")),
       aes(x = short_const, y = V98, group = interaction(case, need_adapt),
           fill = need_adapt)) +
    geom_boxplot() +
    theme(legend.position = 'bottom') + 
    labs(y = 'V98 (%)', x = element_blank())
ggsave(paste0(plots_dir, '/target/','target_V98_weights_need_adapt', '.pdf'),
       width = 16, height = 8, units = "cm")

ggplot(data = filter(dt.stats, struct == 'CTV',
                     method %in% c("Plan", "Weights")),
       aes(x = short_const, y = D2_D98, group = interaction(case, need_adapt),
           fill = need_adapt)) +
    geom_boxplot() +
    theme(legend.position = 'bottom') + 
    labs(y = 'D2-D98 [Gy(RBE)]', x = element_blank())
ggsave(paste0(plots_dir, '/target/','norm_target_D2_D98_weights_need_adapt', '.pdf'),
       width = 16, height = 8, units = "cm")

ctv.averages <- dt.stats %>%
    filter(oar == FALSE) %>%
    group_by(patient, case, method, constraint, short_const) %>%
    # summarise_at(.vars = c("V98", "D2_D98"), .funs = c(mean = "mean", sd = "sd")) %>%
    summarise_at(.vars = c("V98", "D2_D98", "mean"), .funs = c(mean = "mean")) %>%
    group_by(patient) %>%
    mutate(mean_mean_ratio_plan = 100*(mean_mean - mean_mean[method == "Plan"])/mean_mean[method == "Plan"]) %>%
    mutate(V98_mean_ratio_plan = 100*(V98_mean - V98_mean[method == "Plan"])/V98_mean[method == "Plan"]) %>%
    mutate(D2_D98_mean_ratio_plan = 100*(D2_D98_mean - D2_D98_mean[method == "Plan"])/D2_D98_mean[method == "Plan"]) %>%
    mutate(mean_mean_ratio_frac = 100*(mean_mean - mean_mean[method == "None"])/mean_mean[method == "None"]) %>%
    mutate(V98_mean_ratio_frac = 100*(V98_mean - V98_mean[method == "None"])/V98_mean[method == "None"]) %>%
    mutate(D2_D98_mean_ratio_frac = 100*(D2_D98_mean - D2_D98_mean[method == "None"])/D2_D98_mean[method == "None"]) %>%
    mutate(mean_mean_ratio_free = 100*(mean_mean - mean_mean[constraint == "Free"])/mean_mean[constraint == "Free"]) %>%
    mutate(V98_mean_ratio_free = 100*(V98_mean - V98_mean[constraint == "Free"])/V98_mean[constraint == "Free"]) %>%
    mutate(D2_D98_mean_ratio_free = 100*(D2_D98_mean - D2_D98_mean[constraint == "Free"])/D2_D98_mean[constraint == "Free"])

ctv.averages <- dt.stats %>%
    filter(oar == FALSE) %>%
    group_by(patient) %>%
    select(patient, case, D2_D98) %>%
    spread(case, D2_D98)


  
ggplot(data = filter(ctv.averages, method == "Weights"),
       aes(x = short_const, y = D2_D98_mean_ratio_free,
           fill = constraint)) +
    geom_boxplot() +
    theme(legend.position = 'bottom') + 
    labs(y = 'D2-D98 w.r.t. plan [%]', x = element_blank())
ggsave(paste0(plots_dir, '/target/','norm_target_D2_D98_weights_need_adapt', '.pdf'),
       width = 16, height = 8, units = "cm")

ggplot(data = filter(dt.stats, struct == 'CTV',
                     method %in% c("Plan", "Weights")),
       aes(x = short_const, y = D2_D98,
           fill = constraint)) +
    stat_summary(fun.y = "mean", geom = "bar") +
    stat_summary(fun.data = mean_se, geom = "errorbar") +
    labs(y = 'Average D2-D98 [Gy(RBE)]', x = element_blank())
ggsave(paste0(plots_dir, '/target/','target_D2_D98_weights_comparison', '.pdf'),
       width = 16, height = 8, units = "cm")

### Method comparison: OAR -----------------------------

dt.stats <- filter(dt.stats, !(struct %in% c("R. Parotid", "Larynx")))

ggplot(data = filter(dt.stats, main_oar == TRUE,
                     method %in% c("None", "Plan", "Weights")),
       aes(x = short_const, y = mean, group = case,
           fill = constraint)) +
    geom_boxplot() +
    facet_wrap(~ struct, scales = "free") +
    theme(legend.position = 'none') +
    labs(y = 'Mean dose [Gy(RBE)]', x = element_blank())
ggsave(paste0(plots_dir, '/OARs/','oars_mean_weights', '.pdf'), width = 28, height = 14, units = "cm")

ggplot(data = filter(dt.stats, main_oar == TRUE,
                     method %in% c("None", "Plan", "Weights")),
       aes(x = short_const, y = max, group = case,
           fill = constraint)) +
    geom_boxplot() +
    facet_wrap(~ struct, scales = "free") +
    theme(legend.position = 'none') +
    labs(y = 'Max dose [Gy(RBE)]', x = element_blank())
ggsave(paste0(plots_dir, '/OARs/','oars_max_weights', '.pdf'), width = 28, height = 14, units = "cm")


ggplot(data = filter(dt.stats, main_oar == TRUE,
                     method %in% c("None", "Plan", "Weights"),
                     constraint %in% c("None", "Plan", "Free", "Isocenter")),
       aes(x = week.name, y = max, group = interaction(case, week.name),
           fill = constraint)) +
    geom_boxplot() +
    facet_wrap(~ struct, scales = "free") +
    theme(legend.position = 'none') +
    labs(y = 'Max dose [Gy(RBE)]', x = element_blank())
ggsave(paste0(plots_dir, '/OARs/','oars_max_weights', '.pdf'), width = 25, height = 14, units = "cm")

















comp_data <- dt.stats %>%
    filter(main_oar == TRUE, !constraint %in% c('Fixed', 'Isocenter', 'Range shifter')) %>%
    filter(!method %in% c('Robust', 'Geometric')) %>%
    group_by(patient, struct) %>%
    mutate(min = min/min[method == 'Plan']) %>%
    mutate(max = max/max[method == 'Plan']) %>%
    filter(method != 'Plan')

comp_data2 <- comp_data %>%
    ungroup() %>% group_by(patient, struct, week.no) %>%
    select(patient, patient.no, struct, week.no, method, mean) %>%
    spread(method, mean) %>%
    mutate(improvement = Weights - None)

cbPalette3 <- c('#000000')

ggplot(data = comp_data2, aes(x = struct, y = improvement, fill = struct)) +
    geom_boxplot() +
    geom_jitter(data = filter(comp_data2, patient.no == 8),
                aes(color = patient), alpha = 0.6, width = 0.1) +
    scale_fill_manual('Type', values = cbPalette, guide = FALSE) +
    scale_color_manual('Patient', values = cbPalette3) +
    scale_x_discrete(labels = c("L. Submandible gland" = "L. SubMGland",
                                "R. Submandible gland" = "R. SubMGland",
                                "Mandible" = "Mandible",
                                "Larynx" = "Larynx")) +
    labs(title = 'OARs: Mean dose reduction by adaptation', x = element_blank(), y = 'Dose Gy(RBE)') +
    theme(legend.position = c(0.15, 0.15))
ggsave('~/Desktop/methods_OARs.pdf', width = 14.0, height = 10.0, units = "cm")


ggplot(data = comp_data,
       aes(x = struct, y = mean, fill = constraint)) +
    geom_boxplot() +
    geom_boxplot(data = transform(comp_data, patient = 'All Patients')) +
    facet_wrap(~ patient, scales = 'free') +
    scale_x_discrete(labels = c("L. Submandible gland" = "L. SubMGland",
                                "R. Submandible gland" = "R. SubMGland",
                                "Mandible" = "Mandible",
                                "Larynx" = "Larynx")) +
    labs(y = 'Mean dose ratio to plan', x = element_blank())
ggsave(paste0(plots_dir, '/methods_OAR_mean.pdf'), width = 30.0, height = 15.0, units = "cm")

ggplot(data = comp_data,
       aes(x = struct, y = max, fill = constraint)) +
    geom_boxplot() +
    geom_boxplot(data = transform(comp_data, patient = 'All Patients')) +
    facet_wrap(~ patient, scales = 'free') +
    scale_x_discrete(labels = c("L. Submandible gland" = "L. SubMGland",
                                "R. Submandible gland" = "R. SubMGland",
                                "Mandible" = "Mandible",
                                "Larynx" = "Larynx")) +
    labs(y = 'Max dose ratio to plan', x = element_blank())
    ggsave(paste0(plots_dir, '/methods_OAR_max.pdf'), width = 30.0, height = 15.0, units = "cm")
    