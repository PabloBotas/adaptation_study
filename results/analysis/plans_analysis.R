library(ggplot2)
library(plyr)
library(dplyr)
library(data.table)
library(tidyr)

cbPalette <- c('#1f78b4', '#a6cee3', '#b2df8a', '#33a02c', '#fb9a99', '#000000', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#b15928')

setwd('/Users/pbotas/Desktop/gmoc/Desktop/pablo/adaptive_project/results/analysis')
source('plans_parse_data.R')
patients <- c('P01', 'P02', 'P03', 'P15')
location <- '../data/dvh'
target_dose <- 60
dose_fraction <- 60
plan_fractions <- 1
dvhs <- readData(patients, location, plan_fractions, FALSE)

### Get stats, normalize Plan dose and recalculate stats ------------------
source('plans_utils.R')
# print_summaries(dvhs)
dt.stats <- calculateStats(dvhs, dose_fraction)
temp <- normalize_plan(dvhs, dt.stats, 'V95', dose_fraction, 100/dose_fraction)
#temp <- normalize_plan(dvhs, dt.stats, 'V95', dose_fraction, 100/dose_fraction)
dvhs <- temp$dvhs
dt.stats <- temp$stats
rm(temp)

levels(dt.stats$patient) <- c("Pat. 1", "Pat. 2", "Pat. 3", "Pat. 4")
levels(dvhs$patient) <- c("Pat. 1", "Pat. 2", "Pat. 3", "Pat. 4")

### CTV Target Boxplots ------------------------------------------------
# Mean dose per week CTV
ggplot(filter(dt.stats, struct == 'CTV', case != 'Week 7'), aes(x = case, y = mean, fill = case)) +
    geom_boxplot(alpha = 1) +
    geom_smooth(data = filter(dt.stats, struct == 'CTV', case != 'Plan', case != 'Week 7'), method = lm, se = FALSE, aes(group = 1), color = alpha("black", 0.5)) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0), legend.position = 'none') +
    labs(y = "Mean Dose (Gy)")
ggsave('plots_plans/target/Mean.CTV.week.pdf', width = 12, height = 6, units = "cm")
ggsave('plots_plans/target/Mean.CTV.week.png', width = 12, height = 6, units = "cm")
# Max dose per week CTV
ggplot(filter(dt.stats, struct == 'CTV', case != 'Week 7'), aes(x = case, y = max, fill = case)) +
    geom_boxplot(alpha = 1) +
    geom_smooth(data = filter(dt.stats, struct == 'CTV', case != 'Plan', case != 'Week 7'), method = lm, se = FALSE, aes(group = 1), color = alpha("black", 0.5)) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90), legend.position = 'none') +
    labs(y = "Max Dose (Gy(RBE))")
ggsave('plots_plans/target/Max.CTV.week.pdf', width = 12, height = 8, units = "cm")
ggsave('plots_plans/target/Max.CTV.week.png', width = 12, height = 8, units = "cm")
# D98 per week CTV
ggplot(filter(dt.stats, struct == 'CTV', case != 'Week 7'), aes(x = case, y = D98, fill = case)) +
    geom_hline(aes(yintercept = 0.95*dose_fraction, color = 'red', linetype = "dashed")) +
    geom_boxplot(alpha = 1) +
    geom_smooth(data = filter(dt.stats, struct == 'CTV', case != 'Plan', case != 'Week 7'), method = lm, se = FALSE, aes(group = 1), color = alpha("black", 0.5)) +
    geom_text(aes(x = Inf, y = Inf), size = 3, color = 'black', label = 'Line at 95% prescription', hjust = 1, vjust = 1) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90), legend.position = 'none') +
    labs(y = "D98 (Gy)")
ggsave('plots_plans/target/D98.CTV.week.pdf', width = 12, height = 8, units = "cm")
ggsave('plots_plans/target/D98.CTV.week.png', width = 12, height = 8, units = "cm")
# D95 per week CTV
ggplot(filter(dt.stats, struct == 'CTV', case != 'Week 7'), aes(x = case, y = D95, fill = case)) +
    geom_hline(aes(yintercept = 0.95*dose_fraction, color = 'red', linetype = "dashed")) +
    geom_boxplot(alpha = 1) +
    geom_smooth(data = filter(dt.stats, struct == 'CTV', case != 'Plan', case != 'Week 7'), method = lm, se = FALSE, aes(group = 1), color = alpha("black", 0.5)) +
    geom_text(aes(x = Inf, y = Inf), size = 3, color = 'black', label = 'Line at 95% prescription', hjust = 1, vjust = 1) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90), legend.position = 'none') +
    labs(y = "D95 (Gy)")
ggsave('plots_plans/target/D95.CTV.week.pdf', width = 12, height = 8, units = "cm")
ggsave('plots_plans/target/D95.CTV.week.png', width = 12, height = 8, units = "cm")
# V95 per week CTV
ggplot(filter(dt.stats, struct == 'CTV', case != 'Week 7'), aes(x = case, y = V95, fill = case)) +
    geom_hline(aes(yintercept = 98, color = 'red', linetype = "dashed")) +
    geom_boxplot(alpha = 1) +
    geom_text(aes(x = Inf, y = Inf), size = 3, color = 'black', label = 'Line at 98% of target volume', hjust = 1, vjust = 1) +
    geom_smooth(data = filter(dt.stats, struct == 'CTV', case != 'Plan', case != 'Week 7'), method = lm, se = FALSE, aes(group = 1), color = alpha("black", 0.5)) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90), legend.position = 'none') +
    labs(y = "V95 (%)")
ggsave('plots_plans/target/V95.CTV.week.pdf', width = 12, height = 8, units = "cm")
ggsave('plots_plans/target/V95.CTV.week.png', width = 12, height = 8, units = "cm")
# V107 per week CTV
ggplot(filter(dt.stats, struct == 'CTV', case != 'Week 7'), aes(x = case, y = V107, fill = case)) +
    geom_hline(aes(yintercept = 5, color = 'red', linetype = "dashed")) +
    geom_boxplot(alpha = 1) +
    geom_text(aes(x = Inf, y = -Inf), size = 3, color = 'black', label = 'Line at 5% of target volume', hjust = 1, vjust = -1) +
    geom_smooth(data = filter(dt.stats, struct == 'CTV', case != 'Plan', case != 'Week 7'), method = lm, se = FALSE, aes(group = 1), color = alpha("black", 0.5)) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90), legend.position = 'none') +
    labs(y = "V107 (%)")
ggsave('plots_plans/target/V107.CTV.week.pdf', width = 12, height = 8, units = "cm")
ggsave('plots_plans/target/V107.CTV.week.png', width = 12, height = 8, units = "cm")
# D5-D95 per week CTV
ggplot(filter(dt.stats, struct == 'CTV', case != 'Week 7'), aes(x = case, y = D5_D95, fill = case)) +
    geom_hline(aes(yintercept = 5, color = 'red', linetype = "dashed")) +
    geom_boxplot(alpha = 1) +
    geom_text(aes(x = Inf, y = -Inf), size = 3, color = 'black', label = 'Line at 5 Gy', hjust = 1, vjust = -1) +
    geom_smooth(data = filter(dt.stats, struct == 'CTV', case != 'Plan', case != 'Week 7'), method = lm, se = FALSE, aes(group = 1), color = alpha("black", 0.5)) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90), legend.position = 'none') +
    labs(y = "D5-D95 (Gy(RBE))")
ggsave('plots_plans/target/D5_D95.CTV.week.pdf', width = 12, height = 8, units = "cm")
ggsave('plots_plans/target/D5_D95.CTV.week.png', width = 12, height = 8, units = "cm")

### CTV Target per patient plots_plans ------------------------------------------
# Mean dose per week CTV
ggplot(filter(dt.stats, struct == 'CTV', case != 'Week 7'), aes(x = case, y = mean, group = patient, color = patient, shape = patient)) +
    geom_line(alpha = 0.2) + geom_point(alpha = 0.75, size = 2) +
    geom_smooth(method = lm, se = FALSE) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90)) +
    labs(y = "Mean Dose (Gy)")
ggsave('plots_plans/target/pat.Mean.CTV.week.pdf', width = 12, height = 8, units = "cm")
ggsave('plots_plans/target/pat.Mean.CTV.week.png', width = 12, height = 8, units = "cm")
# Max dose per week CTV
ggplot(filter(dt.stats, struct == 'CTV', case != 'Week 7'), aes(x = case, y = max, group = patient, color = patient, shape = patient)) +
    geom_line(alpha = 0.2) + geom_point(alpha = 0.75, size = 2) +
    geom_smooth(method = lm, se = FALSE) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90)) +
    labs(y = "Max Dose (Gy(RBE))")
ggsave('plots_plans/target/pat.Max.CTV.week.pdf', width = 12, height = 8, units = "cm")
ggsave('plots_plans/target/pat.Max.CTV.week.png', width = 12, height = 8, units = "cm")
# D98 per week CTV
ggplot(filter(dt.stats, struct == 'CTV', case != 'Week 7'), aes(x = case, y = D98, group = patient, color = patient, shape = patient)) +
    geom_hline(aes(yintercept = 0.95*dose_fraction), color = 'red', linetype = "dashed") +
    geom_line(alpha = 0.5) + geom_point(size = 2) +
    geom_text(aes(x = Inf, y = Inf), size = 3, color = 'black', label = 'Line at 95% prescription', hjust = 1, vjust = 1) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90)) +
    labs(y = "D98 (Gy)")
ggsave('plots_plans/target/pat.D98.CTV.week.pdf', width = 12, height = 8, units = "cm")
ggsave('plots_plans/target/pat.D98.CTV.week.png', width = 12, height = 8, units = "cm")
# D95 per week CTV
ggplot(filter(dt.stats, struct == 'CTV', case != 'Week 7'), aes(x = case, y = D95, group = patient, color = patient, shape = patient)) +
    geom_hline(aes(yintercept = 0.95*dose_fraction), color = 'red', linetype = "dashed") +
    geom_line(alpha = 0.5) + geom_point(size = 2) +
    geom_text(aes(x = Inf, y = Inf), size = 3, color = 'black', label = 'Line at 95% prescription', hjust = 1, vjust = 1) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90)) +
    labs(y = "D95 (Gy)")
ggsave('plots_plans/target/pat.D95.CTV.week.pdf', width = 12, height = 8, units = "cm")
ggsave('plots_plans/target/pat.D95.CTV.week.png', width = 12, height = 8, units = "cm")
# V95 per week CTV
ggplot(filter(dt.stats, struct == 'CTV', case != 'Week 7'), aes(x = case, y = V95, group = patient, color = patient, shape = patient)) +
    geom_hline(aes(yintercept = 98), color = 'red', linetype = "dashed") +
    geom_line(alpha = 0.5) + geom_point(size = 2) +
    geom_text(aes(x = Inf, y = Inf), size = 3, color = 'black', label = 'Line at 98% of target volume', hjust = 1, vjust = 1) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90)) +
    labs(y = "V95 (%)")
ggsave('plots_plans/target/pat.V95.CTV.week.pdf', width = 12, height = 8, units = "cm")
ggsave('plots_plans/target/pat.V95.CTV.week.png', width = 12, height = 8, units = "cm")
# V107 per week CTV
ggplot(filter(dt.stats, struct == 'CTV', case != 'Week 7'), aes(x = case, y = V107, group = patient, color = patient, shape = patient)) +
    geom_hline(aes(yintercept = 5), color = 'red', linetype = "dashed") +
    geom_line(alpha = 0.5) + geom_point(size = 2) +
    geom_text(aes(x = Inf, y = -Inf), size = 3, color = 'black', label = 'Line at 5% of target volume', hjust = 1, vjust = -1) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90)) +
    labs(y = "V107 (%)")
ggsave('plots_plans/target/pat.V107.CTV.week.pdf', width = 12, height = 8, units = "cm")
ggsave('plots_plans/target/pat.V107.CTV.week.png', width = 12, height = 8, units = "cm")
# D5-D95 per week CTV
ggplot(filter(dt.stats, struct == 'CTV', case != 'Week 7'), aes(x = case, y = D5_D95, group = patient, color = patient, shape = patient)) +
    geom_hline(aes(yintercept = 5), color = 'red', linetype = "dashed") +
    geom_line(alpha = 0.5) + geom_point(size = 2) +
    geom_text(aes(x = Inf, y = -Inf), size = 3, color = 'black', label = 'Line at 5 Gy', hjust = 1, vjust = -1) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90)) +
    labs(y = "D5-D95 (Gy(RBE))")
ggsave('plots_plans/target/pat.D5_D95.CTV.week.pdf', width = 12, height = 8, units = "cm")
ggsave('plots_plans/target/pat.D5_D95.CTV.week.png', width = 12, height = 8, units = "cm")



### OARs Boxplots ------------------------------------------------

# Mean dose per week and structures
p <- ggplot(filter(dt.stats, oar == TRUE, case != 'Week 7'), aes(x = case, y = mean, fill = case)) +
    geom_boxplot(alpha = 1) +
    facet_wrap(~struct, scales = 'free_y') +
    geom_smooth(data = filter(dt.stats, oar == TRUE, case != 'Plan', case != 'Week 7'), method = lm, se = FALSE, aes(group = 1), color = alpha("black", 0.5)) +
    theme(axis.text.x = element_text(angle = 90), legend.position = 'none') +
    labs(x = "Treatment Week", y = "Mean Dose (Gy)")
ggsave('plots_plans/OARs/Mean.OARs.week.pdf', width = 18, height = 12, units = "cm")
ggsave('plots_plans/OARs/Mean.OARs.week.png', width = 18, height = 12, units = "cm")
meanD_stats <- ggplot_build(p)$data[[1]]
# Max dose per week and structures
p <- ggplot(filter(dt.stats, oar == TRUE, case != 'Week 7'), aes(x = case, y = max, fill = case)) +
    geom_boxplot(alpha = 1) +
    facet_wrap(~struct, scales = 'free_y') +
    geom_smooth(data = filter(dt.stats, oar == TRUE, case != 'Plan', case != 'Week 7'), method = lm, se = FALSE, aes(group = 1), color = alpha("black", 0.5)) +
    theme(axis.text.x = element_text(angle = 90), legend.position = 'none') +
    labs(x = "Treatment Week", y = "Max Dose (Gy)")
ggsave('plots_plans/OARs/Max.OARs.week.pdf', width = 18, height = 12, units = "cm")
ggsave('plots_plans/OARs/Max.OARs.week.png', width = 18, height = 12, units = "cm")
maxD_stats <- ggplot_build(p)$data[[1]]
# Mean dose selected structures
ggplot(filter(dt.stats, struct %in% c('Spinal Cord', 'Larynx'), case != 'Week 7'), aes(x = case, y = mean, fill = case)) +
    geom_boxplot(alpha = 1) +
    facet_wrap(~struct, scales = 'free_y') +
    geom_smooth(data = filter(dt.stats, struct %in% c('Spinal Cord', 'Larynx', 'Oral cavity'), case != 'Plan', case != 'Week 7'), method = lm, se = FALSE, aes(group = 1), color = alpha("black", 0.5)) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0), legend.position = 'none') +
    labs(y = "Mean Dose (Gy)")
ggsave('plots_plans/OARs/Mean.OARs.reduced.week.pdf', width = 20, height = 7, units = "cm")
ggsave('plots_plans/OARs/Mean.OARs.reduced.week.png', width = 18, height = 7, units = "cm")
# Max dose selected structures
ggplot(filter(dt.stats, struct %in% c('Spinal Cord', 'Larynx', 'Oral cavity'), case != 'Week 7'), aes(x = case, y = max, fill = case)) +
    geom_boxplot(alpha = 1) +
    facet_wrap(~struct, scales = 'free_y') +
    geom_smooth(data = filter(dt.stats, struct %in% c('Spinal Cord', 'Larynx', 'Oral cavity'), case != 'Plan', case != 'Week 7'), method = lm, se = FALSE, aes(group = 1), color = alpha("black", 0.5)) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90), legend.position = 'none') +
    labs(y = "Max Dose (Gy)")
ggsave('plots_plans/OARs/Max.OARs.reduced.week.pdf', width = 18, height = 7, units = "cm")
ggsave('plots_plans/OARs/Max.OARs.reduced.week.png', width = 18, height = 7, units = "cm")

### OARs per patient ------------------------------------------------
# Mean dose per week and structures
p <- ggplot(filter(dt.stats, oar == TRUE, case != 'Week 7'), aes(x = case, y = mean, group = patient, color = patient, shape = patient)) +
    geom_line(alpha = 0.2) + geom_point(alpha = 0.75, size = 2) +
    facet_wrap(~struct, scales = 'free_y') +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "Treatment Week", y = "Mean Dose (Gy)")
ggsave('plots_plans/OARs/pat.Mean.OARs.week.pdf', width = 18, height = 12, units = "cm")
ggsave('plots_plans/OARs/pat.Mean.OARs.week.png', width = 18, height = 12, units = "cm")
meanD_stats <- ggplot_build(p)$data[[1]]
# Max dose per week and structures
p <- ggplot(filter(dt.stats, oar == TRUE, case != 'Week 7'), aes(x = case, y = max, group = patient, color = patient, shape = patient)) +
    geom_line(alpha = 0.2) + geom_point(alpha = 0.75, size = 2) +
    facet_wrap(~struct, scales = 'free_y') +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "Treatment Week", y = "Max Dose (Gy)")
ggsave('plots_plans/OARs/pat.Max.OARs.week.pdf', width = 18, height = 12, units = "cm")
ggsave('plots_plans/OARs/pat.Max.OARs.week.png', width = 18, height = 12, units = "cm")
maxD_stats <- ggplot_build(p)$data[[1]]
# Mean dose selected structures
ggplot(filter(dt.stats, struct %in% c('Spinal Cord', 'Larynx', 'Oral cavity'), case != 'Week 7'), aes(x = case, y = mean, group = patient, color = patient, shape = patient)) +
    geom_line(alpha = 0.2) + geom_point(alpha = 0.75, size = 2) +
    facet_wrap(~struct, scales = 'free_y') +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90)) +
    labs(y = "Mean Dose (Gy)")
ggsave('plots_plans/OARs/pat.Mean.OARs.reduced.week.pdf', width = 18, height = 7, units = "cm")
ggsave('plots_plans/OARs/pat.Mean.OARs.reduced.week.png', width = 18, height = 7, units = "cm")
# Max dose selected structures
ggplot(filter(dt.stats, struct %in% c('Spinal Cord', 'Larynx', 'Oral cavity'), case != 'Week 7'), aes(x = case, y = max, group = patient, color = patient, shape = patient)) +
    geom_line(alpha = 0.2) + geom_point(alpha = 0.75, size = 2) +
    facet_wrap(~struct, scales = 'free_y') +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90)) +
    labs(y = "Max Dose (Gy)")
ggsave('plots_plans/OARs/pat.Max.OARs.reduced.week.pdf', width = 18, height = 7, units = "cm")
ggsave('plots_plans/OARs/pat.Max.OARs.reduced.week.png', width = 18, height = 7, units = "cm")

### Plot DVHs --------------------------------
ggplot(filter(dvhs, !(checkpoints %in% 'Other')),
       aes(x = dose, y = vol, color = struct, alpha = checkpoints)) +
    geom_path(aes(linetype = checkpoints)) +
    facet_wrap(~ patient) +
    scale_x_continuous(limits = c(0, dose_fraction*1.3), breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(limits = c(0, 100), breaks = scales::pretty_breaks(n = 10)) +
    scale_colour_manual('Contour', values = cbPalette) +
    scale_alpha_manual(values = c(1,0.75,0.75), guide = FALSE) +
    scale_linetype_manual('Fraction', values = c('solid', 'dotted', 'dashed')) +
    theme_light() +
    labs(x = "Dose (Gy)", y = "Contour volume (%)")
ggsave('plots_plans/DVHs.pdf', width = 26, height = 16, units = "cm")


ggplot(filter(dvhs, patient == 'Pat. 3', !(checkpoints %in% 'Other')),
       aes(x = dose, y = vol, color = struct, linetype = case)) +
    geom_path() +
    scale_colour_manual('Contour', values = cbPalette) +
    #scale_linetype_manual('Fraction', values = c('solid', 'dotted', 'dashed')) +
    theme_light() +
    labs(x = "Dose (Gy)", y = "Contour volume (%)")
ggsave('plots_plans/DVH.Pat. 3.pdf', width = 18, height = 10, units = "cm")
ggsave('plots_plans/DVH.Pat. 3.png', width = 18, height = 10, units = "cm")

ggplot(filter(dvhs, patient == 'Pat. 3', (week.no < 2 | week.no == 6)),
       aes(x = dose, y = vol, color = struct, linetype = case)) +
    geom_path() +
    scale_colour_manual('Contour', values = cbPalette) +
    theme_light() +
    labs(x = "Dose (Gy)", y = "Contour volume (%)")

for (pat in levels(dvhs$patient))
{
    summary(dvhs[[pat]])
    ggplot(filter(dvhs, patient == pat, (week.no < 2 | week.no == 6), struct != 'R. Submandible gland', struct != 'Spinal Cord'),
           aes(x = dose, y = vol, color = struct)) +
        geom_path(aes(linetype = case)) +
        scale_colour_manual('Contour', values = cbPalette) +
        scale_linetype_manual('Fraction', values = c('solid', 'dotted', 'dashed')) +
        theme_light() +
        guides(fill=guide_legend(
            keywidth=0.1,
            keyheight=0.01,
            default.unit="inch")
        )+
        labs(x = "Dose (Gy)", y = "Contour volume (%)")
    ggsave(paste('plots_plans/DVH.', pat, '.pdf', sep = ''), width = 18, height = 9, units = "cm")
    ggsave(paste('plots_plans/DVH.', pat, '.png', sep = ''), width = 18, height = 7, units = "cm")
}
