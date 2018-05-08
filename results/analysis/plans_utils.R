library(dplyr)
library(data.table)

print_summaries <- function(l)
{
    for (i in l)
    {
        print(summary(i))
    }
}

maxDose <- function(dose, vol)
{
    return(max(dose[vol != 0]))
}
minDose <- function(dose, vol)
{
    m <- ifelse(vol[length(vol)] == 100, 0, min(dose[vol != 100]))
    return(m)
}
meanDose <- function(dose, vol)
{
    vol.dif <- diff(vol)
    vol.dif <- c(0, vol.dif)
    vol.dif <- vol.dif/sum(vol.dif)
    return(sum(dose*vol.dif))
}
minVolumeAtDose <- function(dose, vol, D, target_dose = NA)
{
    if( ! is.na(target_dose) )
        dose <- 100*dose/target_dose
    dif <- abs(dose - D)
    idx <- which(dif == min(dif))
    return(min(vol[idx]))
}
maxDoseToVolumePctg <- function(dose, vol, V, target_dose = NA)
{
    if( ! is.na(target_dose) )
        dose <- 100*dose/target_dose
    dif <- abs(vol - V)
    idx <- which(dif == min(dif))
    return(max(dose[idx]))
}

calculateStats <- function(dvhs_input, target_dose)
{
    dt.stats <- dvhs_input %>%
        group_by(patient, case, struct, oar) %>%
        summarise(min  = minDose(dose, vol),
                  max  = maxDose(dose, vol),
                  mean = meanDose(dose, vol),
                  D2   = maxDoseToVolumePctg(dose, vol, 2),
                  D5   = maxDoseToVolumePctg(dose, vol, 5),
                  D95  = maxDoseToVolumePctg(dose, vol, 95),
                  D98  = maxDoseToVolumePctg(dose, vol, 98),
                  D99  = maxDoseToVolumePctg(dose, vol, 99),
                  V20  = minVolumeAtDose(dose, vol, 20, target_dose),
                  V90  = minVolumeAtDose(dose, vol, 90, target_dose),
                  V95  = minVolumeAtDose(dose, vol, 95, target_dose),
                  V98  = minVolumeAtDose(dose, vol, 98, target_dose),
                  V100 = minVolumeAtDose(dose, vol, 100, target_dose),
                  V107 = minVolumeAtDose(dose, vol, 107, target_dose),
                  V110 = minVolumeAtDose(dose, vol, 110, target_dose),
                  V115 = minVolumeAtDose(dose, vol, 115, target_dose),
                  V120 = minVolumeAtDose(dose, vol, 120, target_dose)
        ) %>%
        mutate(D5_D95 = D5 - D95, D2_D98 = D2 - D98) %>%
        ungroup() %>% as.data.table()
        

    return(dt.stats)
}

calculateSingleStat <- function(dvhs_input, target_dose, stat)
{
    if (grepl('_', stat)) {
        splt <- strsplit(stat, '_')
        splt1 <- splt[[1]][1]
        splt2 <- splt[[1]][2]
        dt1 <- calculateSingleStat(dvhs_input, target_dose, splt1)
        dt2 <- calculateSingleStat(dvhs_input, target_dose, splt2)
        dt3 <- merge(dt1, dt2)
        dt <- dt3 %>%
            mutate_(stat = paste(splt1, splt2, sep = '-')) %>%
            select(-one_of(c(splt1, splt2)))

    } else if (grepl('V', stat)) {
        val <- as.numeric(gsub("[[:alpha:]]", "", stat))
        dt <- dvhs_input %>% group_by(patient, case, struct, oar) %>%
            summarise(stat = minVolumeAtDose(dose, vol, val, target_dose)) %>%
            ungroup() %>% as.data.table()
        
    } else if (grepl('D', stat)) {
        val <- as.numeric(gsub("[[:alpha:]]", "", stat))
        dt <- dvhs_input %>% group_by(patient, case, struct, oar) %>%
            summarise(stat = maxDoseToVolumePctg(dose, vol, val)) %>%
            ungroup() %>% as.data.table()
        
    } else if (grepl('mean', stat)) {
        dt <- dvhs_input %>% group_by(patient, case, struct, oar) %>%
            summarise(mean = meanDose(dose, vol)) %>%
            ungroup() %>% as.data.table()
        
    } else if (grepl('max', stat)) {
        dt <- dvhs_input %>% group_by(patient, case, struct, oar) %>%
            summarise(max = maxDose(dose, vol)) %>%
            ungroup() %>% as.data.table()
        
    } else if (grepl('min', stat)) {
        dt <- dvhs_input %>% group_by(patient, case, struct, oar) %>%
            summarise(min = minDose(dose, vol)) %>%
            ungroup() %>% as.data.table()
    }
    names(dt)[names(dt) == 'stat'] <- stat

    return(dt)
}

normalize_plan <- function(dvhs_input, dt.stats, col, dose_fraction, val = 1)
{
    dt.norm <- dt.stats %>% ungroup() %>%
        filter(case == 'Plan', struct == 'CTV') %>%
        select_('patient', col) %>%
        mutate_(scale_factor = paste(val*dose_fraction, '/', col)) %>%
        mutate(patient_char = as.character(patient)) %>%
        select(-patient) %>% as.data.table()
    
    dvhs_input <- dvhs_input %>%
        group_by(patient) %>%
        mutate(dose = dose*dt.norm[patient_char == unique(as.character(patient)), scale_factor]) %>%
        ungroup()
    
    dt.stats <- calculateStats(dvhs_input, dose_fraction)
    
    return(list(dvhs = dvhs_input, stats = dt.stats))
}

# normalize_plan <- function(dvhs_input, dt.stats_input, col, dose_fraction, scale = 1)
# {
#     ## This function returns a list containing normalized versions of the dvhs_input and dt.stats_input.
#     ## Parameters explained:
#     ##     - col           : Magnitude used for normalization. It is a column
#     ##                       name (string) of dt.status_input.
#     ##     - dose_fraction : The dose per fraction in dt.stats_input units.
#     ##     - scale         : Value to which we want to normalize.
#     ## 
#     ## col can represent a dose (min, max, mean, D98...) or a volume percentage
#     ## (V90, V95, V100...). The logic of the function will interpret the col 
#     ## type and decide how the normalization should be perfored. This will
#     ## change the meaning of the scale factor.
#     ##
#     ## If col is a dose, then the scaling factor will be the scale parameter
#     ## times the dose per fraction parameter. The scale parameter is then a
#     ## a multiplier of the dose per fraction. Any number (num) can be set as
#     ## reference if scale = num/dose_fraction is set.
#     ##
#     ## If col is a volume percentage, for example V95, then scale is a dose
#     ## column in dt.stats_input to be user as reference. In this case, 95 will
#     ## be the percentage of the dose per fraction to which the aforementioned
#     ## column will be set.
#     ##
#     ## Examples:
#     ##  1) Set D95 to 95% of the dose per fraction
#     ##         normalize_plan(dvhs, dt.stats, 'D95', dose_fraction, 0.95)
#     ##  2) Set V95 to 100% of the volume
#     ##         normalize_plan(dvhs, dt.stats, 'V95', dose_fraction, 'D100')
#     ##         normalize_plan(dvhs, dt.stats, 'V95', dose_fraction, 'min')
#     ##  3) Set V99 to 99% of the volume
#     ##         normalize_plan(dvhs, dt.stats, 'V99', dose_fraction, 'D99')
#     ##  4) Set max dose to 314 (units of dose of dt.stats_input)
#     ##         normalize_plan(dvhs, dt.stats, 'max', dose_fraction, 314/dose_fraction)
#     
#     if (!grepl('V', col)) {
#         ## Because we are calculating the scaling factor with doses, we can use this approach
#         dt.norm <- dt.stats_input %>% ungroup() %>%
#             filter(case == 'Plan', struct == 'CTV') %>%
#             select_('patient', col) %>%
#             mutate_(scale_factor = paste(scale*dose_fraction, '/', col)) %>%
#             mutate(patient_char = as.character(patient)) %>%
#             select(-patient) %>% as.data.table()
#     } else {
#         ## If the stat does not exist...add it!!
#         if ( !(scale %in% names(dt.stats_input)) ) {
#             temp <- calculateSingleStat(dvhs_input, dose_fraction, scale)
#             dt.stats_input <- merge(dt.stats_input, temp)
#         }
#         
#         ## If setting V95 (for example), the scaling factor must be calculated from min and set it to 95% dose
#         prctg <- as.numeric(gsub("[[:alpha:]]", "", col))/100
#         reference <- ifelse(scale == 'D100', 'min', scale)
# 
#         dt.norm <- dt.stats_input %>% ungroup() %>%
#             filter(case == 'Plan', struct == 'CTV') %>%
#             select_('patient', reference) %>%
#             mutate_(scale_factor = paste(prctg*dose_fraction, '/', reference)) %>%
#             mutate(patient_char = as.character(patient)) %>%
#             select(-patient) %>% as.data.table()
#     }
#     
#     dvhs <- dvhs_input %>%
#         group_by(patient) %>%
#         mutate(dose = dose*dt.norm[patient_char == unique(as.character(patient)), scale_factor]) %>%
#         ungroup() %>% as.data.table()
#     
#     dt.stats <- calculateStats(dvhs, dose_fraction)
#     
#     return(list(dvhs = dvhs, stats = dt.stats))
# }
