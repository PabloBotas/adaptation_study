
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
volumeAtDose <- function(dose, vol, D, target_dose = NA)
{
    if ( !is.na(target_dose) )
        dose <- 100*dose/target_dose
    dif <- abs(dose - D)
    idx <- which(dif == min(dif))
    return(min(vol[idx]))
}
doseAtVolume <- function(dose, vol, V, target_dose = NA)
{
    if ( !is.na(target_dose) )
        dose <- 100*dose/target_dose
    dif <- abs(vol - V)
    idx <- which(dif == min(dif))
    return(max(dose[idx]))
}


calculateStats <- function(dvhs_input, prescription, norm_expression)
{
    if (!missing(norm_expression) && norm_expression != '' && !is.null(norm_expression)) {
        # Parse norm expression
        norm_expression = gsub(" ", "", norm_expression)
        norm_expression = strsplit(norm_expression, '=')
        lhs = norm_expression[[1]][1]
        rhs = norm_expression[[1]][2]
        norm_par = tolower(lhs)
        norm_val = as.numeric(gsub('v|d', '', norm_par))
        objective_pct = as.numeric(rhs)
        print(paste('Normalizing', norm_par, 'to', objective_pct, '%'))
        
        # Deduce variables to normalize
        if (startsWith(norm_par, "v")) {
            target_norm = prescription*norm_val/100
            norm.stats <- dvhs_input %>%
                filter(struct == 'CTV') %>%
                group_by(patient.no, week.no, case) %>%
                summarise(orig_val = doseAtVolume(dose, vol, objective_pct)) %>%
                mutate(scale_norm = target_norm/orig_val,
                       case = as.character(case)) %>%
                ungroup()
        } else if (startsWith(norm_par, "d")) {
            ref_vol = norm_val
            ref_dose_pct = objective_pct
            ref_dose = prescription*ref_dose_pct/100
            norm.stats <- dvhs_input %>%
                filter(struct == 'CTV') %>%
                group_by(patient.no, week.no, case) %>%
                summarise(orig_val = 100*doseAtVolume(dose, vol, ref_vol)/prescription) %>%
                mutate(scale_norm = ref_dose_pct/orig_val,
                       case = as.character(case)) %>%
                ungroup()
        }
    } else {
        norm.stats <- dvhs_input %>%
            filter(struct == 'CTV') %>%
            group_by(patient.no, week.no, case) %>%
            summarise(scale_norm = 1) %>%
            ungroup()
        lhs <- 'V98'
    }
    names(norm.stats) <- recode(names(norm.stats),
                                patient.no = 'patno',
                                week.no    = 'wno',
                                case       = 'c',
                                orig_val   = 'orig_val',
                                scale_norm = 'scl')

    ## Apply dose normalization
    dvhs_input <- dvhs_input %>%
        group_by(patient.no, week.no, case) %>%
        mutate(scale_norm = filter(norm.stats,
                                  (patno == unique(patient.no)) &
                                  (wno == unique(week.no)) &
                                  (c == unique(as.character(case))))$scl) %>%
        ungroup() %>%
        mutate(dose_original = dose) %>%
        mutate(dose = dose*scale_norm, dose_pct = dose_pct*scale_norm)
    
    ## Calculate stats
    dt.stats <- dvhs_input %>%
        group_by(patient, patient.no, pat.week, patient_orig, struct, main_oar, oar, week.no,
                 week.name, case, method, constraint, short_const, scale_norm, stage, checkpoints) %>%
        mutate(dose = scale_norm*dose) %>%
        summarise(min  = 100*minDose(dose, vol)/prescription,
                  max  = 100*maxDose(dose, vol)/prescription,
                  mean = 100*meanDose(dose, vol)/prescription,
                  D1   = 100*doseAtVolume(dose, vol, 1)/prescription,
                  D2   = 100*doseAtVolume(dose, vol, 2)/prescription,
                  D5   = 100*doseAtVolume(dose, vol, 5)/prescription,
                  D50  = 100*doseAtVolume(dose, vol, 50)/prescription,
                  D95  = 100*doseAtVolume(dose, vol, 95)/prescription,
                  D98  = 100*doseAtVolume(dose, vol, 98)/prescription,
                  D99  = 100*doseAtVolume(dose, vol, 99)/prescription,
                  V20  = volumeAtDose(dose, vol, 20, target_dose),
                  V90  = volumeAtDose(dose, vol, 90, target_dose),
                  V95  = volumeAtDose(dose, vol, 95, target_dose),
                  V98  = volumeAtDose(dose, vol, 98, target_dose),
                  V99  = volumeAtDose(dose, vol, 99, target_dose),
                  V100 = volumeAtDose(dose, vol, 100, target_dose),
                  V102 = volumeAtDose(dose, vol, 102, target_dose),
                  V105 = volumeAtDose(dose, vol, 105, target_dose),
                  V107 = volumeAtDose(dose, vol, 107, target_dose),
                  V110 = volumeAtDose(dose, vol, 110, target_dose),
                  V115 = volumeAtDose(dose, vol, 115, target_dose),
                  V120 = volumeAtDose(dose, vol, 120, target_dose)
            ) %>%
        mutate(D5_D95 = D5 - D95,
               D2_D98 = D2 - D98,
               D1_D99 = D1 - D99,
               Dmax_Dmin = max - min) %>%
        mutate(auc = DescTools::AUC(x = c(D2, D5, D50, D95, D98), y = c(2, 5, 50, 95, 98),
                         method = "trapezoid")) %>%
        mutate(ref_auc = 98*(D2_D98)) %>%
        mutate(shape = 1 - auc/ref_auc) %>%
        ungroup()
    
    plan.stats <- dt.stats %>%
        filter(week.no == 0) %>%
        rename(PATIENT = patient, STRUCT = struct)
    
    temp <- dt.stats %>%
        group_by(patient, struct) %>%
        mutate_at(.vars = vars(min:Dmax_Dmin),
                  .funs = funs(. - filter(plan.stats, PATIENT == patient, STRUCT == struct)$.)) %>%
        ungroup() %>%
        rename_at(.vars = vars(min:Dmax_Dmin),
                  .funs = funs(paste0(., '_plan')))
    
    dt.stats <- bind_cols(dt.stats, select(temp, min_plan:Dmax_Dmin_plan))
    
    print(select(filter(dt.stats, struct == 'CTV' & week.no == 0),
                 patient, method, struct, scale_norm, `lhs`))

    return(list("stats" = data.table(dt.stats), "mod_dvhs" = data.table(dvhs_input)))
}
