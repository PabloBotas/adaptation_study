library(data.table)

organize_structures <- function(dt)
{
    dt$struct <- recode(dt$struct,
                        gtv               = 'GTV',
                        ctv               = 'CTV',
                        esophagus         = 'Esophagus',
                        spinalcord        = 'Spinal Cord',
                        cord              = 'Spinal Cord',
                        cochlea_r         = 'R. Cochlea',
                        cochlea_l         = 'L. Cochlea',
                        trachea           = 'Trachea',
                        parotid_l         = 'L. Parotid',
                        parotid_r         = 'R. Parotid',
                        submandibgld_l    = 'L. Submandible gland',
                        submandibgld_r    = 'R. Submandible gland',
                        oralcavity        = 'Oral cavity',
                        esophconstrictors = 'Esophagus constrictors',
                        larynx            = 'Larynx',
                        patient.target    = 'Patient',
                        patient           = 'Whole Patient',
                        mandible          = 'Mandible')
    
    main_oars <- list(P01 = c('R. Submandible gland', 'L. Submandible gland'),
                      P02 = c('R. Submandible gland'),
                      P03 = c('L. Submandible gland'),
                      P04 = c('R. Parotid', 'Larynx', 'Mandible', 'Oral cavity'),
                      P05 = c('Larynx', 'R. Submandible gland', 'L. Submandible gland'),
                      P07 = c('Mandible'),
                      P10 = c('Larynx', 'Esophagus constrictors'),
                      P14 = c('Mandible', 'R. Submandible gland'),
                      P15 = c('Larynx'),
                      P16 = c('R. Submandible gland', 'Oral cavity'))

    dt <- dt %>%
        group_by(patient, struct) %>%
        mutate(main_oar = ifelse(struct %in% unlist(main_oars[as.character(patient)]), TRUE, FALSE)) %>%
        ungroup()
    
    # Make GTV then CTV first
    if ('GTV' %in% levels(dt$struct))
        dt$struct <- relevel(dt$struct, 'GTV')
    if ('CTV' %in% levels(dt$struct))
        dt$struct <- relevel(dt$struct, 'CTV')

    return(dt)
}

map_cases <- function(data, inputs) {
    outputs <- NULL
    for (i in inputs) {
        num <- gsub("[^0-9]", "", i)
        case <- gsub(paste0('cbct_', num, '_adapt_'), '', i)
        case <- gsub('_', ' ', case)
        case <- tools::toTitleCase(case)
        case <- gsub(' e ', ' E ', case)
        outputs <- c(outputs, paste(case, num))
    }
    return(plyr::mapvalues(data, warn_missing = FALSE, from = inputs, to = outputs))
}

organize_cases <- function(dt) {
    summary(dt)
    dt$case <- plyr::mapvalues(dt$case, warn_missing = FALSE,
                         from = c('cbct_1_robust', 'cbct_2_robust', 'cbct_3_robust',
                                  'cbct_4_robust', 'cbct_5_robust', 'cbct_6_robust', 'cbct_7_robust'),
                         to = c('Week 1 Robust', 'Week 2 Robust', 'Week 3 Robust',
                                'Week 4 Robust', 'Week 5 Robust', 'Week 6 Robust', 'Week 7 Robust'))
    dt$case <- plyr::mapvalues(dt$case, warn_missing = FALSE,
                         from = c('base', 'cbct_1', 'cbct_2', 'cbct_3', 'cbct_4', 'cbct_5', 'cbct_6', 'cbct_7'),
                         to = c('Plan', 'Week 1', 'Week 2', 'Week 3', 'Week 4', 'Week 5', 'Week 6', 'Week 7'))
    inputs <- c('cbct_1_geometric_free', 'cbct_2_geometric_free', 'cbct_3_geometric_free',
                'cbct_4_geometric_free', 'cbct_5_geometric_free', 'cbct_6_geometric_free',
                'cbct_7_geometric_free')
    dt$case <- map_cases(dt$case, inputs)
    inputs <- c('cbct_1_geometric_iso_shift', 'cbct_2_geometric_iso_shift', 'cbct_3_geometric_iso_shift',
                'cbct_4_geometric_iso_shift', 'cbct_5_geometric_iso_shift', 'cbct_6_geometric_iso_shift',
                'cbct_7_geometric_iso_shift')
    dt$case <- map_cases(dt$case, inputs)
    inputs <- c('cbct_1_cold_spots_free', 'cbct_2_cold_spots_free', 'cbct_3_cold_spots_free',
                'cbct_4_cold_spots_free', 'cbct_5_cold_spots_free', 'cbct_6_cold_spots_free',
                'cbct_7_cold_spots_free')
    dt$case <- map_cases(dt$case, inputs)
    inputs <- c('cbct_1_cold_spots_iso_shift', 'cbct_2_cold_spots_iso_shift', 'cbct_3_cold_spots_iso_shift',
                'cbct_4_cold_spots_iso_shift', 'cbct_5_cold_spots_iso_shift', 'cbct_6_cold_spots_iso_shift',
                'cbct_7_cold_spots_iso_shift')
    dt$case <- map_cases(dt$case, inputs)
    
    return(dt)
}

filter_structures <- function(dt, patient)
{
    # normalize names to lower case
    setnames(dt, names(dt), tolower(names(dt)))
    # replace weird characters
    setnames(dt, names(dt), sub("\ \\(gy\\)", '', names(dt)))
    # keep only structs contained in list
    if ( patient == "P01" ) {
        setnames(dt, names(dt), sub("ctv op", 'ctv', names(dt)))
        dt <- select(dt, dose, ctv, spinalcord, parotid_l, parotid_r,
                     oralcavity, submandibgld_r, submandibgld_l, patient-target)
    } else if ( patient == "P02" ) {
        setnames(dt, names(dt), sub("gtv tongue", 'gtv', names(dt)))
        setnames(dt, names(dt), sub("ctv_tongue", 'ctv', names(dt)))
        dt <- select(dt, dose, ctv, larynx, parotid_r, oralcavity, mandible,
                     submandibgld_r, cord, patient-target)
    } else if ( patient == "P03" ) {
        setnames(dt, names(dt), sub("gtv l tonsil", 'gtv', names(dt)))
        setnames(dt, names(dt), sub("ctv l tonsil", 'ctv', names(dt)))
        dt <- select(dt, dose, ctv, larynx, parotid_r, mandible, submandibgld_r,
                     submandibgld_l, spinalcord, patient-target)
    } else if ( patient == "P04" ) {
        dt <- select(dt, -ctv)
        setnames(dt, names(dt), sub("ctv_neck_r", 'ctv', names(dt)))
        dt <- select(dt, dose, patient, spinalcord, esophconstrictors, larynx,
                     mandible, parotid_r, oralcavity, parotid_l, ctv)
    } else if ( patient == "P05" ) {
        dt <- select(dt, dose, ctv, patient, spinalcord, esophconstrictors,
                     parotid_r, parotid_l, oralcavity, mandible, submandibgld_l,
                     larynx, submandibgld_r)
    } else if ( patient == "P07" ) {
        dt <- select(dt, dose, ctv, patient, spinalcord, parotid_l, mandible,
                     esophconstrictors, larynx)
    } else if ( patient == "P10" ) {
        setnames(dt, names(dt), sub("ctv_larynx", 'ctv', names(dt)))
        dt <- select(dt, dose, ctv, patient, esophconstrictors, spinalcord,
                     parotid_l, parotid_r, oralcavity, mandible, larynx,
                     trachea, esophagus, submandibgld_l, submandibgld_r,
                     cochlea_r, cochlea_l)
    } else if ( patient == "P14" ) {
        setnames(dt, names(dt), sub("ctvtongue", 'ctv', names(dt)))
        dt <- select(dt, dose, ctv, patient, mandible, spinalcord, parotid_l,
                     parotid_r, submandibgld_r, esophconstrictors, larynx,
                     trachea)
    } else if ( patient == "P15" ) {
        setnames(dt, names(dt), sub("ctv tonsil r", 'ctv', names(dt)))
        dt <- select(dt, dose, ctv, patient, spinalcord, parotid_r, mandible,
                     esophconstrictors, larynx, parotid_l, oralcavity,
                     submandibgld_r)
    } else if ( patient == 'P16' ) {
        setnames(dt, names(dt), sub('ctv bot l', 'ctv', names(dt)))
        dt <- select(dt, dose, ctv, patient, spinalcord, parotid_r, mandible,
                     esophconstrictors, larynx, parotid_l, oralcavity,
                     submandibgld_r)
    } else {
        print('WARNING! No explicit selection of structures. Automatically selecting:')
        print(names(dt))
    }

    return(dt)
}

process_single_file <- function(file_name, location, scans, pat_no)
{
    print(file_name)
    # Parse file name
    path <- paste(location, file_name, sep = '/')
    patient <- sub("_.*", "", file_name)
    case <- sub(".*[0-9]{3}_.*?", "", file_name)
    case <- sub(".dvh", "", case)
    case <- sub("_plan", "", case)
    # Fast solution for P04 data, NOT IDEAL!!
    case <- sub("Neck_R_", "", case)
    
    # Read file and remove dummy line
    dt <- as.data.table(fread(path))
    
    # Process DT
    dt <- dt %>%
        filter_structures(patient) %>%
        gather(struct, vol, -dose, factor_key = TRUE) %>%
        filter(vol > 0) %>%
        mutate(dose = dose/plan_fractions, vol = vol*100) %>%
        mutate(dose_pct = 100*dose/target_dose) %>%
        mutate(oar = ifelse(struct %in% c('gtv','ctv', 'GTV', 'CTV'), FALSE, TRUE)) %>%
        mutate(scans = scans) %>%
        mutate(patient = factor(patient)) %>%
        mutate(patient.no = factor(pat_no)) %>%
        mutate(case = case) %>%
        ## Organize and map input
        organize_structures() %>%
        organize_cases() %>%
        mutate(struct = as.factor(struct)) %>%
        ## Timing information
        mutate(case = as.character(case)) %>%
        mutate(week.no = ifelse(case == 'Plan', 0,
                                ifelse(grepl('cumulative', case), 99, as.numeric(gsub("\\D", "", case))))) %>%
        mutate(pat.week = paste0(patient, '.', week.no)) %>%
        mutate(week.name = ifelse(week.no == 0, 'Plan',
                                  ifelse(week.no == 99, 'Cum.', paste('Week', week.no)))) %>%
        mutate(week.name = as.factor(week.name)) %>%
        mutate(checkpoints = ifelse(week.no == 0, 'Plan',
                                    ifelse(week.no == scans, 'Last',
                                           ifelse(week.no == scans %/% 2, 'Mid',
                                                  ifelse(week.no == 99, 'Cum.', 'Other'))))) %>%
        mutate(checkpoints = factor(checkpoints)) %>%
        mutate(stage = ifelse(week.no == 0, 'Plan',
                              ifelse(week.no == 99, 'Cum.', 'Weekly'))) %>%
        mutate(stage = factor(stage))
        ## Adaptation Methods
    dt <- dt %>%
        mutate(method = ifelse(str_detect(case, 'Robust'), 'Robust', case)) %>%
        mutate(method = ifelse(str_detect(case, 'Week') && !str_detect(case, 'Robust'), 'None', method)) %>%
        mutate(method = ifelse(str_detect(method, 'adapt_geometric'), 'Geometric', method)) %>%
        mutate(method = ifelse(str_detect(method, 'adapt_dij'), 'Dij', method)) %>%
        mutate(method = ifelse(str_detect(method, 'adapt_cold_spots'), 'Weights', method)) %>%
        mutate(method = ifelse(str_detect(method, 'cumulative'), 'None', method)) %>%
        mutate(method = as.factor(method)) %>%
        ## Adaptation Methods Constraints
        mutate(constraint = ifelse(str_detect(case, '_free'), 'Free', case)) %>%
        mutate(constraint = ifelse(str_detect(constraint, '_fixed'), 'Fixed', constraint)) %>%
        mutate(constraint = ifelse(str_detect(constraint, '_iso_shift_range_shifter'), 'Isocenter - Range shifter', constraint)) %>%
        mutate(constraint = ifelse(str_detect(constraint, '_range_shifter_iso_shift'), 'Isocenter - Range shifter', constraint)) %>%
        mutate(constraint = ifelse(str_detect(constraint, '_iso_shift'), 'Isocenter', constraint)) %>%
        mutate(constraint = ifelse(str_detect(constraint, '_range_shifter'), 'Range shifter', constraint)) %>%
        mutate(constraint = ifelse(str_detect(constraint, '_|Week|cumulative'), 'None', constraint)) %>%
        mutate(constraint = factor(constraint)) %>%
        mutate(short_const = as.character(constraint)) %>%
        mutate(short_const = ifelse(short_const == 'Isocenter', 'Iso',
                                    ifelse(short_const == 'Range shifter', 'RS',
                                           ifelse(short_const == 'Isocenter - Range shifter', 'Iso-RS', short_const)))) %>%
        mutate(short_const = factor(short_const)) %>%
        ## Case
        mutate(case = paste0(method,'/',constraint)) %>%
        mutate(case = as.factor(case)) %>%
        ## Other factors and vars
        mutate(week.no = factor(week.no),
               week.name = factor(week.name),
               patient_orig = patient) %>%
        ## To output format
        as.data.table()
    
    return(dt)
}

readData <- function(patients, location, plan_fractions, target_dose, verbose = FALSE)
{
    require(parallel)
    # Loop around files
    cohort_dt <- data.table()
    pat_no <- 1
    for (pat in patients) {
        file.names <- dir(location, pattern = pat)
        scans <- max(as.numeric(sub(".*cbct_([1-9]+?)[.,_].*", "\\1", file.names[grep('cbct_', file.names)])))
        pat_dt <- mclapply(file.names, process_single_file, location = location,
                           scans = scans, pat_no = pat_no, mc.cores = 4)
        # pat_dt <- lapply(file.names, process_single_file, location = location,
                         # scans = scans, pat_no = pat_no)
        pat_dt <- rbindlist(pat_dt)

        if (verbose)
            print(summary(pat_dt))
        cohort_dt <- rbindlist(list(cohort_dt, pat_dt))
        pat_no <- pat_no + 1
    }

    ## Reorder factor levels: method
    if ('None' %in% levels(cohort_dt$method))
        cohort_dt$method <- relevel(cohort_dt$method, 'None')
    if ('Plan' %in% levels(cohort_dt$method))
        cohort_dt$method <- relevel(cohort_dt$method, 'Plan') # Will be the first
    if ('Isocenter - Range shifter' %in% levels(cohort_dt$constraint))
        cohort_dt$constraint <- relevel(cohort_dt$constraint, 'Isocenter - Range shifter')
    if ('Range shifter' %in% levels(cohort_dt$constraint))  
        cohort_dt$constraint <- relevel(cohort_dt$constraint, 'Range shifter')
    if ('Isocenter' %in% levels(cohort_dt$constraint))
        cohort_dt$constraint <- relevel(cohort_dt$constraint, 'Isocenter')
    if ('Free' %in% levels(cohort_dt$constraint))
        cohort_dt$constraint <- relevel(cohort_dt$constraint, 'Free')
    if ('Plan' %in% levels(cohort_dt$constraint))
        cohort_dt$constraint <- relevel(cohort_dt$constraint, 'Plan') # Will be the first
    
    ## Reorder factor levels: week.name
    if ('Cum.' %in% levels(cohort_dt$week.name))
        cohort_dt$week.name <- relevel(cohort_dt$week.name, 'Cum.') # Will be the second
    if ('Plan' %in% levels(cohort_dt$week.name))
        cohort_dt$week.name <- relevel(cohort_dt$week.name, 'Plan') # Will be the first
    
    cohort_dt <- cohort_dt %>%
        mutate(short_const)
    
    ## Rename patient files
    levels(cohort_dt$patient) <- sprintf("Patient %s", 1:pat_no)
    
    return(as.data.table(cohort_dt))
}
