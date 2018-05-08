library(data.table)
library(dplyr)

organize_structures <- function(dt)
{
    dt$struct <- revalue(dt$struct, warn_missing = FALSE,
                         c('gtv'               = 'GTV',
                           'ctv'               = 'CTV',
                           'esophagus'         = 'Esophagus',
                           'spinalcord'        = 'Spinal Cord',
                           'cord'              = 'Spinal Cord',
                           'trachea'           = 'Trachea',
                           'parotid_l'         = 'L. Parotid',
                           'parotid_r'         = 'R. Parotid',
                           'submandibgld_l'    = 'L. Submandible gland',
                           'submandibgld_r'    = 'R. Submandible gland',
                           'oralcavity'        = 'Oral cavity',
                           'esophconstrictors' = 'Esophagus constrictors',
                           'larynx'            = 'Larynx',
                           'patient'           = 'Patient',
                           'mandible'          = 'Mandible')
                        )
    # Make GTV then CTV first
    if ("GTV" %in% levels(dt$struct))
        dt$struct <- relevel(dt$struct, "GTV")
    if ("CTV" %in% levels(dt$struct))
        dt$struct <- relevel(dt$struct, "CTV")

    return(dt)
}

organize_cases <- function(dt){
    dt$case <- mapvalues(dt$case, warn_missing = FALSE,
              from = c('base', 'cbct_1', 'cbct_2', 'cbct_3', 'cbct_4', 'cbct_5', 'cbct_6', 'cbct_7'),
              to = c('Plan', 'Week 1', 'Week 2', 'Week 3', 'Week 4', 'Week 5', 'Week 6', 'Week 7'))
    return(dt)
}

filter_structures <- function(dt, patient)
{
    # normalize names to lower case
    names(dt) <- tolower(names(dt))
    # replace weird characters
    names(dt) <- sub("\\..gy.", '', names(dt))
    # keep only structs contained in list
    if ( patient == "P01" ) {
        names(dt) <- sub("ctv.op", 'ctv', names(dt))
        l <- c('dose', 'ctv', 'spinalcord', 'esophconstrictors', 'larynx', 'parotid_l', 'parotid_r', 'oralcavity', 'submandibgld_r', 'submandibgld_l')
    } else if ( patient == "P02" ) {
        names(dt) <- sub("gtv.tongue", 'gtv', names(dt))
        names(dt) <- sub("ctv_tongue", 'ctv', names(dt))
        l <- c('dose', 'ctv', 'esophconstrictors', 'larynx', 'parotid_r', 'oralcavity', 'mandible', 'submandibgld_r', 'submandibgld_l', 'cord') # 'trachea', 'esophagus', 'parotid_l', 'brainstem???', 'cochlea_r???'
    } else if ( patient == "P03" ) {
        names(dt) <- sub("gtv.l.tonsil", 'gtv', names(dt))
        names(dt) <- sub("ctv.l.tonsil", 'ctv', names(dt))
        l <- c('dose', 'ctv', 'esophconstrictors', 'larynx', 'parotid_l', 'parotid_r', 'oralcavity', 'mandible', 'submandibgld_r', 'submandibgld_l', 'spinalcord')
    } else if ( patient == "P15" ) {
        names(dt) <- sub("gtv.tonsil.r", 'gtv', names(dt))
        names(dt) <- sub("ctv.tonsil.r", 'ctv', names(dt))
        l <- c('dose', 'ctv', 'esophconstrictors', 'larynx', 'parotid_r', 'oralcavity', 'mandible', 'submandibgld_l', 'spinalcord')
    }

    dt <- select(dt, one_of(l))

    return(dt)
}

readData <- function(patients, location, plan_fractions, verbose = FALSE)
{
    # Loop around files
    cohort_dt <- data.table()
    for (pat in patients)
    {
        pat_dt <- data.table()
        file.names <- dir(location, pattern = pat)
        file.names <- file.names[!grepl('adapt', file.names)]
        for (file in file.names)
        {
            print(file)
            # Parse file name
            path <- paste(location, file, sep = '/')
            patient <- sub("_.*", "", file)
            case <- sub(".*[0-9]{3}_.*?", "", file)
            case <- sub(".dvh", "", case)
            case <- sub("_plan", "", case)

            # Read file and remove dummy line
            dt <- as.data.table(read.table(path, header = TRUE, sep = ","))

            # Process DT
            dt <- dt %>%
                filter_structures(patient) %>%
                gather(struct, vol, -dose, factor_key = TRUE) %>%
                filter(vol > 0) %>%
                mutate(dose = dose/plan_fractions, vol = vol*100) %>%
                mutate(oar = ifelse(tolower(struct) %in% c('gtv','ctv'), FALSE, TRUE)) %>%
                mutate(fractions = length(file.names) - 1) %>%
                mutate(patient = as.factor(patient)) %>%
                mutate(case = as.factor(case)) %>%
                organize_structures() %>%
                organize_cases() %>%
                mutate(week.no = ifelse(as.character(case) == 'Plan', 0, as.numeric(gsub("\\D", "", as.character(case))))) %>%
                mutate(checkpoints = ifelse(week.no == 0, 'Plan',
                                            ifelse(week.no == fractions, 'Last',
                                                   ifelse(week.no == fractions %/% 2, 'Mid', 'Other'
                                                          )
                                                   )
                                            )
                       ) %>%
                mutate(checkpoints = as.factor(checkpoints)) %>%
                as.data.table()

            pat_dt <- rbind(pat_dt, dt)
        }
        if (verbose)
            print(summary(pat_dt))
        cohort_dt <- rbind(cohort_dt, pat_dt)
    }

    return(as.data.table(cohort_dt))
}
