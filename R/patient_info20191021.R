setwd("patient information/")
all346_urine_info <- 
  readr::read_csv("SmartD_all346urine.csv")

clinic_varibale_summary <-
  readxl::read_xlsx("SmartD_ClinicalVariables_PartiallySummarized.xlsx", 
                    sheet = 1)

###This is batch data set
smartd_rplc_batch1 <- 
  readr::read_csv("E:/project/smartD/data_analysis20190828/prediction/identification_table/smartd_rplc_batch1.csv")

smartd_rplc_batch2 <- 
  readr::read_csv("E:/project/smartD/smartD_batch1_2/RPLC/POS_NEG/ms1_data_rplc_batch2.csv")

###we must get the samples with 'P' coresponding to which sample
all346_urine_info$sample_id_RPLC %>% stringr::str_sort(numeric = TRUE) %>% 
  grep("SFU_B", ., value = TRUE)

##we have 198 sampels in total
colnames(smartd_rplc_batch2) %>% stringr::str_sort(numeric = TRUE) %>% 
  grep("X", ., value = TRUE)

colnames(smartd_rplc_batch2) %>% stringr::str_sort(numeric = TRUE) %>% 
  grep("P", ., value = TRUE)

###so it is very clear that X178 is P1, X179 is P2, and X180 is P3. X198 is P21.

colnames(smartd_rplc_batch2)[grep("P", colnames(smartd_rplc_batch2))] <-
  colnames(smartd_rplc_batch2)[grep("P", colnames(smartd_rplc_batch2))] %>% 
  stringr::str_replace("P", "") %>% 
    as.numeric() %>% 
    `+`(177) %>% 
    paste("X", ., sep = "")

name_batch1 <-
  smartd_rplc_batch1 %>% 
  select(-(name:Database)) %>% 
  colnames()

name_batch2 <- 
  smartd_rplc_batch2 %>% 
  select(-(name:rt)) %>% 
  select(-contains("blk")) %>% 
  select(-contains("QC")) %>% 
  colnames()

sample_id_RPLC <-
  all346_urine_info$sample_id_RPLC

sample_id_RPLC[grep("SFU_B", sample_id_RPLC)] <- 
  sample_id_RPLC[grep("SFU_B", sample_id_RPLC)] %>% 
  gsub(pattern = "SFU_B", replacement = "", x = .) %>% 
  paste("X", ., sep = "")

name_batch1 %in% sample_id_RPLC
name_batch2 %in% sample_id_RPLC


all346_urine_info$sample_id_RPLC <-
  sample_id_RPLC

all346_info <- 
  all346_urine_info

all346_info <- 
  all346_info %>% 
  select(Patient_ID = PTID, Sample_ID = sample_id_RPLC, 
         Visit, Date.Acquired, GA = Visit.GA) %>% 
  arrange(Patient_ID, Visit)

##some sample's GA are NA, should be removed
# all346_info %>% 
#   filter(is.na(GA)) %>% 
#   pull(Sample_ID)
# 
# patient_info <- 
#   readr::read_csv("patient_info.csv")
# 
# patient_info$Sample_ID <- 
#   patient_info$Sample_ID %>% 
#   gsub(pattern = "SFU_B", replacement = "X", x = .)
#   
# patient_info %>% 
#   select(Sample_ID, `Visit GA`) %>% 
#   arrange(Sample_ID)
# 
# all346_info %>% 
#   select(Sample_ID, GA) %>% 
#   filter(stringr::str_starts(Sample_ID, "X")) %>% 
#   arrange(Sample_ID)
# 
# 
# intersect(all346_info$Sample_ID, patient_info$Sample_ID)

###modify GA
GA <- all346_info$GA

GA <- 
  (GA - trunc(GA))*10/7 + trunc(GA)

all346_info$GA <-
  GA

##patient_info
final_info <- 
  readr::read_csv("FINALSMARTDiaphragm2_DATA_2019-07-01_1615.csv")

###change the Patient ID
all346_info$Patient_ID[-grep(pattern = "SF", all346_info$Patient_ID)] <-
  all346_info$Patient_ID[-grep(pattern = "SF", all346_info$Patient_ID)] %>% 
  paste("SF", ., sep = "")
  
participant_id <- 
final_info$participant_id 

participant_id <-
  gsub("sf", "SF", participant_id)

participant_id[-grep('SF', participant_id)] <-
  paste("SF", participant_id[-grep('SF', participant_id)], sep = "")
  
 
  final_info$participant_id <-
    participant_id
  
  
setdiff(all346_info$Patient_ID, final_info$participant_id)
setdiff(final_info$participant_id, all346_info$Patient_ID)

final_info$redcap_event_name <- 
  gsub("study_visit_|_arm_1", "", final_info$redcap_event_name) %>% 
  as.numeric()


sample_info_191021 <- 
  left_join(all346_info, final_info, by = c("Patient_ID" = "participant_id",
                                          "Visit" = "redcap_event_name"))

sample_info_191021 %>% 
  filter(Sample_ID %in% name_batch1) %>% 
  pull(GA)

sample_info_191021 %>% 
  filter(Sample_ID %in% name_batch2) %>% 
  pull(GA)

clinic_varibale_summary$ID <- 
  stringr::str_to_upper(clinic_varibale_summary$ID)

clinic_varibale_summary$ID[-grep("SF", clinic_varibale_summary$ID)] <-
  paste("SF", clinic_varibale_summary$ID[-grep("SF", clinic_varibale_summary$ID)], sep = "")
  
sample_info_191021$Patient_ID %in% clinic_varibale_summary$ID

sample_info_191021 <- 
  sample_info_191021 %>% 
  left_join(clinic_varibale_summary, by = c("Patient_ID" = "ID"))

  

write.csv(sample_info_191021, "sample_info_191021.csv", row.names = FALSE)  






# patient_info$Sample_ID
# 
# info1 <- 
#   sample_info_191021 %>% 
#   filter(stringr::str_starts(Sample_ID, "X")) %>% 
#   select(Patient_ID:GA) %>% 
#   arrange(Sample_ID)
# 
# info2 <- 
#   patient_info %>% 
#   filter(stringr::str_starts(Sample_ID, "SFU_B")) %>%
#   select(Patient_ID:`Visit GA`) %>% 
#   arrange(Sample_ID) %>% 
#   mutate(Sample_ID = stringr::str_replace(Sample_ID, "SFU_B", "X"))
  


#####the clinic information 





