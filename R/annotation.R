

library(metID)


#####RPLC positive
setwd("H:/smartD/RPLC/POS/annotation/NCE25/")
parameter1 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "positive", 
                               ce = "all", 
                               column = "rp", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "msDatabase_rplc0.0.1", 
                               threads = 4)

parameter2 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "positive", 
                               ce = "all", 
                               column = "rp", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "hmdbDatabase0.0.1", 
                               threads = 4)

parameter3 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "positive", 
                               ce = "all", 
                               column = "rp", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "massbankDatabase0.0.1", 
                               threads = 4)

parameter4 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "positive", 
                               ce = "all", 
                               column = "rp", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "metlinDatabase0.0.1", 
                               threads = 4)

parameter5 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "positive", 
                               ce = "all", 
                               column = "rp", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "monaDatabase0.0.1", 
                               threads = 4)

parameter6 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "positive", 
                               ce = "all", 
                               column = "rp", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "nistDatabase0.0.1", 
                               threads = 4)

parameter7 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "positive", 
                               ce = "all", 
                               column = "rp", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "orbitrapDatabase0.0.1", 
                               threads = 4)

parameter8 <- mzIdentifyParam(ms1.match.ppm = 25, 
                              polarity = "positive", 
                              column = "rp", 
                              candidate.num = 3,
                              database = "HMDB.metabolite.data", 
                              threads = 4)
ms1.data <- grep("\\.csv", dir(), value = TRUE)
ms2.data <- grep("mgf", dir(), value = TRUE)
result.pRPLC.nce25 <- metIdentify4all(ms1.data = ms1.data,
                                      ms2.data = ms2.data, 
                                      parameter.list = c(parameter1,
                                                         parameter2,
                                                         parameter3,
                                                         parameter4,
                                                         parameter5,
                                                         parameter6,
                                                         parameter7,
                                                         parameter8),
                                      path = ".")

save(result.pRPLC.nce25, file = "result.pRPLC.nce25")
print(1)


setwd("H:/smartD/RPLC/POS/annotation/NCE50/")

ms1.data <- grep("\\.csv", dir(), value = TRUE)
ms2.data <- grep("mgf", dir(), value = TRUE)
result.pRPLC.nce50 <- metIdentify4all(ms1.data = ms1.data,
                                      ms2.data = ms2.data, 
                                      parameter.list = c(parameter1,
                                                         parameter2,
                                                         parameter3,
                                                         parameter4,
                                                         parameter5,
                                                         parameter6,
                                                         parameter7),
                                      path = ".")

save(result.pRPLC.nce50, file = "result.pRPLC.nce50")
print(1)

#####RPLC negative
setwd("H:/smartD/RPLC/NEG/annotation/NCE25/")
parameter1 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "negative", 
                               ce = "all", 
                               column = "rp", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "msDatabase_rplc0.0.1", 
                               threads = 4)

parameter2 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "negative", 
                               ce = "all", 
                               column = "rp", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "hmdbDatabase0.0.1", 
                               threads = 4)

parameter3 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "negative", 
                               ce = "all", 
                               column = "rp", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "massbankDatabase0.0.1", 
                               threads = 4)

parameter4 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "negative", 
                               ce = "all", 
                               column = "rp", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "metlinDatabase0.0.1", 
                               threads = 4)

parameter5 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "negative", 
                               ce = "all", 
                               column = "rp", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "monaDatabase0.0.1", 
                               threads = 4)

parameter6 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "negative", 
                               ce = "all", 
                               column = "rp", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "nistDatabase0.0.1", 
                               threads = 4)

parameter7 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "negative", 
                               ce = "all", 
                               column = "rp", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "orbitrapDatabase0.0.1", 
                               threads = 4)

parameter8 <- mzIdentifyParam(ms1.match.ppm = 25, 
                              polarity = "negative", 
                              column = "rp", 
                              candidate.num = 3,
                              database = "HMDB.metabolite.data", 
                              threads = 4)

ms1.data <- grep("\\.csv", dir(), value = TRUE)
ms2.data <- grep("mgf", dir(), value = TRUE)
result.nRPLC.nce25 <- metIdentify4all(ms1.data = ms1.data,
                                      ms2.data = ms2.data, 
                                      parameter.list = c(parameter1,
                                                         parameter2,
                                                         parameter3,
                                                         parameter4,
                                                         parameter5,
                                                         parameter6,
                                                         parameter7,
                                                         parameter8),
                                      path = ".")

save(result.nRPLC.nce25, file = "result.nRPLC.nce25")
print(1)

setwd("H:/smartD/RPLC/NEG/annotation/NCE50/")
ms1.data <- grep("\\.csv", dir(), value = TRUE)
ms2.data <- grep("mgf", dir(), value = TRUE)
result.nRPLC.nce50 <- metIdentify4all(ms1.data = ms1.data,
                                      ms2.data = ms2.data, 
                                      parameter.list = c(parameter1,
                                                         parameter2,
                                                         parameter3,
                                                         parameter4,
                                                         parameter5,
                                                         parameter6,
                                                         parameter7),
                                      path = ".")

save(result.nRPLC.nce50, file = "result.nRPLC.nce50")
print(1)

##intergrate different anotation results from different databases
####RPLC positive
setwd("H:/smartD/RPLC/POS/annotation/")
load("NCE25/result.pRPLC.nce25")
load("NCE50/result.pRPLC.nce50")
library(tidyverse)
identification.table1 <- getIdentificationTable(result.pRPLC.nce25[[1]],
                                                result.pRPLC.nce50[[1]],
                                                candidate.num = 1,
                                                type = "old")

identification.table1 <- 
  identification.table1 %>% 
  filter(!is.na(Identification)) %>% 
  mutate(Level = 1)

dim(identification.table1)

identification.table2 <-
  getIdentificationTable(
    result.pRPLC.nce25[[2]],
    result.pRPLC.nce25[[3]],
    result.pRPLC.nce25[[4]],
    result.pRPLC.nce25[[5]],
    result.pRPLC.nce25[[6]],
    result.pRPLC.nce25[[7]],
    result.pRPLC.nce50[[2]],
    result.pRPLC.nce50[[3]],
    result.pRPLC.nce50[[4]],
    result.pRPLC.nce50[[5]],
    result.pRPLC.nce50[[6]],
    result.pRPLC.nce50[[7]],
    candidate.num = 1,
    type = "old"
  )

identification.table2 <-
  identification.table2 %>% 
  filter(!is.na(Identification), !is.element(name, identification.table1$name)) %>% 
  mutate(Level = 2)

dim(identification.table2)


identification.table3 <- 
  getIdentificationTable2(result.pRPLC.nce25[[8]],
                          candidate.num = 3,
                          type = "old")


identification.table3 <-
  identification.table3 %>% 
  filter(!is.na(Identification), !is.element(name, identification.table1$name), !is.element(name, identification.table2$name)) %>% 
  mutate(Level = 3)

dim(identification.table3)





identification.table <- 
  dplyr::full_join(rbind(identification.table1, identification.table2), 
                   identification.table3, 
                   by = colnames(identification.table3))

write.csv(identification.table, "identification.table.csv", row.names = FALSE)


identification.table.new <- trans2newStyle(identification.table = identification.table)
write.csv(identification.table.new, "identification.table.new.csv", row.names = FALSE)


####RPLC negative
setwd("H:/smartD/RPLC/NEG/annotation/")
load("NCE25/result.nRPLC.nce25")
load("NCE50/result.nRPLC.nce50")

identification.table1 <- getIdentificationTable(result.nRPLC.nce25[[1]],
                                                result.nRPLC.nce50[[1]],
                                                candidate.num = 1,
                                                type = "old")

identification.table1 <- 
  identification.table1 %>% 
  filter(!is.na(Identification)) %>% 
  mutate(Level = 1)

dim(identification.table1)

identification.table2 <-
  getIdentificationTable(
    result.nRPLC.nce25[[2]],
    result.nRPLC.nce25[[3]],
    result.nRPLC.nce25[[4]],
    result.nRPLC.nce25[[5]],
    result.nRPLC.nce25[[6]],
    result.nRPLC.nce25[[7]],
    result.nRPLC.nce50[[2]],
    result.nRPLC.nce50[[3]],
    result.nRPLC.nce50[[4]],
    result.nRPLC.nce50[[5]],
    result.nRPLC.nce50[[6]],
    result.nRPLC.nce50[[7]],
    candidate.num = 1,
    type = "old"
  )

identification.table2 <-
  identification.table2 %>% 
  filter(!is.na(Identification), !is.element(name, identification.table1$name)) %>% 
  mutate(Level = 2)

dim(identification.table2)


identification.table3 <- 
  getIdentificationTable2(result.nRPLC.nce25[[8]],
                          candidate.num = 3,
                          type = "old")


identification.table3 <-
  identification.table3 %>% 
  filter(!is.na(Identification), !is.element(name, identification.table1$name), !is.element(name, identification.table2$name)) %>% 
  mutate(Level = 3)

dim(identification.table3)



identification.table <- 
  dplyr::full_join(rbind(identification.table1, identification.table2), 
                   identification.table3, 
                   by = colnames(identification.table3))

write.csv(identification.table, "identification.table.csv", row.names = FALSE)


identification.table.new <- trans2newStyle(identification.table = identification.table)
write.csv(identification.table.new, "identification.table.new.csv", row.names = FALSE)











#####HILIC positive
setwd("/home/shenxt/plasma/HILIC-POS-MS2/CE50")
library(metIdentify)
parameter1 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "positive", 
                               ce = "all", 
                               column = "hilic", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "msDatabase_hilic0.0.1", 
                               threads = 4)

parameter2 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "positive", 
                               ce = "all", 
                               column = "hilic", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "hmdbDatabase0.0.1", 
                               threads = 4)

parameter3 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "positive", 
                               ce = "all", 
                               column = "hilic", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "massbankDatabase0.0.1", 
                               threads = 4)

parameter4 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "positive", 
                               ce = "all", 
                               column = "hilic", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "metlinDatabase0.0.1", 
                               threads = 4)

parameter5 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "positive", 
                               ce = "all", 
                               column = "hilic", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "monaDatabase0.0.1", 
                               threads = 4)

parameter6 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "positive", 
                               ce = "all", 
                               column = "hilic", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "nistDatabase0.0.1", 
                               threads = 4)

parameter7 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "positive", 
                               ce = "all", 
                               column = "hilic", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "orbitrapDatabase0.0.1", 
                               threads = 4)

parameter8 <- mzIdentifyParam(ms1.match.ppm = 25, 
                              polarity = "positive", 
                              column = "hilic", 
                              candidate.num = 3,
                              database = "HMDB.metabolite.data", 
                              threads = 4)
ms1.data <- grep("\\.csv", dir(), value = TRUE)
ms2.data <- grep("mgf", dir(), value = TRUE)
result.pHILIC.nce50 <- metIdentify4all(ms1.data = ms1.data,
                                       ms2.data = ms2.data, 
                                       parameter.list = c(parameter1,
                                                          parameter2,
                                                          parameter3,
                                                          parameter4,
                                                          parameter5,
                                                          parameter6,
                                                          parameter7,
                                                          parameter8),
                                       path = ".")

save(result.pHILIC.nce50, file = "result.pHILIC.nce50")
print(1)














#####HILIC negative
library(metIdentify)
setwd("/home/shenxt/plasma/HILIC-NEG-MS2/CE25")
parameter1 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "negative", 
                               ce = "all", 
                               column = "hilic", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "msDatabase_hilic0.0.1", 
                               threads = 4)

parameter2 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25,
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "negative", 
                               ce = "all", 
                               column = "hilic", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "hmdbDatabase0.0.1", 
                               threads = 4)

parameter3 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "negative", 
                               ce = "all", 
                               column = "hilic", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "massbankDatabase0.0.1", 
                               threads = 4)

parameter4 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "negative", 
                               ce = "all", 
                               column = "hilic", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "metlinDatabase0.0.1", 
                               threads = 4)

parameter5 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "negative", 
                               ce = "all", 
                               column = "hilic", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "monaDatabase0.0.1", 
                               threads = 4)

parameter6 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "negative", 
                               ce = "all", 
                               column = "hilic", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "nistDatabase0.0.1", 
                               threads = 4)

parameter7 <- metIdentifyParam(ms1.ms2.match.mz.tol = 25, 
                               ms1.ms2.match.rt.tol = 60, 
                               ms1.match.ppm = 25, 
                               ms2.match.tol = 0.4, 
                               rt.match.tol = 90,
                               polarity = "negative", 
                               ce = "all", 
                               column = "hilic", 
                               total.score.tol = 0, 
                               candidate.num = 3,
                               database = "orbitrapDatabase0.0.1", 
                               threads = 4)

parameter8 <- mzIdentifyParam(ms1.match.ppm = 25, 
                              polarity = "negative", 
                              column = "hilic", 
                              candidate.num = 3,
                              database = "HMDB.metabolite.data", 
                              threads = 4)
ms1.data <- grep("\\.csv", dir(), value = TRUE)
ms2.data <- grep("mgf", dir(), value = TRUE)
result.nHILIC.nce25 <- metIdentify4all(ms1.data = ms1.data,
                                       ms2.data = ms2.data, 
                                       parameter.list = c(parameter1,
                                                          parameter2,
                                                          parameter3,
                                                          parameter4,
                                                          parameter5,
                                                          parameter6,
                                                          parameter7,
                                                          parameter8),
                                       path = ".")

save(result.nHILIC.nce25, file = "result.nHILIC.nce25")
print(1)









