#load .txt file of state data
statedata <- read.delim("IRK state data.txt")

#load the dplyr package
library(dplyr)

#generate new dataframes, one for each exposure
expo10sec <- filter(statedata, Exposure == 0.167)
expo1min <- filter(statedata, Exposure == 1.000)
expo10min <- filter(statedata, Exposure == 10.000)
expo1hr <- filter(statedata, Exposure > 59 & Exposure < 61)
expo2hrs <- filter(statedata, Exposure > 119 & Exposure < 121)

exposures <- c("expo10sec", "expo1min", "expo10min", "expo1hr", "expo2hrs")

#create dataframe df with all exposures for all peptides
df <- data.frame("Start" = expo10sec["Start"], "End" = expo10sec['End'], expo10sec["BECUptake"], expo1min["BECUptake"], expo10min["BECUptake"], expo1hr["BECUptake"], expo2hrs["BECUptake"])
#rename column names
colnames(df)[colnames(df)=="BECUptake"] <- exposures[1]
colnames(df)[colnames(df)=="BECUptake.1"] <- exposures[2]
colnames(df)[colnames(df)=="BECUptake.2"] <- exposures[3]
colnames(df)[colnames(df)=="BECUptake.3"] <- exposures[4]
colnames(df)[colnames(df)=="BECUptake.4"] <- exposures[5]

#create new dataframe perresiavgs with all residues and all time points
resis <- min(df['Start']):max(df['End'])
perresiavgs <- data.frame("Residues" = resis, "Expo10sec" = NA, "Expo1min" = NA, "Expo10min" = NA, "Expo1hr" = NA, "Expo2hrs" = NA )

for (a in 1:5) {

  for (i in 1:nrow(perresiavgs)) {
    count <- 0
    sum <- 0
    
    for (j in 1:nrow(df)) {
        if (resis[i] >= df[j,'Start'] & resis[i] <= df[j,'End']) {
        sum <- sum + df[j,a+2]
        count <- count + 1
        }
    }
    avg <- sum / count
    perresiavgs[i,a+1] = avg
  }
}

#change the residue numbers so that they align with residue numbers from 1irk.pdb
perresiavgspdb <- data.frame(perresiavgs)
perresiavgspdb[[1]] <- perresiavgspdb[[1]] - 27

#replace all NAs with zero -- these represent regions with no peptide coverage
perresiavgspdb[is.na(perresiavgspdb)] <- 0

#create dataframes for each time point
irk10sec <- data.frame(perresiavgspdb[1], perresiavgspdb[2])
irk1min <- data.frame(perresiavgspdb[1], perresiavgspdb[3])
irk10min <- data.frame(perresiavgspdb[1], perresiavgspdb[4])
irk1hr <- data.frame(perresiavgspdb[1], perresiavgspdb[5])
irk2hrs <- data.frame(perresiavgspdb[1], perresiavgspdb[6])

#generate .txt files in the format for input for the python color as b factor script for each time point
write.table(irk10sec, file = "irk10sec.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(irk1min, file = "irk1min.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(irk10min, file = "irk10min.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(irk1hr, file = "irk1hr.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(irk2hrs, file = "irk2hrs.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

