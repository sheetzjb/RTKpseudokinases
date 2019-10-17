#run avgperresi.R script prior to running this script

perresiavgspdbROR1 <- read.delim("perresiavgspdbROR1.txt")
perresiavgspdbROR2 <- read.delim("perresiavgspdbROR2.txt")
perresiavgspdbRYK <- read.delim("perresiavgspdbRYK.txt")
perresiavgspdbPTK7 <- read.delim("perresiavgspdbPTK7.txt")

perresiavgspdb[is.na(perresiavgspdb)] <- 0
perresiavgspdbROR1[is.na(perresiavgspdbROR1)] <- 0
perresiavgspdbROR2[is.na(perresiavgspdbROR2)] <- 0
perresiavgspdbRYK[is.na(perresiavgspdbRYK)] <- 0
perresiavgspdbPTK7[is.na(perresiavgspdbPTK7)] <- 0

allpseudoavg <- data.frame("Residues" = 983:1276, "expo10sec" = NA, "expo1min" = NA, "expo10min" = NA, "expo1hr" = NA, "expo2hrs" = NA)
allpseudoSD <- data.frame("Residues" = 983:1276, "expo10sec" = NA, "expo1min" = NA, "expo10min" = NA, "expo1hr" = NA, "expo2hrs" = NA)
allpseudoavg[allpseudoavg == 0] <- NA
irkdevpseudo <- data.frame("Residues" = 978:1276, "expo10sec" = 0, "expo1min" = 0, "expo10min" = 0, "expo1hr" = 0, "expo2hrs" = 0)

for (a in 2:6) {

  for (i in 1:nrow(perresiavgspdb)) {
    sum <- 0
    count <- 0
    SD <- 0
    difror1 <- 0
    difror2 <- 0
    difptk7 <- 0
    difryk  <- 0
    
    sum <- sum + perresiavgspdbROR1[i,a] + perresiavgspdbROR2[i,a] + perresiavgspdbRYK[i,a] + perresiavgspdbPTK7[i,a]
    if (perresiavgspdbROR1[i,a] != 0) {
      count <- count + 1
    }
    if (perresiavgspdbROR2[i,a] != 0) {
      count <- count + 1
    }
    if (perresiavgspdbRYK[i,a] != 0) {
      count <- count + 1
    }
    if (perresiavgspdbPTK7[i,a] != 0) {
      count <- count + 1
    }
    if (count != 0) {
      avg <- sum / count
    } else {
      avg <- 0
    }
    allpseudoavg[i,a] <- avg
    
    if (perresiavgspdbROR1[i,a] != 0) {
      difror1 <- perresiavgspdbROR1[i,a] - avg
    }
    if (perresiavgspdbROR2[i,a] != 0) {
      difror2 <- perresiavgspdbROR2[i,a] - avg
    }
    if (perresiavgspdbRYK[i,a] != 0) {
      difptk7 <- perresiavgspdbPTK7[i,a] - avg
    }
    if (perresiavgspdbPTK7[i,a] != 0) {
      difryk <- perresiavgspdbRYK[i,a] - avg
    }
    if (count != 0) {
      SD <- sqrt(((difror1)^2 + (difror2)^2 + (difptk7)^2 + (difryk)^2) / count)
    } else {
      SD <- 0
    }
    allpseudoSD[i,a] <- SD
    
  }
  for (j in 1:nrow(allpseudoavg)) {
    if (allpseudoSD[j,a] != 0 & perresiavgspdb[j,a] != 0) {
      irkdevpseudo[j+5,a] <- (allpseudoavg[j,a] - perresiavgspdb[j,a]) / allpseudoSD[j,a]
    } else {
      irkdevpseudo[j+5,a] <- 0
    }
    if (irkdevpseudo[j+5,a] < 2 & irkdevpseudo[j+5,a] > -2) {
      irkdevpseudo[j+5,a] <- 0
    }
  }
}


#create dataframes for each time point
irkdev10sec <- data.frame(irkdevpseudo[1], irkdevpseudo[2])
irkdev1min <- data.frame(irkdevpseudo[1], irkdevpseudo[3])
irkdev10min <- data.frame(irkdevpseudo[1], irkdevpseudo[4])
irkdev1hr <- data.frame(irkdevpseudo[1], irkdevpseudo[5])
irkdev2hrs <- data.frame(irkdevpseudo[1], irkdevpseudo[6])
irkdevavg <- data.frame(irkdevpseudo[1], "Average Index" = 0)

for (k in 1:nrow(irkdevpseudo)) {
  irkdevavg[k,2] <- sum(irkdevpseudo[k,2], irkdevpseudo[k,3], irkdevpseudo[k,4], irkdevpseudo[k,5], irkdevpseudo[k,6]) / 5
}


write.table(irkdev10sec, file = "index10sec.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(irkdev1min, file = "index1min.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(irkdev10min, file = "index10min.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(irkdev1hr, file = "index1hr.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(irkdev2hrs, file = "index2hrs.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(irkdevavg, file = "indexavg.txt", sep = "\t", row.names = FALSE, col.names = FALSE)