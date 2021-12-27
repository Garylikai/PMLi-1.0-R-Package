## code to prepare `DATASET` dataset goes here
data0 <- read.table(
  "http://www.ams.sunysb.edu/~pfkuan/Teaching/AMS597/Data/d_logret_6stocks.txt", 
   header = TRUE)
pfizer <- unlist(data0$Pfizer)
intel <- unlist(data0$Intel)
matrix1 <- cbind(pfizer, intel)
citigroup <- unlist(data0$Citigroup)
amerexp <- unlist(data0$AmerExp)
matrix2 <- cbind(citigroup, amerexp)
exxon <- unlist(data0$Exxon)
genmotor <- unlist(data0$GenMotor)
matrix3 <- cbind(exxon, genmotor)
data1 <- data.frame(rbind(matrix1, matrix2, matrix3))
names(data1) = c("Sample1", "Sample2")
data1[101:150, 2] <- NA
data1[151:190, 1] <- NA
data1[191:192,] <- NA

pm <- data1

usethis::use_data(pm, overwrite = TRUE)
