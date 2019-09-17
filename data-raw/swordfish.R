


library(MSEtool)

swo <- read.csv("data/datalist/swordfish.csv")


swordfish <- new("Data")
swordfish@Name <- "North Atlantic Swordfish (Source: ASPIC software)"
swordfish@Year <- swo$Year
swordfish@Cat <- matrix(swo$Catch, nrow = 1)
swordfish@Ind <- matrix(swo$Index, nrow = 1)

usethis::use_data(swordfish, overwrite = TRUE)
