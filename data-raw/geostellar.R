## code to prepare `geostellar` dataset

geostellar <- read.csv(file = "inst/extdata/geo_stellar_chrono_scale_2018.9.20.csv",
utils::data(strat2012)
had <- strat2012[rep(116, each=6),]
rownames(had) <- 116:121
had$era <- c("Neohadean", "Mesohadean", "Palaeohadean")
had$period <- c("Prometean", "Acastan", "Procrustean", "Canadian", "Jacobian", "Hephestean")
had$MA <- c(3900, 4000, 4100, 4200, 4300, 4400)
chao <- strat2012[rep(116, each=4),]
rownames(chao) <- 122:125
chao$eon <- rep("Chaotian", 4)
chao$era <- c(rep("Neochaotian", 2), rep("Eochaotian", 2))
chao$period <- c("Titanomachean", "Hyperitian", "Erebrean", "Nephelean")
chao$MA <- c(4500, 4560, 4567, 4730)
usethis::use_data(geostellar)
