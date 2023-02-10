US_tmp <- read.csv("env_covar/US_temperature_2010_18.csv")
EU_tmp <- read.csv("env_covar/EU_temperature_2010_18.csv")
EU_tmp <- subset(EU_tmp, select = -MT)

quantile(c(unlist(US_tmp[, -1]), unlist(EU_tmp[, -1])))

EU_pop_tab <- read.csv("env_covar/EU19_pop_table.csv")
colnames(EU_pop_tab) <- c("code", "pop", "state")
write.table(EU_pop_tab, file = "env_covar/EU_pop.tsv", quote = FALSE, sep='\t', row.names = FALSE)
write.table(EU_pop_tab$code, file = "EU_states.tsv", quote = FALSE, sep='\t', row.names = FALSE, col.names = FALSE)
