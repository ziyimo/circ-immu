args <- commandArgs(trailingOnly=TRUE)
inf <- args[1]
load(inf)
print(ss_fit)
print(seed_sh)
#myrange <- 1:length(ss_fit)
#for (state in myrange) {
#
#    cat(paste(seed_sh,collapse=":"),ss_fit[[state]],names(ss_fit[state]),"\n")
#}
