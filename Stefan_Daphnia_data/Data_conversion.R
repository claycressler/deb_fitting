for (geno in c('113','322','324','413','711','715','725','815','822')) {
    len <- read.csv(paste0(geno,'-A-','GROWTH.csv'), header=F, colClasses='character')
    egg <- read.csv(paste0(geno,'-A-','ASEXUAL.csv'), header=F, colClasses='character')
    age <- as.numeric(len[6:nrow(len),1])
    data <- vector(mode='list', length=(ncol(len)-1))
    for (i in 1:length(data)) {
    	Lobs <- as.numeric(len[6:nrow(len),i+1])
    	end <- which(Lobs < 0)-1
	if (length(end) < 1)
	   end <- min(which(is.na(Lobs)))-1
    	Lobs <- Lobs[1:end]*0.7/2.3
    	Robs <- as.numeric(egg[6:nrow(egg),i+1])
    	Robs <- cumsum(Robs[1:end])
    	data[[i]] <- as.data.frame(cbind(age=age[1:end],cbind(Lobs,Robs)))
    }
    save(data, file=paste0(geno,'-A_all_growth_egg_data.rda'))
}

