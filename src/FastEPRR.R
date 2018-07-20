
library(FastEPRR)

start=1
end = 20000
outf = paste('./rho_',start,'_',end,'fprr.txt',sep='')
FastEPRR_ALN(alnFilePath="/data/EBV/msas/up150k/all150.fas",format = 1, outputFilePath=outf,erStart=as.character(start),erEnd=as.character(end))


start = 1
end = 100000
wl= 20000
FastEPRR_ALN(alnFilePath="/data/EBV/msas/up150k/all150.fas",format = 1, outputFilePath=outf,erStart=as.character(start),erEnd=as.character(end),winLength=as.character(wl))
