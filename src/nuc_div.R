## this Rscript compute nucleotide diversity of a msa splited in regions


library(ape)
library(pegas)
fasta = '/data/EBV/msas/ids_171_ebv/ids_171_ebv_msa_gap1.51_trimmed0.05.fas'
region_file = '~/Projects/Gutierrez/EBV-recomb/recomb/rdp4_results/ids_171_ebv/10Miter/results/result_by_region_10Miters.csv'
gene_file = '~/Projects/Gutierrez/EBV-recomb/recomb/rdp4_results/ids_171_ebv/10Miter/results/result_by_gene_10Miters.csv'



regions = read.csv(region_file)
genes = read.csv(gene_file)
aln = ape::read.FASTA(fasta)
mat = as.matrix(aln)


nt_gc <-  function(m){
    #aln = as.DNAbin(m)
    gc = GC.content(m)
    nd = nuc.div(m,pairwise.deletion = F)
    res = c(gc,nd)
    names(res) = c('GC','NucDiv')
    return(res)
    }



    
# analysis by region
reg = lapply(1:nrow(regions),function(i){
x <<- regions[i,]
print(i)
submat = mat[,x$trimmed_st:x$trimmed_end]
return(nt_gc(submat))
})

reg = do.call('rbind',reg)
regions_with_div = cbind(regions,reg)
out_regions = paste(region_file,'_with_ndiv.csv',sep ='')
write.table(regions_with_div,file = out_regions,row.names = F,sep =',')


# analysis by gene
gen = lapply(1:nrow(genes),function(i){
x <<- regions[i,]
print(i)
submat = mat[,x$trimmed_st:x$trimmed_end]
return(nt_gc(submat))
})

gen = do.call('rbind',gen)
genes_with_div = cbind(genes,gen)
out_genes = paste(gene_file,'_with_ndiv.csv',sep ='')
write.table(genes_with_div,file = out_genes,row.names = F,sep =',')


#sum(is.na(genes_with_div$NucDiv))

    
