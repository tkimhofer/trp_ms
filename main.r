### import tryptophane pw - tarteted ms data from vendor format
### perform unsupervised analysis
### correct for batch and drift effects detected by LTR
# tkimhofer @ Murdoch Uni (V1, March 2021)

source('hlp_fct.r')

# get name and directory of files to read-in
path='/path/to/dir/' # change this to desired directory
fil=list.files(path, full.names = T, pattern='file_pattern') # change file pattern, e.g. to `pattern='barwon'`

# import data (rows = samples, columns=metabolites)
trp=import_trp(fil)

# check number of missing values per metabolite
nas=data.frame(apply(trp, 2, function(x){table(is.na(x))}))
nas$missingValues=rownames(nas)
nas=melt(nas, id.vars = 'missingValues')
ggplot(nas, aes(variable, value, fill=missingValues))+geom_bar(stat='identity')+coord_flip()+theme_bw()+labs(y=NULL, x=NULL)

# remove standards and impute missing values with column mean
idx=grep('SIL', colnames(out), invert = T)
trp=apply(trp[,idx], 2, function(x){x[is.na(x)]=mean(x,na.rm = T); x})

# perform pca 
mod=pca(trp)
g1=plotscores(mod, qc=grep('LTR', rownames(trp)))
g2=plotload_cat(mod)
grid.arrange(g1, g2, ncol=2)




