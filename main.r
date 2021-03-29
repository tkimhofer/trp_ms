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



# get run order
# there is one re-run 46r on plate 02 -> remove sample 46
spl=strsplit(rownames(trp), '_')
ro=rank(as.numeric(gsub('BARWINp', '', sapply(spl, '[[', 6)))*100 +as.numeric(sapply(spl, function(x){x[length(x)]})))

plro=sapply(spl,function(x){x[length(x)]})
ra=as.numeric(gsub('BARWINp', '', sapply(spl, '[[', 6)))*100
idx=grep('[aA-zZ]', plro)
iid=as.numeric(gsub('[aA-zZ]', '', plro[idx]))+ra[idx]

idx_keep=which(!ro %in% iid)

ro=rank(as.numeric(gsub('BARWINp', '', sapply(spl, '[[', 6)))*100 +as.numeric(gsub('[aA-zZ]', '', sapply(spl, function(x){x[length(x)]}))))
idx_keep=which(!ro %in% iid)

ro1=ro[idx_keep]
trp1=trp[idx_keep,]
ra1=ra[idx_keep]

an=data.frame(ro, ltr=grepl('LTR', rownames(trp)), qc=grepl('QC', rownames(trp)), id=rownames(trp),ra ,stringsAsFactors = F)
trp=data.frame(trp, an)
dt=melt(trp, id.vars = colnames(an))


ggplot(dt, aes(ro, value, shape=ltr, color=factor(ltr)))+geom_point(shape=1)+geom_point(data=dt[dt$ltr==T,])+facet_wrap(.~variable)+scale_y_continuous(trans='log10')+theme_bw()


X=trp[, !colnames(trp) %in% colnames(an)]

# batch correction
idx=grep('LTR', rownames(trp))
test=batch_cor(X, idx_qc=idx, idx_batch=an$ra, qc_mean=T)
colnames(test)=colnames(X)

# vis
dt1=melt(data.frame(test, an), id.vars = colnames(an))
ggplot(dt1, aes(ro, value, shape=ltr, color=ltr))+geom_point(shape=1)+geom_point(data=dt1[dt1$ltr==T,])+facet_wrap(.~variable)+scale_y_continuous(trans='log10')+theme_bw()

# drift correction
test1=drift_cor_loess(test, run_order=rank(an$ro), idx_qc=an$ltr, smoothing=0.2)













