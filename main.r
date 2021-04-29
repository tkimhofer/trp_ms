### import targeted LC-MS data (biogenic amines or trp assay) 
### reshape vendor data format to 2D matrix
### QC: check missing values, coeffVar
### NA imputation and PCA
### correct for batch and drift effects using QC/LTR
### check CV
# tkimhofer @ Murdoch Uni (V1, March 2021)

source('hlp_fct.r')

# get name and directory of files to read-in
path=[...] # change to directory containing tsv files
fil=list.files(path, full.names = T, pattern=[...]) # change file naming pattern, e.g. to `pattern='study_id'`

# import tryptophane data (tab format)
trp=import_trp(fil)

# import amino acid data (excel file format)
# trp=import_aa(fil)


# check number of missing values per metabolite
cs=apply(trp, 2, function(x){c('TRUE'=length(which(is.na(x))), 'FALSE'=length(which(!is.na(x))))})
nas=as.data.frame(cs, stringsAsFactors = F)
nas$missingValues=rownames(nas)
nas=melt(nas, id.vars = 'missingValues')
ggplot(nas, aes(variable, value, fill=missingValues))+geom_bar(stat='identity')+coord_flip()+theme_bw()+labs(y=NULL, x=NULL)+scale_y_continuous(sec.axis = sec_axis(~./nrow(trp)*100, name='%'),name='n')

# remove standards and impute missing values with column mean
idx=grep('SIL', colnames(trp), invert = T)
trp=apply(trp[,idx], 2, function(x){x[is.na(x)]=mean(x,na.rm = T); x})

# perform pca 
mod=pca(trp, scale = 'UV', center = T)
g1=plotscores(mod, qc=grep('LTR', rownames(trp)))
g2=plotload_cat(mod)
grid.arrange(g1, g2, ncol=2)



# coefficient of variation in LTR samples
idx_ltr<-grep('LTR', rownames(trp))
cvs=cvar(trp, idx_ltr, cv_thres = 30)
ccv=data.frame(cvs, stringsAsFactors = F)
ccv$id=rownames(ccv)
ggplot(ccv, aes(id, cvs))+geom_bar(stat='identity')+coord_flip()

# outlined below is the procedure for correcting batch and run-order effects
# define run order
spl=strsplit(rownames(trp), '_')
ro=rank(as.numeric(gsub('.*p', '', sapply(spl, '[[', 6)))*100 + as.numeric(sapply(spl, function(x){x[length(x)]})))
rack=as.numeric(gsub('.*p', '', sapply(spl, '[[', 6)))

# visualise feature stability w/o correction
an=data.frame(ro, ltr=grepl('LTR', rownames(trp)), qc=grepl('QC', rownames(trp)), id=rownames(trp),stringsAsFactors = F)
an$ra=rack

dt=melt(data.frame(trp, an), id.vars = colnames(an))
ggplot(dt, aes(ro, value, shape=ltr, color=factor(ltr)))+geom_point(shape=1)+geom_point(data=dt[dt$ltr==T,])+facet_wrap(.~variable, scales='free_y')+scale_y_continuous(trans='log10', name='Conc')+theme_bw()+labs(x='Run Order')

ggplot(dt, aes(ro, value, shape=ltr, color=factor(ra)))+geom_point(shape=1)+geom_point(data=dt[dt$ltr==T,])+facet_wrap(.~variable, scales='free_y')+scale_y_continuous(trans='log10', name='Conc')+theme_bw()+labs(x='Run Order')

# perform batch correction using rack information
trp_bc=batch_cor(trp, idx_qc=idx_ltr, idx_batch=an$ra, qc_mean=T)

# visualise feature stability after batch-correction
dt_bc=melt(data.frame(trp_bc, an), id.vars = colnames(an))
ggplot(dt_bc, aes(ro, value, shape=ltr, color=ltr))+geom_point(shape=1)+geom_point(data=dt_bc[dt_bc$ltr==T,])+facet_wrap(.~variable, scales='free_y')+scale_y_continuous(trans='log10')+theme_bw()

# perform drift correction using run-order information
trp_bdc=drift_cor_loess(trp_bc, run_order=rank(an$ro), idx_qc=an$ltr, smoothing=0.4)

# visualise feature stability after drift-correction
dt_dbc=melt(data.frame(trp_bdc, an), id.vars = colnames(an))
ggplot(dt_dbc, aes(ro, value, shape=ltr, color=ltr))+geom_point(shape=1)+geom_point(data=dt_bc[dt_bc$ltr==T,])+facet_wrap(.~variable, scales='free_y')+scale_y_continuous(trans='log10')+theme_bw()


# coefficient of variation (LTR samples) of corrected data
cvs_bdc=cvar(trp_bdc, idx_ltr)


