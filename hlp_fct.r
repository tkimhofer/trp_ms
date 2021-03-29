
#############
# function to import targeted MS data
#############
import_trp<-function(fil){
  # fil: path to files to be imported
  # output: 2d data matrix of conc. with samples in rows and variables in cols
  # (c) T Kimhofer, V1 (03/2021)
  
  # for each file
  comp=lapply(fil, function(fid){
    print(fid)
    # read and extract compounds
    ds=read.table(fid, skip = 0, fill = T, sep='\t', comment.char = '')
    idx=grep('compound [0-9]', ds$V1, ignore.case = T, value = F)
    comp=gsub('Compound [0-9].?: .?', '', ds$V1[idx])
    
    # collate individual compound information
    comp_info=lapply(seq(length(idx)), function(m){
      #browser()
      i=idx[m]
      iid=grep("^[0-9].?", ds$V1[i:nrow(ds)], invert = T)
      iid1=grep("^$", ds$V3[i+iid[2]:nrow(ds)], invert = F)
      if(i==idx[length(idx)]){;out=ds[seq((i+iid[2]), nrow(ds)),]}else{
        out=ds[seq((i+iid[2]), (i+iid1[1]-1)),]
      }
      #rm na
      #idx_notna=which(apply(out, 1, function(x){length(which(is.na(x)))})<= round(ncol(out)/3))
      #out[idx_notna,]
      colnames(out)=ds[i+iid[2]-1,]
      out$compound=comp[m]
      out
    })
    
    comp_info=do.call(rbind, comp_info)
    comp_info=comp_info[comp_info$Type=='Analyte',] # exclude calibration data
    
    return(comp_info)
  })
  
  comp=do.call(rbind, comp)
  
  comp_un=unique(comp)
  if(nrow(comp_un)!=nrow(comp)){comp=comp_un; message('There are double entries, discarding doublets.')}
  
  #browser()
  # reshape
  dd=dcast(comp, Name~compound, value.var = 'Conc.')
  mat=apply(dd[,-1], 2, as.numeric)
  rownames(mat)=dd$Name
  
  return(mat)
  
}


# batch correction targeted ms

# helper function batch correction 
# input 
# x: feature vector
# idx_qc: index of QC samples in x
# grand: grand median (in case of multiple batches), 
# trim: trimmed mean quantile
# qc_mean: boolean qc step fct (in case of multiple batches)
batch_substr=function(x, idx_qc, grand, trim, qc_mean){
  # browser()
  if(qc_mean){x_qc_new=mean(x[idx_qc], trim = trim, na.rm = T)}else{
    x_qc_new=mean(x, trim = trim, na.rm = T)
  }
  return(x+(grand-x_qc_new))
}

# # helper function batch correction 
# X: feature matrix
grand_median = function(X){
  apply(X, 2, median, na.rm=T)
}

#batch correction fct
# X: feature matrix
# idx_qc: row-index indicating QC samples
# idx_batch: batch definition as numeric array, with length equals nrow(X)
batch_cor=function(X, idx_qc, idx_batch, qc_mean=T){
  
  grand_x=grand_median(X)
  batch=data.frame(batch=idx_batch, idx=1:length(idx_batch), type='Sample')
  batch$type[idx_qc]='QC'
  
  out=dlply(batch, .(batch), function(b){
    idx_qc_batch=which(b$type=='QC')
    b_cor=sapply(seq(ncol(X)), function(i, iid=idx_qc_batch){
      batch_substr(X[b$idx,i], idx_qc=iid, grand_x[i], trim=0.2,  qc_mean)
    })
    rownames(b_cor)=b$idx
    return(list(b, b_cor))
  })
  batch_reord=do.call(rbind, lapply(out, '[[', 1))
  out_reord=do.call(rbind, lapply(out, '[[', 2))
  
  idx=order(batch_reord$idx)
  out=out_reord[idx,]
  return(out)
}

# drift correction function
# X: feature matrix without batch effects
# run_order: run order of samples represented in X
# idx_qc: index of QC samples
# smoothing: granularity of correction (0.1-1)
drift_cor_loess=function(X, run_order, idx_qc, smoothing=0.2){
  
  if(!all(run_order %in% 1:nrow(X))){stop('Check run order')}
  
  X[X==0]=NA
  batch=data.frame(ro=run_order, idx=1:nrow(X), type='Sample')
  batch$type[idx_qc]='QC'
  batch$ro_qc=NA
  batch$ro_qc[idx_qc]=rank(batch$ro[idx_qc])
  #browser()
  b_cor=apply(X, 2, function(x, iid=idx_qc){
    #browser()
    mod=loess(x[iid]~batch$ro[iid], span=smoothing, degree = 2)
    
    #
    # test=predict(mod, newdata=as.numeric(x[iid]))
    # plot(batch$ro_qc[iid], x[iid])
    # points(batch$ro_qc[iid], mod$fitted, col='red', type='l')
    #
    #
    x_pred=predict(mod, batch$ro)
    # browser()
    #points(batch$ro_qc, x_pred, col='cyan')
    gres=mean(x_pred, na.rm = T)-x_pred
    x_new=(x + gres)
    x_new
  })
  
  return(b_cor)
}

