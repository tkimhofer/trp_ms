library(metabom8)
library(plyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(scales)
library(readxl)


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
    ds=read.table(fid, skip = 0, fill = T, sep='\t', comment.char = '', stringsAsFactors = F)
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



#############
# function to import targeted MS data (amino acids)
#############
import_aa<-function(fil){
  require(readxl)
  require(ggplot2)
  require(scales)
  require(reshape2)
  # fil: path to files to be imported
  # output: 2d data matrix of conc. with samples in rows and variables in cols
  # (c) T Kimhofer, V1 (03/2021)
  
  # for each file
  comp=lapply(fil, function(fid){
    print(fid)
    # read and extract compounds
    ds=suppressMessages(data.frame(read_xlsx(fid), stringsAsFactors = F))
    ds=ds[which(apply(ds, 1, function(x){length(which(is.na(x)))})<ncol(ds)),]
    
    # filter for analytes
    idx=grep('[A-Z]+[0-9]+$|[0-9]+[A-Z]+$|Acc.?QTag$', ds$Analyte.Name, value=F, invert = T)
    ds=ds[idx,]
    
    # remove QC and blank
    idx=grep('QC|Blank', ds$Data.Set, ignore.case = T, invert = T)
    ds=ds[idx,]
    
    # get non-numeric values
    idx1=suppressWarnings(which(is.na(as.numeric(ds$Quantity..units.))))
    if(length(idx1)>0){
      message('Removing non-numeric values from column \"Quantity..units\":')
      print(table(ds$Quantity..units.[idx1]))
      cat('\n')
    }
    
    comp_info=ds
    
    return(comp_info)
  })
  
  comp=do.call(rbind, comp)
  
  comp_un=unique(comp)
  if(nrow(comp_un)!=nrow(comp)){comp=comp_un; message('There are double entries, discarding doublets.')}
  
  g1=ggplot(comp, aes(as.numeric(RT..min.), as.numeric(m.z.meas.), colour=Analyte.Name))+
    geom_point(shape=2)+theme_bw()+labs(x='rt (min)', y='m/z')+
    theme(legend.position = 'bottom')+
    scale_x_continuous(sec.axis = sec_axis(~.*60, name='rt (s)', breaks = pretty_breaks()))
  
  suppressWarnings(plot(g1))
  
  dd=dcast(comp, Data.Set~Analyte.Name, value.var = 'Quantity..units.')
  mat=suppressWarnings(apply(dd[,-1], 2, as.numeric))
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
  #browser()
  grand_x=grand_median(X)
  batch=data.frame(batch=idx_batch, idx=1:length(idx_batch), type='Sample', stringsAsFactors = F)
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
  colnames(out)=colnames(X)
  return(out)
}




read_aa_V1=function(fils, filter=T, plotting=T, interactive=T){
  require(ggplot2)
  require(reshape2)
  require(scales)
  #fils=rep('/Users/TKimhofer/Downloads/PAI-05_Barwon_AA_Plate1.TSV', 2)
  dats=lapply(fils, function(x){
    cat(x);cat('\n')
    rfile=read.delim(x, sep = '\t', header = T, fileEncoding="latin1")
    diag=rfile[,colnames(rfile) %in% c('m.z', 'Retention.Time.min.', 'AnalyteName', 'AnalysisName')]
    
    ds=rfile[,colnames(rfile) %in% c('AnalysisName', 'Quantity.w..unit', 'AnalyteName')]
    ds$Quantity.w..unit=as.numeric(gsub(' units', '', ds$Quantity.w..unit))

    if(filter){
      # filter for analytes
      idx1=!grepl('[A-Z]+[0-9]+$|[0-9]+[A-Z]+$|Acc.?QTag$|\\[IS\\]', ds$AnalyteName)
      idx2=!grepl('Cal|QC|Blank', ds$AnalysisName)
      ds=ds[(idx1 & idx2),]
    }
    
    sub=ds[,colnames(ds) %in% c('AnalysisName', 'AnalyteName')]
    if(nrow(sub)!=nrow(unique(sub))){
      iid=paste(sub$AnalysisName, sub$AnalyteName, sep='-')
      idx_rm=which(iid %in% names(which(table(iid)>1)))
      message('The following IDs have more than one entry in the imported data file: ')
      print(iid[idx_rm[-1]])
      message('Removing doublicates for further processing.')
      ds=ds[-idx_rm,]
    }
    
    # remove cystine 
    length(table(ds$AnalyteName))
    idx=grepl('^Cystine', ds$AnalyteName, ignore.case = T)
    if(any(idx)){message('Removing cystein from assay'); ds=ds[!idx,]}
    
    dout=dcast(ds, AnalysisName~AnalyteName, value.var='Quantity.w..unit')
    print(dim(dout))
    cat('\n')
    
    # if(any(is.na(dout))){browser()}
    # print(length(which(is.na(ds$Quantity.w..unit))))
    return(list(dout, diag))
  })

  out=do.call(rbind, lapply(dats, '[[', 1))

  id=out$AnalysisName
  if(length(id)!=length(unique(id))){cat('File names are not unique - appending index.'); id=paste0(id, 1:length(id))}
  
  #browser()
  rownames(out)=id
  idx_ckeep=which(!colnames(out) %in% 'AnalysisName')
  out=out[,idx_ckeep]
  
  
  if(plotting){
    comp=do.call(rbind, lapply(dats, '[[', 2))
    
    if(filter){
      # filter for analytes
      idx=grep('[A-Z]+[0-9]+$|[0-9]+[A-Z]+$|Acc.?QTag$|\\[IS\\]', comp$AnalyteName, value=F, invert = T)
      comp=comp[idx,]
    }
    
    
    
    if(!interactive){
      g1=ggplot(comp, aes(as.numeric(Retention.Time.min.), as.numeric(m.z), colour=AnalyteName))+
        geom_point(shape=2)+theme_bw()+labs(x='rt (min)', y='m/z')+
        theme(legend.position = 'bottom')+
        scale_x_continuous(sec.axis = sec_axis(~.*60, name='rt (s)', breaks = pretty_breaks()))
      
      suppressWarnings(plot(g1))
    }else{
        require(plotly)
      comp$m.z=as.numeric(comp$m.z)
      comp$Retention.Time.min.=as.numeric(comp$Retention.Time.min.)
      cmap=c('mz'='m.z', 'rt'='Retention.Time.min.', 'compound'='AnalyteName', 'fid'='AnalysisName')
      colnames(comp)=names(cmap)[match(colnames(comp), cmap)]

      fig <- suppressWarnings(plot_ly(comp, x = ~rt, y = ~mz, text = ~paste(compound, '<br>Sample:', fid), color = ~compound, marker=list(), type="scatter", mode='markers'))

      (suppressWarnings(fig))
      }
    
      }else{fig=NULL}

  return(list(data=out, figure=fig))
  
}




# drift correction function
# X: feature matrix without batch effects
# run_order: run order of samples represented in X
# idx_qc: index of QC samples
# smoothing: granularity of correction (0.1-1)
drift_cor_loess=function(X, run_order, idx_qc, smoothing=0.2){
  
  if(!all(run_order %in% 1:nrow(X))){stop('Check run order')}
  
  X[X==0]=NA
  batch=data.frame(ro=run_order, idx=1:nrow(X), type='Sample', stringsAsFactors = F)
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


cvar<-function(X, qc_idx, plot=T, cv_thres=20){
  
  cv<-apply(X[qc_idx,], 2, function(x){
    sd(x, na.rm=T)/mean(x, na.rm=T)*100
  })
  if(plot){
    df=data.frame(id=colnames(X), cv, cv_pass=cv<cv_thres, stringsAsFactors = F)
    if(nrow(df)>100){message('Large number of features! Plotting first 100.');df=df[1:100,]}
    g1=ggplot(df, aes(id, cv, fill=cv_pass))+geom_bar(stat='identity', alpha=1)+theme_bw()+labs(y='CV (%)', x='Feature')+geom_hline(yintercept=cv_thres, col='black', linetype=2)+coord_flip()+guides(fill=F)+scale_fill_manual(values=c('TRUE'='#79FFFE', 'FALSE'='#FF8B8B'))+theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
    plot(g1)
  }
  
  return(cv)
  
}

