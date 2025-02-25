#Runs SSF with contrasum including reverse method to get all parameters
#outputs parms and/or model with aic

#Dat_T- data
#NLCDT- table of NLCD levels
#set.colname- colname used for strata in clogit model
#ssf_opt- option to run by nlcd alone or with sl: "nlcd_only" or "nlcd_sl"
#out_opt- output options: "parms", "model", "modelandparms"
#set.colname="animal_p_stp"
#ssf_opt="nlcd_only"
#out_opt="modelandparms"
contrSum_iSSF<-function(Dat_T,NLCDT,set.colname,ssf_opt,out_opt){
  #drop factor levels
  Dat_T<-droplevels(Dat_T)
  
  colnames(Dat_T)[which(colnames(Dat_T)==set.colname)]<-"set"
  
  #get summaries of step/lc counts
  #NLCDT<-as.data.frame(table(Dat_T$nlcd,Dat_T$case_,dnn=c("nlcd","case")))
  #NLCDT<-subset(NLCDT,NLCDT$case=="TRUE"&NLCDT$Freq>0)
  NLCD_ct<-dplyr::count(NLCDT)
  
  #if any LC vals in subset (should always be true, stop if not)
  if(NLCD_ct$n==0) {stop("no land classes in subset?")}
  
    #Get nlcd frequencies/levels to do contra sums
    Dat_T2<-merge(Dat_T,NLCDT, by="nlcd")
    Dat_T2<-droplevels(Dat_T2)
    Dat_T_out<-Dat_T2 #for write out
    
    #use NLCD_ct count to get contra sum for nlcd
    contrasts(Dat_T2$nlcd) =contr.sum(as.numeric(NLCD_ct$n))
  
 
  if(ssf_opt=="nlcd_only"){
    T1<-survival::clogit(case_ ~ nlcd_rev + strata(set), data=Dat_T2)
   
    AICT<-AIC(T1)
    #pull coefficients from clogit output
    sumT<-summary(T1)$coefficients
    #get ordered df of contrasts
    contr=contrasts(Dat_T2$nlcd)
    #same order, needs be rownames of coefficients
    #last one was dropped, remove
    rownames(sumT)=rownames(contr)[1:(nrow(contr)-1)]
    
    
    #Get reverse order NLCD table
    NLCDTr<-as.data.frame(table(Dat_T$nlcd_rev,Dat_T$case_,dnn=c("nlcd_rev","case")))
    NLCDTr<-subset(NLCDTr,NLCDTr$case=="TRUE"&NLCDTr$Freq>0)
    NLCD_ctr<-dplyr::count(NLCDTr)
    #Get frequency count for contrasums
    Dat_T2<-merge(Dat_T,NLCDTr, by="nlcd_rev")  
    #drop levels
    Dat_T2<-droplevels(Dat_T2)
    
    #Get contrasts
    contrasts(Dat_T2$nlcd_rev) =contr.sum(as.numeric(NLCD_ctr$n))
    #Run clogot model, reversed contrasums
    T1<-survival::clogit(case_ ~ nlcd_rev + strata(set), data=Dat_T2)
    #pull coefficients
    sumTrev<-summary(T1)$coefficients
    #only need the first row that was dropped in original model
    first_row <- sumTrev[1,,drop=FALSE]
    #combine together, get all coefs
    Sum2<-rbind(sumT,first_row)
    #get lc names from contrast ordering
    contr=contrasts(Dat_T2$nlcd)
    rownames(Sum2)=rownames(contr)[1:(nrow(contr))]

    if(out_opt=="parms"){
    out=list("params"=Sum2,"aic"=AICT)
    }
    
    if(out_opt=="model"){
      out=list("model"=T1,"dat"=Dat_T_out)
    }
    
    if(out_opt=="modelandparms"){
      out=list("model"=T1,"dat"=Dat_T_out,"params"=Sum2,"aic"=AICT)
    }
    
    }
  
  if(ssf_opt=="nlcd_sl"){
    T1<-survival::clogit(case_ ~ nlcd*sl_ + strata(set), data=Dat_T2)
  
    AICT<-AIC(T1)
    #pull coefficients from clogit output
    sumT<-summary(T1)$coefficients
    #get ordered df of contrasts
    contr=contrasts(Dat_T2$nlcd)
    #same order, needs be rownames of coefficients
    #last one was dropped, remove
    #get correct sequence of rownames   
    sl_name_seq=c(rownames(contr)[1:(nrow(contr)-1)],"sl_",paste0(rownames(contr)[1:(nrow(contr)-1)],":sl_"))
    
    #same order, needs be rownames of coefficients
    #last one was dropped, remove
    rownames(sumT)=sl_name_seq    
    
    #Get reverse order NLCD table
    NLCDTr<-as.data.frame(table(Dat_T$nlcd_rev,Dat_T$case_,dnn=c("nlcd_rev","case")))
    NLCDTr<-subset(NLCDTr,NLCDTr$case=="TRUE"&NLCDTr$Freq>0)
    NLCD_ctr<-dplyr::count(NLCDTr)
    #Get frequency count for contrasums
    Dat_T2<-merge(Dat_T,NLCDTr, by="nlcd_rev")  
    #drop levels
    Dat_T2<-droplevels(Dat_T2)
    #Get contrasts
    contrasts(Dat_T2$nlcd_rev) =contr.sum(as.numeric(NLCD_ctr$n))
    
    
    #Run clogot model, reversed contrasums
    T1<-survival::clogit(case_ ~ nlcd_rev*sl_ + strata(set), data=Dat_T2)
    #pull coefficients
    sumTrev<-summary(T1)$coefficients
    #only need the first row that was dropped in original model
    first_row <- sumTrev[1,,drop=FALSE]
    #combine together, get all coefs
    Sum2<-rbind(sumTrev,first_row)
    #get lc names from contrast ordering
    contr=contrasts(Dat_T2$nlcd_rev)
    
    #get correct sequence of rownames   
    sl_name_seq=c(rownames(contr)[1:nrow(contr)],"sl_",paste0(rownames(contr)[1:(nrow(contr)-1)],":sl_"))
    
    #same order, needs be rownames of coefficients
    #last one was dropped, remove
    rownames(Sum2)=sl_name_seq
    
    if(out_opt=="parms"){
    out=list("params"=Sum2,"aic"=AICT)
    }
    
    if(out_opt=="model"){
      out=list("model"=T1,"dat"=Dat_T_out)
    }
    
    if(out_opt=="modelandparms"){
      out=list("model"=T1,"dat"=Dat_T_out,"params"=Sum2,"aic"=AICT)
    }
    
    }
  
  return(out)
}
  









