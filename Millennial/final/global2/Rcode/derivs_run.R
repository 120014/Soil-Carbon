
simmethod = c("simi", "wrb10","wrb9")[1]
if(clab){
  # fixCp_cold = paste0(fixCp,"_oldc")
  # fixCp_cplmwc = paste0(fixCp,"_plantc2lmwc")
  fixCp_cold = paste0(fixCp,"_oldc")
  fixCp_cplmwc = paste0(fixCp,"_plantc2lmwc")
}else{
  fixCp_cold = NULL
  fixCp_cplmwc = NULL
}
outname = c("year",fixCp,"Tco2",allflux.name,fixCp_c13,"Tco2_c13", fixCp_cold,"Toldco2", fixCp_cplmwc,"Tpi2LMWCco2")
derivnm = paste0(rep(outname,npar),".",rep(1:npar,each=length(outname)))

#提取第一年日尺度或月尺度数据
firstYear <- function(df0){
  if(ncol(df0) %% 12 == 0){
    cols_per_year <- 12
  }else if(ncol(df0) %% 365 == 0) {
    cols_per_year <- 365
  }
  df0[] = df0[,1:cols_per_year]
  return(df0)
}


if(length(which(forc_npp_all<0))>0){
    forc_npp_all[forc_npp_all<0] = 0
}

forc_npp_Oa_all = forc_npp_all * fexu_oa/100   #有机酸通量
forc_npp_Exu_all = forc_npp_all - forc_npp_Oa_all  #非有机酸部分 分泌物

rm(forc_npp_all)

# 所有驱动数据的"baseline2015"
forc_st_all_2015 =  firstYear(forc_st_all)
forc_npp_Oa_all_2015 =  firstYear(forc_npp_Oa_all)
forc_npp_Exu_all_2015 =  firstYear(forc_npp_Exu_all)
forc_sw_all_2015 =  firstYear(forc_sw_all)
forc_litt_all_2015 =  firstYear(forc_litt_all)

#单因子改变
#基准ssp都不变2015
if(TsNPPnm == "ssp2015"){
  forc_st_all =  forc_st_all_2015
  forc_npp_Oa_all =  forc_npp_Oa_all_2015
  forc_npp_Exu_all =  forc_npp_Exu_all_2015
  forc_sw_all =  forc_sw_all_2015
  forc_litt_all =  forc_litt_all_2015
  
}else if(TsNPPnm == "TsChang"){
  forc_npp_Oa_all =  forc_npp_Oa_all_2015
  forc_npp_Exu_all =  forc_npp_Exu_all_2015
  forc_sw_all =  forc_sw_all_2015
  forc_litt_all =  forc_litt_all_2015
  
}else if(TsNPPnm == "OaChang" | TsNPPnm == "OaChangDes0"){
  forc_st_all =  forc_st_all_2015
  forc_npp_Exu_all =  forc_npp_Exu_all_2015
  forc_sw_all =  forc_sw_all_2015
  forc_litt_all =  forc_litt_all_2015
  
}else if(TsNPPnm == "ExuChang"){
  forc_st_all =  forc_st_all_2015
  forc_npp_Oa_all =  forc_npp_Oa_all_2015
  forc_sw_all =  forc_sw_all_2015
  forc_litt_all =  forc_litt_all_2015
  
}else if(TsNPPnm == "OaExuChang"){
  forc_st_all =  forc_st_all_2015
  forc_sw_all =  forc_sw_all_2015
  forc_litt_all =  forc_litt_all_2015
}else if(TsNPPnm == "TsOaChangDes0"){
  forc_npp_Exu_all =  forc_npp_Exu_all_2015
  forc_sw_all =  forc_sw_all_2015
  forc_litt_all =  forc_litt_all_2015
  
}else if(TsNPPnm == "TsEqnoOaChang"){ 
  # 添加等量的其他根系分泌物 + 升温
  forc_npp_Oa_all =  forc_npp_Oa_all_2015
  forc_sw_all =  forc_sw_all_2015
  forc_litt_all =  forc_litt_all_2015
  
}else if(TsNPPnm == "TsOaChang"){
  # 有机酸增加 + 升温
  forc_npp_Exu_all =  forc_npp_Exu_all_2015
  forc_sw_all =  forc_sw_all_2015
  forc_litt_all =  forc_litt_all_2015
  
}else if(TsNPPnm == "EqnoOaChang"){
  # 添加等量的其他根系分泌物 + 升温
  forc_npp_Oa_all =  forc_npp_Oa_all_2015
  forc_sw_all =  forc_sw_all_2015
  forc_litt_all =  forc_litt_all_2015
  forc_st_all =  forc_st_all_2015

}


#基准ssp都变,控制单因子不变
if(TsNPPnm == "Ts"){
  forc_st_all =  forc_st_all_2015
}else if(TsNPPnm == "NPP_Oa"){
  forc_npp_Oa_all =  forc_npp_Oa_all_2015
}else if(TsNPPnm == "NPP_Exud"){
  forc_npp_Exu_all =  forc_npp_Exu_all_2015
}else if(TsNPPnm == "NPP_OaAExud"){
  forc_npp_Oa_all =  forc_npp_Oa_all_2015
  forc_npp_Exu_all =  forc_npp_Exu_all_2015
}

rm(forc_st_all_2015, forc_npp_Oa_all_2015, forc_npp_Exu_all_2015, 
   forc_sw_all_2015 ,forc_litt_all_2015)

stR_time = Sys.time()
cpools1 = data.frame()
if(scase=="ss"){
  outdf = data.frame()
}
outlist = list()

for(s in sites){
  stS_time = Sys.time()
  #  cat("------------ s = ", s, "   at ", stS_time)
  if(simmethod == "simi"){
    simid = soilPsim$sim[s]
  }else if(simmethod == "wrb10"){
    simid = soilPsim$wrb10[s]
  }else if(simmethod == "wrb9"){
    simid = soilPsim$wrb9[s]
  }
  
  Simsite = c("黄棕壤1","黑土","栗钙土","火山灰土","红壤","黄壤",
              "棕壤","草甸土","黄棕壤2","赤红壤")[simid]
  if(dim(forc_st_all)[2]%%365 ==0){
    forc_st0 <- forc_st_all[s,]
  }else{
    forc_st0 <- forc_repday(forc_st_all[s,])
  }
  if(dim(forc_sw_all)[2]%%365 ==0){
    forc_sw0 <- forc_sw_all[s,]
  }else{
    forc_sw0 <- forc_repday(forc_sw_all[s,])
  }
  if(dim(forc_litt_all)[2]%%365 ==0){
    forc_lit0 <- forc_litt_all[s,]
    forc_acid0 <- forc_npp_Oa_all[s,]
    forc_exud0 <- forc_npp_Exu_all[s,]
  }else{
    forc_lit0 <- forc_repday(forc_litt_all[s,])
    forc_acid0 <- forc_repday(as.numeric(forc_npp_Oa_all[s,]))
    forc_exud0 <- forc_repday(forc_npp_Exu_all[s,])
  }
  
  if(sstime%%length(forc_st0)!=0){
    stop(cat("forc_st0 length ERROR! forc_st0 = ",length(forc_st0),"sstime=",sstime,"\n"))
  }
  if(sstime%%length(forc_sw0)!=0){
    stop(cat("forc_sw0 length ERROR! forc_sw0 = ",length(forc_sw0),"sstime=",sstime,"\n"))
  }
  if(sstime%%length(forc_lit0)!=0){
    stop(cat("forc_lit0 length ERROR! forc_lit0 = ",length(forc_lit0),"sstime=",sstime,"\n"))
  }
  if(sstime%%length(forc_acid0)!=0){
    stop(cat("forc_npp0 length ERROR! forc_npp0 = ",length(forc_acid0),"sstime=",sstime,"\n"))
  }
  
  # if(soilPsim$IGBP[s] ==1){
  #   fnpp_oa = 0.2 #0.15 #0.1 #森林
  # }else if(soilPsim$IGBP[s] ==2){
  #   fnpp_oa = 0.2 #0.3 #灌木
  # }else if(soilPsim$IGBP[s] ==3){
  #   fnpp_oa = 0.2 #0.45 #草地
  # } 
  
  # forc_acid0 <-  forc_Oa0 * fnpp_oa *0.5 #有机酸比例
  # forc_exud0 <-  forc_Exud0 * fnpp_oa *0.5
  
  #如果n凋落物为0，结束此循环
  npar = 20
  cpools1p = NULL
  outss = NULL
  if(mean(forc_lit0) <= 1e-5 | is.na(mean(forc_lit0))){
    print("average litter input is 0")
    if(scase=="ss"){
      outss = rep(NA,(7*4+1)*npar)
      # outss = rep(NA,(8+7*3+1)*npar)
    }else{
      outss = as.data.frame(matrix(NA,nrow = nyear,ncol = length(derivnm)))
      cpools1p = outss[1,]
    }
    
    #next
  }else{
    #土壤属性
    fitpars$depth <- 0.2 #m
    fitpars$param_pH = soilPsim[s,"pH"]
    fitpars$param_bulkd <- soilPsim[s,"BD"]
    fitpars$param_clay <- soilPsim[s,"clay"]
    fitpars$param_claysilt <- soilPsim[s,"clay"] + soilPsim[s,"silt"]
    fitpars$porosity <- soilPsim[s,"porosity"]
    fitpars$param_em <- soilPsim[s,"param_em"]
    #    cat("soilPsim pH=",soilPsim[s,"pH"],"fitpars pH=",fitpars$param_pH,"clay=",fitpars$param_clay,"em=",fitpars$param_em,'\n')
    calpars = pars_resl[[Simsite]]$par
    
    for(pari in 1:npar){
      fitpars[colnames(calpars)]=calpars[pari,]
      
      if(scase=="ss"){
        soilCin = c(LIT=1,SOM=4, POM=1, LMWC=1, MIC=1, MAOM=1, MA = 0, MD = 0, Tco2=0)
        soilCin["MA"] = fitpars$r0 * soilCin["MIC"]
        soilCin["MD"] = (1-fitpars$r0) * soilCin["MIC"]
        
        soilCin13 = c(SOM=0, POM=0, LMWC=0, MIC=0, MAOM=0, MA = 0, MD = 0, Tco2=0)
      }else{
        objsite = which(as.numeric(rownames(outdf))==s) #site match
        objp = which(colnames(outdf) %in% paste0(c(fixCp,fixCp_c13,fixCp_cold, fixCp_cplmwc),".",pari))
        soilCin = outdf[objsite,objp[1:7]]
        # soilCin = outdf[objsite,objp[1:8]]
        names(soilCin) = fixCp
        soilCin["Tco2"] = 0
        soilCin["LIT"] = 0
        
        if(clab){
          soilCin13 = outdf[objsite,objp[8:14]]
          # soilCin13 = outdf[objsite,objp[9:15]]
          # names(soilCin13) = fixCp
          names(soilCin13) = fixCp
          soilCin13["Tco2"] = 0

          # label old carbon
          soilCinOld = outdf[objsite,objp[15:21]]
          # soilCinOld = outdf[objsite,objp[16:22]]
          names(soilCinOld) = fixCp
          soilCinOld["Tco2"] = 0

          # label plant carbon directly to lmwc
          soilCinPlantC2toLMWC = outdf[objsite,objp[22:28]]
          # soilCinPlantC2toLMWC = outdf[objsite,objp[23:29]]
          names(soilCinPlantC2toLMWC) = fixCp
          soilCinPlantC2toLMWC["Tco2"] = 0

        }else{
          soilCin13 = c(SOM=0, POM=0, LMWC=0, MIC=0, MAOM=0, MA = 0, MD = 0, Tco2=0)
          soilCinOld = c(SOM=0, POM=0, LMWC=0, MIC=0, MAOM=0, MA = 0, MD = 0, Tco2=0)
          soilCinPlantC2toLMWC = c(SOM=0, POM=0, LMWC=0, MIC=0, MAOM=0, MA = 0, MD = 0, Tco2=0)

        }
        
        if(scase=="hist"){
          ssgx = outdf[objsite, which(colnames(outdf) %in% paste0("ssgx.",pari))]
          ssgx = ifelse(is.na(ssgx),-9999,ssgx)
          if(ssgx!=1){
            soilCin["SOM"]==-9999
          }
        }
      }
      
      if(soilCin["SOM"]==-9999 | is.na(soilCin["SOM"])){
        outcppdf = as.data.frame(matrix(NA,nrow = nyear,ncol = length(derivnm)/npar))
      }else{
        OaChangDes0Cond = ifelse(TsNPPnm == "TsOaChangDes0", TRUE, FALSE)
        outcpp = derivs_V2_MM_AD_CO2(sstime,fitpars,as.list(soilCin),forc_st0,forc_sw0,forc_lit0,
                                     forc_acid0,forc_exud0,outflux,nloop,clab, as.list(soilCin13),OaChangDes0Cond,
                                     as.list(soilCinOld), as.list(soilCinPlantC2toLMWC)
                                    ) #return list,class(outcpp[[1]]) is numeric
        outcppdf = do.call(rbind,outcpp)
      }
      if(scase=="ss"){
        outss = c(outss, unlist(outcpp))
        #        if(unlist(outcpp)[length(outcpp[[1]])]==1){
        #           cat("ssgx=",unlist(outcpp)[length(outcpp[[1]])])
        #        }
      }else{
        
        if(is.null(outss)){
          outss = outcppdf
        }else{
          outss = cbind(outss, outcppdf)
        }
        #把最后一天的也输出保存
        cpools1p = c(cpools1p, outcppdf[nrow(outcppdf),])
        
      }
      
    }#pari ends
    
  }#lit > 0
  cat("---s=",s,"Total time: ", as.numeric(Sys.time()-stS_time, units = "mins"), "\n")
  
  if(scase=="ss"){
    names(outss) = derivnm
    outdf = rbind(outdf,outss)
  }else {
    colnames(outss) = derivnm
    outlist[[as.character(s)]] = outss
    names(cpools1p) = derivnm
    if(ncol(cpools1)!=0){
      colnames(cpools1) = derivnm
    }
  }
  #  print(colnames(cpools1))
  #  print(derivnm)
  #  print(names(cpools1p))
  cpools1 = rbind(cpools1, cpools1p)
  #  print("error")
  
}#s ends

print("foreach ending")
curtime = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
cat("site time:",curtime,"total time: ",as.numeric(Sys.time()-stR_time,units = "mins"),"\n")

#outss = do.call(outlist,rbind)
if(scase!="ss"){
  colnames(cpools1) = derivnm
  main_names <- gsub("\\.\\d+$", "", colnames(cpools1))
  outdf = cpools1[,which(main_names %in% c(fixCp,fixCp_c13,fixCp_cold, fixCp_cplmwc))]
}else{
  colnames(outdf) = paste0(rep(outname,npar),".",rep(1:npar,each=length(outname)))
}
rownames(outdf) = sites
