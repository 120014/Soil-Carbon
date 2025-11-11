#多个模型稳态模拟
# sts =  1
# ends = 1000
# 动态驱动数据当前时期多年平均 #hist multi-year average
#forc_npp_all[forc_npp_all<0] = 0
# forc_npp_all_h[forc_npp_all_h<0] = 0
# forc_npp_allMy = matrix(mutliY_mean(forc_npp_all_h),ncol = ncol(forc_npp_all),nrow=nrow(forc_npp_all))
if(clab){
  fixCp_cold = paste0(fixCp,"_oldc")
  fixCp_cplmwc = paste0(fixCp,"_plantc2lmwc")
  # fixCp_cold = paste0(fixCp[2:8],"_oldc")
  # fixCp_cplmwc = paste0(fixCp[2:8],"_plantc2lmwc")
}else{
  fixCp_cold = NULL
  fixCp_cplmwc = NULL
}
outname = c("year",fixCp,"Tco2",allflux.name,fixCp_c13,"Tco2_c13", fixCp_cold,"Toldco2", fixCp_cplmwc,"Tpi2LMWCco2")
# outname = c("year", fixCp, "Tco2", allflux.name, 
#             paste0(fixCp[2:8],"_c13"), "Tco2_c13",  # 跳过LIT
#             paste0(fixCp[2:8],"_oldc"), "Toldco2",  # 跳过LIT
#             paste0(fixCp[2:8],"_plantc2lmwc"), "Tpi2LMWCco2")  # 跳过LIT
derivnm = paste0(rep(outname,npar),".",rep(1:npar,each=length(outname)))

cat("derivnm 长度:", length(derivnm), "\n")




if(ncol(forc_npp_all) %% 12 == 0){
  cols_per_year <- 12
}else if(ncol(forc_npp_all) %% 365 == 0) {
  cols_per_year <- 365
}
forc_npp_allMy = forc_npp_all
forc_npp_allMy[] = forc_npp_allMy[,1:cols_per_year]

if(TsNPPnm == "Ts"){
  if(ncol(forc_st_all) %% 12 == 0){
    cols_per_year <- 12
  }else if(ncol(forc_st_all) %% 365 == 0) {
    cols_per_year <- 365
  }
  # forc_st_all[] =  mutliY_mean(forc_st_all)
  forc_st_all[] =  forc_st_all[,1:cols_per_year]
}

stR_time = Sys.time()
#sink(file.path(Rgl,paste0(modi,"_",forc_mod,"_",nloop,"_",scase,"_output.txt")))
cpools1 = data.frame()
if(scase=="ss"){
  outdf = data.frame()
}
outlist = list()
#out <- foreach(s = sites, .export = pvars, .packages("Rcpp")) %dopar%  {
for(s in sites){
  stS_time = Sys.time()
  #  cat("------------ s = ", s, "   at ", stS_time)
  
  Simsite = c("黄棕壤1","黑土","栗钙土","火山灰土","红壤","黄壤",
              "棕壤","草甸土","黄棕壤2","赤红壤")[soilPsim$sim[s]]
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
    forc_npp0 <- forc_npp_all[s,]
    forc_nppMy0 <- forc_npp_allMy[s,]
  }else{
    forc_lit0 <- forc_repday(forc_litt_all[s,])
    forc_npp0 <- forc_repday(as.numeric(forc_npp_all[s,]))
    forc_nppMy0 <- forc_repday(forc_npp_allMy[s,])
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
  if(sstime%%length(forc_npp0)!=0){
    stop(cat("forc_npp0 length ERROR! forc_npp0 = ",length(forc_npp0),"sstime=",sstime,"\n"))
  }
  
  if(soilPsim$IGBP[s] ==1){
    fnpp_oa = 1 #0.2 #0.15 #0.1 #森林
  }else if(soilPsim$IGBP[s] ==2){
    fnpp_oa = 1 #0.2 #0.3 #灌木
  }else if(soilPsim$IGBP[s] ==3){
    fnpp_oa = 1 #0.2 #0.45 #草地
  } 

  forc_acid0 <-  forc_npp0 * fnpp_oa *fexu_oa/100 #有机酸比例
  forc_exud0 <-  forc_npp0 - forc_acid0
  
#  if(TsNPPnm == "NPP_Oa"){
#    forc_acid0 =  forc_nppMy0 * fnpp_oa *0.5 #有机酸比例
#  }else if(TsNPPnm == "NPP_Exud"){
#    forc_exud0 =  forc_nppMy0 * fnpp_oa *0.5 #有机酸比例
#  }else if(TsNPPnm == "NPP_OaAExud"){
#    forc_acid0 =  forc_nppMy0 * fnpp_oa *0.5 #有机酸比例
#    forc_exud0 =  forc_nppMy0 * fnpp_oa *0.5 #有机酸比例
#  }
  
  #如果n凋落物为0，结束此循环
  npar = 20
  cpools1p = NULL
  outss = NULL
  if(mean(forc_lit0) <= 1e-5 | is.na(mean(forc_lit0))){
    print("average litter input is 0")
    if(scase=="ss"){
      outss = rep(NA,(7*4+1)*npar)
      # outss = rep(NA,(8+7*4+1)*npar)
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
        # objp = which(colnames(outdf) %in% paste0(c(fixCp,fixCp_c13,fixCp_cold, fixCp_cplmwc),".",pari))
        objp = which(colnames(outdf) %in% paste0(c(fixCp,fixCp_c13,fixCp_cold, fixCp_cplmwc),".",pari))
        soilCin = outdf[objsite,objp[1:7]]
        # soilCin = outdf[objsite,objp[1:8]]
        names(soilCin) = fixCp
        soilCin["Tco2"] = 0
        soilCin["LIT"] = 0
        
        if(clab){
          soilCin13 = outdf[objsite,objp[8:14]]
          # soilCin13 = outdf[objsite,objp[9:15]]
          names(soilCin13) = fixCp
          # names(soilCin13) = fixCp[2:8]
          soilCin13["Tco2"] = 0

          # label old carbon
          # soilCinOld = outdf[objsite,objp[15:21]]
          soilCinOld = outdf[objsite,objp[16:22]]
          # names(soilCinOld) = fixCp
          names(soilCinOld) = fixCp[2:8]
          soilCinOld["Tco2"] = 0
          # label plant carbon directly to lmwc
          soilCinPlantC2toLMWC = outdf[objsite,objp[22:28]]
          # soilCinPlantC2toLMWC = outdf[objsite,objp[23:29]]
          names(soilCinPlantC2toLMWC) = fixCp
          # names(soilCinPlantC2toLMWC) = fixCp[2:8]
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
        outcpp = derivs_V2_MM_AD_CO2(sstime,fitpars,as.list(soilCin),forc_st0,forc_sw0,forc_lit0,
                                     forc_acid0,forc_exud0,outflux,nloop,clab, as.list(soilCin13), FALSE,
                                     as.list(soilCinOld), as.list(soilCinPlantC2toLMWC)
                                     ) #return list,class(outcpp[[1]]) is numeric
        
        cat("C++ output length:", length(outcpp), "expected:", length(derivnm)/npar, "\n")
        cat("outss 列数:", ncol(outss), "\n")
        cat("cpools1 列数:", ncol(cpools1), "\n")
        
        
        
        
        
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
# save(outss,file = file.path(outdir,paste0(modi,"_",forc_mod,sprintf("_ss_%04d_%04d.Rdata",sts,ends))))
# print("outdf save finished")

# 提前计算需要分析的结果？算了
