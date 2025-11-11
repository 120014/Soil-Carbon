#多目标参数率定-clab

#---------
#率定*****只能单个点率定，如果多个点的实验和模拟数据顺序不对应*******
#testoutnm = c("_Vdown10","_rmPOCtoMAOC","_rmPOCtoMAOC_Vdown10_MAOCmm","_rmPOCtoMAOC_Vdown10_MAOCmm_limitfoa60", "_Vdown10_MAOCmm_limitfoa60",#5
#              "_rmPOCtoMAOC_MAOCmm_Vp5","_rmPOCtoMAOC_MAOCmm_Vp5_fixkm","_rmPOCtoMAOC_MAOCmm_Vp1_fixkm","_rmPOCtoMAOC_MAOCmm_Vpdown31",#9
#              "_rmPOCtoMAOC_MAOCmm_Vpdown31_sw",#10
#              "_rmPOCtoMAOC_MAOCmm_Vpdown31_sw_MAdtoPOC",#11
#              "_rmPOCtoMAOC_MAOCmm_Vpdown31_wpMAD_sw", #12
#              "_rmPOCtoMAOC_MAOCmm_Vpdown31_wpMAD_MAdtoPOC_sw" #13
#             )[13]

testoutnm = c("_rmPOCtoMAOC_MAOCmm_Eap66_sw05",#1
              "_rmPOCtoMAOC_MAOCmm_Eap66_sw05_vdown31",  #2
              "_rmPOCtoMAOC_MAOCmm_Eap66_sw05_vdown31_madtopoc",   #3
              "_rmPOCtoMAOC_MAOCmm_Eap66_sw05_madtopoc",
              "_rmPOCtoMAOC_rmMAOCmm_Eap66_sw05_vdown31_madtopoc"    #5
             )[5]
do_mod_fit=1
pars_type = c("alpha","Vmax_tref","alpha_vmd")[2]
naga2condim=ifelse(pars_type=="Vmax_tref",2,3)
c13werror = FALSE
qmax_cue_error = FALSE

ss_cond_obj = c(TRUE,FALSE)[1] #率定参数是否要增加稳态判断惩罚项
print(ss_cond_obj)
if(ss_cond_obj){
  SST = FALSE
  ss.num = 0
  # sstype = c("runsteady","ND")[2]
  nloop = 25*50#20
  if(sstype=="runsteady"){
    library(rootSolve)
    # POM =1 
    # LMWC = 1
    # MAOM =1
    # MIC =1
    # MA=0.1 # par
    # MD=0.9
    qmax_cue_error_ss = FALSE
  }else if(sstype=="ND"){
    rm(POM,LMWC,MAOM,MIC,MA,MD)
  }else if(sstype=="ND_cpp"){
    library(Rcpp)
#    sourceCpp(file.path(indirRexp,"Rmillennial2.cpp"))

    if(testMd == "_VMA_VMD"){
      sourceCpp(file.path(indirRexp,"Rmillennial2_VMA_VMD.cpp"))
    }else{
      sourceCpp(file.path(indirRexp,"Rmillennial2.cpp"))
    }

  }
}else{
  sstype=NULL
}

obj_factor = 100 #目标函数缩放因子
#nobj = 3 #ifelse(testMd != "_objCO2R2MARE",5,3)
if(testMd == "_objCO2R2MARE" | testMD2 | testMd == "_fixbeta_Obj"){
  nobj = 3
}else{
  nobj = 5
}
nobj = ifelse(clab[1]!=FALSE, nobj+1, nobj)
nobj = ifelse(ss_cond_obj, nobj+1, nobj)
cat("----c13lab=",clab,"SS=",ss_cond_obj,"  nobj=",nobj,"  sstype=",sstype," pars_type=",pars_type,'\n')
source(file.path(indirRexp,"pars_cal_choice.R"))

if(do_mod_fit==1){
  
  #目标函数
  
  Objective.run <- function(x, parset = names(x)) {
    
    sink(file.path(outdir,paste0(si,"_",Simsite,"_",ngen,"_",npop,"_vmax_vgcue_cueacid_time",sstype,testMd,"_",testMD2,testoutnm,"_litter",".txt"))) #输出到文件中的信息同时打印在console???, split=TRUE,append=TRUE)
    print(Simsite)
    itensga2 <<- itensga2 +1
    print(paste("迭代次数: ",itensga2))
    if(ss_cond_obj){
      ss.num <<- ifelse(SST, ss.num+1, ss.num)
      cat("the number of steady is ",ss.num,'\n')
    }
    # print(paste("x:",x))
    
    # cat("1.names(x): ",names(x),'\n')
    names(x) = names(parToFit.lower)
    # cat("2.names(x): ",names(x),'\n')
    # cat("parset: ",parset,'\n')
    pars[parset] <- x
    # cat(" pars[parset]: ", pars[parset],'\n')
    # cat(" pars: ", pars,'\n')
    # print(pars[parset])
    
    #fitpars = unlist(pars)
    #pars = fitpars
    
    #达到稳态后求解实验部分
    c13werror <<- FALSE
    qmax_cue_error <<- FALSE
    # cat("out前qmax_cue_error=",qmax_cue_error,'\n')
    out <- Model.pools(pars,inputdata[which(inputdata$site %in% Simsite),],expT=expT) 
    # cat("out后=",qmax_cue_error,'\n')
    #print(out)
    #模拟单个点还是多个点，并计算目标
    #方法1
    Obs.co2_day0 = Obs.co2[which(Obs.co2$site %in% Simsite),"CO2"] #速率
    Obs.pools0 =  Obs.pools[which(Obs.pools$site %in% Simsite),]
    Obs.co2_C13_day = Obs.co2_C13[which(Obs.co2_C13$site %in% Simsite),"C13_CO2"] #C13_CO2速率
    Obs.co2_tt0 = Obs.co2_tt[which(Obs.co2_tt$site %in% Simsite),"CO2"] #累积
    Obs.co2_C13_tt0 = Obs.co2_C13_tt[which(Obs.co2_C13_tt$site %in% Simsite),"C13_CO2"]#C13_CO2累计
    
    
    #方法1 if(Simsite=="黄棕壤1" | Simsite=="黑土" | Simsite=="栗钙土" | Simsite=="火山灰土")
    if(si<=4){
      # out.CO2_day = out[vector_gap(1+CO2_30_r,simd), "CO2_day"]
      # out.CO2_tt = out[vector_gap(1+CO2_30,simd), "CO2"] 
      # out.ML =  out[vector_gap(1+ML30,simd), c("LMWC", "MIC")]
      CO2t = CO2_30
      CO2_rt = c(0,CO2t)
      MLt = ML30
      c13CO2t = c13CO2_30
      c13CO2_rt = c(0,c13CO2t)
    }else{
      # out.CO2_day = out[vector_gap(1+CO2_10_r,simd), "CO2_day"]
      # out.CO2_tt = out[vector_gap(1+CO2_10,simd), "CO2"]
      # out.ML =  out[vector_gap(1+ML10,simd), c("LMWC", "MIC")]
      CO2t = CO2_10
      CO2_rt = c(0,CO2t)
      MLt = ML10
      c13CO2t = c13CO2_10
      c13CO2_rt = c(0,c13CO2t)
      
    }
    
    out[,c("CO2_day")] = calculate_diff(out[,"Tco2",drop=FALSE]) #多个温度也直接减？因为保留了初始值可以的
    
    out.CO2_day = out[which(out$time %in% (1+CO2_rt)), "CO2_day"] 
    out.CO2_tt = out[which(out$time %in% CO2t), "Tco2"] 
    out.ML =  out[which(out$time %in% MLt), c("LMWC", "MIC")]
    colnames(out.ML) = paste0("sim",colnames(out.ML))
    #只取添加有机酸的元素
    
    if(clab[1]!=FALSE){
      out[,c("C13_CO2_day")] = calculate_diff(out[,"Tco2_c13_Cinoa",drop=FALSE])
      c13row = which(out$time %in% (1+c13CO2_rt))
      out.CO2_c13_day = out[c13row[-(1:length(c13row) / 2)], "C13_CO2_day"] 
      c13row = which(out$time %in% c13CO2t)
      out.CO2_c13_tt = out[c13row[-(1:length(c13row) / 2)], "Tco2_c13_Cinoa"]
    }
    
    #目标函数权重
    if(testMd == "_objCO2R2MARE" | testMD2){
      objw = c(2, 1,1)
      names(objw) = c("R.CO2","M.LMWC","M.MIC")
    }else{
      objw = c(2, 1, 2 ,1,1)
      names(objw) = c("R.CO2","R.CO2_day","M.CO2","M.LMWC","M.MIC")
    }
    
    # if(Simsite=="黄棕壤2"){
    #   objw = c(2, 1, 4 ,1,1)
    # }
    
    if(clab[1]!=FALSE){
      objw = c(objw, M.c13 = 2)
    }
    if(ss_cond_obj){
      objw = c(objw, ss.ARE = 1) #R2_CO2, R2_CO2_day, MARE_CO2, LMWC,MIC,SSARE,没有增加C13_CO2
    }
    objw = objw/sum(objw)
    
    #CO2的目标函数为1-R2
    MARE = c()
    R2.CO2 = c()
    
    Wqmax_cue = ifelse(qmax_cue_error,100,1)
#    if(testMd == "_objCO2R2MARE"){ #后面看是否能直接改成公式计算R2
      tempR2 = compute_R2(Obs.co2_tt0, out.CO2_tt)
#    }else{
#      tempR2 = summary(lm(Obs.co2_tt0~out.CO2_tt))$r.squared
#    } 
#    tempR2 = compute_R2(Obs.co2_tt0, out.CO2_tt)
#    tempR2 = summary(lm(Obs.co2_tt0~out.CO2_tt))$r.squared
#    print(diff(out[c(2,simd)+simd*3, "CO2_day"]))
    Condoutco2r = diff(out[c(168,nrow(out)), "CO2_day"]) >= 0 
    wlimtR = ifelse(Condoutco2r, 90, 1) #限制CO2速率一直增加和

    out30obs10r = out[nrow(out), "Tco2"] / Obs.co2_tt0[length(Obs.co2_tt0)]
    if(out30obs10r > 4){
      Condoutco2Cul = 10
    }else if(out30obs10r > 3){
      Condoutco2Cul = 5
    }else if(out30obs10r > 2.5){
      Condoutco2Cul = 2
    }else if(out30obs10r > 2.2){
      Condoutco2Cul = 1.5
    }else{
      Condoutco2Cul = 1
    }   

    R2.CO2["CO2"] = (1-tempR2)*objw["R.CO2"] *Wqmax_cue * wlimtR * Condoutco2Cul
    
    if(testMd != "_objCO2R2MARE" &  testMD2==FALSE){
      R2.CO2["CO2_day"] = (1-summary(lm(Obs.co2_day0~out.CO2_day))$r.squared)*objw["R.CO2_day"]*Wqmax_cue
      #R2.CO2 = 1-summary(lm(Obs.co20~out.CO2))$r.squared
      
      #DOC和MBC的目标函数MARE
      MARE["CO2"] = mean(abs(out.CO2_tt-Obs.co2_tt0)/Obs.co2_tt0)*objw["M.CO2"]*Wqmax_cue
    }

    # MARE["CO2_day"] = mean(abs(out.CO2_day-Obs.co2_day0)/Obs.co2_day0)
    MARE["LMWC"] = mean(abs(out.ML[,"simLMWC"]-Obs.pools0[,"LMWC"])/Obs.pools0[,"LMWC"])*objw["M.LMWC"]*Wqmax_cue
    MARE["MIC"] = mean(abs(out.ML[,"simMIC"]-Obs.pools0[,"MIC"])/Obs.pools0[,"MIC"])*objw["M.MIC"]*Wqmax_cue
    if(clab[1]!=FALSE){
      c13werror = ifelse(c13_error,100,1)
      # MARE["c13CO2_day"] = mean(abs(out.CO2_c13_day-Obs.co2_C13_day)/Obs.co2_C13_day)*objw["M.c13"]
      MARE["c13CO2"] = mean(abs(out.CO2_c13_tt-Obs.co2_C13_tt0)/Obs.co2_C13_tt0)*objw["M.c13"] * c13werror
    }
    
    # cat("out.c13：",out.CO2_c13_day,'\n')
    # cat("obs.c13：",Obs.co2_C13_day,'\n')
    # cat("objw[6]=",objw[6],"mean=",mean(abs(out.CO2_c13_day-Obs.co2_C13_day)/Obs.co2_C13_day),'\n')
    #能否达到稳态解的惩罚项,如果实验的Qmax和CUE已经错误，那就不再求解稳态
    if(ss_cond_obj & qmax_cue_error==FALSE){
      cat("begin to run steady*********************************** ",sstype,'\n')
      # print(pars)
      outps = calpars_SSfun(sstype = sstype,fitpars=pars) #step,runsteady,返回data.frame
#      print(paste("outps = ",outps))
      if(sstype=="ND_cpp"){
        SST <<- ifelse(outps["ssgx"]==1,TRUE,FALSE)
        ssgx = ifelse(SST, 1, 100)#稳态惩罚项
      }else{
        ssgx = 1
      }
      ssc = which(names(outps) %in% fixCp[2:6]) #3:5
#      Wfmaoc = ifelse((outps["MAOM"]/outps["SOM"])>0.8, 10, 1) #乘权重比例
      fmaocsocSS = outps["MAOM"]/outps["SOM"] * 100 #稳态只用MAOC和MAOC占SOC比例, #乘权重比例
      if(fmaocsocSS>=90){
        Wfmaoc = 100
      }else if(fmaocsocSS>=80){
        Wfmaoc = 5
      }else if(fmaocsocSS>=70){
        Wfmaoc = 2
      }else{
        Wfmaoc = 1
      }
      
      
      #增加一个稳态目标函数POM、LMWC、MIC和MAOM，一个???
      Obsp = (inputdata[which(inputdata$site %in% Simsite)[1],names(outps)[ssc],drop=FALSE])
      meanARE = mean(unlist(abs(outps[ssc]-Obsp[1,])/Obsp[1,])) *objw["ss.ARE"] *Wqmax_cue * Wfmaoc#* ssgx
    } else{
      ssgx = 100
      meanARE = objw["ss.ARE"] *Wqmax_cue #* ssgx
    }
    ssgx = ifelse(ss_cond_obj==FALSE,1, ssgx)
    
    print("obj结束")
    
    sink()
    
    # 在函数内部修改fitpars_iter全局变量
    fitpars_iter$par <<- rbind(fitpars_iter$par, x)
    rownames(fitpars_iter$par) <<- c(1:itensga2)
    # fitpars_iter$value <<- rbind(fitpars_iter$value, c(R2.CO2,MARE,meanARE))
    
    out.obj = c(R2.CO2,MARE)
    out.value = out.obj
    if(ss_cond_obj){
      out.obj = c(out.obj, meanARE)
      names(out.obj) = names(objw)
      out.value = c(out.obj, ssgx = ssgx)
    }
    out.obj = out.obj * ssgx
    out.value[1:nobj] = out.value[1:nobj] * ssgx
    
    fitpars_iter$value <<- rbind(fitpars_iter$value, out.value)
    rownames(fitpars_iter$value) <<- c(1:itensga2)
    # cat("obj的长度:",length(c(R2.CO2,MARE)),'\n')
    #提前输出200*160的结果
    if(itensga2==32160){# | itensga2==16160 | itensga2==24160 
      Fit.pools_iter[[Simsite]] <- fitpars_iter #所有迭代结果
      save(Fit.pools_iter,
           file = file.path(outdir,"V6Rdata/200_160",paste0(si,"_",Simsite,"_",(itensga2-160)/160,"_vmax_vgcue_cueacid_time",sstype,testMd,"_",testMD2,"_litter",".RData")))
      print(paste(Simsite,(itensga2-160)/160,"*160",sstype,"_vmax_vgcue_cueacid_time已输出"))
    }
    return(out.obj*obj_factor)
    
  }
  
  compute_R2 <- function(obs, pred) {
    SS_res <- sum((obs - pred)^2, na.rm = T)
    SS_tot <- sum((obs - mean(obs, na.rm = T))^2, na.rm = T)
    R2 <- 1 - (SS_res / SS_tot)
    return(round(R2,4))
  }
}

