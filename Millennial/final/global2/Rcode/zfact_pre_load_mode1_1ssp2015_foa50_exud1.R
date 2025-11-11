#选择运行平台
run_platf = c("mine","hangdian","linux","linux_gzx")[3]
if(run_platf=="mine"){
  tempdir <- function() "D:/code/R/Rtemp"
}else if(run_platf=="hangdian"){
  tempdir <- function() "D:/liucq/Rtemp"
}else if(run_platf=="linux"){
  tempdir <- function() "/datanode05/liucq/python/Rtmp"
}else{
  tempdir <- function() "/datanode05/guozx/liucq/Rtmp"
}

unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir",tempdir, baseenv())
lockBinding("tempdir",baseenv())
print(tempdir())
print(tempfile())

#####
for(i in c(3,5)){#ssp为3:5
  scase = c("ss","hist","ssp126","ssp245","ssp585")[i]
  modi = c(1,2,3)[1]
  forc_mod = c("CMCC-ESM2","NorESM2-MM","CESM2-WACCM")[modi]
  
  baseline = c("sspChange", "ssp2015")[2]
  if(baseline == "sspChange"){
    TsNPPnm = c("","Ts","NPP_Oa","NPP_Exud","NPP_OaAExud")[4]
  }else if(baseline == "ssp2015"){
    TsNPPnm = c(baseline,"TsChang","OaChang","TsOaChangDes0", # 4
                "TsEqnoOaChang", "TsOaChang")[1] #跑前4个就行 # , "ExuChang","OaExuChang"
  }
  
  fexu_oa = c(25, 50,75)[2] #有机酸比例 % 
  modelnm = paste0(modi,"_",forc_mod)
  outfnm = paste0(modelnm,"_",scase)
  chagnm = c("","_Exudatef15","_Exudatef20","_sampling_8_nloop500_ssloop2")[1]
  fexudmod = c("mean", "BCC-CSM2", "NorEMS2")[2]
  
  npar = 20
  interval = 1#8
#  sites = seq(1,7905,interval) #抽样选一部分点
  outflux = ifelse(scase=="ss",FALSE,TRUE)
  clab = c(FALSE,TRUE)[2]
  fixCp = c("SOM","POM","LMWC","MIC","MAOM","MA","MD")
  
  if(outflux){
    allflux.name = c("f_Lit_PO", "f_Lit_LM", "f_Exud_LM", "f_Oa_MA", "f_Oa_LM",
                     "f_PO_LM", "f_PO_MA", "f_MA_LMm", "f_MA_LMe", "acid_des","f_MA_LM","f_LM_MA", 
                     "f_MA_MD", "f_MD_MA", "f_LM_MB",  "f_MA_maint_co2", "f_MA_growth_co2",
                     "f_MA_co2", "f_MD_co2", "f_MA_growth",
                     "f_MB_LM", "f_MB_MA", "f_MB_turn","f_LM_leach")
  }else{
    allflux.name = NULL
  }
  if(clab){
    fixCp_c13 = paste0(fixCp,"_c13")
  }else{
    fixCp_c13 = NULL
  }
  cat("outfnm=",outfnm,"interval=",interval,"outflux=",outflux,"clab=",clab,"\n")
  
  #usgl = "/datanode05/guozx/liucq/Millennial/final"
  usgl = "/datanode05/liucq/python/Millennial/final"
  outdir0 =  file.path(usgl,"global2/output")
  indirgl = file.path(usgl,"global2/input")
  Rgl = file.path(usgl,"global2/Rcode")
  indir1 = file.path(usgl,"global2/input")
  
  load(file.path(indir1,"soilPsim2.Rdata"))
  cat("soilP2 path: ",file.path(indir1,"soilPsim2.Rdata"),'\n')
  
  library(Rcpp)
  sourceCpp(file.path(Rgl,"Rmillennial2.cpp"))
  
  #3.参数
  load(file.path(indirgl,"parsResl_cpp.RData"))
  pars_in_raw <- read.table(file.path(indirgl,"soilpara_in_fit1_clab.txt"))
  pars_in <- as.list(pars_in_raw$V2)
  names(pars_in) <- pars_in_raw$V1
  pars <- unlist(pars_in)
  rm(pars_in,pars_in_raw)
  
  fitpars = pars
  source(file.path(Rgl,"ss_trans_processfun.R"))
  print("preprocess finished!")
  
  #不同情景模拟
  grid = c(0,500,2000,3000,4000,5000,6000,nrow(soilPsim))#[1:4]##[4:8]
  #grid = c(0,2000,4000,6000,nrow(soilPsim))[2:3]
  sites = seq(1,nrow(soilPsim),interval)[(min(grid)+1):max(grid)]
  sitenml = sprintf("_%04d_%04d",min(sites),max(sites))

#不同情景下模拟输出名字, 2025.08.01
  source(file.path(Rgl,"testoutnm.R")) 
  # chagnm = testoutnm
  sink(file.path(Rgl,"outTXT",paste0(outfnm,chagnm,"_",TsNPPnm,sitenml,"_",fexu_oa,"_",fexudmod,"_output.txt")))

#  sink(file.path(Rgl,"outTXT",paste0(outfnm,"_",TsNPPnm,sitenml,"_fexu_oa-",fexu_oa,"_",fexudmod,"_output.txt")))
  cat("outfnm=",outfnm,"interval=",interval,"outflux=",outflux,"clab=",clab,"\n")
  #chagnm = c("","_sampling_8_nloop500_ssloop2")[2]
  
  if(scase=="ss"){
    load(file.path(indir1,paste0(modi,"_forc_all_",forc_mod,"_hist",".Rdata")))
    outdir =  file.path(outdir0,"ss")
    nloop = 500 # 0->no steady-state; >500->force data is constant; others->global 
    nyear = 25
    sstime = 365* nyear * nloop #*25*50
    
    outname = c(fixCp[1:7],paste0(fixCp[1:7],"_c13"),"ssgx")
    derivnm = paste0(rep(outname,npar),".",rep(1:npar,each=length(outname)))
    
    cat("rows of soilPsim:",nrow(soilPsim),"forc_st=",nrow(forc_st_all),"sitenml=",sitenml,'\n')
    source(file.path(Rgl,"derivs_run.R")) #(1).solve steady state
    save(outdf,file = file.path(outdir,paste0(outfnm,chagnm,sitenml,".Rdata")))
    print("run is finished!")
    
  }else{
#    if(TsNPPnm!=""){
#      load(file.path(indir1,paste0(modi,"_forc_all_",forc_mod,"_hist",".Rdata")))
#      forc_npp_all_h = forc_npp_all
#      forc_st_all_h = forc_st_all
#      rm(forc_st_all,forc_sw_all,forc_litt_all,forc_npp_all)
#    }

    load(file.path(indir1,paste0(modi,"_forc_all_",forc_mod,"_",scase,".Rdata")))

        #增加更新后分泌物数据 2025.7.15
    rm(forc_npp_all)
    if(fexudmod == "mean"){
      load(file.path(usgl,"global2/input/forc_exud_allCase.RData")) #变量名为：forc_exud_hist, forc_exud_ssp126, forc_exud_ssp585
      forc_npp_all = as.matrix(switch(scase, "hist" = forc_exud_hist,
                                           "ssp126" = forc_exud_ssp126,
                                           "ssp585" = forc_exud_ssp585)
      )
      rm(forc_exud_hist,forc_exud_ssp126,forc_exud_ssp585)

    }else {
      load(file.path(usgl,"global2/input/forc_exud_allCase_forc_exud_model2.RData")) 
      forc_npp_all = as.matrix(switch(scase, "hist" = forc_exud_model2[[fexudmod]]$forc_exud_hist,
                                           "ssp126" = forc_exud_model2[[fexudmod]]$forc_exud_ssp126,
                                           "ssp585" = forc_exud_model2[[fexudmod]]$forc_exud_ssp585))
      rm(forc_exud_model2)                                     
    }

    nloop = 0
    outname = c("year",fixCp,"Tco2",allflux.name,fixCp_c13,"Tco2_c13")
    derivnm = paste0(rep(outname,npar),".",rep(1:npar,each=length(outname)))
    outdir =  file.path(outdir0,paste0("fact",scase))
    
    if(scase=="hist"){
      nyear = 25
      inout = file.path(outdir0,"ss")
      load(file.path(inout,paste0(modelnm,"_ss",chagnm,sitenml,".Rdata"))) #initial c pools
#      load(file.path(inout,paste0(modelnm,"_ss",chagnm,"_0001_7527",".Rdata"))) #test
    }else{
      nyear = 86
      inout = file.path(outdir0,"hist")  
      load(file.path(inout,paste0(modelnm,"_hist",chagnm,"_",sitenml,"_fexu_oa-",fexu_oa,"_",fexudmod,"_outdf.Rdata"))) #initial c pools   
 #     load(file.path(inout,paste0(modelnm,"_hist",chagnm,"__0001_7527","_fexu_oa-",fexu_oa,"_outdf.Rdata"))) #test
    }
    incpools = outdf
    sstime = 365* nyear
    cat("rows of soilPsim:",nrow(soilPsim),"forc_st=",nrow(forc_st_all),'\n')
    source(file.path(Rgl,"derivs_run.R"))
    save(outdf,file = file.path(outdir,paste0(outfnm,chagnm,"_",TsNPPnm,sitenml,"_fexu_oa-",fexu_oa,"_",fexudmod,"_outdf.Rdata")))
    save(outlist,file = file.path(outdir,paste0(outfnm,chagnm,"_",TsNPPnm,sitenml,"_fexu_oa-",fexu_oa,"_",fexudmod,"_outlist.Rdata"))) #每个点，每年20组参数结果df
    
  }
  print("pre_load finished!")
  sink()
  rm(list=ls())
}









