#多个模型稳态模拟
# sts =  1
# ends = 1000
#testoutnm = c("","rmPOCtoMAOC_MAOCmm_Vp5","rmPOCtoMAOC_MAOCmm_Vpdown","rmPOCtoMAOC_MAOCmm_Vpdown_noQmaxbreak",#4
#              "rmPOCtoMAOC_MAOCmm_Vpdown31",#5
#              "rmPOCtoMAOC_MAOCmm_Vpdown31_sw", #6 
#              "limitfoa80_rmMAOCmm", #7
#              "rmPOCtoMAOC_MAOCmm_Vpdown31_sw_MAdtoPOC", #8
#              "rmPOCtoMAOC_MAOCmm_Vpdown31_wpMAD_forcmean_error", #9
#              "rmPOCtoMAOC_MAOCmm_Vpdown31_wpMAD_sw_forcmean", #10
#              "rmPOCtoMAOC_MAOCmm_Vpdown31_wpMAD_sw", #11
#              "rmPOCtoMAOC_MAOCmm_Vpdown31_wpMAD_MAdtoPOC_sw",#12
#              "rmPOCtoMAOC_wrb10", #13
#              "rmPOCtoMAOC_wrb9", #14
#              "rmPOCtoMAOC_MAOCmm_Vpdown31_wrb10",#15
#              "rmPOCtoMAOC_MAOCmm_Vpdown31_wrb9",#16
#              "rmPOCtoMAOC_MAOCmm_Vpdown31_wrb10_loopth100"#17
#             )[17]
simmethod = c("simi", "wrb10","wrb9", "meanmin10", "meanmin9")[1]

testoutnm = c("rmPOCtoMAOC_MAOCmm_Eap66_sw05",
              "rmPOCtoMAOC_MAOCmm_Eap66_sw05_vdown31",
              "rmPOCtoMAOC_MAOCmm_Eap66_sw05_madtopoc",#不能用
              "rmPOCtoMAOC_MAOCmm_Eap66_sw05_vdown31_madtopoc", #4
              "rmPOCtoMAOC_MAOCmm_Eap66_sw05_vdown31_madtopoc_noQmaxReturn",  #5 
              "rmPOCtoMAOC_rmMAOCmm_Eap66_sw05_vdown31_madtopoc", #6
              "rmPOCtoMAOC_MAOCmm_Eap66_sw05_vdown31_madtopoc_oldc_noQmaxReturn", #7 
              "rmPOCtoMAOC_rmMAOCmm_Eap66_sw05_vdown31_madtopoc_oldc_noQmaxReturn", #8 
              "rmPOCtoMAOC_MAOCmm_Eap66_sw05_vdown31_madtopoc_oldc_plantC2LMWC_noQmaxReturn", #9
              "rmPOCtoMAOC_MAOCmm_Eap66_sw05_vdown31_madtopoc_oldc_plantC2LMWC_noQmaxReturn_modify", #10
              "rmPOCtoMAOC_MAOCmm_Eap66_sw05_vdown31_madtopoc_oldc_plantC2LMWC_noQmaxReturn_litter" #10
              )[11]
testoutnm = paste0(testoutnm,"_", simmethod)

interval = 1
#sites = seq(1,200,interval) #seq(1,nrow(soilPsim),interval) #抽样选一部分点
sites = seq(1,nrow(soilPsim),interval)
outflux = c(FALSE,TRUE)[1]
fixCp = c("SOM","POM","LMWC","MIC","MAOM","MA","MD","Tco2")

#并行计算
library(foreach)
library(doParallel)

stR_time = Sys.time()
pvars = c("forc_repday","derivs_V2_MM_AD_CO2","pars_resl","fitpars","soilPsim","fixCp","outflux",
          "forc_st_all","forc_sw_all","forc_litt_all","forc_npp_all")
# nloop = 500
# nyear = 25
nloop = 500
nyear = 25

sstime = 365* nyear * nloop #*25*50

# 创建一个集群并注册
#nCPUs <- as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE", unset = 1)) # 或者 SLURM_CPUS_PER_TASK
#cat("nCPUs=",nCPUs,'\n')
#cl <- makeCluster(nCPUs) # 40是设置的线程数
#registerDoParallel(cl)
sink(file.path(Rgl,"outTXT",paste0(modi,"_",forc_mod,"_",fexu_oa,"_",fexudmod,"_interval-",interval,"_",testoutnm,"_output.txt")))
print(length(sites))
outdf = data.frame()
#out <- foreach(s = sites, .export = pvars, .packages("Rcpp")) %dopar%  {
for(s in sites){
  stS_time = Sys.time()
  cat("------------ s = ", s,"\n")#, "   at ", stS_time)

  if(simmethod == "simi"){
    simid = soilPsim$sim[s]
  }else if(simmethod == "wrb10"){
    simid = soilPsim$wrb10[s]
  }else if(simmethod == "wrb9"){
    simid = soilPsim$wrb9[s]
  }else if(simmethod == "meanmin10"){
    simid = soilPsim$meanmin[s]
  }else if(simmethod == "meanmin9"){
    simid = soilPsim$meanminrm4[s]
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
    forc_npp0 <- forc_npp_all[s,]
  }else{
    forc_lit0 <- forc_repday(forc_litt_all[s,])
    forc_npp0 <- forc_repday(forc_npp_all[s,])
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
 
  if(length(which(forc_npp0<0))>0){
    forc_npp0[forc_npp0<0] = 0
  }

  if(soilPsim$IGBP[s] ==1){
    fnpp_oa = 1 #森林, 分泌物占npp比例为0.1；但现在forc_npp就是分泌物，所以为1
  }else if(soilPsim$IGBP[s] ==2){
    fnpp_oa = 1 #灌木
  }else if(soilPsim$IGBP[s] ==3){
    fnpp_oa = 1 #草地, 分泌物占npp比例
  }
  
  forc_acid0 <-  forc_npp0 * fnpp_oa * fexu_oa/100  #有机酸比例
  forc_Exud0 <-  forc_npp0 - forc_acid0 #更新后的npp就是分泌物, forc_Exud0是转入derivs_V2_MM_AD_CO2的非有机酸部分的根系分泌物


  #驱动数据改为多年日平均
  if(F){

    nloop = 25*100 # 0->no steady-state; 1250->force data is constant; 50->global 
    nyear = 1
    sstime = 365* nyear * nloop#*25*50

    day = 365
    forc_st0 = rep(mean(forc_st0, na.rm = T),day)
    forc_sw0 = rep(mean(forc_sw0, na.rm = T),day)
    forc_lit0 = rep(mean(forc_lit0, na.rm = T),day)
    forc_acid0 = rep(mean(forc_acid0, na.rm = T),day)
    forc_Exud0 = rep(mean(forc_Exud0, na.rm = T),day)
  }


  #如果n凋落物为0，结束此循环
  npar = 20
  
  if(mean(forc_lit0) <= 1e-5 | is.na(mean(forc_lit0))){
    print("average litter input is 0")
    outss = rep(NA,(7*4+1)*npar)
    # outss = rep(NA,(8+7*3+1)*npar)
    #next
  }else{
#    nloop = 100 # 0->no steady-state; >500->force data is constant; others->global 
#    nyear = length(forc_st0)/365
#    nyear = 25
#    sstime = 365* nyear * nloop#*25*50
    
    #土壤属性
    fitpars$depth <- 0.2 #m
    fitpars$param_pH = soilPsim[s,"pH"]
    fitpars$param_bulkd <- soilPsim[s,"BD"]
    fitpars$param_clay <- soilPsim[s,"clay"]
    fitpars$param_claysilt <- soilPsim[s,"clay"] + soilPsim[s,"silt"]
    fitpars$porosity <- soilPsim[s,"porosity"]
    fitpars$param_em <- soilPsim[s,"param_em"]
    calpars = pars_resl[[Simsite]]$par
    
    soilCin = c(LIT=1, SOM=4, POM=1, LMWC=1, MIC=1, MAOM=1, MA = 0, MD = 0, Tco2=0)
    clab = c(FALSE,TRUE)[2]
    soilCin13 = c(SOM=0, POM=0, LMWC=0, MIC=0, MAOM=0, MA = 0, MD = 0, Tco2=0)
    soilCinOld = soilCin13  # 2025.08.20
    soilCinPlantC2toLMWC = soilCin13 # 2025.09.01
    
    outss = c()
    for(pi in 1:npar){
      fitpars[colnames(calpars)]=calpars[pi,]
      soilCin["MA"] = fitpars$r0 * soilCin["MIC"]
      soilCin["MD"] = (1-fitpars$r0) * soilCin["MIC"]
      
      outcpp = derivs_V2_MM_AD_CO2(sstime,fitpars,as.list(soilCin),forc_st0,forc_sw0,forc_lit0,
                                   forc_acid0,forc_Exud0,outflux,nloop,clab, as.list(soilCin13),FALSE, 
                                   as.list(soilCinOld), as.list(soilCinPlantC2toLMWC)
                                  ) #return list,class(outcpp[[1]]) is numeric
#      print(outcpp)
      outss = c(outss, unlist(outcpp))
    }#pi ends
    
  }#lit > 0
  cat("---s=",s,"Total time: ", as.numeric(Sys.time()-stS_time, units = "mins"), "\n")
  outname = c(fixCp[1:7],
              paste0(fixCp[1:7],"_c13"),
              paste0(fixCp[1:7],"_oldc"), 
              paste0(fixCp[1:7],"_plantc2lmwc"), 
              "ssgx")
  names(outss) = paste0(rep(outname,npar),".",rep(1:npar,each=length(outname)))
  cat("输出数据维度:", dim(outdf), "\n")
  cat("输出列数:", ncol(outdf), "\n")
  outdf = rbind(outdf,outss)
  
  
  #cat(dim(outdf))
#  return(outss)
}#s ends

# 结束任务
#stopCluster(cl)
print("foreach ending")
curtime = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
cat("site time:",curtime,"total time: ",as.numeric(Sys.time()-stR_time,units = "mins"),"\n")

#outss = do.call(out,rbind)
colnames(outdf) = paste0(rep(outname,npar),".",rep(1:npar,each=length(outname)))
rownames(outdf) = sites
# save(outdf,file = file.path(outdir,paste0(modi,"_",forc_mod,"_sampling_",interval,"_nloop",nloop,"_ssloop2.Rdata")))
# save(outss,file = file.path(outdir,paste0(modi,"_",forc_mod,sprintf("_ss_%04d_%04d.Rdata",sts,ends))))
save(outdf,file = file.path(outdir,paste0(modi,"_",forc_mod,"_fexu_oa-",fexu_oa,"_",fexudmod,"_ss_0001_",max(sites),"_",testoutnm,".Rdata")))
print("outdf save finished")


