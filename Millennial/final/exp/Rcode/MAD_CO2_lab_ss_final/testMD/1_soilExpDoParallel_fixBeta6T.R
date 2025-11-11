#选择运行平台
run_platf = c("mine","hangdian","linux")[3]
if(run_platf=="mine"){
  tempdir <- function() "D:/code/R/Rtemp"
}else if(run_platf=="hangdian"){
  tempdir <- function() "D:/liucq/Rtemp"
}else if(run_platf=="linux"){
  tempdir <- function() "/datanode05/liucq/Rtemp"
}

unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir",tempdir, baseenv())
lockBinding("tempdir",baseenv())


#10个点下标
sts = 6
ends = 6
ngen = 300 #迭代次数
npop = 160 #种群数
testMd = c("_parsMax", "_parsMaxbeta","_fixbeta","_objCO2R2MARE","_rmPOCtoMAOC","_VMA_VMD","_rmMAOCmm")[3]
testMD2 = c(FALSE, TRUE)[2] #范围缩小、优化目标去掉、MAOC米氏方程去掉
sstype = c("runsteady","ND","ND_cpp")[3]
print(paste("site=",sts,"generations=",ngen,"popsize=",npop))
#有机质分解实???
library(openxlsx)
# library(deSolve)
library(mco)
library(dplyr)
library(tidyr)
library(doParallel)
library(foreach)

if(run_platf=="mine"){
  outdir = "E:/code/R/Millennial/final/exp/output"
  indir = "E:/code/R/Millennial/final/input"
  indirexp = "E:/code/R/Millennial/final/exp/input"
  #indirR = "E:/code/R/Millennial/final/Rcode/Fld"
  indirRexp = "E:/code/R/Millennial/final/exp/Rcode/MAD_CO2_lab_ss_final"
}else if(run_platf=="hangdian"){
  outdir = "D:/liucq/Millennial/final/exp/output"
  indir = "D:/liucq/Millennial/final/input"
  indirexp = "D:/liucq/Millennial/final/exp/input"
  indirR = "D:/liucq/Millennial/final/Rcode/Fld"
  indirRexp = "D:/liucq/Millennial/final/exp/Rcode/MAD_CO2_lab_c13"
}else if(run_platf=="linux"){
  outdir = "/datanode05/liucq/python/Millennial/final/exp/output"
  indir = "/datanode05/liucq/python/Millennial/final/input"
  indirexp = "/datanode05/liucq/python/Millennial/final/exp/input"
  # indirR = "/data01nfs/user/liucq/python/Millennial/final/Rcode/Fld"
  indirRexp = "/datanode05/liucq/python/Millennial/final/exp/Rcode/MAD_CO2_lab_ss_final"
}

source(file.path(indirRexp,"derivs_V2_MM_AD_CO2_Clab_vmax_cuetime.R"))
source(file.path(indirRexp,"ss_trans_processfun.R"))
source(file.path(indirRexp,"expdata_fun_clab.R"))
source(file.path(indirRexp,"Pars_expcal_clab_ss.R"))



#并行
if(run_platf=="linux"){
  cl <- makePSOCKcluster(length(ends-sts + 1),
                         timeout = 3600,  # 延长超时（单位：秒）
                         outfile = "/dev/null",  # 关闭日志输出（避免IO阻塞）
                         homogeneous = TRUE,  # 强制所有节点环境一致
                         rscript_args = "--max-ppsize=100000")  # 增加参数栈大小) #当使用的核心数太多时，使用此方法构建
  
}else {
  cl = makeCluster(2)
}

registerDoParallel(cl)  # 关键步骤：注册并行后端
# 在每个节点加载包
clusterEvalQ(cl, {
  library(openxlsx)
  # library(deSolve)
  library(mco)
  library(dplyr)
  library(tidyr)
})
clusterExport(cl,  varlist = ls(envir = .GlobalEnv))  # 共享变量到各节点

#长度365*15
foreach(si = sts:ends, .options.multicore = list(preschedule = FALSE)) %do% {
  
  source(file.path(indirRexp,"Pars_expcal_clab_ss.R"))

    # Fit the model
  Fit.pools = list() #所有点最优参数nsga2
  Fit.pools_iter = list() #所有点每次迭代参数nsga2
  
  { #:10 for(si in vi)
    #si = sts
    start_time = Sys.time()
    print(start_time)
    
    print(si)
    itensga2= 0
    fitpars_iter = list()
    
    Simsite = c("黄棕壤1","黑土","栗钙土","火山灰土","红壤","黄壤",
                "棕壤","草甸土","黄棕壤2","赤红壤")[si] #用火山灰土和棕壤测试
    print(Simsite)
    if(Simsite=="黄棕壤2"){
      parToFit.upper["cue_ref"] = 0.9
    }
    nSimsamep = dim(inputdata[which(inputdata$site %in% Simsite),])[1]
    print(parToFit.lower)
    print(parToFit.upper)
    # sink(file.path(outdir,paste0(Simsite,"_",ngen,"_",npop,"_vmax_vgcue_cueacid_time",sstype,"_nsga2Test.txt"))) #输出到文件中的信息同时打印在console???, split=TRUE,append=TRUE)
    
    Fit.pools[[Simsite]] <- nsga2(Objective.run,idim=length(parToFit.lower),
                                  odim=nobj,#mdist=20,mprob=0.5, #cdist=20,
                                  generations = ngen,popsize = npop,#popsize???4的倍数
                                  lower.bounds = parToFit.lower,upper.bounds = parToFit.upper,
                                  constraints = consfun, cdim = naga2condim)
    
    Fit.pools_iter[[Simsite]] <- fitpars_iter #所有迭代结果
    
  }
  
  end_time = Sys.time()
  print(end_time)
  print(end_time-start_time)
  
  #保存率定参数250_160/
  save(Fit.pools,Fit.pools_iter,
       file = file.path(outdir,"V6Rdata",paste0(si,"_",Simsite,"_vmax_vgcue_cueacid_time",sstype,testMd,"_",testMD2,testoutnm,".RData")))
  print(paste(Simsite,"finished!"))
  
  cat("itensga2=: ",itensga2,'\n') #"ss.num=",ss.num,
  tail(fitpars_iter$value)
  
  print("finished!")
}

print("foreach ending")

stopCluster(cl)

#load(file.path(outdir,"clab_Rdata","9_黄棕壤2__vmax_vgcue_cueacid_timeND_cpp.RData"))
#Fit.pools_iter
