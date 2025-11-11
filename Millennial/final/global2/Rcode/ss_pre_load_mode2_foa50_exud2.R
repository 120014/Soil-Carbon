
modi = c(1,2,3)[2]
forc_mod = c("CMCC-ESM2","NorESM2-MM","CESM2-WACCM")[modi]
print(forc_mod)
fexu_oa = c(25, 50,75)[2] #有机酸比例 %
fexudmod = c("mean", "BCC-CSM2", "NorEMS2")[3]

#选择运行平台
run_platf = c("mine","hangdian","linux")[3]
if(run_platf=="mine"){
  tempdir <- function() "D:/code/R/Rtemp"
}else if(run_platf=="hangdian"){
  tempdir <- function() "D:/liucq/Rtemp"
}else if(run_platf=="linux"){
  tempdir <- function() "/datanode05/liucq/python/Rtmp"
}

unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir",tempdir, baseenv())
lockBinding("tempdir",baseenv())
print(tempdir())
print(tempfile())

usgl = "/datanode05/liucq/python/Millennial/final"
outdir =  file.path(usgl,"global2/output/ss")
indirgl = file.path(usgl,"global2/input")
Rgl =file.path(usgl,"global2/Rcode")
inoutss = file.path(usgl,"global2/output/ss") #trans,ssp

load(file.path(usgl,"global2/input",paste0(modi,"_forc_all_",forc_mod,"_hist.Rdata")))
load(file.path(usgl,"global2/input","soilPsim2.Rdata"))
#增加更新后分泌物数据 2025.7.15
rm(forc_npp_all)
if(fexudmod == "mean"){
  load(file.path(usgl,"global2/input/forc_exud_allCase.RData")) #变量名为：forc_exud_hist, forc_exud_ssp126, forc_exud_ssp585
  forc_npp_all = as.matrix(forc_exud_hist)
}else {
  load(file.path(usgl,"global2/input/forc_exud_allCase_forc_exud_model2.RData")) #变量名为：forc_exud_hist, forc_exud_ssp126, forc_exud_ssp585
  forc_npp_all = as.matrix(forc_exud_model2[[fexudmod]]$forc_exud_hist)
}

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

print("preprocess finished!")
sink()
source(file.path(Rgl,"ss_trans_processfun.R"))
#(1).solve steady state
#source(file.path(Rgl,"1_ss_run_pfor.R"))
source(file.path(Rgl,"1_ss_run.R"))
print("ss_run is finished!")



