#标记

#实验数据准备已经运行函数、提取数据函数等
#-------
testDat = c("Without Oxalic acid","Adding Oxalic acid","both")[3]
expT = c("15","25","35")#[2] #实验温度
if(testDat=="Without Oxalic acid"){
  Jacid = FALSE #是否添加有机???
  # xlsxobs = "不添加草???"
}else if(testDat=="Adding Oxalic acid"){
  Jacid = TRUE #是否添加有机???
  #  xlsxobs = "添加草酸"
}else{
  Jacid = 2 #是否添加有机???
}

xlsxobs = c("不添加草酸","添加草酸")
print(xlsxobs)

C13 = c(F,T)[1] #是否跟踪C13
vlidCO2 = c("day","hr",NULL)[1]  #CO2验证数据是每天还是累???
if(is.na(vlidCO2)){
  vlidCO2 = NULL
}

#-------

pool.names = c("LIT", "POM","LMWC","MIC","MAOM")#add_litter_pool
npools = length(pool.names)
#测的时间
simd = 32 #模拟时间
ML10 = c(2,10)
ML30 = c(3,30)
CO2_10 = seq(1,10)
CO2_30 = c(1,2,3,5,7,10,20,30)
c13CO2_10 = c(2,5,7,10)
c13CO2_30 = c(1,3,5,7,20)
nsite = 10 #10个土壤点
#CO2速率
# CO2_10_r = c(0,CO2_10)
# CO2_30_r = c(0,CO2_30)
# c13CO2_10_r = c(0,c13CO2_10)
# c13CO2_30_r = c(0,c13CO2_30)

#--------驱动数据和碳库背景???
Jexp = read.xlsx(file.path(indirexp,"SITES.xlsx"),sheet="SITES")
# kb = Jexp$`kb-cmfd`
# fmmaoc = Jexp$`fmmaoc-cmfd`
# cue = Jexp$CUE
inputdata = data.frame(site =Jexp$Site, LIT=Jexp$Lit, POM=Jexp$POC,LMWC=Jexp$DOC,MIC=Jexp$MBC,MAOM=Jexp$MAOC,#add_litter_pool
                       pH=Jexp$Ph,BD=Jexp$BULK,depth = 0.2,clay=Jexp$CLAY,claysilt=Jexp$claysilt,
                       SoilTexture = Jexp$SoilTexture,porosity =Jexp$porosity,
                       # forc_st=expT[1],forc_sw=Jexp$Swe,forc_lit=0,forc_acid = Jexp$草酸输入量,
                       forc_st=expT[1],forc_sw=Jexp$Swe,forc_lit = 0,forc_acid = Jexp$草酸输入量,
                       pc=Jexp$pc,param_em = Jexp$foa_maom/100) #,Kb=kb,fmmaoc=fmmaoc,cue=cue
#模拟添加有机???/没添???
if(testDat =="Without Oxalic acid"){
  inputdata$forc_acid=0  
}else if(testDat =="both"){
  inputdata = rbind(inputdata,inputdata)
  inputdata$forc_acid[1:nsite] = 0
}
head(inputdata)
kgtom2 = inputdata$BD * inputdata$depth #(gC/kg soil) to (gC/m2)
inputdata[,c(pool.names,"forc_acid")] = inputdata[,c(pool.names,"forc_acid")]*kgtom2
head(inputdata)

#-------1.观测CO2
#单位换算，mgC/kg soil to (gC/m2)
vlidCO2
expT
CO2unit_kgsoil_tom2<-function(data){
  #data = Obs.co20_tt[,c(3:4)]
  for(i in 1:nsite){
    sitenm = inputdata$site[i]
    data[which(data$site %in% sitenm),2] = data[which(data$site %in% sitenm),2] * kgtom2[i] *1e-3
  }
  return(data)
}
Obs.co20.f=Obs.co21.f =  data.frame()
Obs.co20_tt.f=Obs.co20_tt_hr.f = Obs.co21_tt.f=Obs.co21_tt_hr.f = data.frame()
Obs.co2_C13 = Obs.co2_C13_tt =Obs.co2_C13_tt_hr =data.frame()
for(i in 1:length(expT)){
  #速率
  Obs.co20 = read.xlsx(file.path(indirexp,"实验指标.xlsx"),sheet = paste0(expT[i],"-",xlsxobs[1],"CO2",vlidCO2))
  Obs.co20 = CO2unit_kgsoil_tom2(Obs.co20)
  Obs.co21 = read.xlsx(file.path(indirexp,"实验指标.xlsx"),sheet = paste0(expT[i],"-",xlsxobs[2],"CO2",vlidCO2))
  Obs.co2_C13 = rbind(Obs.co2_C13, CO2unit_kgsoil_tom2(Obs.co21[,c(3:4)]))
  if(i == length(expT)){
    Obs.co2_C13 = Obs.co2_C13[-which(is.na(Obs.co2_C13$site)),]
  }
  Obs.co21 = CO2unit_kgsoil_tom2(Obs.co21[,1:2])
  #合并多个温度CO2
  Obs.co20.f = rbind(Obs.co20.f,Obs.co20)
  Obs.co21.f = rbind(Obs.co21.f,Obs.co21)
  
  
  #累积
  Obs.co20_tt = read.xlsx(file.path(indirexp,"实验指标.xlsx"),sheet = paste0(expT[i],"-",xlsxobs[1],"CO2-累积"))
  Obs.co20_tt_hr = CO2unit_kgsoil_tom2(Obs.co20_tt[,c(3:4)])
  Obs.co20_tt = CO2unit_kgsoil_tom2(Obs.co20_tt[which(!is.na(Obs.co20_tt[,1])),c(1:2)]) #不带小时的累积CO2
  
  Obs.co21_tt = read.xlsx(file.path(indirexp,"实验指标.xlsx"),sheet = paste0(expT[i],"-",xlsxobs[2],"CO2-累积"))
  Obs.co21_tt_hr = CO2unit_kgsoil_tom2(Obs.co21_tt[,c(3:4)])
  Obs.co2_C13_tt = rbind(Obs.co2_C13_tt, CO2unit_kgsoil_tom2(Obs.co21_tt[which(!is.na(Obs.co21_tt[,5])),c(5:6)]))
  Obs.co2_C13_tt_hr = rbind(Obs.co2_C13_tt_hr, CO2unit_kgsoil_tom2(Obs.co21_tt[which(!is.na(Obs.co21_tt[,7])),c(7:8)]))
  Obs.co21_tt = CO2unit_kgsoil_tom2(Obs.co21_tt[which(!is.na(Obs.co21_tt[,1])),c(1:2)])
  
  #合并多个温度
  Obs.co20_tt.f = rbind(Obs.co20_tt.f,Obs.co20_tt)
  Obs.co20_tt_hr.f = rbind(Obs.co20_tt_hr.f,Obs.co20_tt_hr)
  Obs.co21_tt.f = rbind(Obs.co21_tt.f,Obs.co21_tt)
  Obs.co21_tt_hr.f = rbind(Obs.co21_tt_hr.f,Obs.co21_tt_hr)
  
}

if(testDat=="Without Oxalic acid"){
  Obs.co2 = Obs.co20.f #CO2速率
  Obs.co2_tt = Obs.co20_tt.f #累积CO2
  Obs.co2_tt_hr = Obs.co20_tt_hr.f #???4???8h的C累积CO2
  #  Obs.co2_C13 = 0
}else if(testDat=="Adding Oxalic acid"){
  Obs.co2 = Obs.co21.f #CO2速率
  Obs.co2_tt = Obs.co21_tt.f #累积CO2
  Obs.co2_tt_hr = Obs.co21_tt_hr.f #???4???8h的C累积CO2
  
}else if(testDat=="both"){
  Obs.co2 =  rbind(Obs.co20.f,Obs.co21.f)
  Obs.co2_tt =  rbind(Obs.co20_tt.f,Obs.co21_tt.f)
  #  Obs.co2_tt = Obs.co2_tt[-which(is.na(Obs.co2_tt$CO2)),]
  Obs.co2_tt_hr =  rbind(Obs.co20_tt_hr.f,Obs.co21_tt_hr.f)
  
}
head(Obs.co2)
max(Obs.co2[,"CO2"])
#由于率定参数数据库的组成是，没添加草酸-c(黄棕1，黑土..)每个土壤*3个温度


#2.观测LMWC和MIC-------
poolunit_kgsoil_tom2<-function(Obs.pools){
  i=p=1 #单位换算，mgC/kg soil to (gC/m2)
  while(p<nrow(Obs.pools)){
    # print(paste("p",p,"i",i))
    Obs.pools[p:(p+1),c("MIC","LMWC")] = Obs.pools[p:(p+1),c("MIC","LMWC")] * kgtom2[i] *1e-3#add_litter_pool
    p=p+2
    i=i+1
  }
  return(Obs.pools)
}

Obs.pools0.f = Obs.pools1.f =data.frame()
for(i in 1:length(expT)){
  Obs.pools0 = read.xlsx(file.path(indirexp,"实验指标.xlsx"),sheet = paste0(expT[i],"-",xlsxobs[1]))
  Obs.pools1 = read.xlsx(file.path(indirexp,"实验指标.xlsx"),sheet = paste0(expT[i],"-",xlsxobs[2]))
  Obs.pools0 = poolunit_kgsoil_tom2(Obs.pools0)
  Obs.pools1 = poolunit_kgsoil_tom2(Obs.pools1)
  
  #合并多个温度
  Obs.pools0.f = rbind(Obs.pools0.f,Obs.pools0)
  Obs.pools1.f = rbind(Obs.pools1.f,Obs.pools1)
}

if(testDat=="Without Oxalic acid"){
  Obs.pools =  Obs.pools0.f
}else if(testDat=="Adding Oxalic acid"){
  Obs.pools =  Obs.pools1.f
}else if(testDat=="both"){
  Obs.pools =  rbind(Obs.pools0.f,Obs.pools1.f)
}

head(Obs.pools)
max(Obs.pools[,c("LMWC","MIC")])#add_litter_pool

#3.参数-------
pars.in <- NULL
pars.in.raw <- read.table(file.path(indir,"soilpara_in_fit1_clab.txt"))
pars.in <- as.list(pars.in.raw$V2)
names(pars.in) <- pars.in.raw$V1
pars <- unlist(pars.in)
if(testMd == "_rmPOCtoMAOC"){
  pars["param_fpl"] = 1
}

#去掉MAOC的米氏方程
#if(testMd == "_rmMAOCmm"  | testMD2  | testMd == "_fixbeta_rmMAOCmm"){
#  pars["Vml0"] = 0
#}

if(testMd == "_fixbeta"){
  pars["beta"] = 0.001
}
print(pars)
# print(pars)
#4.模型参数设置--------
fixCp = c("SOM","LIT","POM","LMWC","MIC","MAOM","MA","MD","Tco2")#add_litter_pool
flux.cul = c(FALSE,TRUE)[2] #通量输出累计还是每天的
out.allflux = c(FALSE,TRUE)[1] #是否输出所有通量
clab = switch(1,FALSE,c("Cin","Cinoa","Cad")[2]) #标记C类型 #设为全局变量
clab.outvar = c("Cpools_Tco2", "all")[2] #标记C输出哪些变量
sav.init =  c(FALSE,TRUE)[2] #输出结果是否要保存初始值
allflux.name = c("f_Lit_PO", "f_Lit_LM", "f_Exud_LM", "f_Oa_MA", "f_Oa_LM",
                 "f_PO_LM", "f_PO_MA", "f_MA_LMm", "f_MA_LMe", "acid_des","f_MA_LM","f_LM_MA", 
                 "f_MA_MD", "f_MD_MA", "f_LM_MB",  "f_MA_maint_co2", "f_MA_growth_co2",
                 "f_MA_co2", "f_MD_co2", "f_MA_growth",
                 "f_MB_LM", "f_MB_MA", "f_MB_turn","f_LM_leach")
#只输出总碳库和Tco2
nfixvar = length(fixCp)
fixvar = fixCp
if(out.allflux){
  Cout.nflux = nfixvar+length(allflux.name)
  colnms = c(fixvar,allflux.name)
}else{
  Cout.nflux = nfixvar
  colnms =  fixvar
}
#输出标记碳库
if(clab[1]!=FALSE){
  MAOM_c130 = 0 #草酸解吸附初始值
  if(clab.outvar=="all"){
    Cout.nlab = length(c(fixCp,allflux.name[-1:-5])) * length(clab)
    colnms.lab = paste0(c(fixCp,allflux.name[-1:-5]),"_c13_",rep(clab,each=length(c(fixCp,allflux.name[-1:-5]))))
  }else if(clab.outvar=="Cpools_Tco2"){
    Cout.nlab = Cout.nflux * length(clab)
    colnms.lab = paste0(fixCp,"_c13_",rep(clab,each=length(fixCp)))
  }
}else{
  Cout.nlab = 0
  colnms.lab = NULL
}
TCout.nflux = Cout.nflux + Cout.nlab
Tcolnms = c(colnms, colnms.lab)
cat("TCout.nflux=",TCout.nflux,"Tcolnms = ",Tcolnms,'\n')

# Cpools.nm = c("POM","LMWC","MIC","MAOM","MA","MD")
#模拟实验
#5.动态模拟-------
Model.pools <- function(parsin, inputs, expT=expT) {
  #输出结果
  
  output_csv <- TRUE  # 默认不输出CSV
  # output_csv <- FALSE  # 默认不输出CSV
  
  state.run = data.frame()
  for(s in 1:dim(inputs)[1]){ # dim(inputs)[1]
    cat(inputs$site[s],'\n')
    #print(s)
    #parsin = pars
    #parsin = fitpars
    #parsin = fitpars
    #parsin = pars.in
    #inputs = inputdata[which(inputdata$site %in% Simsite),]
    #inputs = inputdata[9,] s=1
    #s=2
    #flux = T
    #forc_st <- inputs[s,"forc_st"]
    forc_sw <- inputs[s,"forc_sw"]
    # forc_lit <- inputs[s,"forc_lit"]
    forc_acid0 <- inputs[s,"forc_acid"]
    forc_lit <- 0  #add_litter_pool
    forc_exud <- 0
    
    parsin$param_pH = inputs[s,"pH"]
    parsin$param_bulkd <- inputs[s,"BD"]
    parsin$depth <- inputs[s,"depth"]
    parsin$param_clay <- inputs[s,"clay"]
    parsin$param_claysilt <- inputs[s,"claysilt"]
    parsin$porosity <- inputs[s,"porosity"]
    parsin$param_em <- inputs[s,"param_em"]
    
    #培养实验的凋落物???0，DOC淋溶速率???0
    parsin$rate_leach = 0
    
    #模拟时间
    run.steps <- seq(1,simd) #培养实验30???
    run.steps.forc <- seq(1,simd+1) #插值函数需要大???1???
    
    #Define forcing functions
    forc_sw <- approxfun(run.steps.forc, rep(forc_sw,simd+1),method = "const") #moisture input function
    forc_lit <- approxfun(run.steps.forc,rep(forc_lit,simd+1),method = "const") #litter input function
    forc_exud <- approxfun(run.steps.forc,rep(forc_exud,simd+1),method = "const") #exudation input function
    if(forc_acid0==0){
      forc_acid <- approxfun(run.steps.forc,rep(forc_acid0,simd+1),method = "const")
    }else{
      forc_acid0 = c(forc_acid0,rep(0,simd)) #只有第一天输入有机酸
      forc_acid <- approxfun(run.steps.forc,forc_acid0)
    }
    
    #温度循环
    for(Ti in 1:length(expT)){
      # cat("Ti = ",expT[Ti],'\n')
      #Ti = 1
      soilCin = c(LIT = inputs[s,"LIT"], #add_litter_pool
                  #LIT = inputs[s,"LIT"],
                  POM = inputs[s,"POM"],LMWC = inputs[s,"LMWC"] , 
                  MIC = inputs[s,"MIC"], MAOM=inputs[s,"MAOM"] ,
                  MA = 0, MD = 0, Tco2=0) #
      #初始化激活态MA和休眠态MD
      soilCin["MA"] = parsin$r0 * soilCin["MIC"]
      soilCin["MD"] = (1-parsin$r0) * soilCin["MIC"]
      
      forc_st <- approxfun(run.steps.forc, rep(as.numeric(expT[Ti]),simd+1),method = "const") #temperature input function
      
      # print(paste("温度：",as.numeric(expT[Ti])))
      out.derivs = derivs_Model_step_trans0_clab(simd,soilCin,parsin, forc_st,forc_sw,
                                                 forc_lit,forc_exud,forc_acid)
      # cat("out.derivs列数=",ncol(out.derivs),"state.run列数=",ncol(state.run),'\n')
      
      out.derivs = cbind(site=inputs$site[s],simT=as.numeric(expT[Ti]),out.derivs)
      # print(head(out.derivs,4))
      # out.derivs$site = inputs$site[s]
      # out.derivs$simT = as.numeric(expT[Ti])
      #接受derivs返回???
      if(s==1 & Ti==1){
        state.run = out.derivs
      }else{
        if(!all(colnames(out.derivs)==colnames(state.run))){
          nc1 = which(colnames(out.derivs) != colnames(state.run))
          cat("out.derivs列名:",colnames(out.derivs)[nc1],"state.run列名:",colnames(state.run)[nc1],'\n')
        }
        # cat("modle colnames=",all(colnames(out.derivs)==colnames(state.run)),'\n')
        state.run = rbind(state.run, out.derivs)
      }
    } #Ti end
  } #s end
  # print(state.run[nrow(state.run),])
  # print(head(state.run))
  
  # 检查这些列是否存在
  carbon_pools_columns <- c("time", "SOM", "LIT", "POM", "LMWC", "MIC", "MAOM", "MA", "MD")
  
  # 确保 state.run 包含这些列
  print("state.run 的列名:")
  print(colnames(state.run))
  print("需要的列名:")
  print(carbon_pools_columns)
  
  # 检查哪些列不存在
  missing_cols <- setdiff(carbon_pools_columns, colnames(state.run))
  if (length(missing_cols) > 0) {
    cat("缺少的列:", missing_cols, "\n")
  }
  
  
  
  
  
  
  csv_filename = file.path(outdir,"pools.csv")
  # 输出CSV文件
  if(output_csv && !is.null(csv_filename)) {
    # 确保目录存在
    dir.create(dirname(csv_filename), showWarnings = FALSE, recursive = TRUE)
    
    # 选择需要输出的碳库列
    carbon_pools_columns <- c("time", "SOM", "LIT", "POM", "LMWC", "MIC", "MAOM", "MA", "MD")
    
    # 如果包含变化量，也输出
    if(any(grepl("^d", colnames(state.run)))){
      delta_columns <- grep("^d", colnames(state.run), value = TRUE)
      carbon_pools_columns <- c(carbon_pools_columns, delta_columns)
    }
    
    # 写入CSV
    write.csv(state.run[, c("site", "simT", carbon_pools_columns)], 
              file = csv_filename, row.names = FALSE)
    cat("碳库时间序列已输出到:", csv_filename, "\n")
  }
  
  
  
  
  
  print(tail(state.run,n=2))
  print(dim(state.run))
  # modeled.pools <- as.data.frame(cbind(state.run$site, state.run[,-ncol(state.run)]))
  # names(modeled.pools) <- c("site","day",names(state.run)[c(-1,-ncol(state.run))])
  print("Model.pools is finished !")
  return(state.run)
}
# modelt = Model.pools(pars,inputdata)
# sim.CO2(modelt)
# sim.ML(modelt)
#一个点多种情况
vector_gap<-function(st,gap){
  #st = CO2_10
  #gap = simd
  result = c()
  i=1
  while(i <=(nSimsamep*length(expT))){
    #print(i)
    midvector = st + (i-1)*gap
    result = c(result, midvector)
    i=i+1
  }
  return(result)
}

#求多列相邻行之差
# 计算data.frame中多列相邻行之差
calculate_diff <- function(df, cols=ncol(df)) {
  
  #df = optim.out
  # print(head(df))
  # cat("class(df)=",class(df),"dim(df)=",dim(df),'\n')
  result <- data.frame(matrix(nrow = nrow(df), ncol = cols))
  colnames(result) <-  cols
  
  for (i in 1:ncol(result)) {
    diff_values <- diff(df[,i]) #数据框的某一列
    result[2:nrow(df), i] <- diff_values
  }
  
  #为满足与验证数据相应行
  result[1,] = 0
  result[which(result[,1]<0),]=0
  colnames(result) = colnames(df)
  return(result)
}




