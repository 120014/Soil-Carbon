#标签添加稳态求解

#求每个点的4个驱动数据
#参数率定稳态部分
#-------------
#加载动态tif提取10个样点的动态数据
if(F){
  forc_mod = c("NorESM2","CMCC-ESM2")[2]
  # load(paste0("D:/liucq/Millennial/Millennial3/input/",forc_mod,"_forctif.Rdata"))
  load(paste0("E:/code/R/Millennial/Millennial3/Mill_linux/input/",forc_mod,"_forctif.Rdata"))
  # crs(st_tif) = CRS("+init=epsg:4326") 
  #提取点的25年平均
  library(raster)
  library(openxlsx)
  # Jexp = read.xlsx("D:/liucq/Millennial/Millennial3/input/SITES.xlsx",sheet="SITES")
  Jexp = read.xlsx("E:/code/R/Millennial/final/input/exp/SITES.xlsx",sheet="SITES")
  
  map.points = Jexp[,c("Lon","Lat")]
  forc_st_all = raster::extract(st_tif,map.points)  - 273.15 #nc文件中温度K转为°
  forc_sw_all = raster::extract(sw_tif,map.points) /1000/0.2 #kg/m2 -> 0-1
  forc_npp_all = raster::extract(npp_tif,map.points) *1000*60*60*24 #gC/m2/day
  
  forc_acid_all = forc_npp_all * 0.1 
  forc_acid_all[forc_acid_all<0] = 0
  
  forc_litt_all = raster::extract(flitt_tif,map.points)*1000*60*60*24 #gC/m2/day
  dim(forc_st_all)
  
  #save(forc_st_all,forc_sw_all,forc_acid_all,forc_litt_all,file="E:/code/R/Millennial/final/input/exp/forc_exp.Rdata")
}

#----
#直接加载实验点的Rdata数据
load(file.path(indirexp,"forc_exp.Rdata"))
range(forc_st_all)
range(forc_sw_all)
range(forc_acid_all)
range(forc_litt_all)
#某个土壤的Simsite
forc_mean <- function(rowi){
  
  forc_data = c(st = NA,sw=NA,lit=NA,acid=NA)
  forc_data["st"] <- mean(forc_st_all[rowi,]) #temperature input function
  forc_data["sw"] <- mean(forc_sw_all[rowi,]) #moisture input function
  forc_data["lit"] <- mean(forc_litt_all[rowi,]) #moisture input function
  forc_data["acid"] <- mean(forc_acid_all[rowi,]) #temperature input function
  
  return(forc_data)
}

#参数率定稳态求解函数
#@ 需要调用forc_mean(si)函数求4个动态驱动
#@ derivs_Model_step_exp()求
calpars_SSfun <- function(sstype,fitpars){
  
  inputs = inputdata[which(inputdata$site %in% Simsite)[1],]
  
  #土壤属性参数和实验一直
  fitpars$param_pH = inputs[,"pH"]
  fitpars$param_bulkd <- inputs[,"BD"]
  fitpars$depth <- inputs[,"depth"]
  fitpars$param_clay <- inputs[,"clay"]
  fitpars$param_claysilt <- inputs[,"claysilt"]
  fitpars$porosity <- inputs[,"porosity"]
  fitpars$param_em <- inputs[,"param_em"]
  
  state0 = c(LIT=1, POM=1, LMWC=1, MIC=1, MAOM=1, MA = 0, MD = 0, Tco2=0) #,CO2=0
  state0["MA"] = fitpars$r0 * state0["MIC"]
  state0["MD"] = (1-fitpars$r0) * state0["MIC"]
  
  # print(fitpars)
  forc_stode = forc_mean(si) #多年平均动态驱动
  forc_stode["acid"] = forc_acid_all2[si] #2.更新有机酸,forc_acid_all是有机酸，不是分泌物
  # SStime = 1e8 #300*30*10000 #time at which steady state solution is evaluated
  SStime = 365*nloop #最多loop100年,100是不是太多？
  
  ##稳态求解的动态驱动数据多年平均值forc_acid->forc_acid
  if(sstype == "runsteady"){
    # POM <<- 1 
    # LMWC <<- 1
    # MAOM <<- 1
    # MIC <<- 1
    # MA <<- fitpars$r0 * MIC
    # MD <<- (1-fitpars$r0)* MIC
    # SOC.r = POM + MAOM + LMWC + MIC
    # state0 = SOC.r
    
    ystate = state0
    LM_MI_erro = FALSE
    qmax_cue_error_ss <<- FALSE
    SS.pools = runsteady(y = ystate, func = derivs_runsteady, parms = fitpars, #soil_pars = soil_pars,
                         forc_st=forc_stode["st"], forc_sw=forc_stode["sw"],
                         forc_lit= forc_stode["lit"],forc_exud=forc_stode["acid"],
                         forc_acid = forc_stode["acid"]
                         ,verbose=FALSE, times = c(1, SStime*200))#hmin=2#所有状态变量（例如你提到的 x 和 y）都需要达到稳态。
    
    cat("ss.soc=",SS.pools$y,"ss.attr=",attr(SS.pools, "steady"),'\n')
    outss = as.data.frame(t(as.data.frame(SS.pools$y)))#,POM, LMWC,MIC, MAOM, MA,MD
    colnames(outss) = fixCp[2:8]
    
    if(SS.pools$y["LMWC"]>1500 | SS.pools$y["MIC"]>1000){
      cat("LMWC and MIC is too big!","LMWC=",SS.pools$y["LMWC"],"MIC=",SS.pools$y["MIC"],'\n')
      LM_MI_erro = TRUE
    }
    SST <<- attr(SS.pools, "steady") & SS.pools$y[1]!=1 & ((qmax_cue_error_ss==FALSE & LM_MI_erro==FALSE)) #稳态逻辑值
  }
  if(sstype == "ND"){
    c13werror <<- FALSE
    qmax_cue_error <<- FALSE
    
    #Define forcing functions
    # print(forc_stode["st"])
    forc_st <- approxfun(1:SStime, rep(forc_stode["st"],SStime),method = "const") #temperature input function
    forc_sw <- approxfun(1:SStime, rep(forc_stode["sw"],SStime),method = "const") #moisture input function
    forc_lit <- approxfun(1:SStime,rep(forc_stode["lit"],SStime),method = "const") #litter input function
    forc_acid <- approxfun(1:SStime,rep(forc_stode["acid"],SStime),method = "const")
    forc_exud <- forc_acid
    # print("step")
    
    outss = derivs_Model_step_trans0_clab(SStime,state0,fitpars,forc_st,forc_sw,
                                          forc_lit,forc_acid,forc_exud,ss_cond=TRUE)
    # print("outss")
    outss = outss[1,]#返回向量
    outss[is.na(outss)] = 0
  }
  if(sstype == "ND_cpp"){
    nloop = 25*100 # 0->no steady-state; 1250->force data is constant; 50->global 
    
    nyear = 1
    sstime = 365* nyear * nloop#*25*50
    if(nyear==1){
      day=365
    }else if(nyear>1){
      day=1
    }
    forc_st = rep(forc_stode["st"],day)
    forc_sw = rep(forc_stode["sw"],day)
    forc_lit = rep(forc_stode["lit"],day)
    forc_acid = rep(forc_stode["acid"],day)
    clab = c(FALSE,TRUE)[1]
    soilCin13 = c(POM=0, LMWC=0, MIC=0, MAOM=0, MA = 0, MD = 0, Tco2=0)
    outcpp = derivs_V2_MM_AD_CO2(sstime,fitpars,as.list(state0),forc_st,forc_sw,forc_lit,
                                 forc_acid,FALSE,nloop,clab, soilCin13) #return list,class(outcpp[[1]]) is numeric
    outss = unlist(outcpp)
    names(outss) = c(fixCp[1:8],"ssgx")
  }
  return(outss)
  
}

#
#step.num0为总模拟时长，不考虑有机酸引起的解吸附时forc_acid=0
#@acid_desp:0-无有机酸输入;1-有机酸输入且有机酸解吸附;其他有机酸输入但是有机酸不解吸附
derivs_Model_step_trans0_clab  <- function(step.num0,soilC,parameters,forc_st,
                                           forc_sw,forc_lit,forc_exud,forc_acid=NULL,
                                           ss_cond=FALSE,acid_desp=1) {
  
  step.df = ss.loop = data.frame() #用于稳态判断，多年平均，每列为一个loop
  loopth = 1
  #有机酸输入为0
  # if(acid_desp==0){
  #   forc_acid <- approxfun(1:step.num0,rep(0,step.num0))
  # }
  # #保留第0天的初始值
  inidf0 = c(time = 0, SOM=sum(soilC[2:5]),soilC, rep(0,TCout.nflux-nfixvar))
  names(inidf0) = c("time",Tcolnms)
  
  out.flux = 0
  for(day in 1:step.num0){#
    
    # cat("day=",day,'\n')
    # print(soilC)
    dpools = derivs_V2_MM_AD_CO2_Clab(step.num=day,soilC,parameters=parameters,
                                      forc_st=forc_st, forc_sw=forc_sw, forc_lit= forc_lit,
                                      forc_exud = forc_exud,forc_acid = forc_acid,ss_cond) 
    # print("------")
    # print(dpools)
    
    if(length(which(is.na(dpools)))>0){
      na.df0 = data.frame(matrix(0,nrow = (step.num0-day+1),ncol=TCout.nflux+1))
      colnames(na.df0) = c("time",Tcolnms)
      na.df0[,"time"] = c(day:step.num0)

      if(sav.init){
        step.df = rbind(inidf0,step.df,na.df0) #保留初始值
      }

      if(ss_cond){
        SST <<- FALSE
        return(tail(step.df,n=1)[,fixCp])
      }else{
        return(step.df)
      }
      
    }
    if(length(which(dpools[-c(1:6)]<0))>0){
      print("flux less 0 ")
      print(dpools)
    }
    
    nmi = sapply(strsplit(names(dpools), "_"), function(y) head(y, 1))
    cpoolsCols = which(nmi %in% fixCp[-nfixvar])
    soilC = dpools[cpoolsCols]
    if(flux.cul == TRUE){
      out.flux = out.flux + dpools[-cpoolsCols] #碳库已经更新，排除碳库外其他通量需要更新
      dpools[-cpoolsCols] = out.flux #求累计通量
    }
    
    # cat("nmi=",nmi,"cpoolsCols= ",cpoolsCols,'\n')
    # cat("一天结束输出的state：",dpools,'\n')
    #循环输出返回为data.frame
    step.df = rbind(step.df,c(time=day,dpools))
    colnames(step.df) = c("time",names(dpools))
    
    #稳态判断
    if(ss_cond){
      if(day%%((step.num0)/nloop)==0 | day==step.num0){
        # cat("loop =  ",loopth,"nrow(ss.loop)=",nrow(ss.loop))
        #print(state)
        ss.loop = rbind(ss.loop, colMeans(step.df))
        colnames(ss.loop) = c("time", fixvar)
        #loop5次达到稳态结束预热
        if((loopth>=15)){

          #5个loop都比较稳定则达到稳态
          if(sscondfun(tail(ss.loop,n=10)[,"SOM"])){
            ##根据Guo,Huang判断稳态时的LMWC和MIC
            SST <<- TRUE
            cat("reaching at steady!"," day=",day,"\n")
            return(tail(ss.loop, n=1)[,fixCp]) #只返回最后loop #和loop-1
          }
          
          if(ss.loop[nrow(ss.loop),"LMWC"]>1500 | ss.loop[nrow(ss.loop),"MIC"]>1000){
            cat("No steady over MIC or LMWC content! LMWC=",ss.loop[nrow(ss.loop),"LMWC"],
                "   MIC=",ss.loop[nrow(ss.loop),"MIC"]," day=",day,"\n")
            SST <<- FALSE
            return(tail(ss.loop,n=1)[,fixCp])
          }
        }
        
        loopth = loopth +1
        step.df = data.frame()  #只存储每个loop的碳库
      }
    }#ss ends
  } # day ends
  
  if(ss_cond){
    #稳态返回
    SST <<- FALSE
    cat("nm=",fixCp," day=",day,"\n")
    print(tail(ss.loop,n=1)[,fixCp])
    return(tail(ss.loop,n=1)[,fixCp])
    
  }else{
    
    if(sav.init){
      step.df = rbind(inidf0,step.df) #保留初始值
    }
    colnames(step.df) = c("time",names(dpools))
    return(step.df)
  }

}


#稳态判断函数
sscondfun <-function(ss.soc0){
  #ss.pools0 = ss.pools
  #正常情况的稳态判断
  if(all(ss.soc0)){
   
    #SOC满足abs(loop50-loop49)<1; abs(loop50-loop49)/loop49 * 100 <0.1%
    sscond2 = abs(diff(ss.soc0)) < 1
    sscond3 = abs(diff(ss.soc0))/ss.soc0[-length(ss.soc0)]*100 < 0.1
    
    if(all(sscond2 & sscond3)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
  
}

#runsteady稳态求解
derivs_runsteady <- function(step.num,ystate,parameters,forc_st,forc_sw,forc_lit,forc_exud,forc_acid=NULL) {
  with(as.list(c(ystate,parameters)), {

    # Soil type properties  
    #Equation 10
    kaff_lm = exp(-param_p1 * param_pH - param_p2) * kaff_des
    
    #Equation 11
    param_qmax = depth*param_bulkd * param_claysilt * param_pc
    #约束率定的pc
    # print(paste("param_qmax=",param_qmax,"MAOM=",MAOM,"step.num=",step.num))
    if(param_qmax < MAOM){
      qmax_cue_error_ss <<- TRUE
      print("param_qmax exit")
      return(list(rep(1,length(ystate))))
    }
    
    # Hydrological properties
    forc_sw0 = min(forc_sw,porosity)
    SWP = fSWC2SWP(forc_sw0,SWCsat=porosity)  #MPa,SWP<0,饱和基质势能=0
    matpot = abs(SWP * 1e3)        #kPa
    
    #Equation 4
    scalar_wd = (forc_sw0 / porosity)^0.5
    
    #Equation 15
    scalar_wb = exp(lambda * -matpot) * (kamin + (1 - kamin) * ((porosity - forc_sw0) / porosity)^0.5) * scalar_wd
    
    # Decomposition
    gas_const <- 8.31446
    
    
    T0 = forc_st + 273.15
    Tref = tae_ref + 273.15
    tp_scalar <- function(sCase,Ea=NULL){
      if(is.null(Ea)){
        switch(sCase,
               "MR" = {Ea <- 20e3},  #J/mol
               "LIG"= {Ea <- 53.0e3},
               "CEL"= {Ea <- 36.3e3},
               "POM"= {Ea <- 63909}, 
               "Km"= {Ea <- 30e3},
               "LMWC"= {Ea <- 47e3},
               "MAOM"= {Ea <- 67e3},# {Ea <- 47e3},
               "MD"= {Ea <- 47e3},
               "KM"= {Ea <- 30e3},
               "LIT"= {Ea <- 62e3}         #GPT建议62e3
        ) 
        # print(Ea)
      }
      tp_scalar0 = exp(-Ea/gas_const * (1/T0-1/Tref))
      return(tp_scalar0)
    }
    
    #外源C输入
    
    
    
    
    
    f_Lit_LIT = forc_lit
    k = k_lit * tp_scalar("LIT") * scalar_wd
    
    # kaff_lit = 1 - f_resp
    
    if (LIT > 0) {
      #f_LIT = max(LIT * vmax_lit * kaff_lit,0)
      f_LIT = max(LIT * k,0)
      f_LIT = min(LIT, f_LIT)
      
    }else{
      f_LIT = 0
    }
    
    f_LIT_co2 = f_LIT * f_resp;
    f_LIT_SOC = f_LIT * (1 - f_resp);
    
    
    
    
    
    f_Lit_PO = f_LIT_SOC * param_pi
    f_Lit_LM = f_LIT_SOC * (1. - param_pi)
    f_Exud_LM = forc_exud
    f_Oa_MA = forc_acid * param_em
    f_Oa_LM = forc_acid * (1. - param_em)
    
    ### Equation 2
    # POM -> LMWC
    # print(paste("POM= ",POM,"MIC= ",MIC))
    # vmax_pl = alpha_pl * exp(-eact_pl / (gas_const * (forc_st + 273.15)))
    vmax_pl = Vpl0 * tp_scalar("POM") * scalar_wd
    kaff_pl = kaff_pl * tp_scalar("KM")
    if(POM>0 && MA>0){
      f_PO = max(vmax_pl * POM * MA / (kaff_pl + MA),0)
      f_PO = min(POM, f_PO)
    }else{
      f_PO=0
    }
    f_PO_LM = f_PO * param_fpl
    f_PO_MA = f_PO * (1. - param_fpl)
    
    #输入有机酸：因为有机酸反应特别快，所以要在加入后立即更新 #添加草酸立马吸附到MAOM和LMWC
    # LMWC = LMWC + f_Oa_LM
    # MAOM = MAOM + f_Oa_MA
    #立马发生解吸附（几分钟就完成了）
    unit_kgTom2 = param_bulkd * depth *1e-3 #(mgC/kg soil) to (gC/m2)
    if( MAOM > 0 ){
      acid_des = ((MAOM/unit_kgTom2 *1e-3)^0.04 *(forc_acid/unit_kgTom2)^0.86
                  * param_clay^(-0.24) * param_pH^(-2.14) * exp(3.52))* unit_kgTom2
      acid_des = min(max(acid_des,0),MAOM)
    }else{
      acid_des = 0
    }
    
    LMWC = LMWC  + f_Oa_LM + acid_des
    MAOM = MAOM + f_Oa_MA - acid_des
    
    ###MAOM的酶促分解反应MAOM->LMWC
    # vmax_ml = max(alpha_ml * exp(-eact_ml / (gas_const * (forc_st + 273.15))),0)
    vmax_ml = Vml0 * tp_scalar("MAOM") * scalar_wd
    kaff_ml = kaff_ml * tp_scalar("KM")
    if(MAOM>0 && MA>0){
      f_MA_LMe = max(vmax_ml * MAOM * MA / (kaff_ml + MA),0)
      f_MA_LMe = min(MAOM, f_MA_LMe)
    }else{
      f_MA_LMe=0
    }
    
    #Equation 12
    # MAOM -> LMWC
    if(MAOM>0){
      f_MA_LMm = max(kaff_des * MAOM / param_qmax,0)
      f_MA_LMm = min(MAOM, f_MA_LMm)
    }else{
      f_MA_LMm = 0
    }
    
    
    #Equation 14
    # vmax_lb = alpha_lb * exp(-eact_lb / (gas_const * (forc_st + 273.15)))
    Vg0 = Vlb0 #假设MEND中的的Vg就是Vd
    Vg = Vg0 * tp_scalar("LMWC") # * scalar_wb ?MA的生长呼吸速率，类似MEND中的Vd
    
    #MA维持呼吸速率：α/(1-α) 占MA生长呼吸的比例。微生物激活MA和休眠MD周转
    Vm = Vg0 * Ma/(1-Ma) * tp_scalar("MR")
    VmA2D = Vm # * tp_scalar("MR") #* scalar_wd #* wp_scalar_low
    VmD2A = Vm # * tp_scalar("MR") #* scalar_wd #* wp_scalar
    
    #MD维持呼吸速率: 是MA维持呼吸的β倍
    VmD = Vm * beta
    
    #Equation 22
    # microbial growth flux, but is not used in mass balance
    CUE = cue_ref - cue_t * (forc_st - tae_ref) 
    CUE = CUE - acue * forc_acid    #加入草酸后CUE降低
    #    print(paste("step.num=",step.num, "CUE=",CUE,"acue=",acue))
    CUE_max = 0.95 #0.9还是0.95？
    CUE_min = 0.01
    if(CUE<CUE_min | CUE>CUE_max){
      qmax_cue_error_ss <<- TRUE
      print("CUE NA exit")
      print(paste("step.num=",step.num, "CUE=",CUE,"acue=",acue))
      return(list(rep(1,length(ystate))))
    }
    
    # LMWC -> MA
    vmax_lb = (Vg + Vm) / CUE
    kaff_lb = kaff_lb * tp_scalar("KM")
    
    #Equation 13
    # LMWC -> MIC
    if(LMWC>0 && MA>0){
      f_LM_MB = max(vmax_lb * MA * LMWC / (kaff_lb + LMWC),0)
      f_LM_MB = min(LMWC, f_LM_MB)
    }else{
      f_LM_MB=0
    }
    # print(paste("LMWC=",LMWC,"LMWC被微生物利用",f_LM_MB))
    
    #Equation 21
    # MIC -> atmosphere,#MA呼吸Flb*(1-CUE),就是现在的f_MB_atm
    # f_MA_co2 = f_LM_MB * (1 - CUE) #method1
    f_MA_growth_co2 = max( Vg * (1/CUE - 1) * MA * LMWC / (kaff_lb + LMWC),0)
    f_MA_maint_co2 = max( Vm * (1/CUE - 1) * MA * LMWC / (kaff_lb + LMWC),0)
    f_MA_co2 = f_MA_growth_co2 + f_MA_maint_co2
    # MA = MA + f_MA_growth
    
    LMWC = LMWC - f_LM_MB #立刻更新，用于吸附-解吸附计算
    f_MA_growth = f_LM_MB * CUE
    # cat("f_LM_MB=",f_LM_MB,"f_MA_growth_co2=",f_MA_growth_co2,"f_MA_maint_co2=",f_MA_maint_co2,
    #     "f_MA_co2=",f_MA_co2,"f_MA_growth=",f_MA_growth,"0相等=",f_LM_MB-(f_MA_growth+f_MA_co2),'\n')
    
    #MD维持呼吸
    # cat("VmD = ",VmD,'\n')
    if(MD>0){
      f_MD_co2 = max(VmD * MD,0) #0.001 
      f_MD_co2 = min(MD, f_MD_co2)
    }else{
      f_MD_co2 = 0
    }
    Tco2 = f_MA_co2 + f_MD_co2 + f_LIT_co2
    # cat("MA=",MA,"MD=",MD,"f_MD_co2=",f_MD_co2)
    
    #Equation 8
    # LMWC -> out of system leaching
    if(LMWC>0){
      f_LM_leach = max(rate_leach * scalar_wd * LMWC,0)
      f_LM_leach = min(LMWC, f_LM_leach)
    }else{
      f_LM_leach=0
    }
    
    #Equation 9
    # LMWC -> MAOM
    if(LMWC>0 && MAOM>0){
      f_LM_MA = max(scalar_wd * kaff_lm * LMWC * (1 - MAOM / param_qmax),0)
      f_LM_MA = min(LMWC, f_LM_MA)
    }else{
      f_LM_MA=0
    }
    
    #Equation 16
    # MB -> MAOM/LMWC #微生物死亡周转
    rate_bd = max(0,rate_bd + rate_Kbd*(forc_st-tae_ref))
    if(MA>0){
      # f_MB_turn = rate_bd * MIC^2.0
      f_MB_turn = max(rate_bd * MA,0) #Liucq
      f_MB_turn = min(MA, f_MB_turn)
    }else{
      f_MB_turn=0
    }
    f_MB_LM = f_MB_turn * (1. - param_pb)
    f_MB_MA = f_MB_turn * param_pb
    
    #Equation 14.1
    #微生物激活MA和休眠MD周转,
    VmA2D = Vm # * tp_scalar("MR") #* scalar_wd #* wp_scalar_low
    VmD2A = Vm # * tp_scalar("MR") #* scalar_wd #* wp_scalar
    sat = LMWC/(kaff_lb + LMWC) #Wang 2015
    # sat = (f_MA_growth/(f_MA_growth + f_MA_maint_co2 + f_MB_turn))^2#Huang
    if(LMWC >0 & MA>0){
      f_MA_MD = max((1-sat) * VmA2D * MA,0)
      f_MA_MD = min(MA, f_MA_MD)
    }else{
      f_MA_MD = 0
    }
    if(LMWC >0 & MD>0){
      f_MD_MA = max(sat * VmD2A * MD,0)
      f_MD_MA = min(MD, f_MD_MA)
    }else{
      f_MD_MA = 0
    }
    # MA = MA + f_MD_MA - f_MA_MD #MA和MD立刻更新，因为有机酸加入之后反应非常快！
    # MD = MD + f_MA_MD - f_MD_MA
    
    #碳库损失通量之和<碳库
    if((f_MA_LMe + f_MA_LMm) > MAOM){
      f_MA_LMe = f_MA_LMm = MAOM/2
    }
    f_MA_LM = acid_des + f_MA_LMe + f_MA_LMm
    
    if((f_LM_leach + f_LM_MA)>LMWC){
      f_LM_leach = f_LM_MA = LMWC/2
    }
    
    if((f_MD_MA + f_MD_co2)>MD){
      f_MD_MA = f_MD_co2 = MD/2
    }
    
    #所有通量
    # Cinput = c(f_Lit_PO, f_Lit_LM, f_Exud_LM, f_Oa_MA, f_Oa_LM)
    # names(Cinput) = allflux.name[1:5]
    # 
    # flux = c(f_PO_LM, f_PO_MA, f_MA_LMm, f_MA_LMe, acid_des, f_MA_LM, f_LM_MA, 
    #          f_MA_MD, f_MD_MA, f_LM_MB,  f_MA_maint_co2, f_MA_growth_co2, f_MA_co2, 
    #          f_MD_co2,f_MA_growth, f_MB_LM, f_MB_MA, f_MB_turn,  f_LM_leach 
    # )
    # names(flux) = allflux.name[-1:-5]
    
    
    
    dLIT = f_Lit_LIT - f_LIT
    
    #Equation 1
    dPOM = f_Lit_PO - f_PO
    
    #Equation 7
    dLMWC = f_Lit_LM + f_Exud_LM + f_PO_LM + f_MB_LM + (f_MA_LMe + f_MA_LMm) - f_LM_leach  - f_LM_MA 
    
    #Equation 20
    dMA =  f_MA_growth + f_MD_MA - f_MA_MD - f_MB_turn # + f_MD_MA  - f_MA_MD #在上面已更新
    dMD =  f_MA_MD - f_MD_MA - f_MD_co2  # f_MA_MD - f_MD_MA #在上面已更新
    dMIC = dMA + dMD
    # cat("dMA=",dMA,"dMD=",dMD ,'\n') #,"dMIC = ",dMIC
    
    #Equation 19
    # dMAOM = C_rate*param_em + f_LM_MA - f_MA_LM + f_MB_turn * param_pb
    #加草酸立马就吸附到MAOM???
    dMAOM = f_LM_MA - (f_MA_LMe + f_MA_LMm) + f_MB_MA + f_PO_MA 
    
    return(list(c(dLIT, dPOM,dLMWC,dMIC,dMAOM,dMA,dMD)))
    # Update state variables
    # POM <<- POM + dPOM
    # LMWC <<- LMWC + dLMWC
    # MA <<- MA + dMA
    # MD <<- MD + dMD
    # # MIC = MIC + dMIC
    # MIC <<- MA + MD 
    # MAOM <<- MAOM + dMAOM
    
    # dSOC0 = dPOM + dMA + dMD + dMAOM +dLMWC
    # return(list(dSOC0))
  })
}


#求多列相邻行之差
# 计算data.frame中多列相邻行之差
#@cols为df的所有列
calculate_diff <- function(df, cols=ncol(df)) {
  
  #df = optim.out
  result <- data.frame(matrix(nrow = nrow(df), ncol = cols))
  colnames(result) <-  cols
  
  for (i in 1:ncol(result)) {
    diff_values <- diff(df[,i]) #数据框的某一列
    result[2:nrow(df), i] <- diff_values
  }
  
  #为满足与验证数据相应行
  result[1,] = 0
  result[which(result[,1]<0),] = 0
  colnames(result) = colnames(df)
  return(result)
}

#动态数据在大月重复31天，其他月份重复30天，总共365天
forc_repday <- function(v){
  # 初始化一个空向量来存储重复后的结果
  result <- c()
  
  # 循环处理每个月份
  for (i in 1:length(v)) {
    if (i %% 12 %in% c(1, 3, 5, 7, 8, 10, 0)) {
      result <- c(result, rep(v[i], 31))
    }else if (i %% 12 == 2) {
      result <- c(result, rep(v[i], 28))
    }else {
      result <- c(result, rep(v[i], 30))
    }
  }
  
  # 输出结果
  return(result)
}

#针对每年累计通量输出结果求多年平均,并最后输出为data.frame
ayearto_myMean <- function(data_list){
  # 计算每个数据框每列的平均值
  column_means <- lapply(data_list, function(x) {
    # 计算每列的平均值
    colMeans(x, na.rm = TRUE)
  })
  
  # 将所有数据框的列平均值合并为一个数据框
  result_df <- do.call(rbind, column_means)
  result_df[result_df==0] = NA #方便后续分位数计算
  # 显示结果
  # print(result_df)
  return(result_df)
}

# 定义要计算的百分位数
percfun <- function(df){
  percentiles = c(0.1, 0.5, 0.9) 
  # percentiles = c(0.025, 0.5, 0.975) 
  # 计算每行数据的指定百分位数
  row_percentiles <- t(apply(df, 1, function(x) quantile(x, percentiles , na.rm=T))) #1为对行操作
  # 转换结果为数据框并添加行名
  row_percentiles_df <- as.data.frame(row_percentiles)
  # print(row_percentiles_df)
  
  return(row_percentiles_df)
}

#求某个数据库多列分位数
#@colby 
#ncolp
perfunPars <- function(df1,colby,ncolp){
  out = NULL
  colnm = c() #列名
  for(stcol in 1:ncolp){
    # print(stcol)
    tpool = percfun(df1[,seq(stcol,ncol(df1),colby),drop = FALSE]) #
    
    if(is.null(out)){
      out = tpool
    }else{
      out = cbind(out,tpool)
    }
    colnm = c(colnm,paste0(colnames(df1)[stcol],"_",colnames(tpool)))
  }
  # print(head(out))
  print(colnames(df1)[1:stcol])
  print(paste("tpool",colnames(tpool)))
  colnames(out) = colnm
  return(out)
}

