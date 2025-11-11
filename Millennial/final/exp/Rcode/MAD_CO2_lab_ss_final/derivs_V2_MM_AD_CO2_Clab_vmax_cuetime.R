#区分MA和MD的生长呼吸速率和维持速率版本2
#标记碳

subVGparameter <- function(SoilTexture) {
  # !! Soil Water Characteristic Curve
  # !! van-Genuchten parameters
  # !! Tuller M, Or D (2004) Retention of water in soil and the soil water characteristic curve. 
  # !! Encyclopedia of soils in the environment, 4, 278-289.
  # 初始化参???
  SWCres <- 0
  SWCsat <- 0
  alpha <- 0
  n <- 0
  
  # 根据SoilTexture选择参数
  switch(SoilTexture,
         "Sand"            = {SWCres <- 0.058; SWCsat <- 0.37; alpha <- 0.035; n <- 3.19},
         "Loamy-Sand"      = {SWCres <- 0.074; SWCsat <- 0.39; alpha <- 0.035; n <- 2.39},
         "Sandy-Loam"      = {SWCres <- 0.067; SWCsat <- 0.37; alpha <- 0.021; n <- 1.61},
         "Loam"            = {SWCres <- 0.083; SWCsat <- 0.46; alpha <- 0.025; n <- 1.31},
         "Silt"            = {SWCres <- 0.123; SWCsat <- 0.48; alpha <- 0.006; n <- 1.53},
         "Silt-Loam"       = {SWCres <- 0.061; SWCsat <- 0.43; alpha <- 0.012; n <- 1.39},
         "Sandy-Clay-Loam" = {SWCres <- 0.086; SWCsat <- 0.40; alpha <- 0.033; n <- 1.49},
         "Clay-Loam"       = {SWCres <- 0.129; SWCsat <- 0.47; alpha <- 0.030; n <- 1.37},
         "Silty-Clay-Loam" = {SWCres <- 0.098; SWCsat <- 0.55; alpha <- 0.027; n <- 1.41},
         "Silty-Clay"      = {SWCres <- 0.163; SWCsat <- 0.47; alpha <- 0.023; n <- 1.39},
         "Clay"            = {SWCres <- 0.102; SWCsat <- 0.51; alpha <- 0.021; n <- 1.20},
         "DEFAULT"         = {SWCres <- 0.095; SWCsat <- 0.45; alpha <- 0.024; n <- 1.66}
  )
  
  return(c(SWCres = SWCres, SWCsat = SWCsat, alpha = alpha, n = n))
}
# 调用函数并传入SoilTexture参数
fSWC2SWP <- function(SWC0,SoilTexture=NULL,SWCsat){
  # van-Genuchten equation, SWP in units of [cm], SWC if fraction (0-1)
  # this function converts SWP [cm] to SWP [MPa]
  # return fSWC2SWP), actually fSWC2SWP<0
  # convert SWC(0-1) to SWP(MPa)
  # SWC needs be converted to fraction first
  #
  # SWPmin = -13.86e0  !!for Microbial Mortality
  const_cm2MPa = 98e-6           #1cm water column
  rlim = 1.01
  
  #求VG参数
  # VGpars = subVGparameter(SoilTexture)
  # SWCres = VGpars["SWCres"]
  # SWCsat = SWCsat #VGpars["SWCsat"]
  # alpha = VGpars["alpha"]
  # n = VGpars["n"]
  
  SWCres = 0.067
  SWCsat = SWCsat #VGpars["SWCsat"] #0.37
  alpha = 0.021
  n = 1.61
  
  
  m = 1-1/n
  #SWC0 = 0.5
  if (SWC0 <= SWCres*rlim){
    #fSWC2SWP = SWPmin
    SWC = SWCres*rlim
  }else{
    SWC = SWC0
  }
  
  if (SWC < SWCsat){
    eff_sat = (SWC - SWCres)/(SWCsat - SWCres) #effective saturation
    fSWC2SWP = (1/(eff_sat**(1/m)) - 1)**(1/n)/alpha
    fSWC2SWP = -1*fSWC2SWP*const_cm2MPa
  }else{
    fSWC2SWP = 0
  }
  
  return(unname(fSWC2SWP))
}

#返回更新后的碳库，C通量
CalLabC <- function(state_c13,tmp.soilC,state0,flux0,c13lab,CUE,param_pb,Ama_lm=1){
  with(as.list(c(state_c13,tmp.soilC,state0,flux0)),{
    
    # cat("state_c13=",state_c13,'\n')
    # cat("tmp.soilC=",tmp.soilC,'\n')
    # cat("state0=",state0,'\n')
    #初始标记碳库
    # POM_c13 = ifelse(clab == Cin, f_Lit_PO)
    # LMWC_c13 = switch(clab, Cin = f_Lit_LM + f_Exud_LM, Cinoa = f_Oa_LM, Cad = acid_des) #外源C输入初始化标记碳库
    # MAOM_c13 = ifelse(clab == Cinoa, f_Oa_MA)
    
    #标记C周转
    if(POM0 >0 ){
      f_PO_LM_c13 = f_PO_LM * POM_c13/POM0
      f_PO_MA_c13 = f_PO_MA * POM_c13/POM0
    }else{
      f_PO_LM_c13 = 0 
      f_PO_MA_c13 = 0
    }
    
    # Ama_lm = 1
    #有机酸输入到MAOM和LMWC碳库
    if(MAOM0 >0 ){
      acid_des_c13 = acid_des * MAOM_c13/MAOM0 * Ama_lm #此MAOM为经过解吸附后的MAOM，
      acid_des_c13 = min(acid_des_c13, MAOM_c13)
    }else{
      acid_des_c13 = 0 
    }
    MAOM_c13 = MAOM_c13 - acid_des_c13
    
    if(c13lab == "Cad"){
      LMWC_c13 = LMWC_c13  + acid_des
    }else{
      LMWC_c13 = LMWC_c13  + acid_des_c13
    }
    
    if(c13lab == "Cinoa"){
      LMWC_c13 = LMWC_c13 + f_Oa_LM 
      MAOM_c13 = MAOM_c13 + f_Oa_MA
    }
    
    ###MAOM的酶促分解反应MAOM->LMWC
    if(MAOM > 0 ){
      f_MA_LMe_c13  = min(f_MA_LMe * MAOM_c13/MAOM  * Ama_lm,MAOM_c13)
      f_MA_LMm_c13 =  min(f_MA_LMm  * MAOM_c13/MAOM * Ama_lm, MAOM_c13)
    }else{
      f_MA_LMe_c13 = 0 
      f_MA_LMm_c13 = 0 
    }

    #MA <-> MD
    if(MA0 > 0 ){
      f_MA_MD_c13 = min(f_MA_MD * MA_c13/MA0,MA_c13) #MA是更新后的，不是初始输入的state！！！MA0是初始
    }else{
      f_MA_MD_c13 = 0 
    }
    
    if(MD0 > 0 ){
      f_MD_MA_c13 = min(f_MD_MA * MD_c13/MD0,MD_c13)
    }else{
      f_MD_MA_c13 = 0 
    }
    
    MA_c13 = MA_c13 + f_MD_MA_c13 - f_MA_MD_c13 #MA和MD立刻更新，因为有机酸加入之后反应非常快！
    MD_c13 = MD_c13 + f_MA_MD_c13 - f_MD_MA_c13
    
    #Equation 13
    # LMWC -> MIC
    # 微生物优先利用草酸Alma_mb = ifelse(,)
    Woa_lm = ifelse(c13lab == "Cinoa",2,1)
    if(LMWC0 > 0 & LMWC_c13>0){
      # f_LM_MB_c13 = min(f_LM_MB * LMWC_c13/(LMWC0 + f_Oa_LM + acid_des),LMWC_c13)
      # f_MA_growth_c13 = f_MA_growth * LMWC_c13/(LMWC0 + f_Oa_LM + acid_des)
      # f_MA_co2_c13 = f_MA_co2  * LMWC_c13/(LMWC0 + f_Oa_LM + acid_des)
      # f_MA_growth_co2_c13 = f_MA_growth_co2 * LMWC_c13/(LMWC0 + f_Oa_LM + acid_des)
      # f_MA_maint_co2_c13 = f_MA_maint_co2 * LMWC_c13/(LMWC0 + f_Oa_LM + acid_des)
      
      WTlm = (LMWC0 + f_Oa_LM + acid_des) + LMWC_c13 * (Woa_lm-1)
      f_LM_MB_c13 = min(f_LM_MB * LMWC_c13*Woa_lm/WTlm, LMWC_c13) #-0.001*LMWC_c13

      f_MA_growth_c13 = f_LM_MB_c13 * CUE
      f_MA_co2_c13 = f_LM_MB_c13 * (1-CUE)
      f_MA_growth_co2_c13 = 0
      f_MA_maint_co2_c13 = 0
      
    }else{
      f_LM_MB_c13 = 0 
      f_MA_growth_c13  = 0 
      f_MA_co2_c13 = 0 
      f_MA_growth_co2_c13  = 0
      f_MA_maint_co2_c13 = 0
    }
    # cat("LMWC_c13=",LMWC_c13,"LMWC0=",LMWC0,"更想LMWC=",(LMWC0 + f_Oa_LM + acid_des),"f_LM_MB_c13",f_LM_MB_c13,"f_LM_MB=",f_LM_MB,'\n')
    LMWC_c13 = LMWC_c13 - f_LM_MB_c13 #立刻更新，用于吸附-解吸附计算
    # MA_c13 = MA_c13 + f_MA_growth_c13
    
    if(MD0 > 0 ){
      f_MD_co2_c13 = min(f_MD_co2  * MD_c13/MD,MD_c13)
    }else{
      f_MD_co2_c13 = 0 
    }
    
    if(LMWC > 0 & LMWC_c13 >0){
      #Equation 8
      # LMWC -> out of system leaching
      f_LM_leach_c13 = min(f_LM_leach  * LMWC_c13/LMWC,LMWC_c13)
      
      #Equation 9
      # LMWC -> MAOM
      f_LM_MA_c13 =  min(f_LM_MA  * LMWC_c13/LMWC,LMWC_c13)
      # cat("f_LM_leach = ",f_LM_leach,"f_LM_leach_c13=",f_LM_leach_c13,"LMWC_c13=",
      #     LMWC_c13,"LMWC=",LMWC,'\n')
    }else{
      f_LM_leach_c13 = 0 
      f_LM_MA_c13 = 0 
    }
    
    
    #Equation 16
    # MB -> MAOM/LMWC #微生物死亡周转
    if(MA > 0 ){
      f_MB_turn_c13 = min(f_MB_turn * MA_c13/MA, MA_c13) #这是总周转，但是进入LMWC和MAOM的呢？
      f_MB_LM_c13 = f_MB_turn_c13* (1. - param_pb)
      f_MB_MA_c13 = f_MB_turn_c13 * param_pb
    }else{
      f_MB_turn_c13 = 0 
      f_MB_LM_c13 = 0 
      f_MB_MA_c13 = 0 
    }
    
    #碳库损失通量之和<碳库
    if((f_MA_LMe_c13 + f_MA_LMm_c13) > MAOM_c13){
      f_MA_LMe_c13 = f_MA_LMm_c13 = MAOM_c13/2
    }
    f_MA_LM_c13 = acid_des_c13 + f_MA_LMe_c13 + f_MA_LMm_c13
    
    if(length(which((c(f_LM_leach_c13,f_LM_MA_c13,LMWC_c13)<0)))>0){
      print("LMWC碳库周转NA或0")
      cat("f_LM_MA=",f_LM_MA,"LMWC=",LMWC,'\n')
      cat("f_LM_leach_c13=",f_LM_leach_c13,"f_LM_MA_c13=",f_LM_MA_c13,"LMWC_c13=",LMWC_c13,'\n')
    }
    if((f_LM_leach_c13 + f_LM_MA_c13)>LMWC_c13){
      f_LM_leach_c13 = f_LM_MA_c13 = LMWC_c13/2
    }
    # if((f_MD_MA_c13 + f_MD_co2_c13)>MD_c13){
    #   f_MD_MA_c13 = f_MD_co2_c13 = MD_c13/2
    # }
    #汇总CO2和解吸附
    Tco2_c13 = f_MA_co2_c13 + f_MD_co2_c13
    
    #Equation 1
    dPOM_c13 = f_Lit_PO - (f_PO_LM_c13 + f_PO_MA_c13)
    
    #Equation 7
    dLMWC_c13 = f_Lit_LM + f_Exud_LM + f_PO_LM_c13 + f_MB_LM_c13 + (f_MA_LMe_c13 + f_MA_LMm_c13) - f_LM_leach_c13  - f_LM_MA_c13 
    
    #Equation 20
    dMA_c13 = f_MA_growth_c13 - f_MB_turn_c13 # + (f_MD_MA_c13 - f_MA_MD_c13)
    dMD_c13 = - f_MD_co2_c13 #+ f_MA_MD_c13 - f_MD_MA_c13 
    # dMIC_c13 = dMA_c13 + dMD_c13
    
    dMAOM_c13 = f_LM_MA_c13 - (f_MA_LMe_c13 + f_MA_LMm_c13) + f_MB_MA_c13 + f_PO_MA_c13 
    
    # Update state variables
    POM_c13 = POM_c13 + dPOM_c13
    LMWC_c13 = LMWC_c13 + dLMWC_c13
    # cat("LMWC_c13=",LMWC_c13,"dLMWC_c13=",dLMWC_c13,'\n')
    MA_c13 = MA_c13 + dMA_c13
    MD_c13 = MD_c13 + dMD_c13
    # MIC_c13 = MIC_c13 + dMIC
    MIC_c13 = MA_c13 + MD_c13
    MAOM_c13 = MAOM_c13 + dMAOM_c13
    
    out.Cp = c(POM_c13, LMWC_c13, MIC_c13, MAOM_c13, MA_c13, MD_c13)
    names(out.Cp) = names(state_c13)
    
    if(clab.outvar == "Cpools_Tco2"){
      out.vect =  c(SOM_c13=sum(out.Cp[1:4]),out.Cp, Tco2_c13=Tco2_c13) #只返回标记的各碳库和总CO2
    }else if(clab.outvar == "all"){
      
      flux = c(f_PO_LM_c13, f_PO_MA_c13, f_MA_LMm_c13, f_MA_LMe_c13, acid_des_c13, f_MA_LM_c13, f_LM_MA_c13, 
               f_MA_MD_c13, f_MD_MA_c13, f_LM_MB_c13,  f_MA_maint_co2_c13, f_MA_growth_co2_c13, f_MA_co2_c13, 
               f_MD_co2_c13,f_MA_growth_c13, f_MB_LM_c13, f_MB_MA_c13, f_MB_turn_c13,  f_LM_leach_c13 
      )
      names(flux) = paste0(names(flux0)[-(1:5)],"_c13")
      # cat("Tco2_c13=",Tco2_c13,'\n')
      out.vect = c(SOM_c13=sum(out.Cp[1:4]),out.Cp, Tco2_c13=Tco2_c13, flux)
    }
    # cat("CalLabC输出的C13",out.vect,'\n')
    return(out.vect)
  })
  
}

derivs_V2_MM_AD_CO2_Clab <- function(step.num,state,parameters,forc_st,forc_sw,
                                     forc_lit,forc_exud,forc_acid=NULL,ss_cond=FALSE) {
  with(as.list(c(state,parameters)), {
    # cat("MD=",MD,'\n')
    # step.num = 1
    # state = soilCin
    # parameters = parsin
    # list2env(as.list(c(state,parameters)), envir = .GlobalEnv)
    
    # cat("最开始MA=",MA,"MD=",MD)
    # 初始化激活态MA和休眠态MD
    # if(step.num == 1){
    #   state[,"MA"] = r0 * MIC
    #   state[,"MD"] = (1. - r0) * MIC
    #   cat("day=1的state MA和MD",state[,c("MA","MD")], "MA=",MA,"MD=",MD,'\n')
    # }
    
    # Soil type properties  
    #Equation 10
    kaff_lm = exp(-param_p1 * param_pH - param_p2) * kaff_des
    
    #Equation 11
    param_qmax = depth*param_bulkd * param_claysilt * param_pc
    #约束率定的pc
    # print(paste("param_qmax=",param_qmax,"MAOM=",MAOM,"step.num=",step.num))
    if(param_qmax < MAOM){
      qmax_cue_error <<- TRUE
      print("param_qmax exit")
      return(rep(NA, length(state)))
    }
    
    # Hydrological properties
    forc_sw0 = min(forc_sw(step.num),porosity)
    SWP = fSWC2SWP(forc_sw0,SWCsat=porosity)  #MPa,SWP<0,饱和基质势能=0
    matpot = abs(SWP * 1e3)        #kPa
    
    #Equation 4
    scalar_wd = (forc_sw0 / porosity)^0.5
    
    #Equation 15
    scalar_wb = exp(lambda * -matpot) * (kamin + (1 - kamin) * ((porosity - forc_sw0) / porosity)^0.5) * scalar_wd
    
    # Decomposition
    gas_const <- 8.31446
    
    T0 = forc_st(step.num) + 273.15
    Tref = tae_ref + 273.15
    tp_scalar <- function(sCase,Ea=NULL){
      if(is.null(Ea)){
        switch(sCase,
               "MR" = {Ea <- 20e3},  #J/mol
               "LIG"= {Ea <- 53.0e3},
               "CEL"= {Ea <- 36.3e3},
#               "POM"= {Ea <- 63909},
               "POM"= {Ea <- 66e3},
               "LMWC"= {Ea <- 57865},
               "MAOM"= {Ea <- 67e3},# {Ea <- 47e3},
               "MD"= {Ea <- 47e3},
               "KM"= {Ea <- 30e3},
               "LIT"= {Ea <- 62e3}
        ) 
        # print(Ea)
      }
      tp_scalar0 = exp(-Ea/gas_const * (1/T0-1/Tref))
      return(tp_scalar0)
    }
    
    #外源C输入
    
    f_Lit_LIT = forc_lit(step.num)
    
    # cat("step.num =", step.num, "forc_lit(step.num) =", f_Lit_LIT, "\n")
    
    k = k_lit * tp_scalar("LIT") * scalar_wd
    
    if (LIT > 0) {
      #f_LIT = max(LIT * vmax_lit * kaff_lit,0)
      f_LIT = max(LIT * k,0)
      f_LIT = min(LIT, f_LIT)
      
    }else{
      f_LIT = 0
    }
    
    f_LIT_co2 = f_LIT * f_resp;
    f_LIT_SOC = f_LIT * (1. - f_resp);
    
    # # # 在LIT分解计算后添加
    # cat("LIT分解调试 - step:", step.num,
    #     # " LIT=", LIT,
    #     " k_lit=", k_lit,
    #     # " f_LIT_co2=", f_LIT_co2,
    #     # " f_LIT_SOC=", f_LIT_SOC,
    #     " f_resp=", f_resp,
    #     # " tp_scalar=", tp_scalar("LIT"),
    #     # " scalar_wd=", scalar_wd,
    #     " k=", k,
    #     # " f_LIT=", f_LIT,
    #     '\n')
    
    
    f_Lit_PO = f_LIT_SOC * param_pi
    f_Lit_LM = f_LIT_SOC * (1. - param_pi)
    f_Exud_LM = forc_exud(step.num)
    f_Oa_MA = forc_acid(step.num) * param_em
    f_Oa_LM = forc_acid(step.num) * (1. - param_em)
    
    ### Equation 2
    # POM -> LMWC
    # print(paste("POM= ",POM,"MIC= ",MIC))
    # vmax_pl = alpha_pl * exp(-eact_pl / (gas_const * (forc_st(step.num) + 273.15)))
    vmax_pl = Vpl0 * tp_scalar("POM") * scalar_wd
#    vmax_pl = Vpl0 * tp_scalar("POM") * (forc_sw0 / porosity)
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
      acid_des = ((MAOM/unit_kgTom2 *1e-3)^0.04 *(forc_acid(step.num)/unit_kgTom2)^0.87
                  * param_clay^(-0.24) * param_pH^(-2.15) * 34)* unit_kgTom2
      acid_des = min(max(acid_des,0),MAOM)
    }else{
      acid_des = 0
    }
    
    LMWC = LMWC  + f_Oa_LM + acid_des
    MAOM = MAOM + f_Oa_MA - acid_des
    
    ###MAOM的酶促分解反应MAOM->LMWC
    # vmax_ml = max(alpha_ml * exp(-eact_ml / (gas_const * (forc_st(step.num) + 273.15))),0)
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
    # vmax_lb = alpha_lb * exp(-eact_lb / (gas_const * (forc_st(step.num) + 273.15)))
    Vg0 = Vlb0 #假设MEND中的的Vg就是Vd
    Vg = Vg0 * tp_scalar("LMWC")  * scalar_wb #MA的生长呼吸速率，类似MEND中的Vd
    
    #MA维持呼吸速率：α/(1-α) 占MA生长呼吸的比例。微生物激活MA和休眠MD周转
    Vm = Vg0 * Ma/(1-Ma) * tp_scalar("MR") * scalar_wb
    
    #MD维持呼吸速率: 是MA维持呼吸的β倍
    VmD = Vm * beta
    
    #Equation 22
    # microbial growth flux, but is not used in mass balance
    CUE = cue_ref - cue_t * (forc_st(step.num) - tae_ref) 
    CUE = CUE - acue * forc_acid(1)#/max(step.num/2,1)  #step.num   #加入草酸后CUE降低,维持30天
    #    print(paste("step.num=",step.num, "CUE=",CUE,"acue=",acue))
    CUE_max = 0.95 #0.9还是0.95？
    CUE_min = 0.01
    if(CUE<CUE_min | CUE>CUE_max){
      qmax_cue_error <<- TRUE
      print("CUE NA exit")
      print(paste("step.num=",step.num, "CUE=",CUE,"acue=",acue))
      return(rep(NA, length(state)))
    }
  
    #Equation 14.1
    #微生物激活MA和休眠MD周转,
     VmA2D = Vm # * tp_scalar("MR") #* scalar_wd #* wp_scalar_low
     VmD2A = Vm # * tp_scalar("MR") #* scalar_wd #* wp_scalar

#水势影响不一样 Wanggs2018
#    w = 1
#    SWPa2d = -0.4 #MPa
#    fwA2D = abs(SWP)^w/(abs(SWP)^w + abs(SWPa2d)^w)
#    VmA2D = Vm /scalar_wb * fwA2D# * tp_scalar("MR") #* scalar_wd #* wp_scalar_low

#    SWPd2a = SWPa2d * tauda #MPa
#    fwD2A = abs(SWPd2a)^w/(abs(SWP)^w + abs(SWPd2a)^w)
#    VmD2A = Vm /scalar_wb * fwD2A# * tp_scalar("MR") #* scalar_wd #* wp_scalar

    vmax_lb = (Vg + Vm) / CUE
    kaff_lb = kaff_lb * tp_scalar("KM")
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
    MA = MA + f_MD_MA - f_MA_MD #MA和MD立刻更新，因为有机酸加入之后反应非常快！
    MD = MD + f_MA_MD - f_MD_MA
    
    # LMWC -> MA
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
    f_MA_co2 = f_LM_MB * (1 - CUE) #method1
    f_MA_growth_co2 = 0#max( Vg * (1/CUE - 1) * MA * LMWC / (kaff_lb + LMWC),0)
    f_MA_maint_co2 = 0#max( Vm * (1/CUE - 1) * MA * LMWC / (kaff_lb + LMWC),0)
    # f_MA_co2 = f_MA_growth_co2 + f_MA_maint_co2
    # MA = MA + f_MA_growth
    
    LMWC = LMWC - f_LM_MB #立刻更新，用于吸附-解吸附计算
    f_MA_growth = f_LM_MB * CUE
    # cat("f_LM_MB=",f_LM_MB,"f_MA_growth_co2=",f_MA_growth_co2,"f_MA_maint_co2=",f_MA_maint_co2,
    #     "f_MA_co2=",f_MA_co2,"f_MA_growth=",f_MA_growth,"0相等=",f_LM_MB-(f_MA_growth+f_MA_co2),'\n')
    
    #MD维持呼吸
    # cat("VmD = ",VmD,'\n')
    if(MD>0){
      if(testMd == "_VMA_VMD"){
        f_MD_co2 = max(VmD * MD  * LMWC / (kaff_lb + LMWC),0) #0.001
      }else{
        f_MD_co2 = max(VmD * MD,0) #0.001
      }
#      f_MD_co2 = max(VmD * MD,0) #0.001 
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
    rate_bd = max(0,rate_bd + rate_Kbd*(forc_st(step.num)-tae_ref))
    if(MA>0){
      # f_MB_turn = rate_bd * MIC^2.0
      f_MB_turn = max(rate_bd * MA,0) #Liucq
      f_MB_turn = min(MA, f_MB_turn)
    }else{
      f_MB_turn=0
    }
#    f_MB_LM = f_MB_turn * (1. - param_pb)
#    f_MB_MA = f_MB_turn * param_pb
#MA死亡进入POC
    f_MB_LM = f_MB_turn * param_pbl
    f_MB_MA = f_MB_turn * param_pb
    f_MB_PO = f_MB_turn * (1. - param_pb - param_pbl)
    
    #碳库损失通量之和<碳库
    if((f_MA_LMe + f_MA_LMm) > MAOM){
      f_MA_LMe = f_MA_LMm = MAOM/2
    }
    f_MA_LM = acid_des + f_MA_LMe + f_MA_LMm
    
    if((f_LM_leach + f_LM_MA)>LMWC){
      f_LM_leach = f_LM_MA = LMWC/2
    }
    
    # if((f_MD_MA + f_MD_co2)>MD){
    #   f_MD_MA = f_MD_co2 = MD/2
    # }
    
    #所有通量
    Cinput = c(f_Lit_PO, f_Lit_LM, f_Exud_LM, f_Oa_MA, f_Oa_LM)
    names(Cinput) = allflux.name[1:5]
    
    flux = c(f_PO_LM, f_PO_MA, f_MA_LMm, f_MA_LMe, acid_des, f_MA_LM, f_LM_MA, 
             f_MA_MD, f_MD_MA, f_LM_MB,  f_MA_maint_co2, f_MA_growth_co2, f_MA_co2, 
             f_MD_co2,f_MA_growth, f_MB_LM, f_MB_MA, f_MB_turn,  f_LM_leach 
    )
    names(flux) = allflux.name[-1:-5]
    
    out.c13 = c()
    tmp.Cpools = c(LIT,POM, LMWC, MIC, MAOM, MA, MD) #立刻更新后的碳库
    names(tmp.Cpools) = fixCp[2:8]
    if(clab[1] != FALSE & ss_cond==FALSE & length(which(is.na(tmp.Cpools)))==0){ #稳态求解不需要运行标记C
      state0 = state
      names(state0) = paste0(names(state),"0") #输入的初始碳库
      
      for(labi in 1:length(clab)){
        if(step.num == 1){
          state_c13 = rep(0,6)
        }else{
          nmi = sapply(strsplit(names(state), "_"), function(y) tail(y, 1))
          state_c13 = state[which(nmi==clab[labi])[2:7]]
        }
        names(state_c13) = c("POM_c13", "LMWC_c13", "MIC_c13", "MAOM_c13", "MA_c13", "MD_c13")
        
        #有机酸MAOC解吸附比例调整
        Ama_lm = 1 #默认为1
        if(clab[labi]=="Cinoa"){
          if(step.num == 1){
            Ama_lm =  Ama_lm0
          }else if(state_c13["MAOM_c13"]>0){
            # cat("step.num=",step.num,'\n')
            if(forc_acid(step.num-1)!=0 & forc_acid(step.num)==0){ #有机酸输入从有到无
              MAOM_c130 <<- state_c13["MAOM_c13"]
              cat("Oacid input = 0, MAOM_c130 = ",MAOM_c130,'\n')
            }else if(forc_acid(step.num)!=0){
              #如果有机酸一直在输入，则默认MAOM_草酸的解吸附就是率定的Ama_oa,
              # Ama_lm =  Ama_lm0
              MAOM_c130 <<- state_c13["MAOM_c13"]
              cat("Oacid input !=0, MAOM_c130 = ",MAOM_c130,'\n')
            }
            Ama_lm =  Ama_lm0 * (state_c13["MAOM_c13"]/(MAOM_c130))^Mma_oa
          }
        }
        #Ama_lm受温度影响
        # Ama_lm = Ama_lm + Kma_oa*(forc_st(step.num) - tae_ref)
        # cat("Ama_lm = ",Ama_lm,'\n')
        if(Ama_lm * state_c13["MAOM_c13"] >= MAOM){
          c13_error <<- TRUE
          cat("Ama_lm is too big! ", "Ama_lm = ",Ama_lm,"MAOM_c13 = ",state_c13["MAOM_c13"],"MAOM = ",MAOM,'\n')
          out.c130 = rep(0,Cout.nlab/length(clab))
          
        }else{
          out.c130 = CalLabC(state_c13, tmp.soilC = tmp.Cpools, state0 = state0 ,
                             flux0 = c(Cinput,flux), c13lab = clab[labi],CUE = CUE,param_pb,Ama_lm)
          # names(out.c130) = paste0(names(out.c130),"_",clab[labi])
        }
        names(out.c130) = colnms.lab
        out.c13 = c(out.c13, out.c130) #多个标记C的输出按行拼接
        # cat("derivs输出的C13",out.c13)
      }
      
    }
    
    dLIT = f_Lit_LIT - f_LIT
    
    #Equation 1
    dPOM = f_Lit_PO - f_PO + f_MB_PO
    
    #Equation 7
    dLMWC = f_Lit_LM + f_Exud_LM + f_PO_LM + f_MB_LM + (f_MA_LMe + f_MA_LMm) - f_LM_leach  - f_LM_MA 
    
    #Equation 20
    dMA =  f_MA_growth  - f_MB_turn #+ f_MD_MA - f_MA_MD #在上面已更新 
    dMD =  - f_MD_co2  # + f_MA_MD - f_MD_MA#在上面已更新
    # dMIC = dMA + dMD
    # cat("dMA=",dMA,"dMD=",dMD ,'\n') #,"dMIC = ",dMIC
    
    #Equation 19
    # dMAOM = C_rate*param_em + f_LM_MA - f_MA_LM + f_MB_turn * param_pb
    #加草酸立马就吸附到MAOM???
    dMAOM = f_LM_MA + f_MB_MA + f_PO_MA - (f_MA_LMe + f_MA_LMm)
    
    # Update state variables
    
    LIT = LIT + dLIT
    POM = POM + dPOM
    LMWC = LMWC + dLMWC
    MA = MA + dMA
    MD = MD + dMD
    # MIC = MIC + dMIC
    MIC = MA + MD 
    MAOM = MAOM + dMAOM
    
    out.Cpools = c(LIT, POM, LMWC, MIC, MAOM, MA, MD)
    names(out.Cpools) = fixCp[2:8]
    
    if(out.allflux ==TRUE){
      out.vect = c(SOM=sum(out.Cpools[2:5]),out.Cpools,Tco2 = Tco2,Cinput,flux,out.c13)
    }else{
      out.vect = c(SOM=sum(out.Cpools[2:5]),out.Cpools, Tco2 = Tco2, out.c13)
      # names(out.vect) = c(names(out.Cpools), "Tco2")
    }
    return(out.vect)
  })
}







