#区分MA和MD的生长呼吸速率和维持速率版本2
#标记碳
#率定参数
#
consfun <- function(x){
  names(x) = names(parToFit.lower)
  constraint =c()
 
  #  DOC8 < POC2 < MAOM5
  if(pars_type=="alpha" | pars_type=="alpha_vmd"){
    constraint[1] <- x["kaff_ml"] - x["kaff_lb"]
    constraint[2] <- x["kaff_pl"] - x["kaff_lb"]
  }

  constraint[length(constraint)+1] <- x["rate_bd"] - x["rate_Kbd"]
  # 新增约束：param_pb + param_pbl < 1
 constraint[length(constraint) + 1] <- 0.99 - (x["param_pb"] + x["param_pbl"])
  
  return(matrix(constraint, nrow = length(constraint), ncol = 1))
}

#设置需要率定的参数
#吸附DOC用参数p2，增加（1-param_fpl）为POM分解到MAOM的比例，beta=0.25,
if(pars_type=="alpha"){
  parToFit.pools <- c(
    
    kaff_pl = pars.in$kaff_pl,       #2
    # Vpl0 = 1.5,
    alpha_pl = pars.in$alpha_pl,     #3
    # eact_pl = pars.in$eact_pl,       #4
    
    kaff_ml = pars.in$kaff_ml,       #5
    # Vml0 = 1,
    alpha_ml = pars.in$alpha_ml,     #6
    # eact_ml = pars.in$eact_ml,       #7
    
    kaff_lb = pars.in$kaff_lb,       #8
    # Vlb0 = 3.5,
    alpha_lb = pars.in$alpha_lb,     #9
    # eact_lb = pars.in$eact_lb,       #10
    
    param_fpl = 0.5,
    
    # param_em = pars.in$param_em,     #15
    #      param_ne = pars.in$param_ne,      #16
    #       param_p1 = pars.in$param_p1,     #19
    param_p2 = pars.in$param_p2,     #20
    
    #      rate_leach = pars.in$rate_leach, #17
    kaff_des = pars.in$kaff_des,     #18
    
    cue_t = pars.in$cue_t          #17
    # cue_th = pars.in$cue_t
    # tae_ref = pars.in$tae_ref,       #18
    # matpot = pars.in$matpot,         #19
    # lambda = pars.in$lambda,         #20
    # kamin = pars.in$kamin,           #21
  ) 
  
  parToFit.lower <- c(
    kaff_pl = 5e3, #5e3,  #
    # Vpl0 = 0.1,
    alpha_pl = 1e11, #1e10, 
    # eact_pl = 63e3,
    
    kaff_ml = 5e3, #5e3,  #
    # Vml0 = 0.01,
    alpha_ml = 1e10,  #1e10,
    # eact_ml = 65e3,
    
    kaff_lb = 100, #10e3, #
    # Vlb0 = 0.1,
    alpha_lb = 5e10,  #6e9, 
    # eact_lb = 55e3,
    
    param_fpl = 0.1, #POM?纸?直?咏???DOC????
    # param_em = 0.5,
    #param_ne = 0,
    # param_p1 = 0.01,#0, #
    param_p2 = 0.1, #????通量???0.5-1.5???????2
    kaff_des = 0.02, #0 #
    #                        rate_leach = 0,
    # param_p2 = 0
    cue_t = 0.001#-Inf,
    # cue_th = 0.001
    # tae_ref = -Inf,
    # matpot = 0,
    # lambda = 0,
    # kamin = 0,
  )
  parToFit.upper <- c(
    kaff_pl = 5e4,  #1e5, 5e4,  #
    # Vpl0 = 5,
    alpha_pl = 5e12, #5e10,
    # eact_pl = 65e3,
    
    kaff_ml = 7e4,
    # Vml0 = 5,
    alpha_ml = 5e12, #5e12,
    # eact_ml = 70e3,
    
    kaff_lb = 1e3, #5e4,  #
    # Vlb0 = 150,
    alpha_lb =  1e12, #9e10,
    # eact_lb = 60e3,
    
    param_fpl = 0.9,
    # param_em = 0.99,
    #                        param_ne = 0.1,
    # param_p1 = 1.5, #0.5,#
    param_p2 = 5,
    kaff_des = 2, #2
    #                        rate_leach = 0.002
    #param_p2 = 0.5 
    cue_t = 0.016 #Inf,
    # cue_th = 0.017
    # tae_ref = Inf,
    # matpot = Inf,
    # lambda = Inf,
    # kamin = Inf,
  )
  
  parToFit.pools <- c(parToFit.pools,
                      rate_bd = 0.003,       #11
                      rate_Kbd = 3e-4,
                      cue_ref = 0.7,       #12
                      param_pb = 0.5     #13
  )
  parToFit.lower <- c(parToFit.lower,rate_bd = 0.00001, rate_Kbd = 1e-4,
                      cue_ref = 0.1
                      ,param_pb = 0.01)
  parToFit.upper <- c(parToFit.upper,rate_bd = 0.005, rate_Kbd = 5e-4,  
                      cue_ref = 0.8, #0.6,???????2???0.8???统?????????值为0.9???????????0.6
                      param_pb = 0.9)
  
  
  parToFit.pools <- c(parToFit.pools, param_pc = pars.in$param_pc     #14
                      #      param_ne = parameters[[1]]$param_ne,      #16
                      #      rate_leach = parameters[[1]]$rate_leach, #17
  ) 
  parToFit.lower <- c(parToFit.lower,param_pc = 0.3
                      #                         ,param_ne = 0,rate_leach = 0,
  )
  parToFit.upper <- c(parToFit.upper,param_pc = 0.99
                      #,param_ne = 0.1,rate_leach = 0.002
  )
  
  parToFit.pools <- c(parToFit.pools, r0 = 0.1,Ma = 0.05,beta=0.005
                      #,wdorm = 4,SWP_A2D=0.4, tau=0.25
                      # , gamma=0.01#,Vlb_mod = 1
                      ,acue = 0.005
  )
  parToFit.lower <- c(parToFit.lower, r0 = 0.01,Ma = 0.0001,beta=5.e-4#1.0e-7
                      #,wdorm = 1,SWP_A2D=0.01, tau=0.1
                      # , gamma=0.0001#,Vlb_mod = 0.01
                      ,acue = 0.00001
  )
  parToFit.upper <- c(parToFit.upper, r0 = 1,Ma = 0.5,beta=0.25#5e-4
                      #,wdorm = 6,SWP_A2D=6, tau=0.95
                      # , gamma=5#,Vlb_mod = 20
                      ,acue = 0.5
  )
  
  #???????矢???Vmax_lb????
  # rmbd = which(names(parToFit.pools)=="rate_bd")
  # parToFit.pools = parToFit.pools[-rmbd]
  # parToFit.lower = parToFit.lower[-rmbd]
  # parToFit.upper = parToFit.upper[-rmbd]
  
}

#设置需要率定的参数
#吸附DOC用参数p2，增加（1-param_fpl）为POM分解到MAOM的比例，beta=0.25,
if(pars_type=="Vmax_tref"){
#  if(testMd == "_parsMax" | testMd == "_parsMaxbeta"){
    parid = 2
#  }else{
#    parid = 1 #
#  }

  parToFit.pools <- c(
    
    k_lit = pars.in$k_lit,    #add_litter_pool
    f_resp = pars.in$f_resp,   #add_litter_pool
    
    
    kaff_pl = pars.in$kaff_pl,       #2
    Vpl0 = 10,
    
    kaff_ml = pars.in$kaff_ml,       #5
    Vml0 = 0.4,
    
    kaff_lb = pars.in$kaff_lb,       #8
    Vlb0 = 0.35,
    
    param_fpl = 0.5,
    param_p2 = pars.in$param_p2,     #20
    # rate_leach = pars.in$rate_leach, #17
    kaff_des = pars.in$kaff_des,     #18
    cue_t = pars.in$cue_t,          #17
    
    rate_bd = 0.003,       #11
    rate_Kbd = 3e-4,
    cue_ref = 0.7,       #12
    param_pb = 0.5,     #13
    
    param_pc = pars.in$param_pc,
    r0 = 0.1,Ma = 0.05,beta=0.005,acue = 0.005
    # ,Ama_lm0=100,Mma_oa = 2,Kma_oa=1
#    ,param_pbl = 0.5  
#  ,tauda = 0.2
  ) 
  
  parToFit.lower <- c(
    
    k_lit = 0.01,    #add_litter_pool
    f_resp = 0.1,   #add_litter_pool
    
    
    
    kaff_pl = 5e3, 
    Vpl0 = c(0.1, 0.001)[parid],
    
    kaff_ml = 5e3, 
    Vml0 = c(0.005, 0.0001)[parid],
    
    kaff_lb = 10, #10e3, #
    Vlb0 = c(0.01, 0.0001)[parid],
    
    param_fpl = c(0.1,0.8)[2], #POM?纸?直?咏???DOC????
    param_p2 = 0.1, #????pH?????缀?力Klm?木???系??
    kaff_des = 0.02, #0 #
    # rate_leach = 0,
    cue_t = 0.001,
    
    rate_bd = 0.00001, 
    rate_Kbd = 1e-4,
    cue_ref = 0.1,
    param_pb = 0.01,
    
    param_pc = 0.3,
    r0 = 0.01,Ma = 0.0001,beta=5.e-4,acue = 1e-5
    # ,Ama_lm0=1,Mma_oa = 0, Kma_oa=-10
#    ,param_pbl = 0.01
#    ,tauda = 0.1
  )
  parToFit.upper <- c(
    
    k_lit = 0.1,    #add_litter_pool
    f_resp = 0.4,   #add_litter_pool
    
    
    kaff_pl = 20e3,  #1e5, 5e4,  #
    Vpl0 = c(30, 10, 5, 1,3)[5],#20,
    
    kaff_ml = c(70e3,20e3)[2],
    Vml0 = c(10, 1, 0.5, 0.1)[2],#1,
    
    kaff_lb = 1e3, 
    Vlb0 = c(50, 0.5, 0.2)[3],#0.5,10,2
    
    param_fpl = c(0.9,0.95)[2],
    param_p2 = 5,
    kaff_des = 2, #2
    # rate_leach = 0.002
    cue_t = 0.016,
    
    rate_bd = 0.1,#0.005, 
    rate_Kbd = 5e-4,  
    cue_ref = 0.8, #0.6,???????2???0.8???统?????????值为0.9???????????0.6
    param_pb = 0.9,
    
    param_pc = c(0.99,0.9)[1],
    r0 = 1,Ma = 0.5,beta=1,#0.25,
    acue = 0.5
    # ,Ama_lm0=300,Mma_oa = 10, Kma_oa = 10
#    ,param_pbl = 0.9 
#    ,tauda = 0.95
  )

  if(testMd == "_parsMaxbeta"  | testMD2){
    parToFit.upper["beta"] = 0.01
  }

  if(testMd == "_fixbeta"){
    rmpars = -which(names(parToFit.pools) %in% "beta")
    parToFit.pools = parToFit.pools[rmpars]
    parToFit.lower = parToFit.lower[rmpars]
    parToFit.upper = parToFit.upper[rmpars]
  }
  
#  if(testMd == "_rmPOCtoMAOC"  | testMD2){
#    rmpars = -which(names(parToFit.pools) %in% "param_fpl")
#    parToFit.pools = parToFit.pools[rmpars]
#    parToFit.lower = parToFit.lower[rmpars]
#    parToFit.upper = parToFit.upper[rmpars]
#  }
  
#  if(testMd == "_rmMAOCmm"  | testMD2  | testMd == "_fixbeta_rmMAOCmm"){
#    rmpars = -which(names(parToFit.pools) %in% c("kaff_ml","Vml0"))
#    parToFit.pools = parToFit.pools[rmpars]
#    parToFit.lower = parToFit.lower[rmpars]
#    parToFit.upper = parToFit.upper[rmpars]
#  }

  #固定VP半饱和常数
  if(F){
    rmpars = -which(names(parToFit.pools) %in% c("kaff_ml","kaff_pl"))
    parToFit.pools = parToFit.pools[rmpars]
    parToFit.lower = parToFit.lower[rmpars]
    parToFit.upper = parToFit.upper[rmpars]
  }  
  
}

if(pars_type=="alpha_vmd"){
  parToFit.lower <- c(
    kaff_pl = 100,#5e3, 
    # Vpl0 = 0.1,
    alpha_pl = 1e11,
    
    kaff_ml = 100,#5e3, 
    # Vml0 = 0.005,
    alpha_ml = 1e10, 
    
    kaff_lb = 1, #100, #10e3, #
    # Vlb0 = 0.0001,
    alpha_lb = 1e7,#5e10,
    
    param_fpl = 0.1, #POM?纸?直?咏???DOC????
    param_p2 = 0.1, #????pH?????缀?力Klm?木???系??
    kaff_des = 0.02, #0 #
    # rate_leach = 0,
    cue_t = 0.001,
    
    rate_bd = 0.00001, 
    rate_Kbd = 1e-4,
    cue_ref = 0.1,
    param_pb = 0.01,
    
    param_pc = 0.3,
    r0 = 0.01,Ma = 0.0001,beta=5.e-4,acue = 0.00001
    # ,Ama_lm0=1,Mma_oa = 0, Kma_oa=-10
  )
  parToFit.upper <- c(
    kaff_pl = 5e4,  #1e5, 5e4,  #
    # Vpl0 = 5,
    alpha_pl = 5e12, #5e10,
    
    kaff_ml = 7e4,
    # Vml0 = 5,
    alpha_ml = 5e12, 
    
    kaff_lb = 1e3, #5e4,  #
    # Vlb0 = 150,
    alpha_lb =  1e10,#1e12,
    
    param_fpl = 0.9,
    param_p2 = 5,
    kaff_des = 2, #2
    # rate_leach = 0.002
    cue_t = 0.016,
    
    rate_bd = 0.1,#0.005, 
    rate_Kbd = 5e-4,  
    cue_ref = 0.8, #0.6,???????2???0.8???统?????????值为0.9???????????0.6
    param_pb = 0.9,
    
    param_pc = 0.99,
    r0 = 1,Ma = 0.5,beta=1,#0.25,
    acue = 0.5
    # ,Ama_lm0=300,Mma_oa = 10, Kma_oa = 10
    
  )
  pars = pars[!names(pars)%in%c("Vpl0","Vml0","Vlb0")]
  
}

if(ss_cond_obj){
  parToFit.pools["param_pi"] = pars.in$param_pi
  parToFit.lower["param_pi"] = c(0.1,0.8)[2]
  parToFit.upper["param_pi"] = c(0.9,0.95)[2]
}

