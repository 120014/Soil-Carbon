#ѡ10????????ͼ
#ѡ10????????ͼ
#??????ͼ------
# ?????????İ?
library(ggplot2)
library(gridExtra)
# ʹ??ggplot????????ͼ??
intervalggplot <-function(df,expT0="20",nline=2,cust.ylab){
  colnames(df) = colnames(plotdf)[-nc13]
  p<-ggplot(df, aes(x = x)) + 
    geom_line(aes(y = ysim1.50), color = "#B9181A") +
    geom_line(aes(y = ysim1.10), color = "#a7413c", alpha = 0.5) +
    geom_line(aes(y = ysim1.90), color = "#a7413c", alpha = 0.5) + 
    geom_ribbon(aes(ymin = ysim1.10, ymax = ysim1.90), fill = "#a7413c", alpha = 0.5)+
    geom_point(aes(y = yobs1),col="#19647E")
  
  if(nline ==2){
    p <- p + 
      geom_line(aes(y = ysim0.50), color = "#eba59b") +
      geom_line(aes(y = ysim0.10), color = "#eba59b", alpha = 0.5) +
      geom_line(aes(y = ysim0.90), color = "#eba59b", alpha = 0.5) + 
      geom_ribbon(aes(ymin = ysim0.10, ymax = ysim0.90), fill = "#eba59b", alpha = 0.5)+
      geom_point(aes(y = yobs0),col="#84A3A9")
  }
  
  p <- p +
    labs(#title=paste(Simsite,expT0,"??"),
      x="day",y = sitnm[i],size=0.3)+
    # ylab(cust.ylab) + 
    theme_bw()+ #??ʾ??ɫ????
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black")
    )
  #theme_classic() #????ʾ???񣬵?ֻ??ʾx????y???߿򣬲???ʾ4???߿?
  #theme_bw() #??ʾ??ɫ????
}

#ggplot????״ͼ
ebarggplot <-function(data,var,error=F){
  #data = df
  data$x = factor(data$x,unique(dfMIC$x))
  data$group = factor(data$group,c("Obs0", "Sim0", "Obs1", "Sim1"))
  if(error==TRUE){
    ymax_plot = max(data[,c(-1,-2)])
  }else{
    ymax_plot = max(data$value)
  }
  
  p1 <- ggplot(data,  aes(x, weight = value, fill = group)) + #ʹ?? weight ????��ͳ?Ʒ????ڱ?��ֵ֮??
    #    geom_hline(yintercept = seq(10, (ymax_plot+1), (ymax_plot+1)%/%2), color = 'gray') +
    geom_bar(color = "#404040", width = .7, position = 'dodge') + #, alpha = 0.65
    
    labs( title=paste(Simsite,expT[Ti],"??"),  size=0.5,x="day" )+
    ylab(bquote(.(var)~"(gC m"^-2*")")) +
    #  scale_fill_brewer(palette = "Set3")+ #n??????ɫ
    scale_fill_manual(values = c("#84A3A950","#eba59b50","#19647ED9","#a7413999")
                      #,limits = c("Obs-without","Sim-without","Obs-adding","Sim-adding")
    )+
    scale_y_continuous(expand = c(0,0),limits = c(0, ymax_plot+1),breaks = seq(0,ymax_plot+1,(ymax_plot+1)%/%5)) +
    theme(plot.title = element_text(size=12,hjust=0.5))+ #,legend.position = "top",legend.justification = c(1,1)
    theme_classic()+
    guides(fill = guide_legend(title = NULL)) #ȥ??ͼ???е?group????????
  
  if(error==TRUE){
    p1 <- p1 + geom_errorbar(aes(ymin = ylow, ymax = yup), width = 0.25, linewidth = 0.3, position = position_dodge(0.7)) #mean + se
  }
  return(p1)
}
#???ݴ??���??
dfbargg <- function(dfperc,p.pl,errbar=c("10%","50%","90%")){
  simdrow = Obs.pools.t$day[1:2]
  initip = inputdata[which(inputdata$site==Simsite)[1],p.pl]
  df = as.data.frame(cbind(Obs0 = c(initip,Obs.pools.t[1:2,p.pl]),
                           Obs1 = c(0,Obs.pools.t[3:4,p.pl]),
                           Sim0 = c(0,dfperc[simdrow,errbar[2]]),
                           Sim1 = c(0,dfperc[simdrow+30,errbar[2]]),
                           x = c(0,simdrow)))
  df = pivot_longer(df,c("Obs0","Sim0","Obs1","Sim1"),names_to = "group",values_to = "value")
  
  #???????ޣ?10%
  errdflow = as.data.frame(cbind(Obs0 = rep(0,3),Obs1 = rep(0,3),
                                 Sim0 = c(0,dfperc[simdrow,errbar[1]]),
                                 Sim1 = c(0,dfperc[simdrow+30,errbar[1]]),
                                 x = c(0,simdrow)))
  df[,"ylow"]  = pivot_longer(errdflow,c("Obs0","Sim0","Obs1","Sim1"),names_to = "group",values_to = "ylow")[,"ylow"]
  #???????ޣ?90%
  errdfup = as.data.frame(cbind(Obs0 = rep(0,3),Obs1 = rep(0,3),
                                Sim0 = c(0,dfperc[simdrow,errbar[3]]),
                                Sim1 = c(0,dfperc[simdrow+30,errbar[3]]),
                                x = c(0,simdrow)))
  df[,"yup"]  = pivot_longer(errdfup,c("Obs0","Sim0","Obs1","Sim1"),names_to = "group",values_to = "yup")[,"yup"]
  
  #???û?ͼ????
  df$x = as.character(df$x)
  
  # if(p.pl=="MIC"){
  #   pmi[[Ti]] <<- ebarggplot(df,var = p.pl,error=TRUE)
  # }else{
  #   plm[[Ti]] <<- ebarggplot(df,var = p.pl,error=TRUE)
  # }
  return(df)
}
# ???ƺ???ͼ
boxgg <- function(df){
  #df = dfnpar[,c(1,2)]
  df[,2] = as.numeric(df[,2])
  ylab2 = colnames(df)[ncol(df)]
  
  colnames(df) = c("group","value")
  df$group = factor(df$group,levels = unique(df$group))
  
  ggplot(df, aes(x = group, y = value)) +
    geom_boxplot()+
    ylab(ylab2)+
    #theme()+
    stat_summary(fun = "mean", geom = "point", shape = 1, size = 1, color = "#404040")+
    theme_bw()+ #??ʾ??ɫ????
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1,size=7),
          axis.title.x = element_blank() #????ʾgroup
    )
}
#???֣?3??5?У?30??ֻ??ʾһ??ͼ????????????
layshow <-function(p){
  print(ggarrange(p[[1]],p[[4]],p[[7]], p[[10]],p[[13]],
                  p[[2]],p[[5]],p[[8]], p[[11]],p[[14]],
                  p[[3]],p[[6]],p[[9]], p[[12]],p[[15]],
                  #?˺?????ֱ??????ֻ??ʾһ??ͼ????#?????޷?????????,ȡ??˳????һ???Ϳ???
                  nrow=3,ncol=5,common.legend = TRUE, legend="bottom"))
  print(ggarrange(p[[16]],p[[19]],p[[22]], p[[25]],p[[28]],
                  p[[17]],p[[20]],p[[23]], p[[26]],p[[29]],
                  p[[18]],p[[21]],p[[24]], p[[27]],p[[30]],
                  #?˺?????ֱ??????ֻ??ʾһ??ͼ????#?????޷?????????,ȡ??˳????һ???Ϳ???
                  nrow=3,ncol=5,common.legend = TRUE, legend="bottom"))
  # print(ggarrange(p[[16]],p[[19]],p[[22]], #p[[25]],#p[[28]],
  #                 p[[17]],p[[20]],p[[23]], #p[[26]],#p[[29]],
  #                 p[[18]],p[[21]],p[[24]], #p[[27]],#p[[30]],
  #                 #?˺?????ֱ??????ֻ??ʾһ??ͼ????#?????޷?????????,ȡ??˳????һ???Ϳ???
  #                 nrow=3,ncol=3,common.legend = TRUE, legend="bottom"))
  
}
layshow3 <-function(p){
  print(ggarrange(p[[1]],p[[2]],p[[3]], 
                  p[[4]],p[[5]],p[[6]], 
                  p[[7]],p[[8]],p[[9]], 
                  p[[12]],p[[15]],
                  #?˺?????ֱ??????ֻ??ʾһ??ͼ????#?????޷?????????,ȡ??˳????һ???Ϳ???
                  nrow=3,ncol=5,common.legend = TRUE, legend="bottom"))
  print(ggarrange(p[[16]],p[[19]],p[[22]], p[[25]],p[[28]],
                  p[[17]],p[[20]],p[[23]], p[[26]],p[[29]],
                  p[[18]],p[[21]],p[[24]], p[[27]],p[[30]],
                  #?˺?????ֱ??????ֻ??ʾһ??ͼ????#?????޷?????????,ȡ??˳????һ???Ϳ???
                  nrow=3,ncol=5,common.legend = TRUE, legend="bottom"))
  # print(ggarrange(p[[16]],p[[19]],p[[22]], #p[[25]],#p[[28]],
  #                 p[[17]],p[[20]],p[[23]], #p[[26]],#p[[29]],
  #                 p[[18]],p[[21]],p[[24]], #p[[27]],#p[[30]],
  #                 #?˺?????ֱ??????ֻ??ʾһ??ͼ????#?????޷?????????,ȡ??˳????һ???Ϳ???
  #                 nrow=3,ncol=3,common.legend = TRUE, legend="bottom"))
  
}
# ?޸??Ѿ?????ͼ??ymax
ymaxggp <- function(p,pi,ymax){
  # ybreak =  c(seq(0,ymax+1,10),round(ymax,0)) #?̶ȼ???10????Щ̫ϸ????
  ybreak = seq(0,ymax+1,(ymax)%/%5) 
  
  p[[pi-2]] <- p[[pi-2]] + scale_y_continuous(expand = c(0,0),limits = c(0, ymax+1),breaks = ybreak)
  p[[pi-1]] <- p[[pi-1]] + scale_y_continuous(expand = c(0,0),limits = c(0, ymax+1),breaks = ybreak)
  p[[pi]] <- p[[pi]] + scale_y_continuous(expand = c(0,0),limits = c(0, ymax+1),breaks = ybreak)
  
  # p[[pi-2]] <- p[[pi-2]] + coord_cartesian(ylim = c(0, ymax))
  # p[[pi-1]] <- p[[pi-1]] + coord_cartesian(ylim = c(0, ymax))
  # p[[pi]] <- p[[pi]] + coord_cartesian(ylim = c(0, ymax))
  return(p)
}


# ????Ҫ?????İٷ?λ??------------
percfun <- function(df){
  percentiles = c(0.1, 0.5, 0.9) 
  # percentiles = c(0.025, 0.5, 0.975) 
  # ????ÿ?????ݵ?ָ???ٷ?λ??
  row_percentiles <- t(apply(df, 1, function(x) quantile(x, percentiles))) #1Ϊ???в???
  
  # ת??????Ϊ???ݿ???????????
  row_percentiles_df <- as.data.frame(row_percentiles)
  return(row_percentiles_df)
}
mpars_meansd <- function(df0, var0, outsd=TRUE,onlyMean_sd = FALSE){
  colnames(df0) = gsub("\\.\\d+$", "", colnames(df0))
  outdf = NULL
  for(vi in 1:length(var0)){
    dfvar = df0[,which(colnames(df0)==var0[vi])]
    custM = rowMeans(dfvar,na.rm = TRUE)
    custSD = apply(dfvar,1,sd,na.rm = TRUE)
    if(outsd==FALSE){
      nvarrep = 1
      out1var = data.frame(mean = custM)
    }else{
      nvarrep = 3
      out1var = as.data.frame(cbind(sd0 = custM - custSD, mean=custM , sd1 = custM + custSD))
    }
    
    if(onlyMean_sd){
      nvarrep = 2
      out1var = as.data.frame(cbind(mean=custM, sd = custSD))
    }
    
    if(is.null(outdf)){
      outdf = out1var
    }else{
      outdf = cbind(outdf, out1var)
    }
  }
  colnames(outdf) = paste0(rep(var0,each=nvarrep),"_",colnames(outdf))
  return(outdf)
}
# source("D:/liucq/Millennial/Millennial3/derivs_V2_MM_MAD.R")
# Rdatadir = "D:/liucq/Millennial/Millennial3/input/finalRdata/250_160"
# Rdatadir = "D:/liucq/Millennial/Millennial3/input/finalRdata"

# source("E:/E/code/R/Millennial/Millennial3/derivs_V2_MM_MAD - ????.R")
library(tidyr)
library(openxlsx)
library(dplyr)
# outdir = "E:/E/code/R/Millennial/final/exp/output"
# indir = "E:/E/code/R/Millennial/final/input"
# indirexp = "E:/E/code/R/Millennial/final/exp/input"
# indirR = "E:/E/code/R/Millennial/final/Rcode" #/Fld
# indirRexp = "E:/E/code/R/Millennial/final/exp/Rcode"  #/Fld
# source("E:/E/code/R/Millennial/final/exp/Rcode/expdata_fun.R")

# if(data=="first"){
#   source("E:/E/code/R/Millennial/final/exp/output/Rdata/200_160/test/derivs_V2_MM_MAD.R")
#   Rdatadir = "E:/E/code/R/Millennial/final/exp/output/Rdata/200_160/test"
#   
# }
# if(data=="Fld"){
#   source("E:/E/code/R/Millennial/final/Rcode/Fld/derivs_V2_MM_MAD.R")
#   Rdatadir = "E:/E/code/R/Millennial/final/exp/output/Fld_Rdata"
# }
# if(data=="v0"){
#   # source("E:/E/code/R/Millennial/final/Rcode/Fld/derivs_V2_MM_MAD.R")
#   Rdatadir = "E:/E/code/R/Millennial/final/exp/output/clab_Rdata/300noss"
# }
testMd = c("_parsMax", "_parsMaxbeta","_fixbeta","_fixbeta_Obj","_fixbeta_rmMAOCmm",
           "_objCO2R2MARE","_rmPOCtoMAOC","_VMA_VMD","_rmMAOCmm","rmPOCtoMAOC_MAOCmm_Vpdown31_MAdtoPOC_sw")[3]
print(testMd)
testMD2 = c(FALSE, TRUE)[2] #??Χ??С???Ż?Ŀ??ȥ????MAOC???Ϸ???ȥ??

# testoutnm = c("rmPOCtoMAOC_Vdown10_MAOCmm",
#               "rmPOCtoMAOC_Vdown10_MAOCmm_limitfoa60",
#               "Vdown10_MAOCmm_limitfoa60",
#               "rmPOCtoMAOC_MAOCmm_Vp5",
#               "rmPOCtoMAOC_MAOCmm_Vp5_fixkm",
#               "rmPOCtoMAOC_MAOCmm_Vp1_fixkm",
#               "rmPOCtoMAOC_MAOCmm_Vpdown",     #7
#               "rmPOCtoMAOC_MAOCmm_Vpdown31",   #8
#               "rmPOCtoMAOC_MAOCmm_Vpdown31_sw", #9
#               "rmPOCtoMAOC_MAOCmm_Vpdown31_sw_MAdtoPOC", #10
#               "rmPOCtoMAOC_MAOCmm_Vpdown31_wpMAD_sw", #11
#               "rmPOCtoMAOC_MAOCmm_Vpdown31_wpMAD_MAdtoPOC_sw")[8]

testoutnm = c("rmPOCtoMAOC_MAOCmm_Eap66_sw05",
              "rmPOCtoMAOC_MAOCmm_Eap66_sw05_vdown31",
              "rmPOCtoMAOC_MAOCmm_Eap66_sw05_madtopoc",
              "rmPOCtoMAOC_MAOCmm_Eap66_sw05_vdown31_madtopoc", #4
              "rmPOCtoMAOC_rmMAOCmm_Eap66_sw05_vdown31_madtopoc"
              )[5]
print(testoutnm)

{
  library(openxlsx)
  indir = "E:/E/code/R/Millennial/final/input"
  indirexp = "E:/E/code/R/Millennial/final/exp/input"
  indirRexp = "E:/E/code/R/Millennial/final/exp/Rcode/MAD_CO2_lab_ss_final"
  source(file.path(indirRexp,"derivs_V2_MM_AD_CO2_Clab_vmax_cuetime.R"))
  source(file.path(indirRexp,"ss_trans_processfun.R"))
  source(file.path(indirRexp,"expdata_fun_clab.R"))
  # Rdatadir = "E:/E/code/R/Millennial/final/exp/output/cppRdata/300" #??һ??????
  Rdatadir = file.path("E:/E/code/R/Millennial/final/exp/output/V6Rdata",testoutnm) #????
  fitpars = pars
  
  #ͬʱģ????̬
  {
    sstype = c("runsteady","ND","ND_cpp")[3]
    ss_cond_obj = c(TRUE,FALSE)[1] #?ʶ??????Ƿ?Ҫ??????̬?жϳͷ???
    print(ss_cond_obj)
    nloop = 25*50#20
    library(Rcpp)
    sourceCpp(file.path(indirRexp,"Rmillennial2.cpp"))
    alloutss = data.frame(matrix(NA, nrow = 10,ncol = 20*7))
    colnames(alloutss) = rep(fixCp[1:7],times = 20)
  }
}



#############
string = grep(testMd, list.files(Rdatadir,pattern = ".RData"), value = T)
string = grep(testMD2, string, value = T)
print(string)
sitei = c()
Fit.pools_iter0 = list()
for(ri in 1:(length(string))){
  sitename = unlist(strsplit(string[ri],split="_"))[2]
  sitei = c(sitei,as.numeric(unlist(strsplit(string[ri],split="_"))[1])) # ??ȡ"1"
  
  load(file.path(Rdatadir,string[ri]))
  Fit.pools_iter0[[sitename]] = Fit.pools_iter[[sitename]]  #???е?????
  # Fit.pools_iter0[[sitename]] = Fit.pools[[sitename]] #160??
  print(paste("sitename=",sitename,"sitei=",sitei[ri]))
}
Fit.pools_iter = Fit.pools_iter0
sitei = sitei[1:ri]#sitei[1:10]

fitparnam = colnames(Fit.pools_iter[[1]]$par)
npar = 20 #50 #length(res.par)#
cond_stode = c(FALSE,TRUE)[2]
tss.df = pmi= plm = pco2=pco2.c13=list()
dfnpar=data.frame()
par(mfrow=c(2,5)) #??̬̼?Ⲽ??
mpars_type=c("quatile","mean")[2]
# layout(matrix(c(1:15),3))#,byrow=F
# layout.show(15)#:nsite
pars_resl = list() #ѡ???Ĳ?????????ȫ??ģ??
sitnm = c("TM","EEGN","CF","WDLC","TC","SR","BJ","HB","WH","GL")
for(i in c(1:length(string))){
  
  # par(mfrow = c(2,4))
  # Simsite = c("??????1","????","??????","??ɽ????","????","????",
  #             "????","?ݵ???","??????2","??????")[i]
  # i = 2
  si = sort(sitei)[i]
  
  Simsite = c("??????1","????","??????","??ɽ????","????","????",
              "????","?ݵ???","??????2","??????")[si]
  print(Simsite)
  nSimsamep = 2
  Fit.pools0 = list()
  Obs.pools$SimT = rep(rep(expT,each=nrow(Obs.pools)/(nSimsamep*length(expT))),nSimsamep)
  Obs.co2_tt$SimT = rep(rep(expT,each=nrow(Obs.co2_tt)/(nSimsamep*length(expT))),nSimsamep)
  Obs.co2_C13_tt$SimT = rep(rep(expT,each=nrow(Obs.co2_C13_tt)/(nSimsamep*length(expT))),nSimsamep)
  #ѡ????????
  #
  
  optim.out10 = NULL#data.frame(day=rep(1:simd,2))
  fitpars = unlist(pars)
  
  Fit.pp = Fit.pools_iter[[Simsite]] #Fit.pools[[Simsite]] #
  #ֻ??200??160??ѡ????32160
  # Fit.pp = lapply(Fit.pp, function(df) {
  #   return(df[1:32160,])
  # })
  
  # Fit.pp$value = Fit.pp$value/100
  #??ѡ??????valueС??1??
  # Fit.pp$value[,1] = Fit.pp$value[,1]/10
  res.par=c()
  for(nr in 1:nrow(Fit.pp$value)){
    if(length(which(Fit.pp$value[nr, ]<=1))==ncol(Fit.pp$value)){
      res.par = c(res.par,nr)
    }
  }
  #ֻɸѡ?ﵽ??̬??  res.par=c()
  if(F){
    res.par = which(Fit.pp$value[,"ssgx"]==1)
    min(Fit.pp$value[res.par,1])
    min(Fit.pp$value[res.par,2])
    min(Fit.pp$value[res.par,3])
  }


  if(length(res.par)==1){
    Fit.pools0$par = t(as.data.frame(Fit.pp$par[res.par,]))
    Fit.pools0$value = t(as.data.frame(Fit.pp$value[res.par,]))
  }else{
    Fit.pools0$par = Fit.pp$par[res.par,]
    Fit.pools0$value = Fit.pp$value[res.par,]
  }
  
  if(is.null(colnames(Fit.pools0$par))){
    colnames(Fit.pools0$par) = fitparnam
  }
  
  #length(which(Fit.pools0$value[, 7]!=1)) #????Inf?Ĳ????ж?????
  #res.par = c(1:nrow(Fit.pools0$value))
  
  if(testMd == "_objCO2R2MARE" | testMd == "_fixbeta"| testMD2){
    R2.co2 = Fit.pools0$value[,1]
    # R2.co2.day = Fit.pools0$value[,2]
    # MARE.co2 = Fit.pools0$value[,3]
    MARE.lmwc = Fit.pools0$value[,2]
    MARE.mic = Fit.pools0$value[,3]
  }else{
    R2.co2 = Fit.pools0$value[,1]
    R2.co2.day = Fit.pools0$value[,2]
    MARE.co2 = Fit.pools0$value[,3]
    MARE.lmwc = Fit.pools0$value[,4]
    MARE.mic = Fit.pools0$value[,5]
  }

  
  if(cond_stode){
    # meanARE.SS = Fit.pools0$value[res.par,6]
    if(testMd == "_objCO2R2MARE" | testMd == "_fixbeta"| testMD2){
      meanARE.SS = Fit.pools0$value[,4]
      wtpar = c(2,1,1,1) #c(2,1,2,1,1,1)
      wtpar = switch(Simsite, 
                     "??????1"= c(20,1,2,1),#c(2,1,3,0.5),
                     # "????"= c(2,1,1,1),
                     "??????"= c(2,1,1,0.5),
                     "??ɽ????"= c(20,1,1,1),
                     "????"= c(20,1,1,1),#c(2,1,1,0.5), #
                     "????"= c(5,1,1,0.5),
                     "????"= c(5,1,1,1),
                     "?ݵ???"= c(10,1,1,1),
                     "??????2"= c(10,1,1,0.5),
                     "??????"= c(10,1,1,0.5),
                     wtpar)
      
      eff = sqrt(((1-R2.co2)^2 *wtpar[1] + (1-MARE.lmwc)^2*wtpar[2]+ 
                    (1-MARE.mic)^2*wtpar[3]+(1-meanARE.SS)^2 *wtpar[4])/sum(wtpar))
      
    }else{
      meanARE.SS = Fit.pools0$value[,6]
      wtpar = c(2,1,2,1,1,1) #c(2,1,2,1,1,1)
      wtpar = switch(Simsite, 
                     "??????1"= c(2,1,2,1,3,0.5),
                     "??????"= c(2,1,5,1,1,0.5),
                     "????"= c(2,1,5,1,1,0.5),
                     "????"= c(2,2,1,1,1,1),
                     wtpar)
      eff = sqrt(((1-R2.co2)^2 *wtpar[1] + (1-R2.co2.day)^2 *wtpar[2] 
                  + (1-MARE.co2)^2 *wtpar[3] 
                  + (1-MARE.lmwc)^2*wtpar[4]+ (1-MARE.mic)^2*wtpar[5]+(1-meanARE.SS)^2 *wtpar[6])/sum(wtpar))
      
    }
    
  }
  # else if(clab[1]!=FALSE){
  #   # eff = sqrt(((1-R2.co2)^2 *2 + (1-R2.co2.day)^2 *2 + (1-MARE.co2)^2 *2  + (1-MARE.lmwc)^2+ (1-MARE.mic)^2)/8)
  #   MARE.co2.c13 = Fit.pools0$value[,6]
  #   wtpar = c(1,1,2,1,1,2) #c(2,1,2,1,1,2)
  #   eff = sqrt(((1-R2.co2)^2 *wtpar[1] + (1-R2.co2.day)^2 *wtpar[2]
  #               + (1-MARE.co2)^2 *wtpar[3]
  #               + (1-MARE.lmwc)^2*wtpar[4]+ (1-MARE.mic)^2*wtpar[5]+(1-MARE.co2.c13)^2 *wtpar[6])/sum(wtpar))
  # }
  else if(ncol(Fit.pp$value)==5){
    wtpar = c(2,1,2,1,1) #
    eff = sqrt(((1-R2.co2)^2 *wtpar[1] + (1-R2.co2.day)^2 *wtpar[2]
                + (1-MARE.co2)^2 *wtpar[3]
                + (1-MARE.lmwc)^2*wtpar[4]+ (1-MARE.mic)^2*wtpar[5])/sum(wtpar))
  }
  
  parsgrp = order(eff,decreasing = T)[1:npar] #eff?????????±?
  
  #????ѡ???Ĳ?????????ȫ??ģ??
  pars_resl[[Simsite]]$value = Fit.pools0$value[parsgrp,]
  pars_resl[[Simsite]]$par = Fit.pools0$par[parsgrp,]
 
  sort(eff,decreasing=T)[1:npar]
  res.par[parsgrp] #ѡ???ĵڼ??ε???
  Fit.pools0$value[parsgrp,]
  #range(Fit.pools0$value[,6])
  Fit.pools0$par[parsgrp,]
  #ѡ???Ĳ?????????ͼ
  dfnpar = rbind(dfnpar,as.data.frame(cbind(group=rep(Simsite,npar),
                                            Fit.pools0$par[parsgrp,],
                                            Fit.pools0$value[parsgrp,])))
  #}   
  if(i==nsite){
    colnames(dfnpar) = c("group",colnames(Fit.pools0$par))
  }
  
  
  #ģ?ⲻͬ?¶ȣ???ͬ????
  for(Ti in 1:length(expT)){
    #?۲?????
    Obs.pools.t =  Obs.pools[which((Obs.pools$site %in% Simsite) &(Obs.pools$SimT %in% expT[Ti])),]
    # Obs.co2_C130 = Obs.co2_C13[which(Obs.co2_C13$site %in% Simsite),]
    Obs.co2_tt0 = Obs.co2_tt[which((Obs.co2_tt$site %in% Simsite) &(Obs.co2_tt$SimT %in% expT[Ti])),] #?ۻ?
    Obs.co2_C13_tt0 = Obs.co2_C13_tt[which(Obs.co2_C13_tt$site %in% Simsite & (Obs.co2_C13_tt$SimT %in% expT[Ti])),"C13_CO2"]#C13_CO2?ۼ?
    
    plotdf = data.frame(x=1:30,yobs0=NA,yobs1=NA,yobs0_c13=NA,yobs1_c13=NA,
                        ysim0.10=NA,ysim0.50=NA,ysim0.90=NA
                        ,ysim1.10=NA,ysim1.50=NA,ysim1.90=NA)
    
    if(Simsite=="??????1" | Simsite=="????" | Simsite=="??????" | Simsite=="??ɽ????"){
      plotdf$yobs0[CO2_30] = Obs.co2_tt0[c(1:(length(CO2_30))),"CO2"]
      plotdf$yobs1[CO2_30] = Obs.co2_tt0[-c(1:(length(CO2_30))),"CO2"]
      plotdf$yobs1_c13[c13CO2_30] = Obs.co2_C13_tt0
      
      Obs.pools.t$day = rep(ML30,2)
    }else{
      plotdf$yobs0[CO2_10] = Obs.co2_tt0[c(1:(length(CO2_10))),"CO2"]
      plotdf$yobs1[CO2_10] = Obs.co2_tt0[-c(1:(length(CO2_10))),"CO2"]
      plotdf$yobs1_c13[c13CO2_10] = Obs.co2_C13_tt0
      Obs.pools.t$day = rep(ML10,2)
    }
    
    #parsgrp = res.multipar
    for(pg in 1:length(which(!is.na(parsgrp)))){
      if(pg==1){
        optim.out10 = NULL
      }
      fitpars[colnames(Fit.pools0$par)] = Fit.pools0$par[parsgrp[pg],]
    
      optim.out0 = Model.pools(fitpars,inputdata[which(inputdata$site %in% Simsite),],expT = expT[Ti])  #???ڷ???MIC??LMWC??CO2
      optim.out = optim.out0[which(optim.out0$simT ==expT[Ti]),]
      obj_df = optim.out[-which(optim.out$time %in% c(0,31,32)),which(colnames(optim.out) %in% c("Tco2","LMWC","MIC","Tco2_c13_Cinoa"))]#?Ѿ??ѱ???ֵ??31????ɾ????
      
      # ??̬
      outps = calpars_SSfun(sstype = sstype,fitpars=fitpars) #step,runsteady,????data.frame
      sscid = which(names(outps) %in% fixCp[1:7]) #??̬ȥ??POC
      ssc = outps[sscid]
      
      #??10??ģ?????????кϲ?
      if (is.null(optim.out10)) {
        optim.out10 <- obj_df
        outss = ssc
      } else {
        optim.out10 <- cbind(optim.out10, obj_df) #optim.out10???ݿ????кϲ?????Ϊ0??0??
        outss = c(outss, ssc)
      }
      
    } #????ѭ??????
    
    #ѡ??10?????????н?????0.1??0.9??λ??
    if(mpars_type=="quatile"){
      LMWCp = percfun(optim.out10[,seq(1,ncol(optim.out10),length(unique(colnames(optim.out10))))]) #LMWC
      MICp = percfun(optim.out10[,seq(2,ncol(optim.out10),length(unique(colnames(optim.out10))))]) #MIC
      CO2p = percfun(optim.out10[,seq(3,ncol(optim.out10),length(unique(colnames(optim.out10))))]) #CO2
      CO2p.c13 = percfun(optim.out10[,seq(4,ncol(optim.out10),length(unique(colnames(optim.out10))))]) #CO2
      
    }else if(mpars_type=="mean"){
      if(pg==1){
        LMWCp = MICp =CO2p =CO2p.c13 = data.frame(matrix(0, nrow = nrow(optim.out10), ncol = 3))
        #??ͬ????Ϊһ?β???
        LMWCp[,2] = optim.out10[,"LMWC"]
        MICp[,2] = optim.out10[,"MIC"]
        CO2p[,2] = optim.out10[,"Tco2"]
        #CO2p.c13[,2] = optim.out10[,"Tco2_c13_Cinoa"]
        colnames(LMWCp) = paste0("LMWC", c("_sd0","_mean","_sd1"))
        colnames(MICp) = paste0("MIC", c("_sd0","_mean","_sd1"))
        colnames(CO2p) = paste0("Tco2", c("_sd0","_mean","_sd1"))
        colnames(CO2p.c13) = paste0("Tco2_c13_Cinoa", c("_sd0","_mean","_sd1"))
        
      }else{
        #??ͬ????Ϊһ?β???
        LMWCp = mpars_meansd(optim.out10,"LMWC")
        MICp = mpars_meansd(optim.out10,"MIC")
        CO2p = mpars_meansd(optim.out10,"Tco2")
        CO2p.c13 = mpars_meansd(optim.out10,"Tco2_c13_Cinoa")
      }

      
    }
    
    plotdf[,c("ysim0.10","ysim0.50","ysim0.90")] = CO2p[1:30,]
    plotdf[,c("ysim1.10","ysim1.50","ysim1.90")] = CO2p[-c(1:30),]
    
    plotdf[,c("ysim0.10_c13","ysim0.50_c13","ysim0.90_c13")] = NA
    plotdf[,c("ysim1.10_c13","ysim1.50_c13","ysim1.90_c13")] = CO2p.c13[-c(1:30),]
    
    
    # print(intervalggplot(plotdf[1:30,],expT0=expT[Ti]))
    
    #??CO2????ͼ
    # pi = Ti #????3???¶?
    pi = ((i-1)*3 + Ti) #10????*3???¶?
    print(pi)
    #??c13??CO2??ѡ??��
    nc13 = which(sapply(strsplit(colnames(plotdf), "_"), function(y) tail(y, 1)) == "c13")
    pco2[[pi]] <- intervalggplot(plotdf[1:30,-nc13],expT0=expT[Ti],nline=2,
                                 cust.ylab = bquote("Cumulative CO"[2]~"(gC m"^-2*")"))
    pco2.c13[[pi]] <- intervalggplot(plotdf[1:30,c(1,nc13)],expT0=expT[Ti],nline=1,
                                     cust.ylab = bquote("Cumulative C13-CO"[2]~"(gC m"^-2*")"))
    #??DOC??MIC??״??????ͼ
    dfMIC = dfbargg(MICp,"MIC",colnames(MICp))
    pmi[[pi]] <- ebarggplot(dfMIC,var = "MIC",error=TRUE)
    dfLM = dfbargg(LMWCp,"LMWC",colnames(LMWCp))
    plm[[pi]] <- ebarggplot(dfLM,var = "LMWC",error=TRUE)
    
    #????ÿ???¶??ǵ??��????ģ??´?ѭ????û???ˣ?
    #ͬһ???��?ͬ????y?᷶Χδͳһ,ֻ?ܱ?????��ÿ???¶ȵ?????ֵ???????޸?ymax
    if(Ti==1){
      ymax_co2 = 0
      ymax_co2_c13 = 0
      ymax_mi = 0
      ymax_lm = 0
    }
    ymax_co2 = max(ymax_co2, max(plotdf[,c(-1,-nc13)],na.rm = T))
    ymax_co2_c13 = max(ymax_co2_c13, max(plotdf[,nc13],na.rm = T))
    ymax_mi = max(ymax_mi, max(dfMIC[,3:5]))
    ymax_lm = max(ymax_lm, max(dfLM[,3:5]))
  } #Ti ends
  #ͬһ??3???¶?y??һ??
  pco2 = ymaxggp(p=pco2,pi = pi, ymax = ymax_co2)
  # pco2.c13 = ymaxggp(p=pco2.c13,pi = pi, ymax = ymax_co2_c13)
  pmi = ymaxggp(p=pmi,pi = pi, ymax = ymax_mi)
  plm = ymaxggp(p=plm,pi = pi, ymax = ymax_lm)
  
  #??̬
  alloutss[si,1:length(outss)] = outss
}
#????̬????ƽ??ֵ
colnames(alloutss) = rep(fixCp[1:7],times = 20)
print(testoutnm)
round(rowMeans(alloutss[,grep("MAOM",colnames(alloutss))]/alloutss[,grep("SOM",colnames(alloutss))], na.rm=T) * 100,0)
mean(round(rowMeans(alloutss[,grep("MAOM",colnames(alloutss))]/alloutss[,grep("SOM",colnames(alloutss))]) * 100,0), na.rm=T)
c(round(mean(colMeans(alloutss[,grep("SOM",colnames(alloutss))], na.rm=T))),
  round(mean(colMeans(alloutss[,grep("MAOM",colnames(alloutss))], na.rm=T))),
  round(mean(colMeans(alloutss[,grep("POM",colnames(alloutss))], na.rm=T)))
)
round(alloutss[1:10,1:7])
round(alloutss[,grep("SOM|POM|MAOM",colnames(alloutss))])


#????POC????
{
  fitpars0 = fitpars
  soccP = data.frame()
  for(vpi in c(0.1,0.5,1,5,10)){
    fitpars0["Vpl0"] = fitpars["Vpl0"] *  vpi
    outps = calpars_SSfun(sstype = sstype,fitpars=fitpars0)
    soccP = rbind(soccP, c(fitpars0["Vpl0"], outps[c("SOM","POM","MAOM")]))
  }
  colnames(soccP) = c("vpl","SOM","POM","MAOM")
  pltdf = pivot_longer(soccP, cols = c("SOM","POM","MAOM"), values_to = "cfraction",names_to = "group") 
  ggplot(pltdf, aes(x = vpl, y = cfraction, colour = group)) + 
    geom_line(aes(colour = group), linewidth = 2)+
    geom_point(aes(colour = group), size = 6) + 
    geom_hline(yintercept = 1500) + 
    theme_bw()
}
#??co2????
{
  co2rate = apply(optim.out10[,grep("Tco2",colnames(optim.out10))],2 , function(col0){diff(col0)})
  
  out = Model.pools(fitpars,inputdata[which(inputdata$site %in% Simsite),],expT = c("15","25","35"))
  out[,c("CO2_day")] = calculate_diff(out[,"Tco2",drop=FALSE]) #?????¶?Ҳֱ?Ӽ?????Ϊ?????˳?ʼֵ???Ե?
  Condoutco2r = diff(out[c(168,nrow(out)), "CO2_day"]) >= 0
  
  Obs.co2_tt0 = Obs.co2_tt[which(Obs.co2_tt$site %in% Simsite),"CO2"] #?ۻ?
  Condoutco2Cul = out[nrow(out), "Tco2"] > (Obs.co2_tt0[length(Obs.co2_tt0)]*2.2)
  wlimtR = ifelse(Condoutco2r | Condoutco2Cul, 100, 1) #????CO2????һֱ???Ӻ?
  
  out[nrow(out), "Tco2"]/Obs.co2_tt[tail(which(Obs.co2_tt$site %in% Simsite),1),"CO2"] 
                  
}

#????ˮ??Ӱ?????Ӵη?ֵ?????õ?һ??????
{
  soccP = data.frame()
  
  #1??ָ????????0.5??
  sourceCpp(file.path(indirRexp,"Rmillennial2.cpp"))
  outps = calpars_SSfun(sstype = sstype,fitpars=fitpars)
  soccP = c(exp0 = "0.5", outps)
  
  #2????????1??
  sourceCpp(file.path(indirRexp,"Rmillennial2.cpp"))
  outps = calpars_SSfun(sstype = sstype,fitpars=fitpars)
  soccP = rbind(soccP, c(exp0 = "1", outps))
  
  #3????????poc????1??maoc????0.5
  sourceCpp(file.path(indirRexp,"Rmillennial2.cpp"))
  outps = calpars_SSfun(sstype = sstype,fitpars=fitpars)
  soccP = rbind(soccP, c(exp0 = "1_0.5", outps))

  soccP$exp0   = c("0.5", "1", "1_0.5")
  round(soccP)
}

#??MA?仯
{
  plot(optim.out[,"MA"])
}

#MA??MDת??ˮ??
{
  inputdata[1:10,]
  SWP = c()
  for(pi in 1:10){
    SWP = c(SWP, fSWC2SWP(inputdata[pi,"forc_sw"],SWCsat=inputdata[pi,"porosity"]))  #MPa,SWP<0,???ͻ???????=0
  }
 
  round(SWP,3)
  
  w = 1
  SWPa2d = -0.4 #MPa
  fwA2D = round(abs(SWP)^w/(abs(SWP)^w + abs(SWPa2d)^w),3)
  
  tauda = dfnpar[,1:22] %>% group_by(group) %>% summarise(mean(tauda))
  
  SWPd2a = SWPa2d * tauda #MPa
  fwD2A = round(abs(SWPd2a)^w/(abs(SWP)^w + abs(SWPd2a)^w),3)
  
  VmA2D = Vm /scalar_wb * fwA2D# * tp_scalar("MR") #* scalar_wd #* wp_scalar_low
  VmD2A = Vm /scalar_wb * fwD2A# * tp_scalar("MR") #* scalar_wd #* wp_scalar
  
}
# save(pars_resl,file = paste0(indirRexp,"/new/parsResl_cpp_",testoutnm,".RData")) #????ȷ????????????ѡ??????
#???㻭ͼ
multiplot(plotlist = pco2[16:30], cols = 5,layout = matrix(c(1:15),nrow=3,byrow=F))

library(Rmisc)
multiplot(plotlist = c(pco2[(pi-2):pi],pmi[(pi-2):pi],plm[(pi-2):pi]), 
          cols = 3,layout = matrix(c(1:9),nrow=3,byrow=F))

multiplot(plotlist = pco2[(pi-2):pi], cols = 3,layout = matrix(c(1:15),nrow=3,byrow=F))
multiplot(plotlist = pco2.c13[(pi-2):pi], cols = 3,layout = matrix(c(1:15),nrow=3,byrow=F))
multiplot(plotlist = pmi[(pi-2):pi], cols = 3,layout = matrix(c(1:9),nrow=3,byrow=F))
multiplot(plotlist = plm[(pi-2):pi], cols = 3,layout = matrix(c(1:9),nrow=3,byrow=F))



#????,??CO2??MBC??DOC
library(ggpubr)
print(string)

layshow(pco2)
layshow(pmi)
layshow(plm)

multiplot(plotlist = pco2[16:30],cols = 5,layout = matrix(c(1:15),nrow=3,byrow=F))
multiplot(plotlist = pmi[16:30],cols = 5,layout = matrix(c(1:15),nrow=3,byrow=F))
multiplot(plotlist = plm[16:30],cols = 5,layout = matrix(c(1:15),nrow=3,byrow=F))

#????????????ͼ
ppars = list()
for(np in 1:ncol(Fit.pools0$par)){
  ppars[[np]] <- boxgg(dfnpar[,c("group",colnames(Fit.pools0$par)[np])])
}
#???԰??????򣬵???ͼ??û?취ֻ??ʾһ??
library(Rmisc)
p=ppars
multiplot(plotlist = p[1:10], cols = 3,layout = matrix(c(1:10),nrow=2,byrow=F))
multiplot(plotlist = p[11:20], cols = 3,layout = matrix(c(1:10),nrow=2,byrow=F))

#????value
perval = data.frame()
for(vi in 1:10){
  st = (vi-1)*50 +1
  end = vi*npar
  perval = rbind(perval,round(apply(dfnpar[st:end,22:27], 2, function(x) quantile(as.numeric(x), seq(0,1,1/4))),3))
}
plot(perval[,1],type="o",col="orange")


######ģ?Ͳ???????########
######???????е?ģ??ģ??????#############
{

  testT = c("15","25","35")
  (as.numeric(testT[-1])-as.numeric(testT[1]))/(simd+1)#????????
  load(paste0(indirRexp,"/new/parsResl_cpp_",testoutnm,".RData"))#pars_resl[[Simsite]]$par[pi,]
  flux.cul = c(FALSE,TRUE)[2]
  names(pars_resl)

  simVar = c("Tco2","MIC","LMWC","MAOM","POM","MA","MD")[1:7]
  Rco2mean10 = allsoil.out = NULL
  npar = 20
  pco2r = list()
  for(si in 1:10){
    Simsite = c("??????1","????","??????","??ɽ????","????","????",
                "????","?ݵ???","??????2","??????")[si]
    Fit.pools0 = pars_resl[[Simsite]]
    
    for(pi in 1:nrow( Fit.pools0$par)){
      if(pi==1){
        optim.out10 = NULL
      }
      fitpars[colnames(Fit.pools0$par)] = Fit.pools0$par[pi,]
      #ֻģ?????Ӳ?????
      # optim.out = Model.pools(fitpars,inputdata[which(inputdata$site %in% Simsite)[2],],expT = testT) 
      #ģ??��??
      optim.out = Model.pools(fitpars,inputdata[which(inputdata$site %in% Simsite),],expT = testT) 
      
      # ??̬
      outps = calpars_SSfun(sstype = sstype,fitpars=fitpars) #step,runsteady,????data.frame
      sscid = which(names(outps) %in% fixCp[1:7]) #??̬ȥ??POC
      ssc = outps[sscid]
      
      if (is.null(optim.out10)) {
        optim.out10 <- optim.out
        outss = ssc
      } else {
        optim.out10 <- cbind(optim.out10, optim.out) 
        outss = c(outss, ssc)
      }
    }
    
    #20??????????????̼ͨ��??̼???仯????
    # timedel = -which(optim.out10$time %in% c(31,32))
    # meanout = mpars_meansd(optim.out10[timedel,which(colnames(optim.out10)%in% Tcolnms)],
    #                        Tcolnms,outsd = FALSE)
    colnames(optim.out10) = gsub("\\.\\d+$", "", colnames(optim.out10))
    if(is.null(allsoil.out)){
      allsoil.out=optim.out10[,which(colnames(optim.out10)%in% simVar)]
    }else{
      allsoil.out = cbind(allsoil.out, optim.out10[,which(colnames(optim.out10)%in% simVar)])
    }
    alloutss[si,] = outss
    
  }
  
  ##################????20??????ƽ??co2??MBC??DOC
  outlist = list()
  
  sitnm = c("TM","EEGN","CF","WDLC","TC","SR","BJ","HB","WH","GL")
  for(vi in simVar){
    allsoil.out0 = allsoil.out[,grep(paste0(vi,"*"),colnames(allsoil.out))]
    colnames(allsoil.out0) = rep(sitnm,each=20)
    
    outlist[[vi]] = cbind(Simday = rep(0:32,times = 6), mpars_meansd(allsoil.out0,sitnm,onlyMean_sd=TRUE))
    
  }
  #????̬????ƽ??ֵ
  colnames(alloutss) = rep(fixCp[1:7],times = 20)
  round(rowMeans(alloutss[,grep("MAOM",colnames(alloutss))]/alloutss[,grep("SOM",colnames(alloutss))]) * 100,0)
  
  writexl::write_xlsx(outlist, path = paste0("E:/E/code/R/Millennial/final/exp/output/","allsoil.out_ModelSimOUT_",testoutnm,".xlsx"))
}

#????̬ռ??
{
  plot(allsoil.out[,"MA"]/allsoil.out[,"MIC"])
  plot(allsoil.out[,"LMWC"])
}

###  ?????õ?   #############
#10??????Q10ֵ#################
testT = c("15","16","17","20","25") #??15?濪ʼ????
testT = c("25","26","27","30","35") #??25?濪
#
testT = c("15","25","35")
(as.numeric(testT[-1])-as.numeric(testT[1]))/(simd+1)#????????
load(file.path(indirRexp,"parsResl_cpp.RData"))#pars_resl[[Simsite]]$par[pi,]
flux.cul = c(FALSE,TRUE)[2]
names(pars_resl)

simVar = c("Tco2","MIC","LMWC")[1:3]
Rco2mean10 = allsoil.out = NULL
npar = 20
pco2r = list()
for(si in 1:10){
  Simsite = c("??????1","????","??????","??ɽ????","????","????",
              "????","?ݵ???","??????2","??????")[si]
  Fit.pools0 = pars_resl[[Simsite]]
  
  for(pi in 1:npar){
    if(pi==1){
      optim.out10 = NULL
    }
    fitpars[colnames(Fit.pools0$par)] = Fit.pools0$par[pi,]
    #ֻģ?????Ӳ?????
    # optim.out = Model.pools(fitpars,inputdata[which(inputdata$site %in% Simsite)[2],],expT = testT) 
    #ģ??��??
    optim.out = Model.pools(fitpars,inputdata[which(inputdata$site %in% Simsite),],expT = testT) 
    
    
    if (is.null(optim.out10)) {
      optim.out10 <- optim.out
    } else {
      optim.out10 <- cbind(optim.out10, optim.out) 
    }
  }
  
  #20??????????????̼ͨ��??̼???仯????
  # timedel = -which(optim.out10$time %in% c(31,32))
  # meanout = mpars_meansd(optim.out10[timedel,which(colnames(optim.out10)%in% Tcolnms)],
  #                        Tcolnms,outsd = FALSE)
  colnames(optim.out10) = gsub("\\.\\d+$", "", colnames(optim.out10))
  if(is.null(allsoil.out)){
    allsoil.out=optim.out10[,which(colnames(optim.out10)%in% simVar)]
  }else{
    allsoil.out = cbind(allsoil.out, optim.out10[,which(colnames(optim.out10)%in% simVar)])
  }
  
}

##################????20??????ƽ??co2??MBC??DOC
outlist = list()

sitnm = c("TM","EEGN","CF","WDLC","TC","SR","BJ","HB","WH","GL")
for(vi in c("Tco2","MIC","LMWC")){
  allsoil.out0 = allsoil.out[,grep(paste0(vi,"*"),colnames(allsoil.out))]
  colnames(allsoil.out0) = rep(sitnm,each=20)
  
  outlist[[vi]] = mpars_meansd(allsoil.out0,sitnm,onlyMean_sd=TRUE)
 
}
writexl::write_xlsx(outlist, path = file.path(indirexp,"allsoil.out_ModelSimOUT.xlsx"))

#co2????,20*10
co2Ratio = data.frame()
for(bsr in 2:(length(testT))){
  sts = (bsr-1)*33 + 1
  ends = bsr*33
  co2Rt = allsoil.out[sts:ends,]/allsoil.out[which(optim.out$simT==testT[1]),]
  co2Ratio = rbind(co2Ratio,co2Rt)
  #ÿ???¶ȵ?ƽ??ֵ
  
}
colnames(co2Ratio) = rep(names(pars_resl),each=npar)
#ÿ?????��?20????????ƽ??ֵ
RatinM20 = mpars_meansd(co2Ratio,names(pars_resl),outsd = FALSE)
tst = which(is.na(RatinM20[,1]))
co2Rm20 = data.frame()
for(ti in 1:(length(testT)-1)){
  ends = ifelse(ti==4, nrow(RatinM20),tst[ti+1])
  co2Rm20 = rbind(co2Rm20, colMeans(RatinM20[tst[ti]:ends,],na.rm = T)) #ģ??simd??ƽ???仯
}
rowMeans(co2Rm20)
RatinM20[tst[4]+1,] #??һ??
RatinM20[tst[4]+simd,] #????һ??
allsoil.out
dim(allsoil.out)

#ʵ??10?????��????л?????????̼????----------
load(file.path(indirRexp,"parsResl_cpp.RData"))#pars_resl[[Simsite]]$par[pi,]
names(pars_resl)

#??????��????
clab = switch(1,FALSE,c("Cin","Cinoa","Cad")[3]) #????C???? #??Ϊȫ?ֱ?��
clab.outvar = c("Cpools_Tco2", "all")[1] #????C??????Щ??��
#????????̼??
{
  if(clab[1]!=FALSE){
    MAOM_c130 = 0 #????????????ʼֵ
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
  
}
Rco2mean10 = allsoil.out = NULL
npar = 20
pco2r = list()
for(si in 1:10){
  Simsite = c("??????1","????","??????","??ɽ????","????","????",
              "????","?ݵ???","??????2","??????")[si]
  maxSimsite = envsite[,"maxSimNM"][si] #?????Ƶ??Ĳ???
  Fit.pools0 = pars_resl[[c(Simsite,maxSimsite)[2]]] 
  
  for(pi in 1:npar){
    if(pi==1){
      optim.out10 = NULL
    }
    fitpars[colnames(Fit.pools0$par)] = Fit.pools0$par[pi,]
    optim.out = Model.pools(fitpars,inputdata[which(inputdata$site %in% Simsite)[1:2],],expT = expT) 
    
    if (is.null(optim.out10)) {
      optim.out10 <- optim.out
    } else {
      optim.out10 <- cbind(optim.out10, optim.out) 
    }
  }
  #1.?л?????????CO2ռ??CO2????
  if(clab!=FALSE){
    timedel = -which(optim.out10$time %in% c(0,31,32))
    Rco2 = optim.out10[timedel,which(colnames(optim.out10)=="Tco2_c13_Cad")]/optim.out10[timedel,which(colnames(optim.out10)=="Tco2")]*100
    Rco2mean = mpars_meansd(Rco2,colnames(Rco2)[1])
    Rco2meanp = cbind(day=1:30,T15 = Rco2mean[1:30,], T25=Rco2mean[31:60,],T35=Rco2mean[61:90,])
    pco2r[[si]] = PointLineBar(Rco2meanp)
    if(is.null(Rco2mean10)){
      Rco2mean10 = Rco2mean
    }else{
      Rco2mean10 = cbind(Rco2mean10,Rco2mean)
    }
  }

  
  #20??????????????̼ͨ��??̼???仯????
  # timedel = -which(optim.out10$time %in% c(31,32))
  # meanout = mpars_meansd(optim.out10[timedel,which(colnames(optim.out10)%in% Tcolnms)],
  #                        Tcolnms,outsd = FALSE)
  if(is.null(allsoil.out)){
    allsoil.out=optim.out10#meanout
  }else{
    allsoil.out = cbind(allsoil.out, optim.out10)#meanout)
  }
  
}


library(Rmisc)
multiplot(plotlist = pco2r[1:10], cols = 5,layout = matrix(c(1:10),nrow=2,byrow=T))
library(ggpubr)
# pco2r[[1]] = pco2r[[1]] + theme(legend.text = element_text(size = 14),
#                                 legend.spacing.x = unit(10, 'cm'))
print(ggarrange(pco2r[[1]],pco2r[[2]],pco2r[[3]], pco2r[[4]],pco2r[[5]],
                pco2r[[6]],pco2r[[7]],pco2r[[8]], pco2r[[9]],pco2r[[10]],
                nrow=2,ncol=5,common.legend = TRUE, legend="top",
                legend.grob = get_legend(pco2r[[1]])))

  #10??????ƽ????????̼????ͨ��----------
#1.?л?????????ռ??CO2????
CO2s10 = Rco2mean10[,which(colnames(Rco2mean10)=="Tco2_c13_Cad_mean")]
CO2s10mean = mpars_meansd(CO2s10,"Tco2_c13_Cad_mean")
dim(CO2s10)
CO2s10meanp = cbind(day=1:30,T15 = CO2s10mean[1:30,], T25=CO2s10mean[31:60,],T35=CO2s10mean[61:90,])
PointLineBar(CO2s10meanp)
mean(CO2s10mean[,"Tco2_c13_Cad_mean_mean"])#3???¶?ƽ??
write.xlsx(CO2s10meanp,file="E:/E/carbon/mechanism/root/ģ?ͽ???????????/?л?????????CO2ռ?????????????ٷֱ?.xlsx")

#2.̼?⡢ͨ��
s10mean = mpars_meansd(allsoil.out, paste0(Tcolnms,"_mean"),outsd = FALSE)
dim(s10mean)
length(which(allsoil.out[ ,which(colnames(allsoil.out)=="f_LM_leach_c13_Cad_mean")]>0))
round(colMeans(s10mean[-seq(1,nrow(s10mean),31),]),2) #????ͨ��ƽ??
round(colMeans(s10mean[seq(31,nrow(s10mean),31),]),2) #3???¶??ۼƵ???30????ƽ??ֵ!!!!!!
round(colMeans(s10mean[seq(1,nrow(s10mean),31),]),0) #̼????ʼֵ
s10mean[31,]

  #3.?л?????????ͨ��/????DOC??ͨ�� = Facid/(DOC??????ͨ��+DOC0)--------
TLM_in = c("f_Lit_LM_mean","f_Exud_LM_mean","f_Oa_LM_mean","f_PO_LM_mean","f_MA_LM_mean","f_MB_LM_mean")
Facid_doc = allsoil.out[-seq(1,nrow(s10mean),31),which(colnames(allsoil.out)=="acid_des_mean")]

Tlmcol = which(colnames(allsoil.out) %in% TLM_in)
Tin_doc = NULL
for(i in 1:10){
  stl = (i-1)*length(TLM_in)+1
  endl = i*length(TLM_in)
  colx = Tlmcol[stl:endl]
  if(is.null(Tin_doc)){
    Tin_doc = rowSums(allsoil.out[,colx])
  }else{
    Tin_doc = cbind(Tin_doc, rowSums(allsoil.out[,colx]))
  }
}
Tdoc = Tin_doc[-seq(1,nrow(s10mean),31),] + allsoil.out[-seq(31,nrow(s10mean),31),which(colnames(allsoil.out)=="LMWC_mean")]

rco2f = Facid_doc/Tdoc * 100

CO2s10mean = mpars_meansd(rco2f,colnames(rco2f)[1])
dim(CO2s10)
CO2s10meanp = cbind(day=1:30,T15 = CO2s10mean[1:30,], T25=CO2s10mean[31:60,],T35=CO2s10mean[61:90,])
PointLineBar(CO2s10meanp)
mean(CO2s10mean[,"acid_des_mean_mean"])

  # ?????????????ĵ???ͼ---------
library(ggplot2)
PointLineBar <- function(data){
  p <- ggplot(data, aes(x = day)) +
    geom_errorbar(aes(ymin = data[[2]], ymax = data[[4]]), width = 0.1,color = "DimGrey") +  # ??????????
    geom_line(aes(y = data[[3]], color = "15??", group = 1),linewidth = 1) +  
    geom_point(aes(y = data[[3]], color = "15??"), size = 2) +  # ???????????ӵ?
    
    geom_errorbar(aes(ymin = data[[5]], ymax = data[[7]]), width = 0.1,color = "DimGrey") + 
    geom_line(aes(y = data[[6]], color = "25??", group = 1),linewidth = 1) +  
    geom_point(aes(y = data[[6]], color = "25??"), size = 2) +  #
    
    geom_errorbar(aes(ymin = data[[8]], ymax = data[[10]]), width = 0.1,color = "DimGrey") +
    geom_line(aes(y = data[[9]], color = "35??", group = 1),linewidth = 1) +  
    geom_point(aes(y = data[[9]], color = "35??"), size = 2) 
  
  
  p <- p + 
    geom_hline(yintercept = 5, linetype = "dashed", color = "DimGrey") +
    scale_color_manual(name = "", 
                       #values = c("15??"= "#19647E","25??"= "#D7B09E","35??" ="#B9181A"), 
                       values = c("15??"= "#2892C7","25??"= "#D7B09E","35??" ="#FF7F7F")
                       #breaks = c("?л???????", "?л?????????")  
    ) +
    
    labs(title = Simsite, x = "day", y = "????(%)") +  # ??????????ǩ
    theme_bw() +  # ???ü?????????
    theme(
      legend.position = c(0.95, 0.95),     # ͼ???????????Ͻǿ?????
      legend.justification = c(1, 1),      # ͼ??ê??Ϊ???Ͻ?
      legend.text = element_text(size = 14),
      legend.background = element_blank(), # ȡ??ͼ??????
      legend.box.background = element_blank(), # ȡ??ͼ???߿?
      legend.title = element_blank()       # ????ͼ??????
    )+ guides(color = guide_legend(nrow = 1, byrow = TRUE,reverse = F))
  return(p)
}




#####
#????̬ʱ??̼?⺬��--------
library(rootSolve)
SS.analyzefun <- function(st_npar,end_npar,si0=si){
  
  ssp.df = data.frame()
  do_forc_mean=1
  do_forc_dynam =0
  notss = 0 #ͳ?ƴﲻ????̬????Ŀ
  for(pgi in st_npar:end_npar){ 
    fitpars0 = pars#
    # print(Fit.pools0$par[parsgrp[pg],])
    print(paste("pg=",pgi))#,"parsgrp[pg]=",parsgrp[pg]))
    # print(Fit.pools0$value[parsgrp[pg],])
    fitpars0[colnames(Fit.pools0$par)] = Fit.pools0$par[parsgrp[pgi],]
    
    state0=c(POM=1,LMWC=1, MIC=1, MAOM=1) #,CO2=0
    MA0 = fitpars0["r0"] * state0["MIC"]
    MD0 = (1-fitpars0["r0"]) * state0["MIC"]
    state0 = c(state0,MA = MA0, MD =MD0)
    state0 = as.numeric(state0)
    names(state0) = c("POM","LMWC","MIC","MAOM","MA","MD")#"CO2",
    
    if(do_forc_mean){
      forc_stode0 = forc_mean(si0) #????ƽ??ֵ
      # print(forc_stode0)
      SStime0 = 1e8
      SS.p = runsteady(y = state0, func = derivs_Model_SS, parms = fitpars0, #soil_pars = soil_pars,
                       forc_st=forc_stode0["st"], forc_sw=forc_stode0["sw"], 
                       forc_lit= forc_stode0["lit"],forc_acid = forc_stode0["acid"]
                       ,verbose=FALSE, times = c(1, SStime0))
      # print(SS.p)
    }else if(do_forc_dynam){
      num.years0=5
      forc_stode0 = forc_dynam(Simsite,num.years=num.years0) #????ƽ??ֵ
      SStime0 = num.years0*ncol(forc_st_all)*30 #1e8
      
      start_time = Sys.time()
      SS.p = runsteady(y = state0, func = derivs_Model_SS_dynam, parms = fitpars, #soil_pars = soil_pars,
                       forc_st=forc_stode0$st, forc_sw=forc_stode0$sw, 
                       forc_lit= forc_stode0$lit,forc_acid = forc_stode0$acid
                       ,verbose=FALSE, times = c(1, SStime0-1))
      end_time = Sys.time()
      print(end_time-start_time)
    }
    
    # print(SS.p$y["MAOM"]/SS.p$y["POM"])
    if(attr(SS.p, "steady")==FALSE){
      print(paste('Not steady!',"pgi = ",pgi))
      print(SS.p$y)
      print(Fit.pools0$value[parsgrp[pgi],])
      notss <<- notss + 1 #?ﲻ????̬??
    }
    ssp.df = rbind(ssp.df,cbind(t(as.data.frame(SS.p$y)),attr(SS.p, "steady")))
    
  }
  print(head(ssp.df))
  
  colnames(ssp.df) = c(names(SS.p$y),"steady")
  rownames(ssp.df) = seq(st_npar,end_npar)
  # par(mfrow=c(1,1))
  xn = end_npar-st_npar+1
  plot(1:xn,ssp.df[,"MAOM"],col="blue",type="o",main=Simsite,
       ylim=c(min(ssp.df[,c("POM","MAOM")]), max(ssp.df[,c("POM","MAOM")])))+
    points(1:xn,ssp.df[,"POM"],col="orange",type="o")
  
  # print(quantile(ssp.df[,"MAOM"]/ssp.df[,"POM"],seq(0,1,1/10)))
  return(ssp.df)
}

#ʵ??̼??
obspp = inputdata[which(inputdata$site %in% Simsite)[1],c("POM","LMWC","MIC","MAOM")] #ʵ??̼??
inputdata[which(inputdata$site %in% Simsite)[1],"MAOM"]/inputdata[which(inputdata$site %in% Simsite)[1],"POM"]
t(t(round(inputdata[1:10,"MAOM"]/inputdata[1:10,"POM"],0)))

#??̬Ŀ??????
paran = seq(35,50)
print(Fit.pools0$value[parsgrp[paran],])
SS.analyzefun(35,50,3)
print(paste("the number of no steady : ",notss))

SS.analyzefun(38,38,3)
print(Fit.pools0$value[parsgrp[38],])

paran = 37
for(tt1 in paran ){# 1:npar
  print(tss.df[[Simsite]][tt1,1:4])
  abs(tss.df[[Simsite]][tt1,1:4]-obspp)/obspp
  # print(abs(tss.df[[Simsite]][tt1,1:4]-obspp)/obspp)
  
  par.ss = round(mean(as.numeric(abs(tss.df[[Simsite]][tt1,1:4]-obspp)/obspp))*2/10,3)
  cal.ss = round(Fit.pools0$value[parsgrp[tt1],6],3)
  if(par.ss!=cal.ss){
    print(paste("tt1=",tt1,"par.ss=",par.ss,"cal.ss=",cal.ss))
    print(paste("par.ss*ssgx=",par.ss*5, (par.ss*5) == cal.ss))
    
    #CO2??DOC??Ŀ??ֵһ??????
    x=Obs.co2_tt[which(Obs.co2_tt$site %in% Simsite),"CO2"]
    
    fitpars22 = pars#
    fitpars22[colnames(Fit.pools0$par)] = Fit.pools0$par[parsgrp[tt1],]
    tout = Model.pools(fitpars22,inputdata[which(inputdata$site %in% Simsite),],expT = expT,pc=pc,Pmic=Pmic)  #???ڷ???MIC??LMWC??CO2
    y=tout[vector_gap(1+CO2_30,simd), "CO2"]
    
    print(paste("R2_CO2=",(1-summary(lm(x~y))$r.squared)*2/10*5))
    
    #̼??
    xp = Obs.pools[which(Obs.pools$site %in% Simsite),c("LMWC","MIC")]
    yp = tout[vector_gap(1+ML30,simd), c("LMWC", "MIC")]
    print(colMeans(abs(xp-yp)/xp)*1/10*5)
    #?????ʶ?values
    print(Fit.pools0$value[parsgrp[tt1],])
  }
  
}


#ͨ��????----
#1??????20????????һ????20??????3???¶?????λ??,???????????кϲ?
site.flux = list() #Ԫ??Ϊ10???㣬ÿ????2??data.frame??????0??????1????Ϊ?????¶?*30,??Ϊ?????ϲ?

#????????
{
  flux.cul = c(FALSE,TRUE)[2] #ͨ��?????ۼƻ???ÿ????
  out.allflux = c(FALSE,TRUE)[2] #?Ƿ?????????ͨ��
  clab = switch(1,FALSE,c("Cin","Cinoa","Cad")[2]) #????C???? #??Ϊȫ?ֱ?��
  clab.outvar = c("Cpools_Tco2", "all")[2] #????C??????Щ??��
  sav.init =  c(FALSE,TRUE)[2] #?????????Ƿ?Ҫ??????ʼֵ
  allflux.name = c("f_Lit_PO", "f_Lit_LM", "f_Exud_LM", "f_Oa_MA", "f_Oa_LM",
                   "f_PO_LM", "f_PO_MA", "f_MA_LMm", "f_MA_LMe", "acid_des","f_MA_LM","f_LM_MA", 
                   "f_MA_MD", "f_MD_MA", "f_LM_MB",  "f_MA_maint_co2", "f_MA_growth_co2",
                   "f_MA_co2", "f_MD_co2", "f_MA_growth",
                   "f_MB_LM", "f_MB_MA", "f_MB_turn","f_LM_leach")
  #ֻ??????̼????Tco2
  nfixvar = length(fixCp)
  fixvar = fixCp
  if(out.allflux){
    Cout.nflux = nfixvar+length(allflux.name)
    colnms = c(fixvar,allflux.name)
  }else{
    Cout.nflux = nfixvar
    colnms =  fixvar
  }
  #????????̼??
  if(clab[1]!=FALSE){
    MAOM_c130 = 0 #????????????ʼֵ
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
}
cond_stode = TRUE
for(i in 1:length(string)){
  # i=3
  Fit.pools0 = list()
  si = sort(sitei)[i]
  Simsite = c("??????1","????","??????","??ɽ????","????","????",
              "????","?ݵ???","??????2","??????")[si]
  print(Simsite)
  nSimsamep = 2
  #ѡ????????
  fitpars = unlist(pars)
  
  Fit.pp = Fit.pools_iter[[Simsite]] #Fit.pools[[Simsite]] #
  
  # Fit.pp$value = Fit.pp$value/100
  #??ѡ??????valueС??1??
  res.par=c()
  for(nr in 1:nrow(Fit.pp$value)){
    if(length(which(Fit.pp$value[nr, ]<=1))==ncol(Fit.pp$value)){
      res.par = c(res.par,nr)
    }
  }
  
  Fit.pools0$par = Fit.pp$par[res.par,]
  if(is.null(colnames(Fit.pools0$par))){
    colnames(Fit.pools0$par) = fitparnam
  }
  
  Fit.pools0$value = Fit.pp$value[res.par,]
  #length(which(Fit.pools0$value[, 7]!=1)) #????Inf?Ĳ????ж?????
  #res.par = c(1:nrow(Fit.pools0$value))
  
  if(testMd == "_objCO2R2MARE" | testMd == "_fixbeta"| testMD2){
    R2.co2 = Fit.pools0$value[,1]
    # R2.co2.day = Fit.pools0$value[,2]
    # MARE.co2 = Fit.pools0$value[,3]
    MARE.lmwc = Fit.pools0$value[,2]
    MARE.mic = Fit.pools0$value[,3]
  }else{
    R2.co2 = Fit.pools0$value[,1]
    R2.co2.day = Fit.pools0$value[,2]
    MARE.co2 = Fit.pools0$value[,3]
    MARE.lmwc = Fit.pools0$value[,4]
    MARE.mic = Fit.pools0$value[,5]
  }
  
  
  if(cond_stode){
    # meanARE.SS = Fit.pools0$value[res.par,6]
    if(testMd == "_objCO2R2MARE" | testMd == "_fixbeta"| testMD2){
      meanARE.SS = Fit.pools0$value[,4]
      wtpar = c(2,1,1,1) #c(2,1,2,1,1,1)
      
      # wtpar = switch(Simsite, 
      #                "??????1"= c(2,1,3,0.5),
      #                "??????"= c(2,1,1,0.5),
      #                "????"= c(2,1,1,0.5),
      #                "????"= c(2,1,1,1),
      #                wtpar)
      wtpar = switch(Simsite, 
                     "??????1"= c(20,1,2,1),#c(2,1,3,0.5),
                     "??????"= c(2,1,1,0.5),
                     "????"= c(10,1,1,0.5),#c(2,1,1,0.5), #
                     "????"= c(10,1,1,0.5),
                     "????"= c(10,1,1,1),
                     "?ݵ???"= c(10,1,1,1),
                     "??????"= c(10,1,1,1),
                     wtpar)
      eff = sqrt(((1-R2.co2)^2 *wtpar[1] + (1-MARE.lmwc)^2*wtpar[2]+ 
                    (1-MARE.mic)^2*wtpar[3]+(1-meanARE.SS)^2 *wtpar[4])/sum(wtpar))
      
    }else{
      meanARE.SS = Fit.pools0$value[,6]
      wtpar = c(2,1,2,1,1,1) #c(2,1,2,1,1,1)
      wtpar = switch(Simsite, 
                     "??????1"= c(2,1,2,1,3,0.5),
                     "??????"= c(2,1,5,1,1,0.5),
                     "????"= c(2,1,5,1,1,0.5),
                     "????"= c(2,2,1,1,1,1),
                     wtpar)
      eff = sqrt(((1-R2.co2)^2 *wtpar[1] + (1-R2.co2.day)^2 *wtpar[2] 
                  + (1-MARE.co2)^2 *wtpar[3] 
                  + (1-MARE.lmwc)^2*wtpar[4]+ (1-MARE.mic)^2*wtpar[5]+(1-meanARE.SS)^2 *wtpar[6])/sum(wtpar))
      
    }
    
  }
  
  parsgrp = order(eff,decreasing = T)[1:npar] #eff?????????±?
  sort(eff,decreasing=T)[1:npar]
  #ģ?ⲻͬ?¶ȣ???ͬ????
  for(Ti in 1:length(expT)){
    
    #?????¶????ݿ?
    pardf0 = pardf1 = data.frame()
    for(pg in 1:npar){
      if(pg==1){
        optim.out10 = NULL
      }
      
      fitpars[colnames(Fit.pools0$par)] = Fit.pools0$par[parsgrp[pg],]
      
      #һ???¶Ȳ???0??????1
      optim.out0 = Model.pools(fitpars,inputdata[which(inputdata$site %in% Simsite),],expT = expT[Ti])#,pc=pc,Pmic=Pmic,flux = TRUE)  #???ڷ???MIC??LMWC??CO2
      #?ϲ?30??
      rowbreaks = which(optim.out0$time %in% c(1,30))
      
      if(pg==1){
        pardf0 = optim.out0[rowbreaks[1]:rowbreaks[2],]
        pardf1 = optim.out0[rowbreaks[3]:rowbreaks[4],]
      }else{
        pardf0 = cbind(pardf0,optim.out0[rowbreaks[1]:rowbreaks[2],])
        pardf1 = cbind(pardf1,optim.out0[rowbreaks[3]:rowbreaks[4],])
      }
      
    } #????ѭ??????
    
    #3???¶Ȱ??кϲ?
    site.flux[[Simsite]]$Oacid0 = rbind(site.flux[[Simsite]]$Oacid0,pardf0)
    site.flux[[Simsite]]$Oacid1 = rbind(site.flux[[Simsite]]$Oacid1,pardf1)
    
  }#?¶?ѭ??????
}

round(colMeans(pardf1[,-1]))

length(site.flux)
colnames(site.flux[[Simsite]]$Oacid0)

#ͨ��??????ƽ??ֵ
{
  dim(site.flux[[Simsite]]$Oacid0 )
  colnames(site.flux[[Simsite]]$Oacid0)
  
  #30???ۼ?ͨ��??ѡ????30??
  CulFlux30d_Mean = list()
  
  for(i in 1:length(string)){
    # i=3
    Fit.pools0 = list()
    si = sort(sitei)[i]
    Simsite = c("??????1","????","??????","??ɽ????","????","????",
                "????","?ݵ???","??????2","??????")[si]
    for(Oani in c("Oacid0","Oacid1")){
      CulFlux30d = site.flux[[Simsite]][[Oani]][nrow(site.flux[[Simsite]][[Oani]])/c(3,2,1),]
      dim(CulFlux30d)
      
      #??20???Ż?????ֵ
      unique(colnames(site.flux[[Simsite]][[Oani]]))
      CulFlux30d_MeanSD = mpars_meansd(CulFlux30d, var0 = grep("site",unique(colnames(site.flux[[Simsite]][[Oani]])),invert = T,value = T),
                                       outsd=FALSE ,onlyMean_sd = TRUE)
      #??Mean
      CulFlux30d_Mean[[Oani]] = rbind(CulFlux30d_Mean[[Oani]],round(CulFlux30d_MeanSD[,grep("_mean", colnames(CulFlux30d_MeanSD))],2))
      
    }
    
  }
  
}
CulFlux30d_Mean[["Oacid0"]][,c("MA_mean","MD_mean","f_MA_co2_mean","f_MD_co2_mean")]
CulFlux30d_Mean[["Oacid1"]][,c("MA_mean","MD_mean","f_MA_co2_mean","f_MD_co2_mean")]

oaRmad0 = round(CulFlux30d_Mean[["Oacid0"]][,c("f_MA_co2_mean")]/CulFlux30d_Mean[["Oacid0"]][,c("f_MD_co2_mean")])
oaRmad1 = round(CulFlux30d_Mean[["Oacid1"]][,c("f_MA_co2_mean")]/CulFlux30d_Mean[["Oacid1"]][,c("f_MD_co2_mean")])

mean(oaRmad0)
mean(oaRmad1)

#?ܽ????????ۼ?10?죩,??10???㰴?кϲ?Ϊ???ݿ?
out = list()
sday = 10 #ʮ?컹??30??
tdescols0 = which(colnames(site.flux[[Simsite]]$Oacid0) %in% "MA_LM") #????��?л????????????ܽ?????
oacidcols = which(colnames(site.flux[[Simsite]]$Oacid0) %in% "MA_LMd") #?л???????????
enycols = which(colnames(site.flux[[Simsite]]$Oacid0) %in% "MA_LMe") 
phycols = which(colnames(site.flux[[Simsite]]$Oacid0) %in% "MA_LM.M") 
co2cols = which(colnames(site.flux[[Simsite]]$Oacid0) %in% "CO2") 

for(i in 1:10){
  Simsite = c("??????1","????","??????","??ɽ????","????","????",
              "????","?ݵ???","??????2","??????")[i]
  print(Simsite)
  
  if(i==1){
    out$tdesp0 = site.flux[[Simsite]]$Oacid0[sday,tdescols0] + site.flux[[Simsite]]$Oacid0[sday,oacidcols]
    out$tdesp1 = site.flux[[Simsite]]$Oacid1[sday,tdescols0] + site.flux[[Simsite]]$Oacid1[sday,oacidcols]
    
    out$eny0 = site.flux[[Simsite]]$Oacid0[sday,enycols]
    out$phy0 = site.flux[[Simsite]]$Oacid0[sday,phycols]
    out$eny1 = site.flux[[Simsite]]$Oacid1[sday,enycols]
    out$phy1 = site.flux[[Simsite]]$Oacid1[sday,phycols]
    out$acid1 = site.flux[[Simsite]]$Oacid1[sday,oacidcols]
    
    out$co20 = site.flux[[Simsite]]$Oacid0[sday,co2cols]
    out$co21 = site.flux[[Simsite]]$Oacid1[sday,co2cols]
  }else{
    out$tdesp0 = rbind(out$tdesp0,site.flux[[Simsite]]$Oacid0[sday,tdescols0] + site.flux[[Simsite]]$Oacid0[sday,oacidcols])
    out$tdesp1 = rbind(out$tdesp1,site.flux[[Simsite]]$Oacid1[sday,tdescols0] + site.flux[[Simsite]]$Oacid1[sday,oacidcols])
    
    out$eny0 = rbind(out$eny0,site.flux[[Simsite]]$Oacid0[sday,enycols])
    out$phy0 = rbind(out$phy0,site.flux[[Simsite]]$Oacid0[sday,phycols])
    out$eny1 = rbind(out$eny1,site.flux[[Simsite]]$Oacid1[sday,enycols])
    out$phy1 = rbind(out$phy1,site.flux[[Simsite]]$Oacid1[sday,phycols])
    out$acid1 = rbind(out$acid1,site.flux[[Simsite]]$Oacid1[sday,oacidcols])
    
    out$co20 = rbind(out$co20,site.flux[[Simsite]]$Oacid0[sday,co2cols])
    out$co21 = rbind(out$co21,site.flux[[Simsite]]$Oacid1[sday,co2cols])
  }
  
}
out$tdesm = out$tdesp1-out$tdesp0      #?ܽ?????֮??
out$tco2m = out$co21 - out$co20

percfun(out$tdesp0) #?ܽ?????
percfun(out$tdesp1) 
percfun(out$tdesm)
#??????????ռ?ܽ?????????
percfun(out$eny0)
percfun(out$eny1)
percfun(out$phy0)
percfun(out$phy1)
percfun(out$acid1)
#ռ??
percfun(out$eny0/out$tdesp0)
percfun(out$eny1/out$tdesp1)
percfun(out$acid1/out$tdesp1)

#co2
percfun(out$tco2m)
percfun(out$co20)
percfun(out$co21)
