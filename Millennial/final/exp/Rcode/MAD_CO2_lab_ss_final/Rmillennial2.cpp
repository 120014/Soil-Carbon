#include <iostream>
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
using namespace std;

struct StruParsSoil{
  //Calibrated Parameters
  float k_lit, f_resp;   //add_litter_pool
  float kaff_pl,  Vpl0;
  float kaff_ml,  Vml0;
  float kaff_lb,  Vlb0;
  float param_pi, param_fpl,  param_pb, r0, Ma, beta;
  float param_p2, kaff_des,   param_pc;
  float cue_t,    cue_ref,    acue;
  float rate_bd,  rate_Kbd;
  //Input soil property parameters
  float depth, param_pH, param_bulkd, param_clay, param_claysilt, porosity, param_em, rate_leach;

  // new add
  float param_pbl;
};
struct StrsoilC{
  double SOM, LIT,  POM, LMWC, MIC, MAOM, MA, MD, Tco2;//add_litter_pool
};
struct Strflux{
  float f_Lit_PO, f_Lit_LM, f_Exud_LM, f_Oa_MA, f_Oa_LM;
  float f_PO_LM, f_PO_MA, f_MA_LMm, f_MA_LMe, acid_des,f_MA_LM,f_LM_MA; 
  float f_MA_MD, f_MD_MA, f_LM_MB,  f_MA_maint_co2, f_MA_growth_co2;
  float f_MA_co2, f_MD_co2, f_MA_growth;
  float f_MB_LM, f_MB_MA, f_MB_turn,f_LM_leach;
};

// Initialize struct parameter
StruParsSoil parsInitial(List parsSoiList) {
  StruParsSoil parsSoil;
  
  parsSoil.k_lit = as<float>(parsSoiList["k_lit"]);
  parsSoil.f_resp = as<float>(parsSoiList["f_resp"]);
  
  parsSoil.kaff_pl = as<float>(parsSoiList["kaff_pl"]);
  parsSoil.Vpl0 = as<float>(parsSoiList["Vpl0"]);
  
  parsSoil.kaff_ml = as<float>(parsSoiList["kaff_ml"]);
  parsSoil.Vml0 = as<float>(parsSoiList["Vml0"]);
  
  parsSoil.kaff_lb = as<float>(parsSoiList["kaff_lb"]);
  parsSoil.Vlb0 = as<float>(parsSoiList["Vlb0"]);
  
  parsSoil.param_pi = as<float>(parsSoiList["param_pi"]);
  parsSoil.param_fpl = as<float>(parsSoiList["param_fpl"]);
  parsSoil.param_pb = as<float>(parsSoiList["param_pb"]);
  parsSoil.r0 = as<float>(parsSoiList["r0"]);
  parsSoil.Ma = as<float>(parsSoiList["Ma"]);
  parsSoil.beta = as<float>(parsSoiList["beta"]);
  
  parsSoil.param_p2 = as<float>(parsSoiList["param_p2"]);
  parsSoil.kaff_des = as<float>(parsSoiList["kaff_des"]);
  parsSoil.param_pc = as<float>(parsSoiList["param_pc"]);
  parsSoil.cue_t = as<float>(parsSoiList["cue_t"]);
  parsSoil.cue_ref = as<float>(parsSoiList["cue_ref"]);
  parsSoil.acue = as<float>(parsSoiList["acue"]);
  
  parsSoil.rate_bd = as<float>(parsSoiList["rate_bd"]);
  parsSoil.rate_Kbd = as<float>(parsSoiList["rate_Kbd"]);
  
  //site input data
  parsSoil.depth = as<float>(parsSoiList["depth"]); 
  parsSoil.param_pH = as<float>(parsSoiList["param_pH"]); 
  parsSoil.param_bulkd = as<float>(parsSoiList["param_bulkd"]); 
  parsSoil.param_clay = as<float>(parsSoiList["param_clay"]); 
  parsSoil.param_claysilt = as<float>(parsSoiList["param_claysilt"]); 
  parsSoil.porosity = as<float>(parsSoiList["porosity"]); 
  parsSoil.param_em = as<float>(parsSoiList["param_em"]); 
  parsSoil.rate_leach = as<float>(parsSoiList["rate_leach"]); 

  // new add
  parsSoil.param_pbl = as<float>(parsSoiList["param_pbl"]);  

  return parsSoil;
}


//Compute effect factor of soil temperature with Arrhenius 
double tp_scalar(string sCase, double T0, double Tref) {
  T0 = T0 + 273.15; //(K)
  Tref = Tref + 273.15; //(K)
  double gas_const = 8.31446;
  double Ea = 0;
  
  if (sCase == "MR") Ea = 20e3;
  else if (sCase == "LIG") Ea = 53e3;
  else if (sCase == "CEL") Ea = 36.3e3;
  else if (sCase == "POM") Ea = 66e3;//63909;
  else if (sCase == "LMWC") Ea = 57865;
  else if (sCase == "MAOM") Ea = 67e3;
  else if (sCase == "MD") Ea = 47e3;
  else if (sCase == "KM") Ea = 30e3;
  else if (sCase == "LIT") Ea = 62e3;          //GPT建议62e3
  
  return exp(-Ea / gas_const * (1 / T0 - 1 / Tref));
}

//Compute soil matrix potential with van-Genuchten equation, return SWP unit is [MPa]
double fSWC2SWP(double SWC0, double SWCsat) {
  const double const_cm2MPa = 98e-6;  // 1cm water column
  const float rlim = 1.01;
  
  // VG parameters
  float SWCres = 0.067;
  // SWCsat input
  float alpha = 0.021;
  float n = 1.61;
  float m = 1 - 1 / n;
  
  float SWC = (SWC0 <= SWCres * rlim) ? SWCres * rlim : SWC0;
  float fSWC2SWP;
  
  if (SWC < SWCsat) {
    double eff_sat = (SWC - SWCres) / (SWCsat - SWCres);  // Effective saturation
    fSWC2SWP = pow(pow(1 / eff_sat, 1 / m) - 1, 1 / n) / alpha;
    fSWC2SWP = -1 * fSWC2SWP * const_cm2MPa;  // to MPa
  } else {
    fSWC2SWP = 0;
  }
  return fSWC2SWP;
}

//Steady-state judgement function // [[Rcpp::export]]
bool sscondfun(NumericVector ss_soc0) {
  int n = ss_soc0.size();
  if (is_true(all(ss_soc0 > 0))) {
    //Calculate the difference between adjacent elements
    NumericVector diff_soc0(n - 1);
    for (int i = 0; i < n - 1; ++i) {
      diff_soc0[i] = std::abs(ss_soc0[i + 1] - ss_soc0[i]);
    }
    
    // abs(diff) < 1
    LogicalVector sscond2 = diff_soc0 < 1;
    
    // abs(diff) / ss_soc0[i] * 100 < 0.1%
    LogicalVector sscond3 = (diff_soc0 / ss_soc0[Range(0, n - 2)] * 100) < 0.1;
    
    // 如果所有条件都满足，则返回 true
    if (is_true(all(sscond2 & sscond3))) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

//Trace CO2 from organic acid result in desorption
void CalLabC(StrsoilC* c13,StrsoilC c0, StrsoilC c1, Strflux flux,float CUE,float param_pb){
  
  Strflux lab;
  
  c13->LMWC =  c13->LMWC + flux.acid_des;
  if(c0.MAOM > 0. ){
    lab.acid_des = fmin(flux.acid_des * c13->MAOM/c0.MAOM, c13->MAOM);
  }else{
    lab.acid_des = 0;
  }
  c13->MAOM =  c13->MAOM - lab.acid_des;
  
  if(c1.MAOM > 0. ){
    lab.f_MA_LMe  = fmin(flux.f_MA_LMe  * c13->MAOM/c1.MAOM, c13->MAOM);
    lab.f_MA_LMm =  fmin(flux.f_MA_LMm  * c13->MAOM/c1.MAOM, c13->MAOM);
  }else{
    lab.f_MA_LMe = 0.f; 
    lab.f_MA_LMm = 0.f; 
  }
  
  if(c0.MA > 0. ){
    lab.f_MA_MD = fmin(flux.f_MA_MD * c13->MA/c0.MA,  c13->MA);
  }else{
    lab.f_MA_MD = 0.; 
  }
  
  if(c0.MD > 0. ){
    lab.f_MD_MA = fmin(flux.f_MD_MA * c13->MD/c0.MD, c13->MD);
  }else{
    lab.f_MD_MA = 0.; 
  }
  c13->MA = c13->MA + lab.f_MD_MA - lab.f_MA_MD;
  c13->MD = c13->MD + lab.f_MA_MD - lab.f_MD_MA;
  
  double LMWC = (c1.LMWC + flux.f_LM_MB);
  if(LMWC > 0. && c13->LMWC > 0.){
    lab.f_LM_MB = fmin(flux.f_LM_MB * c13->LMWC/LMWC, c13->LMWC);
    lab.f_MA_growth = lab.f_LM_MB * CUE; //flux.f_MA_growth * c13->LMWC/LMWC;
    lab.f_MA_co2 = lab.f_LM_MB * (1. - CUE); //flux.f_MA_co2  * c13->LMWC/LMWC;
    // f_MA_growth_co2 = f_MA_growth_co2 * LMWC_c13/LMWC;
    // f_MA_maint_co2 = f_MA_maint_co2 * LMWC_c13/LMWC;
  }else{
    lab.f_LM_MB = 0.; 
    lab.f_MA_growth  = 0. ;
    lab.f_MA_co2 = 0. ;
    // f_MA_growth_co2_c13  = 0.;
    // f_MA_maint_co2_c13 = 0.;
  }
  c13->LMWC = c13->LMWC - lab.f_LM_MB;
  
  
  if(c1.MD > 0. ){
    lab.f_MD_co2 = fmin(flux.f_MD_co2  * c13->MD/c1.MD, c13->MD);
  }else{
    lab.f_MD_co2 = 0.; 
  }
  
  if(c1.LMWC > 0. && c13->LMWC >0.){
    lab.f_LM_MA =  fmin(flux.f_LM_MA  * c13->LMWC/c1.LMWC, c13->LMWC);
    lab.f_LM_leach =  fmin(flux.f_LM_leach  * c13->LMWC/c1.LMWC, c13->LMWC);
  }else{
    lab.f_LM_MA = 0.;
    lab.f_LM_leach = 0.; 
  }
  
  if(c1.MA > 0. ){
    lab.f_MB_turn = fmin(flux.f_MB_turn * c13->MA/c1.MA, c13->MA); 
    lab.f_MB_LM = lab.f_MB_turn * (1. - param_pb);//flux.f_MB_LM * c13->MA/c1.MA;
    lab.f_MB_MA = lab.f_MB_turn * param_pb;//flux.f_MB_MA * c13->MA/c1.MA;
  }else{
    lab.f_MB_turn = 0.; 
    lab.f_MB_LM = 0.; 
    lab.f_MB_MA = 0.; 
  }
  
  if((lab.f_MA_LMe + lab.f_MA_LMm) > c13->MAOM){
    lab.f_MA_LMe = lab.f_MA_LMm = c13->MAOM/2;
  }
  lab.f_MA_LM = lab.acid_des + lab.f_MA_LMe + lab.f_MA_LMm;
  
  if((lab.f_LM_leach + lab.f_LM_MA)>c13->LMWC){
    lab.f_LM_leach = lab.f_LM_MA = c13->LMWC/2;
  }
  
  float Tco2_c13 = lab.f_MA_co2 + lab.f_MD_co2;
  
  float dLMWC = lab.f_MB_LM + (lab.f_MA_LMe + lab.f_MA_LMm) - (lab.f_LM_leach + lab.f_LM_MA);
  float dMA =  lab.f_MA_growth - lab.f_MB_turn ; //+ f_MD_MA - f_MA_MD
  float dMD =  - lab.f_MD_co2;   //f_MA_MD - f_MD_MA
  float dMAOM = lab.f_LM_MA + lab.f_MB_MA  - (lab.f_MA_LMe + lab.f_MA_LMm);
  
  // Update C pools
  c13->POM = 0.;
  c13->LMWC = c13->LMWC + dLMWC;
  c13->MA = c13->MA + dMA;
  c13->MD = c13->MD + dMD;
  c13->MAOM = c13->MAOM + dMAOM;
  c13->MIC = c13->MA + c13->MD;
  c13->SOM = c13->POM + c13->LMWC + c13->MIC + c13->MAOM;
  c13->Tco2 = c13->Tco2 + Tco2_c13;
  return;
  
}

// Calculate the average for the specified range
vector<vector<float>> ssCpools_loops_mean(const vector<vector<float>>& data, int start_row, int end_row, int start_col, int end_col) {
  int nrows = end_row - start_row + 1;  
  int ncols = end_col - start_col + 1;
  vector<vector<float>> averages(1, vector<float>(ncols)); //the average of each column
  
  for (int j = start_col; j <= end_col; j++) {
    float sum = 0.0;
    for (int i = start_row; i <= end_row; i++) {
      sum += data[i][j];
    }
    averages[0][j - start_col] = sum / nrows;  
  }
  return averages;  
}



//soil C turnover function
// [[Rcpp::export]] //must!
vector<vector<double>> derivs_V2_MM_AD_CO2(int sstime, List parsSoiList, List CpoolsList,NumericVector forc_st,
                                           NumericVector forc_sw,NumericVector forc_lit,NumericVector forc_acid, 
                                           bool outflux, int nloop,bool clab,List C13poolsList) {
  //output
  int ncpools = 8, nTco2 = 1, curtime = 0, nyear,nyear_row, nrep,nflux=0;//SOM,Cpools,Tco2
  if(nloop==0 || nloop>500){//1constant dynamic input,test,100
    nyear = sstime/365;
    nrep = 365;
  }else{
    nyear = sstime/365/nloop;
    nrep = nyear*365;
  }

  if(outflux){
    nflux = 24;// no Tco2
  }
  // cout <<"nyear="<<nyear<<endl;
  if(clab ){//exclude parameters calibrate of experiment && nrep != 365
    ncpools = ncpools*2 ; 
    nTco2 = nTco2 + 1;// adding Tco2 col
  }
  nyear_row = fmin(nyear,100);
  vector<vector<double>> out(nyear_row, vector<double>(ncpools+nTco2+1+nflux));//year,co2
  vector<vector<double>> ayear(365, vector<double>(ncpools+nTco2+nflux));//co2,buffer variable
  
  //constant parameters
  float param_p1=0.12f;//fexud = 0.1f, fexud_oa = 0.5f, 
  float lambda = 0.00021, kamin = 0.2;
  int Tref = 20, na_value = -9999, maRco2=0; 
  StruParsSoil pars;
  pars = parsInitial(parsSoiList);
  double LIT = as<double>(CpoolsList["LIT"]);   //add_litter_pool
  double POM = as<double>(CpoolsList["POM"]); 
  double LMWC = as<double>(CpoolsList["LMWC"]); 
  double MIC = as<double>(CpoolsList["MIC"]); 
  double MAOM = as<double>(CpoolsList["MAOM"]); 
  double MA = as<double>(CpoolsList["MA"]); 
  double MD = as<double>(CpoolsList["MD"]); 
  NumericVector forc_exud = forc_acid; // root exudates exclude organic acid
  
  //parameters about steady state solving
  float sumSOM = 0.f;
  NumericVector ssSOM(10);//Initial default is 0. 
  //the change of SOM is less than 0.1gC/m2 for 10 years
  int ssi = 0, loopth = 1,ssy_aloop;
  bool boolss;
  if(nrep==365){
    ssy_aloop=1;
  }else{
    ssy_aloop=nyear;
  }
  vector<vector<double>> outss(1, vector<double>(ncpools+1));//ssgx!
  vector<vector<double>> ssyears(ssy_aloop, vector<double>(ncpools));// length is years of a loop
  // cout <<"begin for t ="<<nyear<<endl;
  
  //label CO2 from organic acid desorption
  StrsoilC* c13 = new StrsoilC;
  StrsoilC c0, c1;
  
  for(int t=0; t<sstime; ++t){
    c0 = {0, LIT, POM, LMWC, MIC, MAOM, MA, MD, 0};
    
    int fdi = t%nrep;
    // Equation 10  DOC adsorption rate (/d)
    float kaff_lm = exp(-param_p1 * pars.param_pH - pars.param_p2) * pars.kaff_des;
 
    // Equation 11 
    float param_qmax = pars.depth * pars.param_bulkd * pars.param_claysilt * pars.param_pc;
    // if(t>ttime){
    //   cout << "param_qmax= " << param_qmax <<"  MAOM=" <<MAOM<<endl;
    // }
    if (param_qmax < MAOM) {
      cout << "ERROR EXIT: Qmax is less than MAOM!" << endl;
      break;
    }

    float forc_sw0 = fmin(forc_sw[fdi], pars.porosity); // Hydrological properties
    float scalar_wd = pow(forc_sw0 / pars.porosity, 0.5);// Equation 4: scalar_wd
    
    // Equation 15: scalar_wb
    float SWP = fSWC2SWP(forc_sw0, pars.porosity); // Placeholder for actual SWP calculation
    float matpot = abs(SWP * 1e3);
    float scalar_wb = exp(lambda * -matpot) * (kamin + (1 - kamin) * pow((pars.porosity - forc_sw0) / pars.porosity, 0.5)) * scalar_wd;
    // cout << "scalar_wb = " << scalar_wb << endl;
    
    // C input
    float f_Lit_LIT = forc_lit[fdi];

    float k = pars.k_lit * tp_scalar("LIT", forc_st[fdi], Tref) *scalar_wd;

    float f_LIT; //add_litter_pool

    if (LIT > 0.f) {

      f_LIT = fmax(LIT * k, 0.f);
      f_LIT = fmin(LIT, f_LIT);
    }else{
      f_LIT = 0.f;
    }
    
    float f_LIT_co2 = f_LIT * pars.f_resp;               // respired immediately to atmosphere
    float f_LIT_SOC = f_LIT * (1.0f - pars.f_resp);
    
    // // 添加调试输出
    // std::cout << "C++ - LIT分解调试 - day: " << t
    //           << " LIT: " << LIT
    //           << " f_LIT: " << f_LIT
    //           << " forc_lit=" << forc_lit[fdi]
    //           << std::endl;
    
    
    
    float f_Lit_PO = f_LIT_SOC * pars.param_pi;
    float f_Lit_LM = f_LIT_SOC * (1. - pars.param_pi);
    float f_Exud_LM = forc_exud[fdi];

    float f_Oa_MA = forc_acid[fdi] * pars.param_em;
    float f_Oa_LM = forc_acid[fdi] * (1. - pars.param_em) ;
    // cout << "f_Oa_MA = " << f_Oa_MA << endl;
    
    // Decomposition rates
    float km_TPscalar = tp_scalar("KM", forc_st[fdi], Tref);
    float vmax_pl = pars.Vpl0 * tp_scalar("POM", forc_st[fdi], Tref) * scalar_wd;
//    float vmax_pl = pars.Vpl0 * tp_scalar("POM", forc_st[fdi], Tref) *  pow(forc_sw0 / pars.porosity, 1.0);
    float kaff_pl = pars.kaff_pl * km_TPscalar;
    // cout << "km_TPscalar = " << km_TPscalar << "  vmax_pl = " << vmax_pl << "  kaff_pl = " << kaff_pl << endl;
    
    // Example decomposition of POM to LMWC
    float f_PO;
    if (POM > 0.f && MA > 0.f) {
      f_PO = fmax(vmax_pl * POM * MA / (kaff_pl + MA), 0.f);
      f_PO = fmin(POM, f_PO);
    }else{
      f_PO = 0.f;
    }
    float f_PO_LM = f_PO * pars.param_fpl;
    float f_PO_MA = f_PO * (1. - pars.param_fpl);
    // cout << "f_PO_LM = " << f_PO_LM << "  f_PO_MA = " << f_PO_MA << endl;
    
    //MAOC finish desorption within a few minutes after organic acid input
    float unit_kgTom2 = pars.param_bulkd * pars.depth *1e-3; //(mgC/kg soil) to (gC/m2)
    float acid_des;
    if( MAOM > 0.f ){
      acid_des = fmax((pow(MAOM/unit_kgTom2 *1e-3,0.04f) *pow(forc_acid[fdi]/unit_kgTom2, 0.87f)
                    * pow(pars.param_clay,-0.24f) * pow(pars.param_pH,-2.15f) * 34.0) * unit_kgTom2,0.f);
      acid_des = fmin(acid_des,MAOM);
    }else{
      acid_des = 0.f;
    }
    
    LMWC = LMWC  + f_Oa_LM + acid_des;
    MAOM = MAOM  + f_Oa_MA - acid_des;
    // if(fdi>ttime){
    //   cout << "acid_des = " << acid_des << "  LMWC = " << LMWC << "  MAOM = " << MAOM <<endl;
    // }
    
    //the enzymatic decompsition of MAOM to LMWC
    // vmax_ml = max(alpha_ml * exp(-eact_ml / (gas_const * (forc_st(step.num) + 273.15))),0)
    float vmax_ml = pars.Vml0 * tp_scalar("MAOM", forc_st[fdi], Tref) * scalar_wd;
    float kaff_ml = pars.kaff_ml * km_TPscalar;
    float f_MA_LMe;
    if(MAOM>0.f && MA>0.f){
      f_MA_LMe = fmax(vmax_ml * MAOM * MA / (kaff_ml + MA),0.f);
      f_MA_LMe = fmin(MAOM, f_MA_LMe);
    }else{
      f_MA_LMe=0.f;
    }
    // cout << "f_MA_LMe = " << f_MA_LMe << "  MA = " << MA << endl;
    
    //Equation 12: physical desorption of MAOM to LMWC
    float f_MA_LMm;
    if(MAOM > 0.f){
      f_MA_LMm = fmax(pars.kaff_des * MAOM / param_qmax,0.f);
      f_MA_LMm = fmin(MAOM, f_MA_LMm);
    }else{
      f_MA_LMm = 0.f;
    }
    // cout << "f_MA_LMm = " << f_MA_LMm  << endl;
    
    //Equation 14: respiration of MA and MD and CUE process
    float Vg = pars.Vlb0 * tp_scalar("LMWC", forc_st[fdi], Tref)* scalar_wb; // * scalar_wb //the growth respiration of MA is similar to Vd in MEND model
    float Vm = pars.Vlb0 * pars.Ma/(1-pars.Ma) * tp_scalar("MR", forc_st[fdi], Tref)* scalar_wb; //the maintenance respiration of MA：α/(1-α) means the propotion of maintenance to growth
    float VmD = Vm * pars.beta; //the maintenance respiration of MD
    
    float CUE_max = 0.95f , CUE_min = 0.01f;
    float CUE = pars.cue_ref - pars.cue_t * (forc_st[fdi] - Tref) ;
    CUE = CUE - pars.acue * forc_acid[fdi];     //CUE decline after organic acid input
    // if(t>ttime){
    //   cout << "t: "<<t<<"CUE = " << CUE<<"POM="<<POM<<" LMWC="<<LMWC<<" MIC="<<(MA + MD)<<" MAOM="<<MAOM<<" MA="<<MA<<" MD="<<MD<<endl;
    // }
    if(CUE<CUE_min || CUE>CUE_max){
      cout << "CUE ERROR EXIT: CUE = " << CUE << "t = " << t << "fdi = " << fdi << endl;
      break;
    }
    

    // Equation 14.1: MA <--> MD,
    float VmA2D = Vm; // * tp_scalar("MR") #* scalar_wd #* wp_scalar_low
    float VmD2A = Vm; // * tp_scalar("MR") #* scalar_wd #* wp_scalar
    float vmax_lb = (Vg + Vm) / CUE;
    float kaff_lb = pars.kaff_lb * km_TPscalar;
    float sat = LMWC/(kaff_lb + LMWC); //Wang 2015
    // float sat = pow(f_MA_growth/(f_MA_growth + f_MA_maint_co2 + f_MB_turn), 2); //Huang2022
    float f_MA_MD, f_MD_MA;
    if(LMWC >0.f && MA > 0.f){
      f_MA_MD = fmax((1-sat) * VmA2D * MA,0.f);
      f_MA_MD = fmin(MA, f_MA_MD);
    }else{
      f_MA_MD = 0.f;
    }
    if(LMWC >0.f && MD>0.f){
      f_MD_MA = fmax(sat * VmD2A * MD,0.f);
      f_MD_MA = fmin(MD, f_MD_MA);
    }else{
      f_MD_MA = 0;
    }
    // if(t>ttime){
    //   cout << "f_MA_MD = " << f_MA_MD <<" MA = "<< MA<< "  f_MD_MA = "<< f_MD_MA <<" MD = "<< MD << endl;
    //   cout << "POM="<<POM<<" LMWC="<<LMWC<<" MIC="<<(MA + MD)<<" MAOM="<<MAOM<<" MA="<<MA<<" MD="<<MD<<endl;
    // }
    MA = MA + f_MD_MA - f_MA_MD; //MA and MD update immediately，because of reacting so quickly！
    MD = MD + f_MA_MD - f_MD_MA;
    
    //Equation 13: MA utilizes LMWC for growth and respiration
    float f_LM_MB;
    if(LMWC > 0.f && MA >0.f){
      f_LM_MB = fmax(vmax_lb * MA * LMWC / (kaff_lb + LMWC),0.f);
      f_LM_MB = fmin(LMWC, f_LM_MB);
    }else{
      f_LM_MB=0.f;
    }
    // cout << "f_LM_MB = " << f_LM_MB  << endl;
    
    //Equation 21: MIC -> atmosphere
    float f_MA_co2, f_MA_growth_co2 = 0, f_MA_maint_co2=0;
    if(maRco2 == 0){
      f_MA_co2 = f_LM_MB * (1.f - CUE); //method1
    }else{
      //Equation 22: microbial growth flux, but is not used in mass balance
      f_MA_growth_co2 = fmax( Vg * (1/CUE - 1) * MA * LMWC / (kaff_lb + LMWC),0.f);
      f_MA_maint_co2 = fmax( Vm * (1/CUE - 1) * MA * LMWC / (kaff_lb + LMWC),0.f);
      f_MA_co2 = f_MA_growth_co2 + f_MA_maint_co2;
    }
    // MA = MA + f_MA_growth
    LMWC = LMWC - f_LM_MB; //Update LMWC immediately for adsorption calculations
    float f_MA_growth = f_LM_MB * CUE; 
    // cout << "f_MA_co2 = " << f_MA_co2 << "  f_MA_growth = " << f_MA_growth << "  LMWC = " << LMWC << endl;
    
    //the maintenance respiration of MD
    float f_MD_co2;
    if(MD > 0.f){
      f_MD_co2 = fmax(VmD * MD,0.f); //0.001 
      f_MD_co2 = fmin(MD, f_MD_co2);
    }else{
      f_MD_co2 = 0.f;
    }
    float Tco2 = f_MA_co2 + f_MD_co2 + f_LIT_co2;
    // cout << "f_MD_co2 = " << f_MD_co2 << "  Tco2 = " << Tco2 << endl;
    
    // Equation 8: LMWC -> out of system leaching
    float f_LM_leach;
    if(LMWC > 0.f){
      f_LM_leach = fmax(pars.rate_leach * scalar_wd * LMWC,0.f);
      f_LM_leach = fmin(LMWC, f_LM_leach);
    }else{
      f_LM_leach=0.f;
    }
    //Equation 9: LMWC -> MAOM
    float f_LM_MA;
    if(LMWC > 0.f && MAOM > 0.f){
      f_LM_MA = fmax(scalar_wd * kaff_lm * LMWC * (1 - MAOM / param_qmax),0.f);
      f_LM_MA = fmin(LMWC, f_LM_MA);
    }else{
      f_LM_MA=0.f;
    }
    
    // Equation 16: MA -> MAOM/LMWC 
    float rate_bd = fmax(0.f, pars.rate_bd + pars.rate_Kbd*(forc_st[fdi] - Tref));
    float f_MB_turn;
    if(MA > 0.f){
      //f_MB_turn = rate_bd * MIC^2.0; // Millennial
      f_MB_turn = fmax(rate_bd * MA,0.f); // Liucq
      f_MB_turn = fmin(MA, f_MB_turn);
    }else{
      f_MB_turn=0.f;
    }
//    float f_MB_LM = f_MB_turn * (1. - pars.param_pb);
//    float f_MB_MA = f_MB_turn * pars.param_pb;
    // cout << "f_LM_leach = " << f_LM_leach << "  f_LM_MA = " << f_LM_MA << "  f_MB_turn = " << f_MB_turn <<"f_MB_LM = " <<f_MB_LM<< endl;
    //MA death to POC
    float f_MB_LM = f_MB_turn * pars.param_pbl;
    float f_MB_MA = f_MB_turn * pars.param_pb;
    float f_MB_PO = f_MB_turn * (1. - pars.param_pb - pars.param_pbl);
    
    // flux balance
    if((f_MA_LMe + f_MA_LMm) > MAOM){
      f_MA_LMe = f_MA_LMm = MAOM/2;
    }
    float f_MA_LM = acid_des + f_MA_LMe + f_MA_LMm;
    
    if((f_LM_leach + f_LM_MA)>LMWC){
      f_LM_leach = f_LM_MA = LMWC/2;
    }
    
    //Update cpools
    //Flux each C pool
    float dLIT = f_Lit_LIT - f_LIT;  //add_litter_pool
    float dPOM = f_Lit_PO - f_PO + f_MB_PO;
    float dLMWC = f_Lit_LM + f_Exud_LM + f_PO_LM + f_MB_LM + (f_MA_LMe + f_MA_LMm) - f_LM_leach  - f_LM_MA;
    float dMA =  f_MA_growth - f_MB_turn ; //+ f_MD_MA - f_MA_MD
    float dMD =  - f_MD_co2;   //f_MA_MD - f_MD_MA
    float dMAOM = f_LM_MA + f_MB_MA + f_PO_MA - (f_MA_LMe + f_MA_LMm);
    // if(t>ttime){
    //   cout <<"t: "<<t<< "dPOM="<<dPOM<<"dLMWC="<<dLMWC<<"dMA="<<dMA<<"dMD="<<dMD<<"dMAOM="<<dMAOM<<endl;
    // }
    
    // Update C pools
    LIT = LIT + dLIT;//add_litter_pool
    POM = POM + dPOM;
    LMWC = LMWC + dLMWC;
    MA = MA + dMA;
    MD = MD + dMD;
    MIC = MA + MD;
    MAOM = MAOM + dMAOM;
    
    //Output 
    int ti = t%365;
    ayear[ti][0] = POM + LMWC + MAOM + (MA + MD);//SOM
    ayear[ti][1] = LIT;
    ayear[ti][2] = POM;
    ayear[ti][3] = LMWC;
    ayear[ti][4] = MIC;
    ayear[ti][5] = MAOM;
    ayear[ti][6] = MA;
    ayear[ti][7] = MD;
    if(ti==0){
      ayear[ti][8] = Tco2;
    }else{
      ayear[ti][8] = ayear[ti-1][8]+ Tco2;
    }
    if(outflux){
      if(ti==0){
        ayear[ti][9] = f_Lit_PO, ayear[ti][10] = f_Lit_LM;
        ayear[ti][11] = f_Exud_LM, ayear[ti][12] = f_Oa_MA;
        ayear[ti][13] = f_Oa_LM, ayear[ti][14] = f_PO_LM;
        ayear[ti][15] = f_PO_MA, ayear[ti][16] = f_MA_LMm;
        ayear[ti][17] = f_MA_LMe, ayear[ti][18] = acid_des;
        ayear[ti][19] = f_MA_LM, ayear[ti][20] = f_LM_MA;
        ayear[ti][21] = f_MA_MD, ayear[ti][22] = f_MD_MA;
        ayear[ti][23] = f_LM_MB, ayear[ti][24] = f_MA_maint_co2;
        ayear[ti][25] = f_MA_growth_co2, ayear[ti][26] = f_MA_co2;
        ayear[ti][27] = f_MD_co2, ayear[ti][28] = f_MA_growth;
        ayear[ti][29] = f_MB_LM, ayear[ti][30] = f_MB_MA;
        ayear[ti][31] = f_MB_turn, ayear[ti][32] = f_LM_leach;
      }else{
        ayear[ti][9] = ayear[ti-1][9]+ f_Lit_PO;
        ayear[ti][10] = ayear[ti-1][10]+ f_Lit_LM;
        ayear[ti][11] = ayear[ti-1][11]+ f_Exud_LM;
        ayear[ti][12] = ayear[ti-1][12]+ f_Oa_MA;
        ayear[ti][13] = ayear[ti-1][13]+ f_Oa_LM;
        ayear[ti][14] = ayear[ti-1][14]+ f_PO_LM;
        ayear[ti][15] = ayear[ti-1][15]+ f_PO_MA;
        ayear[ti][16] = ayear[ti-1][16]+ f_MA_LMm;
        ayear[ti][17] = ayear[ti-1][17]+ f_MA_LMe;
        ayear[ti][18] = ayear[ti-1][18]+ acid_des;
        ayear[ti][19] = ayear[ti-1][19]+ f_MA_LM;
        ayear[ti][20] = ayear[ti-1][20]+ f_LM_MA;
        ayear[ti][21] = ayear[ti-1][21]+ f_MA_MD;
        ayear[ti][22] = ayear[ti-1][22]+ f_MD_MA;
        ayear[ti][23] = ayear[ti-1][23]+ f_LM_MB;
        ayear[ti][24] = ayear[ti-1][24]+ f_MA_maint_co2;
        ayear[ti][25] = ayear[ti-1][25]+ f_MA_growth_co2;
        ayear[ti][26] = ayear[ti-1][26]+ f_MA_co2;
        ayear[ti][27] = ayear[ti-1][27]+ f_MD_co2;
        ayear[ti][28] = ayear[ti-1][28]+ f_MA_growth;
        ayear[ti][29] = ayear[ti-1][29]+ f_MB_LM;
        ayear[ti][30] = ayear[ti-1][30]+ f_MB_MA;
        ayear[ti][31] = ayear[ti-1][31]+ f_MB_turn;
        ayear[ti][32] = ayear[ti-1][32]+ f_LM_leach;
      }
    }
    curtime = curtime + 1;
    
    if(clab ){ //&& nrep != 365
      c1 = {0,LIT, POM, LMWC, MIC, MAOM, MA, MD, 0};
      Strflux flux0 = {f_Lit_PO, f_Lit_LM, f_Exud_LM, f_Oa_MA, f_Oa_LM,
                       f_PO_LM, f_PO_MA, f_MA_LMm, f_MA_LMe, acid_des,f_MA_LM,f_LM_MA, 
                       f_MA_MD, f_MD_MA, f_LM_MB,  f_MA_maint_co2, f_MA_growth_co2,
                       f_MA_co2, f_MD_co2, f_MA_growth,
                       f_MB_LM, f_MB_MA, f_MB_turn,f_LM_leach
      };
        
 
      if(t==0){
        c13->LIT = 0.;
        c13->POM = 0.;//as<double>(C13poolsList["POM"]); 
        c13->LMWC = as<double>(C13poolsList["LMWC"]); 
        c13->MIC = as<double>(C13poolsList["MIC"]); 
        c13->MAOM = as<double>(C13poolsList["MAOM"]); 
        c13->MA = as<double>(C13poolsList["MA"]); 
        c13->MD = as<double>(C13poolsList["MD"]); 
        c13->Tco2 = 0.; 
        c13->SOM = c13->POM + c13->LMWC + c13->MIC + c13->MAOM;
      }

      CalLabC(c13, c0, c1, flux0,CUE, pars.param_pb);
      
      int colth;
      if(outflux){
        colth = 32;
      }else{
        colth = 8;
      }
      ayear[ti][colth] = c13->SOM;
      ayear[ti][colth+1] = c13->LIT;
      ayear[ti][colth+2] = c13->POM;
      ayear[ti][colth+3] = c13->LMWC;
      ayear[ti][colth+4] = c13->MIC;
      ayear[ti][colth+5] = c13->MAOM;
      ayear[ti][colth+6] = c13->MA;
      ayear[ti][colth+7] = c13->MD;
      ayear[ti][colth+8] = c13->Tco2;
    }
    
    int year = t/365;
    // cout << "---1111---error---"<<"curtime="<<curtime<<"year="<<year<<endl;
    if(curtime%365 == 0){
      // cout << "------error---"<<"curtime="<<curtime<<"year="<<year+1<<endl;
      if(nloop==0){
        out[year][0] = year+1;
        for(unsigned int col = 0; col < ayear[0].size(); ++col){
          out[year][col+1] = ayear[364][col];//c pools is the last day,flux is accumulate
        }
        c13->Tco2 = 0.;
      }else{
        for (int ci = 0; ci < ayear[0].size(); ++ci) {
          float sum = 0.f;
          if(ci==8){ //Tco2 exclude
            continue;
          }
          for (int day = 0; day < 365; ++day) {
            sum += ayear[day][ci];
          }
          
          int col;
          if(ci>8){
            col = ci-1;
          }else {
            col = ci;
          }
          if(ssy_aloop==1){
            ssyears[0][col] = sum;  //solve steady-state: C pools,include SOM
          }else{
            ssyears[year % nyear][col] = sum;  //solve steady-state: C pools,include SOM
          }
        }
      }
      ayear.assign(365, std::vector<double>(ayear[0].size(), 0.0));
    }

    // ==========================
    // solve steady state
    // ==========================
    if(nloop!=0){
      sumSOM = sumSOM + ayear[t%365][0] + ayear[t%365][1];
      // cout << "----error--loop-"<<"curtime="<<curtime<<"nrep="<<nrep<<endl;
      if(curtime%nrep==0 || curtime==sstime){
        if (ssi == ssSOM.size() ) {
          ssSOM.erase(ssSOM.begin()); // delete the first element
          ssSOM.push_back(sumSOM / nrep); // Adds the latest value to the last position
        } else {
          ssSOM[ssi] = sumSOM / nrep;
          ssi = ssi + 1;
        }
        
        if(loopth >= 40){
          boolss = sscondfun(ssSOM); //sscondfun(ssSOM[Range(0, 9)])
          // cout << "----error-----nrep="<<nrep<<"  ssyears.size()"<< ssyears.size()<<
          //   " year="<<year<<" boolss="<<boolss<<endl;
          // cout << "POM="<<POM<<" LMWC="<<LMWC<<" MIC="<<(MA + MD)<<" MAOM="<<MAOM<<" MA="<<MA<<" MD="<<MD<<endl;
          // cout << "POM="<<ssyears[year][0]<<" LMWC="<<ssyears[year][1]<<" MIC="<<ssyears[year][2]<<" MAOM="<<ssyears[year][3]
          // <<" MA="<<ssyears[year][4]<<" MD="<<ssyears[year][5]<<endl;
          
          //compute each loop c pools content
          if(nrep == 365){ //if dynamic input is constant;
            for(int ci=0; ci < outss[0].size(); ++ci){
              outss[0][ci] = ssyears[0][ci]/nrep;
            }
          }else{
            for (int ci = 0; ci < outss[0].size(); ++ci) {
              float sum = 0.0f;
              for (unsigned int yi = 0; yi < ssyears.size(); ++yi) {
                sum += ssyears[yi][ci];
              }
              outss[0][ci] = sum / nrep;  
            }
          }
          ssyears.assign(ssyears.size(), std::vector<double>(ssyears[0].size(), 0.0));
          
          // cout << "ssi="<<ssi<<"  loopth="<<loopth<<" boolss="<<boolss<<endl;
          if(boolss || curtime==sstime || (outss[0][1]>1500 || outss[0][2]>1000)){
            if(curtime==sstime){
              boolss = sscondfun(ssSOM[Range(8, 9)]); //the latest two loop
            }
            // add ssgx
            if(boolss){
              outss[0][outss[0].size()-1] = 1;
            }else{
              outss[0][outss[0].size()-1] = 100;
            }

            //cout << "curtime="<<curtime<<endl;
            cout << "curtime="<<curtime
                 << "outss=" << outss[0][outss[0].size()-1]<<endl;
            
            
            delete c13;
            return outss;
          }
          
        }
        sumSOM = 0;
        loopth = loopth + 1;
      }
    }
  }//time ends
  
  //cue or qmax exit tackle
  if(curtime < sstime){
    if(nloop!=0){
      std::fill(outss[0].begin(), outss[0].end(), na_value);//solve steady state
    }else{
      for(unsigned int y=curtime/365; y < out.size(); ++y){
        fill(out[y].begin(), out[y].end(), na_value);
      }
    }
  }
  //return output
  delete c13;
  if(nloop!=0){
    // cout <<"output ------"<<endl;
    return outss;
  }else if(sstime<365){
    return ayear;
  }else{
    return out;
  }
}//derivs_V2_MM_AD_CO2()







