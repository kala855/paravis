#include "cell.cuh"
#include <stdlib.h>
#include <stdio.h>

__host__ __device__ Cell::Cell(){
  // constants
  testChange = 0.0;
  //R = 8.3143;         // gas constant [J/K.mmol];
  T = 310.0;          // temperature [K];
  F = 96.4867;        // faraday constant [C/mmol] ;
  RTF = (R*T)/F;      // J/C
  invRTF = 1.0/RTF;

  Cap = 100.0;        // membrane capacitance [pF]

  Vi = 13668.0;       // intracellular volumen [um^3]
  Vup = 1109.52;      // SR uptake compartment volume [um^3]
  Vrel = 96.48;       // SR release compartment volume [um^3]

  // Cell Geometry
  l = 0.01;           // length of the cell [um]
  a = 0.0008;         // radius of the cell [um]
  pi = 2*acos(0.0);

  Ri = 200.0;         // 200 Ohm x cm = 0.2 K Ohms x cm, Tesis Catalina, page 99 y 115 , Resistividad no resistencia
  Rix = 0.2;
  Riy = 0.2;

  // External concentration
  Ko = 5.4;           // extracellular K concentration [mM]
  Nao = 140.0;        // extracellular Na concentration [mM]
  Coa = 1.8;          // extracellular Ca concentration [mM]

  // Maximal  conductances  [nS/pF]
  GNa = 7.8;
  GK1 =  0.09;
  Gto = 0.1652;
  GKr = 0.0294;
  GKs = 0.129;
  GCaL = 0.1238;
  GbCa = 0.00113;
  GbNa = 0.000674;

  // Maximal currents
  INaK_max = 0.6;     // Maximal INaK [pA/pF]
  INaCa_max = 1600.0; // Maximal INaCa [pA/pF]
  IpCa_max = 0.275;   // Maximal IpCa [pA/pF]
  Kq10 = 3.0;         // Temperature scaling factor for IKur and Ito kinetics
  gamma = 0.35;       // Voltage dependance parameter for INaCa

  // Half-saturation constant for currents
  KmNai = 10.0;       // Nai half-saturation constant of INaK [mM]
  KmKo = 1.5;         // Ko half-saturation constant of INaK [mM]
  KmNa = 87.5;        // Nao half-saturation constant of INaCa [mM]
  KmCa = 1.38;        // Cao half-saturation constant of INaCa

  ksat = 0.1;         // Saturation factor for INaCa


  // Ion Valences
  zna = 1.0;          // Na valence
  zk = 1.0;           // K valence
  zca = 2.0;          // Ca valence

  // Myoplasmic Ca Ion Concentration Changes
  Csqn_max = 10.0;    // Total calsequestrin concentration in SR release compartment [mM]
  Km_csqn = 0.8;      // Ca_rel half-saturation constant of Iup [mM]
  Cmdn_max = 0.050;   // Total calmodulin concentration in myoplasm [mM]
  Trpn_max = 0.070;   // Total troponin concentration in myoplasm [mM]
  kmcmdn = 0.00238;   // Cai half-saturation constant for calmodulin [mM]
  Kmtrpn = 0.0005;    // Cai half-saturation constant for troponin [mM]
  Iup_max = 0.005;    // Maximal Iup [mM/mS]

  // future function "initial conditions"
  V = -8.12e1;          // mV
  h = 9.65e-1;
  d = 1.37e-4;
  xr = 3.29e-5;
  Nai = 1.12e1;         // Initial Intracellular Na (mM)
  Ki = 1.39e2;          // Initial Intracellular Ki (mM)
  Ca_rel = 1.49;
  oi = 9.99e-1;
  ui = 9.99e-1;
  /*
    [Cmdn-Ca2+]i=2.05e-3
    [Csqn-Ca2+]i=6.51
  */

  v = 1.0;             // Activation gate v of Ca release from jsr
  m = 2.91e-3;
  j = 9.78e-1;
  f = 9.99e-1;
  xs = 1.87e-2;
  Cai = 1.02e-4;       // Initial Intracellular Ca
  Ca_up = 1.49;
  oa = 3.04e-2;        /* Paralpha_meters Transient Outward Current ito */
  ua = 4.96e-3;        /* Paralpha_meters Ultra-Rapidly activation K Current ikur */
  fca = 7.75e-1;
  /*
    [Trpn-Ca2+]i=1.18e-2
  */
  u = 0.0;          // Activation gate u of Ca release from jsr // Gates Irel
  w = 9.99e-1;      // Inactivation gate w of Ca release from jsr// Gates Irel
  Itot = 0.0;       // mA Current Total
}

__device__ __host__
db Cell::getItot(db dt){
  compute_currents();
  compute_concentrations(dt);
  compute_gates(dt);
  return Itot;
}

/* Calculates All Currents */
__device__ __host__
void Cell::compute_currents(){
  ECa = ((R*T)/(zca*F)) * log(Coa/Cai);
  ENa = ((R*T)/(zna*F)) * log(Nao/Nai);
  EK = ((R*T)/(zk*F)) * log(Ko/Ki);
  ENC = (F*V) / (R*T);

  comp_ical ();     // Calculates Currents through L-Type Ca Channel
  comp_inaca ();    // Calculates Na-Ca Exchanger Current
  comp_ibna ();     // Calculates Na Background Current
  comp_ibca ();     // Calculates Ca Background Current
  comp_ina ();      // Calculates Fast Na Current
  comp_ikr ();      // Calculates Rapidly Activating K Current
  comp_ipca ();     // Calculates Sarcolemmal Ca Pump Current
  comp_iks ();      // Calculates Slowly Activating K Current
  comp_inak ();     // Calculates Na-K Pump Current
  comp_ik1 ();      // Calculates Time-Independant K Current
  comp_itr();
  comp_ito ();      // Calculates Transient Outward Current
  comp_ikur ();     // Calculates Ultra-Rapidly activation K Current
  comp_itot();      // Calulates Total Current
}


__device__ __host__
void Cell::compute_concentrations(db dt){

  //////////DUDAS SOBRE USO///////////////////////
  // Calsequestrin concentration
  //db Ca_csqn = Csqn_max*(Ca_rel/(Ca_rel+Km_csqn));    // Equation 75, Uso?
  //db Ca_Trpn = Trpn_max*(Cai/(Cai+Kmtrpn));           // Equation 74, no se usa
  //db Ca_Cmdn = Cmdn_max*(Cai/(Cai+kmcmdn));           // Equation 73, no se usa
  //////////////////////////////////////////////////

  comp_iupleak();      // Ca leak current by the NSR
  comp_iup();          // Ca uptake current by the NSR
  comp_irel();         // Ca release current from JSR

  // Intracellular ion concentrations
  conc_nai(dt);        // Equation 21
  conc_ki(dt);         // Equation 22
  conc_cai(dt);        // Equation 23
  conc_ca_up(dt);      // Equation 26
  conc_ca_rel(dt);     // Equation 27
}

__device__ __host__
void Cell::conc_nai(db dt){
  // Compute Intracellular Nai Concentration
  db invViF  = 1.0/(Vi*F);
  db totINa = INa+IbNa+3.0*(INaK+INaca);
  db dNai = dt*(-totINa*invViF);                       // Equation 21
  Nai = dNai + Nai;
}

__device__ __host__
void Cell::conc_ki(db dt){
  // Compute Intracellular Ki Concentration
  // En el paper aparece en IbK, pero no esta.
  db invViF  = 1.0/(Vi*F);
  db totIK = 2.0*INaK-IK1-Ito-IKur-IKr-IKs;//-IbK;
  db dKi = dt*(totIK*invViF);                          // Equation 22
  Ki = dKi + Ki;
}

__device__ __host__
void Cell::conc_cai(db dt){
  // Compute Intracellular Cai Concentration
  db invViF2 = 1.0 / (2.0*Vi*F);
  db b1_left = ((2.0*INaca -IpCa-ICal-IbCa)*invViF2);         // left Ecuation 24
  db b1_right = ((Vup*(Iup_leak-Iup)+(Irel*Vrel))/Vi);        // right Ecuation 24
  db b1cai = b1_left + b1_right;                              // Ecuation 24
  db b2_left=1.0+((Trpn_max*Kmtrpn)/pow((Cai + Kmtrpn),2.0)); // left Ecuation 25
  db b2_right= (Cmdn_max*kmcmdn)/pow((Cai + kmcmdn),2.0);     // right Ecuation 25
  db b2cai = b2_left + b2_right;                              // Ecuation 25
  db dcai = dt*(b1cai/b2cai);                                 // Equation 23

  Cai = dcai + Cai;
}

__device__ __host__
void Cell::conc_ca_up(db dt){
  // Compute Ca2+ concentration in uptake compartment Ca_up //nsr
  db dCa_up = dt*(Iup - Iup_leak - Itr*(Vrel/Vup));            // Equation 26
  Ca_up = dCa_up + Ca_up;
}

__device__ __host__
void Cell::conc_ca_rel(db dt){
  // Compute Ca2+ concentration release compartment Ca_rel //jsr
  db dCa_rel = dt*(Itr-Irel)/(1.0+(Csqn_max*Km_csqn)/pow((Ca_rel+Km_csqn),2.0));  // Equation 27
  Ca_rel = dCa_rel + Ca_rel;
}

/* Calculates Fast Na Current  INa*/
__device__ __host__
void Cell::comp_ina(){
  // Probable explicación de multiplicacion por Cap = 100.
  // Las unidades de la conductancia G son Siemens, pero en el paper
  // de CRN, las unidades son nS/pF, eso muestra que la conductancia que
  // nos estan dando ya fue divida por la capacitancia, y si luego lo volvemos a
  // dividir por Cap en el calculo de B, las unidades quedarian nS/pF^2. e
  // Y la corrient tambien quedaria en pA/pF^2.
  // Al multiplicar por Cap, estamos dejando solo en Siemens, para luego si dividir
  // en el calculo de B por Cap.
  // Cap quedan en Siemnes, par

  INa = Cap*GNa*pow(m,3.0)*h*j*(V-ENa);                         // Equation 29
}

/* Calculates Time-Independant K Current IK1*/
__device__ __host__
void Cell::comp_ik1 (){
  IK1 = Cap*(GK1*(V-EK)) / (1.0+exp(0.07*(V+80.0)));            // Equation 35
}

/* Calculates Transient Outward Current  Ito*/
__device__ __host__
void Cell::comp_ito (){
  Ito = Cap*Gto*pow(oa,3.0)*oi*(V-EK);                          //Equation 36
}

/* Calculates Ultra-Rapidly activation K Current IKur*/
__device__ __host__
void Cell::comp_ikur (){
  db GKur = 0.005+(0.05/(1.0+exp(-(V-15.0)/13.0)));             // Equation 42
  IKur = Cap*GKur*pow(ua,3.0)*ui*(V-EK);                        // Equation 41
}

/* Calculates Rapidly Activating K Current Ikr*/
__device__ __host__
void Cell::comp_ikr (){
  db r = 1.0/(1.0+exp((V+15.0)/22.4));
  IKr = Cap*GKr*xr*r*(V-EK);                                    // Equation 47
}

/* Calculates Slowly Activating K Current  IKs*/
__device__ __host__
void Cell::comp_iks (){
  IKs = Cap*GKs*pow(xs,2.0)*(V-EK);                             // Equation 50
}

/* Calculates Currents through L-Type Ca Channel */
__device__ __host__
void Cell::comp_ical (){
  ICal = Cap*GCaL*d*f*fca*(V-65.0);  // ICal  Equation 53
}

/* Calculates Na-K Pump Current */
__device__ __host__
void Cell::comp_inak (){
  db sigma = (exp(Nao/67.3)-1.0)/7.0;                                     // Equation 59
  db fNaK= 1.0/(1.0+0.1245*exp(-0.1*ENC)+0.0365*sigma*exp(-ENC));         // Equation 58
  INaK = Cap*INaK_max*fNaK*(1.0/(1.0+pow((KmNai/Nai),1.5)))*(Ko/(Ko+KmKo));   // Equation 57
}

/* Calculates Na-Ca Exchanger Current */
__device__ __host__
void Cell::comp_inaca (){
  db phif = exp(gamma*ENC);
  db phir = exp((gamma-1.0)*ENC);
  db nmr  = (phif*pow(Nai,3.0)*Coa)-(phir*pow(Nao,3.0)*Cai);
  db dnm  = (pow(KmNa,3.0)+pow(Nao,3.0))*(KmCa+Coa)*(1.0+(ksat*phir));
  INaca = Cap*INaCa_max*(nmr/dnm);                                             // Equation 60
}

/* Calculates Sarcolemmal Ca Pump Current */
__device__ __host__
void Cell::comp_ipca (){
  IpCa = Cap*(IpCa_max*Cai)/(0.0005+Cai);  // IpCa Equation 63
}

/* Calculates Ca Background Current */
__device__ __host__
void Cell::comp_ibca (){
  IbCa = Cap*GbCa*(V-ECa);                // IbCa  Equation 61
}

/* Calculates Na Background Current ibna */
__device__ __host__
void Cell::comp_ibna (){
  IbNa = Cap*GbNa*(V-ENa);                // IbNa  Equation 62
}

// Compute Ca2+ Release Current From JSR Irel
__device__ __host__
void Cell::comp_irel(){
  db krel = 30.0;  // Rate constant of Ca release from JSR due to overload (ms^-1)
  Irel = krel*pow(u,2.0)*v*w*(Ca_rel-Cai);   // Equation 64
}

// Compute Transfer Current From NSR to JSR Itr
__device__ __host__
void Cell::comp_itr(){
  db tautr = 180.0;               // Time constant of Ca transfer from NSR to JSR(ms) ecu 69
  Itr = (Ca_up - Ca_rel)/tautr;   // Equation 69 for dCa_rel, dCa_up
}

// Compute Ca2+ Uptake Current by NSR Iup
__device__ __host__
void Cell::comp_iup(){
  db Kup= 0.00092;                   // Half-saturation concentration of iup (mM)
  Iup = Iup_max / (1.0+(Kup/Cai));   // Equation 71
}

// Compute Ca2+ Leak Current by the NSR Iup_leak
__device__ __host__
void Cell::comp_iupleak(){
  db Ca_up_max = 15.0;                      //  Max. [Ca] in NSR (m)M
  Iup_leak = (Ca_up/Ca_up_max)*Iup_max;     // Equation 72
}

__device__ __host__
void Cell::comp_itot(){
  db IK,INat,ICa;
  IK = IKr + IKs + IK1 + IKur;
  INat = INa + IbNa + INaK + INaca;
  ICa = ICal + IbCa + IpCa;
  Itot = IK + INat + ICa + Ito;
}

__device__ __host__
void Cell::compute_gates(db dt){
  // Compute gates

  testChange = testChange + 1.0;
  gates_irel(dt);   //u,v,w
  gates_ical(dt);   // d,f,fca
  gates_ina(dt);    // h,j,m
  gates_ikr(dt);    // xr
  gates_iks(dt);    // xs
  gates_ito(dt);    // oa, oi
  gates_ikur(dt);   // ua, ui
}

__device__ __host__
void Cell::gates_irel(db dt){
  // Gates for Irel Current
  db fn = (Vrel * (10e-12) * Irel) -((5.0e-13/F) * (0.5*ICal-0.2*INaca));   // Equation 68
  db tauu = 8.0;                                                            // Equation 65
  db u_inf = 1.0/(1.0+exp(-(fn-3.4175e-13)/13.67e-16));
  db tauv = 1.91+(2.09/(1.0+exp(-(fn-3.4175e-13)/13.67e-16)));              // Equation 66
  db v_inf = 1.0-(1.0/(1.0+exp(-(fn-6.835e-14)/13.67e-16)));
  db tauw = 6.0*(1.0-exp(-(V-7.9)/5.0))/((1.0+0.3*exp(-(V-7.9)/5.0))*(V-7.9));
  db w_inf = 1.0-(1.0/(1.0+exp(-(V-40.0)/17.0)));                           // Equation 67
  //Compute Gates
  u = u_inf+(u-u_inf)*exp(-dt/tauu);   // Activation gate u of Ca release from jsr
  v = v_inf+(v-v_inf)*exp(-dt/tauv);   // Activation gate v of Ca release from jsr
  w = w_inf+(w-w_inf)*exp(-dt/tauw);   // Inactivation gate w of Ca release from jsr
}

__device__ __host__
void Cell::gates_ina(db dt){
  // Gates: m,h,j.
  db alpha_m,beta_m,alpha_h,beta_h,alpha_j,beta_j,tau_m, m_inf, tau_h;
  db h_inf, tau_j, j_inf;

  alpha_m = ((V == -47.13)? 3.2 : 0.32*(V+47.13)/(1.0-exp(-0.1*(V+47.13))));  // Equation 30
  beta_m = 0.08 * exp(-V/11.0);

  if (V < -40.0){ // Equation 31,32,33
    alpha_h = 0.135 * exp(-(80.0+V)/6.8);
    beta_h = 3.56 * exp(0.079*V) + 3.1e5 *exp(0.35*V);
    alpha_j = ( (-127140 * exp(0.2444*V)) - (3.474e-5 * exp(-0.04391*V))) * ((V+37.78)/(1.0+exp(0.311*(V+79.23))));
    beta_j = (0.1212 * exp(-0.01052*V))/(1 + exp(-0.1378 * (V+40.14)));
  } else {
    alpha_h = 0.0;
    beta_h = 1.0 / (0.13 * (1.0+exp(-(V+10.66)/11.1)));
    alpha_j = 0.0;
    beta_j = (0.3 * exp(-2.535e-7*V))/(1.0+exp(-0.1*(V+32.0)));
  }

  tau_m = (1.0 / (alpha_m+beta_m));          // Equation 34
  m_inf = alpha_m * tau_m;
  tau_h = (1.0 / (alpha_h+beta_h));
  h_inf = alpha_h * tau_h;
  tau_j = (1.0 / (alpha_j+beta_j));
  j_inf= alpha_j*tau_j;

  // Update gates
  m = m_inf +(m-m_inf)*exp(-dt/tau_m);       // Equation 77
  h = h_inf +(h-h_inf)*exp(-dt/tau_h);
  j = j_inf+(j-j_inf)*exp(-dt/tau_j);
}

__device__ __host__
void Cell::gates_ito(db dt){
  //ACTUALIZO COMPUERTAS
  db alpha_oa, beta_oa,tau_oa,oa_inf,alpha_oi,beta_oi,tau_oi, oi_inf;

  // Gates: oa,oi.
  alpha_oa = 0.65/(exp(-(V+10.0)/8.5)+exp(-(V-30.0)/59.0));   // Equation 37
  beta_oa = 0.65/(2.5+exp((V+82.0)/17.0));
  tau_oa = 1.0/((alpha_oa+beta_oa)*Kq10);                     // Equation 38
  oa_inf = 1.0/(1.0+exp(-(V+20.47)/17.54));

  alpha_oi= 1.0/(18.53+exp((V+113.7)/10.95));                 // Equation 39
  beta_oi = 1.0/(35.56+exp(-(V+1.26)/7.44));
  tau_oi = 1.0/((alpha_oi+beta_oi)*Kq10);                     // Equation 40
  oi_inf = 1.0/(1.0+exp((V+43.1)/5.3));

  // Updates gates
  oa = oa_inf+(oa-oa_inf)*exp(-dt/tau_oa);                    // Equation 77
  oi = oi_inf+(oi-oi_inf)*exp(-dt/tau_oi);
}

__device__ __host__
void Cell::gates_ikur(db dt){
  //ACTUALIZO COMPUERTAS
  db alpha_ua, beta_ua,tau_ua,ua_inf,alpha_ui,beta_ui, tau_ui;
  db ui_inf;

  // Gates: uo,ui.
  alpha_ua = 0.65/(exp(-(V+10.0)/8.5)+exp(-(V-30.0)/59.0));   // Equation 43
  beta_ua = 0.65/(2.5+exp((V+82.0)/17.0));
  tau_ua = 1.0/((alpha_ua+beta_ua)*Kq10);                     // Equation 44
  ua_inf = 1.0/(1.0+exp(-(V+30.3)/9.6));

  alpha_ui = 1.0/(21.0+exp(-(V-185.0)/28.0));                 // Equation 45
  beta_ui = exp((V-158.0)/16.0);
  tau_ui = 1.0/((alpha_ui+beta_ui)*Kq10);                     // Equation 46
  ui_inf = 1.0/(1.0+exp((V-99.45)/27.48));

  // Updates gates
  ua = ua_inf+(ua-ua_inf)*exp(-dt/tau_ua);
  ui = ui_inf+(ui-ui_inf)*exp(-dt/tau_ui);
}

__device__ __host__
void Cell::gates_ikr(db dt){
  //ACTUALIZO COMPUERTAS
  db alpha_xr, beta_xr, tau_xr, xr_inf;

  alpha_xr = 0.0003 * (( V + 14.1)/(1.0-exp(-(V + 14.1)/5.0)));
  beta_xr = 7.3898e-5*((V -3.3328)/(exp((V-3.3328)/5.1237)-1.0));
  tau_xr = 1.0 / (alpha_xr + beta_xr);
  xr_inf = 1.0 / (1.0 + exp(-(V + 14.1) / 6.5) );

  xr = xr_inf + (xr-xr_inf)*exp(-dt/tau_xr);

}

__device__ __host__
void Cell::gates_iks(db dt){
  //ACTUALIZO COMPUERTAS
  db alpha_xs,beta_xs,tau_xs,xs_inf;

  // Gate: xs
  alpha_xs = 4.0e-5 * ((V-19.9) / (1.0-exp(-(V-19.9)/17.0)));    // Equation 51
  beta_xs = 3.5e-5 * ((V-19.9) / (exp((V-19.9)/9.0)-1.0));
  tau_xs = 0.5 / (alpha_xs+beta_xs);                             // Equation 52
  xs_inf = pow(1.0 + (exp(-(V-19.9)/12.7)),-0.5);

  // Update gate
  xs = xs_inf+(xs-xs_inf)*exp(-dt/tau_xs);                       // Equation 77
}

__device__ __host__
void Cell::gates_ical(db dt){
  //ACTUALIZO COMPUERTAS
  db d_inf, tau_d, f_inf, tau_f, fca_inf, tau_fca;

  tau_d = (1.0 - exp((V+10.0)/-6.24))/ (0.035*(V+10.0)*(1.0+exp((V+10.0)/-6.24))); //Equation 54
  d_inf = 1.0/(1.0+exp((V+10.0)/-8.0));                        // Equation 54
  tau_f = 9.0/(0.0197*exp((-1.0)*pow(0.0337,2.0)*pow((V+10.0),2.0))+0.02);
  f_inf = 1.0/(1.0+exp((V+28.0)/6.9));                         // Equation 55
  fca_inf = 1.0/(1.0+(Cai/0.00035));                           // Equation 56
  tau_fca = 2.0;
  d = d_inf + (d - d_inf)*exp(-dt/tau_d);
  f = f_inf + (f - f_inf)*exp(-dt/tau_f);
  fca = fca_inf + (fca - fca_inf)*exp(-dt/tau_fca);
}
