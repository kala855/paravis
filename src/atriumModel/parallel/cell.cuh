//header cell class
#include <math.h>
#define db double

#define R 8.3143
#define TEMP 310.0
#define F 96.4861


#define RTF (R*TEMP)/F      // J/C
#define INVRTF 1.0/RTF

#define CAP 100.0        // membrane capacitance [pF]

#define VI 13668.0
#define VUP 1109.52      // SR uptake compartment volume [um^3]
#define VREL 96.48       // SR release compartment volume [um^3]

// Cell Geometry
#define  LENGTHCELL 0.01           // length of the cell [um]
#define RADIUSCELL 0.0008         // radius of the cell [um]
#define PINUM 2*acos(0.0)

#define RI 200.0         // 200 Ohm x cm = 0.2 K Ohms x cm, Tesis Catalina, page 99 y 115 , Resistividad no resistencia
#define RIX 0.2
#define RIY 0.2

// External concentration
#define KO 5.4           // extracellular K concentration [mM]
#define NAO 140.0        // extracellular Na concentration [mM]
#define COA 1.8          // extracellular Ca concentration [mM]

#define GNA 7.8
#define GK1 0.09
#define GTO 0.1652
#define GKR 0.0294
#define GKS 0.129
#define GCAL 0.1238
#define GBCA 0.00113
#define GBNA 0.000674

#define INAK_MAX 0.6     // Maximal INaK [pA/pF]
#define INACA_MAX 1600.0 // Maximal INaCa [pA/pF]
#define IPCA_MAX 0.275   // Maximal IpCa [pA/pF]
#define KQ10 3.0        // Temperature scaling factor for IKur and Ito kinetics
#define GAMMA 0.35      // Voltage dependance parameter for INaCa

#define KmNai 10.0       // Nai half-saturation constant of INaK [mM]
#define KmKo 1.5         // Ko half-saturation constant of INaK [mM]
#define KmNa 87.5        // Nao half-saturation constant of INaCa [mM]
#define KmCa 1.38        // Cao half-saturation constant of INaCa

#define ksat 0.1         // Saturation factor for INaCa

#define zna 1.0          // Na valence
#define zk 1.0           // K valence
#define zca 2.0          // Ca valence

#define CSQN_MAX 10.0    // Total calsequestrin concentration in SR release compartment [mM]
#define KM_CSQN 0.8      // Ca_rel half-saturation constant of Iup [mM]
#define CMDN_MAX 0.050   // Total calmodulin concentration in myoplasm [mM]
#define TRPN_MAX 0.070   // Total troponin concentration in myoplasm [mM]
#define KMCMDN 0.00238   // Cai half-saturation constant for calmodulin [mM]
#define KMTRPN 0.0005    // Cai half-saturation constant for troponin [mM]
#define IUP_MAX 0.005    // Maximal Iup [mM/mS]*/

//Some new definitions

#define PHIF exp(GAMMA*ENC)
#define PHIR exp((GAMMA-1.0)*ENC)
#define DNM (pow(KmNa,3.0)+pow(NAO,3.0))*(KmCa+COA)*(1.0+(ksat*PHIR))
#define NMR (PHIF*pow(Nai,3.0)*COA)-(PHIR*pow(NAO,3.0)*Cai)
#define SIGMA (exp(NAO/67.3)-1.0)/7.0
#define FNAK 1.0/(1.0+0.1245*exp(-0.1*ENC)+0.0365*SIGMA*exp(-ENC))
#define INVVIF 1.0/(VI*F)
#define TOTINA INa+IbNa+3.0*(INaK+INaca)
#define DNAI dt*(-TOTINA*INVVIF)                       // Equation 21
#define TOTIK 2.0*INaK-IK1-Ito-IKur-IKr-IKs
#define DKI dt*(TOTIK*INVVIF)
#define invViF2 1.0 / (2.0*VI*F)
#define b1_left ((2.0*INaca -IpCa-ICal-IbCa)*invViF2)         // left Ecuation 24
#define b1_right ((VUP*(Iup_leak-Iup)+(Irel*VREL))/VI)        // right Ecuation 24
#define b1cai  b1_left + b1_right                              // Ecuation 24
#define b2_left 1.0+((TRPN_MAX*KMTRPN)/pow((Cai + KMTRPN),2.0)) // left Ecuation 25
#define b2_right (CMDN_MAX*KMCMDN)/pow((Cai + KMCMDN),2.0)     // right Ecuation 25
#define b2cai b2_left + b2_right                              // Ecuation 25
#define dcai dt*(b1cai/b2cai)                                 // Equation 23
#define dCa_up dt*(Iup - Iup_leak - Itr*(VREL/VUP))            // Equation 26
#define dCa_rel dt*(Itr-Irel)/(1.0+(CSQN_MAX*KM_CSQN)/pow((Ca_rel+KM_CSQN),2.0))  // Equation 27
#define GKur 0.005+(0.05/(1.0+exp(-(V-15.0)/13.0)))             // Equation 42
#define FN (VREL * (10e-12) * Irel) -((5.0e-13/F) * (0.5*ICal-0.2*INaca))   // Equation 68
#define tauu 8.0                                                            // Equation 65
#define u_inf 1.0/(1.0+exp(-(FN-3.4175e-13)/13.67e-16))
#define tauv 1.91+(2.09/(1.0+exp(-(FN-3.4175e-13)/13.67e-16)))              // Equation 66
#define v_inf 1.0-(1.0/(1.0+exp(-(FN-6.835e-14)/13.67e-16)))
#define tauw 6.0*(1.0-exp(-(V-7.9)/5.0))/((1.0+0.3*exp(-(V-7.9)/5.0))*(V-7.9))
#define w_inf 1.0-(1.0/(1.0+exp(-(V-40.0)/17.0)))                           // Equation 67


#define beta_m 0.08 * exp(-V/11.0)
#define tau_m (1.0 / (alpha_m+beta_m))          // Equation 34
#define m_inf alpha_m * tau_m
#define  tau_h (1.0 / (alpha_h+beta_h))
#define  h_inf alpha_h * tau_h
#define  tau_j (1.0 / (alpha_j+beta_j))
#define  j_inf alpha_j*tau_j

#define alpha_oa 0.65/(exp(-(V+10.0)/8.5)+exp(-(V-30.0)/59.0))   // Equation 37
#define beta_oa 0.65/(2.5+exp((V+82.0)/17.0))
#define tau_oa 1.0/((alpha_oa+beta_oa)*KQ10)               // Equation 38
#define oa_inf 1.0/(1.0+exp(-(V+20.47)/17.54))

#define alpha_oi 1.0/(18.53+exp((V+113.7)/10.95))                 // Equation 39
#define beta_oi 1.0/(35.56+exp(-(V+1.26)/7.44))
#define tau_oi 1.0/((alpha_oi+beta_oi)*KQ10)                     // Equation 40
#define oi_inf 1.0/(1.0+exp((V+43.1)/5.3))

#define alpha_ua 0.65/(exp(-(V+10.0)/8.5)+exp(-(V-30.0)/59.0))   // Equation 43
#define beta_ua 0.65/(2.5+exp((V+82.0)/17.0))
#define tau_ua 1.0/((alpha_ua+beta_ua)*KQ10)                     // Equation 44
#define ua_inf 1.0/(1.0+exp(-(V+30.3)/9.6))

#define alpha_ui 1.0/(21.0+exp(-(V-185.0)/28.0))                 // Equation 45
#define beta_ui exp((V-158.0)/16.0)
#define tau_ui 1.0/((alpha_ui+beta_ui)*KQ10)                     // Equation 46
#define ui_inf 1.0/(1.0+exp((V-99.45)/27.48))

#define alpha_xr 0.0003 * (( V + 14.1)/(1.0-exp(-(V + 14.1)/5.0)))
#define beta_xr 7.3898e-5*((V -3.3328)/(exp((V-3.3328)/5.1237)-1.0))
#define tau_xr 1.0 / (alpha_xr + beta_xr)
#define xr_inf 1.0 / (1.0 + exp(-(V + 14.1) / 6.5) )

#define alpha_xs 4.0e-5 * ((V-19.9) / (1.0-exp(-(V-19.9)/17.0)))    // Equation 51
#define beta_xs 3.5e-5 * ((V-19.9) / (exp((V-19.9)/9.0)-1.0))
#define tau_xs 0.5 / (alpha_xs+beta_xs)                             // Equation 52
#define xs_inf pow(1.0 + (exp(-(V-19.9)/12.7)),-0.5)

#define tau_d (1.0 - exp((V+10.0)/-6.24))/ (0.035*(V+10.0)*(1.0+exp((V+10.0)/-6.24))) //Equation 54
#define d_inf 1.0/(1.0+exp((V+10.0)/-8.0))                        // Equation 54
#define tau_f 9.0/(0.0197*exp((-1.0)*pow(0.0337,2.0)*pow((V+10.0),2.0))+0.02)
#define f_inf 1.0/(1.0+exp((V+28.0)/6.9))                         // Equation 55
#define fca_inf 1.0/(1.0+(Cai/0.00035))                          // Equation 56
#define tau_fca 2.0


class Cell{
public:
    //------------    parameters   ----------------
    // future function "initial conditions"
    db V;           // mV
    db Cai;         // Initial Intracellular Ca
    db Nai;         // Initial Intracellular Na (mM)
    db Ki;          // Initial Intracellular Ki (mM)
    db m;
    db h;
    db j;
    db xr ;
    db xs ;
    db d;
    db f;
    db fca;
    db yach;

    /* Paraers Ultra-Rapidly activation K Current ikur */
    db ua;
    db ui;

    /* Paraers Transient Outward Current ito */
    db oa;
    db oi;
    db Ca_up;       // (nsr) Ca2+ concentration in uptake compartment
    db Ca_rel;      // (jsr) Ca2+ concentration release compartment
    db u;           // Activation gate u of Ca release from jsr
    db v;           // Activation gate v of Ca release from jsr
    db w;           // Inactivation gate w of Ca release from jsr
    db Itot;        // mA


    db ibarpca;     // max. ca current through sarcolemmal ca pump (ua/uf);
    db kmpca;       // half-saturation concentration of sarcolemmal ca pump (mm);

    //
    // parameters for ikach -> revisar
    //
    db gkach;
    db ach;

    //
    //external concentration
    //

    //db Cao;        //  initial extracellular ca [mm);

    db kmcmdn;
    db kmtrpn;

    //currents
    db Ito;      /* Transient Outward Potassium Current */
    db IKr;      /* Rapidly Activating Potassium Current */
    db IK1;      /* Time-Independent Potassium Current */
    db IKur;     /* Ultra-Rapidly Activating Potassium Current */
    db ikach;    /* Acetylcholine-Activated Potassium Current */
    db INa;      /* Fast Sodium Current (time dependant) */
    db IbNa;     /* Na Background Current */
    db INaca;    /* Sodium-Calcium Exchanger */
    db ICal;     /* Current through L-type Ca Channel */
    db IKs;      /* Slowly Activating Potassium Current */
    db INaK;     /* Sodium-Potassium Pump */
    db IpCa;     /* Sarcolemmal Ca Pump */
    db IbCa;     /* Ca Background Current */
    db Irel;
    db Itr;
    db Iup;
    db Iup_leak;

    //others
    db mo;
    db ho;
    db jo;

    db ireljsr;

    /* Equilibrium Potencial - Nerst Equation */
    db ENa;
    db EK;
    db ENC;
    db ECa;

    /* Parameters Transient Outward Current ito */
    db Gito;
    db ato;
    db iito;

    /* Parameters Na-Ca Exchanger Current inaca */
    db knaca;     // pa/pf;

    db kmnancx;   // na saturation constant for naca exchanger;
    db kmcancx;   // ca saturation factor for naca exchanger;
    db ksatncx;   // saturation factor for naca exchanger;

    /* Parameters Na-K Pump Current inak */
    db INaKmax;  // max. current through na-k pump (ua/uf)  / 1.0933
    db kmnai;    // half-saturation concentration of nak pump (mm);
    db kmko;     // half-saturation concentration of nak pump (mm);


    __device__ __host__ Cell();

    // main functions
    __device__ __host__ void init();
    __device__ __host__ void compute_currents();
    __device__ __host__ void compute_concentrations(db dt);
    __device__ __host__ void compute_gates(db dt);

    /* Ion Current Functions */
    __device__ __host__ void comp_ina ();        /* Calculates Fast Na Current */
    __device__ __host__ void comp_ical ();       /* Calculates Currents through L-Type Ca Channel */
    __device__ __host__ void comp_ikr ();        /* Calculates Rapidly Activating K Current */
    __device__ __host__ void comp_iks ();        /* Calculates Slowly Activating K Current */
    __device__ __host__ void comp_ik1 ();        /* Calculates Time-Independant K Current */
    __device__ __host__ void comp_ikach();       /* Calculates Acetylcholine-sensitive potassium*/
    __device__ __host__ void comp_ikur ();       /* Calculates Ultra-Rapidly activation K Current*/
    __device__ __host__ void comp_ito ();        /* Calculates Transient Outward Current */
    __device__ __host__ void comp_inaca ();      /* Calculates Na-Ca Exchanger Current */
    __device__ __host__ void comp_inak ();       /* Calculates Na-K Pump Current */
    __device__ __host__ void comp_ipca ();       /* Calculates Sarcolemmal Ca Pump Current */
    __device__ __host__ void comp_ibca ();       /* Calculates Ca Background Current */
    __device__ __host__ void comp_ibna ();       /* Calculates Na Background Current */
    __device__ __host__ void comp_irel();        // Compute Ca2+ Release Current From JSR Irel
    __device__ __host__ void comp_itr();         // Time constant of Ca transfer from NSR to JSR(ms)
    __device__ __host__ void comp_iup();         // Compute Ca2+ Uptake Current by NSR Iup
    __device__ __host__ void comp_iupleak();     // Compute Ca2+ Leak Current by the NSR Iup_leak
    __device__ __host__ void comp_itot ();       /* Calculates Total Current */
    __device__ __host__ db getItot(db dt);       /* Return Itot */
    __device__ __host__ void conc_nai(db dt);    /* Calculates new myoplasmic Na ion concentration */
    __device__ __host__ void conc_ki(db dt);     /* Calculates new myoplasmic K ion concentration */
    __device__ __host__ void conc_cai(db dt);    /* Calculates new myoplasmic Ca ion concentration */
    __device__ __host__ void conc_ca_up(db dt);
    __device__ __host__ void conc_ca_rel(db dt);
    __device__ __host__ void gates_ina(db dt);
    __device__ __host__ void gates_ito(db dt);
    __device__ __host__ void gates_ikur(db dt);
    __device__ __host__ void gates_ikr(db dt);
    __device__ __host__ void gates_iks(db dt);
    __device__ __host__ void gates_ical(db dt);
    __device__ __host__ void gates_irel(db dt);
    __device__ __host__ void comp_itot2(db Istim);

    // device main functions
   __device__  __host__ void d_compute_currents();
    // device Ion Current functions


};
