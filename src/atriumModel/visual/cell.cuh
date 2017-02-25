//header cell class
#include <math.h>
#define db double
//----------------------------parameters---------------------------------------
#define R 8.3143                  // universal gas constant [j/kmol*k];
#define TEMP 310.0                // temperature [k];
#define F 96.4867                 // faraday"s constant [c/mol]
#define RTF (R*TEMP)/F            // j/c
#define invRTF 1.0/RTF
#define CAP 100.0                 // capacitancia valor para courtemanche? 0.1
//-----------------------------intracellular volume--------------------------
#define Vi 13668.0                // micro_m^3
#define Vup 1109.52               // micro_m^3
#define Vrel 96.48                // micro_m^3
//----------------------------cell geometry--------------------------------
#define LENGTHCELL 0.01           // cm
#define RADIUSCELL 0.0008         // cm

#define Ri 200.0                  // Intracell Resistance - 200 Ohm x cm = 0.2 K Ohms x cm, Tesis Catalina, page 99 y 115 Resistivida no resistencia
#define Rix 0.2                   // Specific resistance of intracell liquid
#define Riy 0.2                   // Specific resistance of intracell liquid

//-----------------------------External Concentration-----------------------
#define Ko 5.4                    // Extracellular K concentration [mm]
#define Nao 140.0                 // Extracellular Na concentration [mm]
#define Cao 1.8                   // Extracellular Ca concentration [mm]

//------------------------------Maximal Conductances------------------------
#define GNa 7.8
#define GK1 0.09
#define Gto 0.1652
#define GKr 0.0294
#define GKs 0.129
#define GCaL 0.1238
#define GbCa 0.00113
#define GbNa 0.000674

//-------------------------------Maximal Currents------------------------------
#define INaK_max 0.6                // Max Current through Na-K pump [pA/pF]
#define INaCa_max 1600.0            // Max. INaCa [pA/pF]
#define IpCa_max 0.275              // Max. ca current through sarcolemmal ca pump [pA/pF]
#define Kq10 3.0                    // Temperature Scaling factor for IKur and Ito kinetics
#define gamma 0.35                  // Voltage dependance for INaCa

//-------------------------Half-saturation constant for currents-----------------
#define KmNai 10.0                  // Nai half-saturation constant of INaK[mm]
#define KmKo 1.5                    // Ko half-saturation constant of
#define KmNa 87.5                   // Nao half-saturation constant of
#define KmCa 1.38                   // Cao half-saturation constant of
#define ksat 0.1                    // Saturation factor for INaCa

//-------------------------Ion Valences------------------------------------
#define zna 1.0                     // Na valence
#define zk 1.0                      // K valence
#define zca 2.0                     // Ca valence

//-------------------------Myoplasmic Ca Ion Concentration Changes--------------
#define Csqn_max 10.0               // Total calsequestrin concentration in SR release compartment [mm]
#define Km_csqn 0.8                 // Ca_rel half-saturation constant of Iup [mm]
#define Cmdn_max 0.050              // Total calmodulin concentration in myoplasm [mm]
#define Trpn_max 0.070              // Total troponin concentration in myoplasm [mm]
#define kmcmdn 0.00238              // Cai Half-saturation constant for calmodulin [mm]
#define Kmtrpn 0.0005               // Cai half-saturation constant for troponin [mm]
#define Iup_max 0.005                      // Maximal Iup [mm/ms]

#define gf 0.075
#define ef -22

class Cell{
public:
    db pi;

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
    db ifunny;

    /* Equilibrium Potencial - Nerst Equation */
    db ENa;
    db EK;
    db ENC;
    db ECa;


    __device__ __host__ Cell();

    // main functions
    __device__ __host__ void init();
    __device__ __host__ void compute_currents(db dt);
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
    __device__ __host__ void comp_if(db dt);

    // device main functions
   __device__  __host__ void d_compute_currents();
    // device Ion Current functions


};
