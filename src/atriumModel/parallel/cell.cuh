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



class Cell{
public:
    //------------    parameters   ----------------

    //--- constants -----
   // db R;           // universal gas constant [j/kmol*k);
    //db T;           // temperature [k);
    //db F;           // faraday"s constant [c/mol) ;
   // db RTF;         // j/c
   // db invRTF;

    //--- capacitance ---
   // db Cap;         // valor para courtemanche? 0.1 ;

    //intracellular volume
   // db Vi;          // micro_m^3;
   // db Vup;         // micro_m^3;
   // db Vrel;        // micro_m^3;

    //cell geometry
   // db l;           // length of the cell (cm);
   // db a;           // radius of the cell (cm);
   // db pi;

   /* db Ri;
    db Rix;         // Specific resistance of intracell liquid
    db Riy;         // Specific resistance of intracell liquid

    // External concentration
    db Ko;          // initial extracellular k [mm);
    db Noa;         // initial extracellular na [mm);
    db Coa; */        // initial extracellular ca [mm);

    /* Maximal  conductances  nS/pF;*/
    /*db GNa;
    db GK1;
    db Gto;
    db GKr;
    db GKs;
    db GCaL;
    db GbCa;
    db GbNa;*/

    /* Maximal currents */
    /*db INaK_max;      // Max. current through Na-K pump (pA/uF)
    db INaCa_max;     // pA/pF;
    db IpCa_max;      // max. ca current through sarcolemmal ca pump (ua/uf);
    db Kq10;          // Temperature scaling factor for IKur and Ito kinetics
    db gamma;*/         // Voltage dependance parameter for INaCa;

    /* Half-saturation constant for currents */
    /*db KmNai;       // half-saturation concentration of nak pump (mm);
    db KmKo;        // half-saturation concentration of nak pump (mm);
    db KmNa;        // na saturation constant for naca exchanger;
    db KmCa;        // ca saturation factor for naca exchanger;

    db ksat;*/        // saturation factor for naca exchanger;

    // Ion Valences
   /* db zna;         // Na valence
    db zk;          // K valence
    db zca; */        // Ca valence

    // Myoplasmic Ca Ion Concentration Changes
    /*db Csqn_max;     // Max. [Ca] buffered in CSQN (mM)
    db Km_csqn;      // Equalibrium constant of buffering for CSQN (mM)
    db Cmdn_max;     // Max. [Ca] buffered in CMDN (mM)
    db Trpn_max;     // Max. [Ca] buffered in TRPN (mM)
    db kmcmd;        // Equalibrium constant of buffering for CMDN (mM)
    db Kmtrpn;       // Equalibrium constant of buffering for TRPN (mM)
    db Iup_max;      // Max. current through iup channel (mM/ms)
    db caupm;*/        //  Max. [Ca] in NSR (mM)
    //db tautr;      // Time constant of Ca transfer from NSR to JSR (ms)

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

    db Nao;        //  initial extracellular na [mm);
    db Cao;        //  initial extracellular ca [mm);

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
