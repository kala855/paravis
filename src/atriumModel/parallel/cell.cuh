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
//intracellular volume
#define Vi 13668.0                // micro_m^3
#define Vup 1109.52               // micro_m^3
#define Vrel 96.48                // micro_m^3
//----------------------------cell geometry--------------------------------
#define LENGTHCELL 0.01           // cm
#define RADIUSCELL 0.0008         // cm

#define Ri 200.0                  // Intracell Resistance - 200 Ohm x cm = 0.2 K Ohms x cm, Tesis Catalina, page 99 y 115 Resistivida no resistencia
#define Rix 0.2                   // Specific resistance of intracell liquid
#define Riy 0.2                   // Specific resistance of intracell liquid

class Cell{
public:
    db pi;

    // External concentration
    db Ko;          // initial extracellular k [mm);
    db Noa;         // initial extracellular na [mm);
    db Coa;         // initial extracellular ca [mm);

    /* Maximal  conductances  nS/pF;*/
    db GNa;
    db GK1;
    db Gto;
    db GKr;
    db GKs;
    db GCaL;
    db GbCa;
    db GbNa;

    /* Maximal currents */
    db INaK_max;      // Max. current through Na-K pump (pA/uF)
    db INaCa_max;     // pA/pF;
    db IpCa_max;      // max. ca current through sarcolemmal ca pump (ua/uf);
    db Kq10;          // Temperature scaling factor for IKur and Ito kinetics
    db gamma;         // Voltage dependance parameter for INaCa;

    /* Half-saturation constant for currents */
    db KmNai;       // half-saturation concentration of nak pump (mm);
    db KmKo;        // half-saturation concentration of nak pump (mm);
    db KmNa;        // na saturation constant for naca exchanger;
    db KmCa;        // ca saturation factor for naca exchanger;

    db ksat;        // saturation factor for naca exchanger;

    // Ion Valences
    db zna;         // Na valence
    db zk;          // K valence
    db zca;         // Ca valence

    // Myoplasmic Ca Ion Concentration Changes
    db Csqn_max;     // Max. [Ca] buffered in CSQN (mM)
    db Km_csqn;      // Equalibrium constant of buffering for CSQN (mM)
    db Cmdn_max;     // Max. [Ca] buffered in CMDN (mM)
    db Trpn_max;     // Max. [Ca] buffered in TRPN (mM)
    db kmcmd;        // Equalibrium constant of buffering for CMDN (mM)
    db Kmtrpn;       // Equalibrium constant of buffering for TRPN (mM)
    db Iup_max;      // Max. current through iup channel (mM/ms)
    db caupm;        //  Max. [Ca] in NSR (mM)
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
