
#include "syspara.h"

void val_consts(FILE *fp1)
{
	int i,w;
	double v_old,dvdt,dvdt_new;

	// Cell Geometry */
	// unknown from ten Tusscher, Noble, Noble, Panfilov AJP paper. 
	//	var.length = 0.01;       // Length of the cell (cm)
	//	var.a = 0.0008;     // Radius of the cell (cm)
	//var.vcell = 1000*M_PI*var.a*var.a*var.length; // Cell Volume:3.801e-5 (uL)
	//var.vcell = 20100; // Cell Volume: (um**3)
	//var.ageo = 2*M_PI*var.a*var.a+2.0*M_PI*var.a*var.length;  // eometric membrane area: 7.671e-5 (cm^2)
	//var.acap = var.ageo*2;     // Capacitive membrane area = Capacitance: 108.57e-6 (uF) --> 108.57 (pF)
	var.acap = 185E-12;          // Capacitance: (pF)
	var.vmyo = 16.404E-12;   // Cytoplasmic volume (pico L)
	var.vsr  = 1.094E-12;        // Sarcoplasmic volume (pico L)
	var.vss  = 0.05468E-12;    	//  Subspace volume (pico L)
	var.vr1 = var.acap/(F*var.vmyo);
	var.vr2 = var.acap/(2.0*F*var.vmyo);
	var.vr3 = var.vsr/var.vmyo;
	var.vr4 = var.vsr/var.vss;
	var.vr5 = var.vmyo/var.vss;

	// Q10
	var.K_Q10 = 3.0;

	// Ion Valences
	var.zna = 1.0;  // Na valence
	var.zk = 1.0;   // K valence
	var.zca = 2.0;  // Ca valence

	// invariant constant
	var.RTonF = R*T/F;
	var.RTon2F = R*T/(var.zca*F);
	var.RTon4F = R*T/(var.zca*var.zca*F*F);

	// Extracellular Concentrations
	var.nao = 140.0;     // Initial Bulk Medium Na (mM)
	var.ko = 5.4;      // Initial Bulk Medium K (mM)
	var.cao = 2.0;     // Initial Bulk Medium Ca (mM)

	// Fast sodium current
	ina.Gna = 14.838;	// (nS/pF).
		
	// Inward rectifier K current: Ik1
	ik1.Gk1 = 5.405; // (nS/pF)
	ik1.rategk1 = sqrt(var.ko/5.4);

	// Transient outward current
	if(var.celltype==0){				// Endo
		ito.Gto = 0.073;	// (nS/pF).
	} else if(var.celltype==1){			// Mid
		ito.Gto = 0.294;	// (nS/pF).
	} else if(var.celltype==2){			// Epi
		ito.Gto = 0.294;	// (nS/pF).
	}

	// Rapid delayed rectifier potassium current (Ikr)
		ikr.Gkr = 0.153;  //(nS/pF)
		ikr.rategkr = sqrt(var.ko/5.4);

	// Slow delayed rectifier potassium current (Iks)
	if(var.celltype==0){				//Endo
		iks.Gks = 0.392;  //(nS/pF)
	} else if(var.celltype==1){			// Mid
		iks.Gks = 0.098;  //(nS/pF)
	} else if(var.celltype==2){			// Epi
		iks.Gks = 0.392;  //(nS/pF)
	}
	// L-type calcium current
		ical.gca = 3.980E-5;	// (cm/ms/microF)

	// Na/K Pump
		inak.kmna = 40;		// (mM)
		inak.kmk = 1.0;		// (mM)
		inak.max = 2.724;	// (pA/pF)

	// Sodium-Calcium Exchanger V-S
		ncx.max = 1000.0;	// (pA/pF)
		ncx.kmna = 87.5;	// (mM)
		ncx.kmca = 1.38;	// (mM)
		ncx.ksat = 0.1;
		ncx.gamma = 0.35;
		ncx.a = 2.5;

	// Sarcolemmal Ca Pump
		ipca.G = 0.1238;		// (nS/pF)
		ipca.km = 0.0005;		// Half-saturation concentration of sarcolemmal Ca pump (mM)

	// plateau K Current 
		ikp.G = 0.0146;		// Max. conductance of plateau K current(nS/pF)
	// Ca Background Current 
		icab.G = 0.000592; // (nS/pF)
	// Na Background Current 
		inab.G = 0.000290; //(nS/pF)

	// calcium uptake via SERCA pump (Jup)
		jup.vmaxup = 0.006375;	// Maximal Iup conductance (mM)
		jup.kup =    0.00025;	// (mM)

	// SR calcium release flux, via RyR (Jrel)
		jrel.vrel = 40.8;	// (mM/ms)
		jrel.k1_prime = 0.15;	// R to O and RI to I Irel transition rate (1/mM^2/ms)
		jrel.k2_prime = 0.045;	// O to I and R to RI Irel transition rate (1/mM/ms)
		jrel.k3 = 0.060;		// O to R and I to RI Irel transition rate (1/ms)
		jrel.k4 = 0.000015;		// I to O and RI to I Irel transition rate (1/ms)
		jrel.EC = 1.5;	// Ca_sr half-saturation constant of kcasr (mM)
		jrel.maxsr = 2.5;
		jrel.minsr = 1.0;
	// leakage of Ca from SR (I_leak)
		jleak.vleak = 0.00036;	// Maximal Ileak conductance (mM/ms)
		jxfer.vxfer = 0.0038;	// Maximal Ixfer conductance (mM/ms)

	// Myoplasmic Ca Ion Concentration Changes 
		buf.bufc = 0.2;   		// Max. [Ca] buffer concentration (mM)
		buf.kcai = 0.001;   	// Cai half-saturation constant for cytoplasmic buffer (mM)
		buf.bufsr = 10.0;     	// Total sarcoplasmic buffer concentration (mM)
		buf.kcasr = 0.3;  		// Ca_sr half-saturation constant for sarcoplasmic buffer (mM)
		buf.bufss = 0.4;   		// Total subspace buffer concentration (mM)
		buf.kcass = 0.00025;	// Ca_ss half-saturation constant for subspace buffer (mM)

}

