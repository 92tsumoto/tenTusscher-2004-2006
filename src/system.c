#include "syspara.h"

void function(double x[],double f[],double t)
{

	int i;
	comp_reversal_potential(x);
	comp_ina(x);
	comp_ik1(x);
	comp_ito(x);
	comp_ikr(x);
	comp_iks(x);
	comp_ical(x);
	comp_inak(x);
	comp_ikp(x);
	comp_icab(x);
	comp_inab(x);
	comp_inaca(x);
	comp_ipca(x);
	comp_jup(x);
	comp_jleak(x);
	comp_jrel(x);
	comp_jxfer(x);
	comp_concentration(x);

		if(var.deb==1){ 
            printf("%lf %lf %f %f %e\n", x[0],x[1],x[2],x[3],x[4]);
            printf("%lf %lf %f %f %e\n", x[5],x[6],x[7],x[8],x[9]);
            printf("%lf %e %lf %lf %lf\n", x[10],x[11],x[12],x[13],x[14]);
            printf("%lf %e %lf %lf\n", x[15],x[16],x[17],x[18]);
            printf("time=%lf,Istim=%lf\n",t,var.Istim);
            printf("ENa=%lf, EK=%lf ECa=%lf\n", var.Ena,var.Ek,var.Eca);
        }

	var.Ina_total = ina.fast + inab.na + 3.0*inak.inak + 3.0*ncx.j;
	var.Ik_total = ito.ik + ikr.ik + iks.ik + ik1.ik + ikp.k - 2.0*inak.inak + var.Istim;
	var.Ica_total = ipca.ca + icab.ca - 2.0*ncx.j;

	var.Itotal = ina.fast + ito.ik + + ikr.ik + iks.ik + ik1.ik + ical.ica
					+ ncx.j + inak.inak + inab.na + icab.ca + ikp.k + ipca.ca + var.Istim;

	f[0] = -var.Itotal;
	//Fast sodium current
	f[1] = (ina.mss - x[1])/ina.taum; // m
	f[2] = (ina.hss - x[2])/ina.tauh; // h
	f[3] = (ina.jss - x[3])/ina.tauj; // j
	//Transient outward current
	f[4] = (ito.rss - x[4])/ito.taur;
	f[5] = (ito.sss - x[5])/ito.taus;
	// Ikr
	f[6] = (ikr.xr1ss - x[6])/ikr.tauxr1;
	f[7] = (ikr.xr2ss - x[7])/ikr.tauxr2;
	// Iks
	f[8] = (iks.xsss - x[8])/iks.tauxs;
	// LTCC
	f[9] = (ical.dss - x[9])/ical.taud;
	f[10] = (ical.fss - x[10])/ical.tauf;
	f[11] = (ical.f2ss - x[11])/ical.tauf2;
	f[12] = (ical.fcass - x[12])/ical.taufcass;
	// Jrel
	f[13] = -jrel.k2*x[16]*x[13] +jrel.k4*(1.0-x[13]); 
	// [Ca]i
	f[14] = buf.cai*(-var.Ica_total*var.vr2 + (jleak.ca-jup.ca)*var.vr3 + jxfer.ca);
	// [Ca]sr
	f[15] = buf.casr*(jup.ca - jleak.ca - jrel.ca);
	// [Ca]ss
	f[16] = buf.cass*(-ical.ica*var.acap/(var.zca*F*var.vss)+jrel.ca*var.vr4-jxfer.ca*var.vr5);
	// [Na]i
	f[17] = -var.Ina_total*var.vr1;
	// [K]i
	f[18] = -var.Ik_total*var.vr1;

	//for(i=0;i<NN;i++){
	//	printf("f[%d]=%e\n",i,f[i]);
	//}
}
