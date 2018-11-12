#include "syspara.h"

void comp_ina(double x[])
{
	//MKL_INT iV=0;
	int iV=0;
	double V1,V2,d1,d2;
	
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;
	//printf("iV=%d,V1=%f,V2=%f,d1=%f,d2=%f\n",iV,V1,V2,d1,d2);

	ina.mss = ina.Tmss[iV]*d2 + ina.Tmss[iV+1]*d1;
	ina.taum = ina.Ttaum[iV]*d2 + ina.Ttaum[iV+1]*d1;
	ina.hss = ina.Thss[iV]*d2 + ina.Thss[iV+1]*d1;
	ina.tauh = ina.Ttauh[iV]*d2 + ina.Ttauh[iV+1]*d1;
	ina.jss = ina.Tjss[iV]*d2 + ina.Tjss[iV+1]*d1;
	ina.tauj = ina.Ttauj[iV]*d2 + ina.Ttauj[iV+1]*d1;
	ina.fast = ina.Gna*(x[0]-var.Ena)*x[1]*x[1]*x[1]*x[2]*x[3];
}

// Inward rectifier potassium current (Ik1)
void comp_ik1 (double x[])
{
       
	ik1.ak1 = 0.1/(1.0+exp(0.06*(x[0]-var.Ek-200.0)));
	ik1.bk1 = (3.0*exp(0.0002*(x[0]-var.Ek+100.0))+exp(0.1*(x[0]-var.Ek-10.0)))/(1.0+exp(-0.5*(x[0]-var.Ek)));

	ik1.k1ss = ik1.ak1/(ik1.ak1 + ik1.bk1);

	ik1.ik = ik1.Gk1*ik1.rategk1*ik1.k1ss*(x[0]-var.Ek);

}

// Ito Transient Outward Current
void comp_ito (double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ito.rss = ito.Trss[iV]*d2 + ito.Trss[iV+1]*d1;
	ito.taur = ito.Ttaur[iV]*d2 + ito.Ttaur[iV+1]*d1;

	ito.sss = ito.Tsss[iV]*d2 + ito.Tsss[iV+1]*d1;
	ito.taus = ito.Ttaus[iV]*d2 + ito.Ttaus[iV+1]*d1;

	ito.ik = ito.Gto*(x[0]-var.Ek)*x[4]*x[5];

}

// Rapidly Activating Potassium Current 
void comp_ikr (double x[])
{
	MKL_INT iV=0;	
	double V1,V2,d1,d2;
	
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ikr.xr1ss = ikr.Txr1ss[iV]*d2 + ikr.Txr1ss[iV+1]*d1;
	ikr.tauxr1 = ikr.Ttauxr1[iV]*d2 + ikr.Ttauxr1[iV+1]*d1;
	ikr.xr2ss = ikr.Txr2ss[iV]*d2 + ikr.Txr2ss[iV+1]*d1;
	ikr.tauxr2 = ikr.Ttauxr2[iV]*d2 + ikr.Ttauxr2[iV+1]*d1;

	ikr.ik = ikr.Gkr*ikr.rategkr*(x[0]-var.Ek)*x[6]*x[7];

}

// Slowly Activating Potassium Current 
void comp_iks (double x[])
{
	
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	iks.xsss = iks.Txsss[iV]*d2 + iks.Txsss[iV+1]*d1;
	iks.tauxs = iks.Ttauxs[iV]*d2 + iks.Ttauxs[iV+1]*d1;
	iks.ik = iks.Gks*x[8]*x[8]*(x[0]-var.Eks);

}


// L-type calcium current
void comp_ical(double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	// VDA
	ical.dss = ical.Tdss[iV]*d2 + ical.Tdss[iV+1]*d1;
	ical.taud = ical.Ttaud[iV]*d2 + ical.Ttaud[iV+1]*d1;
	// VDI 
	ical.fss = ical.Tfss[iV]*d2 + ical.Tfss[iV+1]*d1;
	ical.tauf = ical.Ttauf[iV]*d2 + ical.Ttauf[iV+1]*d1;
	// VDI2 
	ical.f2ss = ical.Tf2ss[iV]*d2 + ical.Tf2ss[iV+1]*d1;
	ical.tauf2 = ical.Ttauf2[iV]*d2 + ical.Ttauf2[iV+1]*d1;
	// CDI 
	ical.fcass = 0.4+0.6/(1.0+(x[16]/0.05)*(x[16]/0.05));
	ical.taufcass = 2.0+80.0/(1.0+(x[16]/0.05)*(x[16]/0.05));
	//ical.taufca = 2.0;

	ical.ica =ical.gca*x[9]*x[10]*x[11]*x[12]*(x[0]-15.0)/var.RTon4F*(0.25*x[16]*exp((x[0]-15.0)/var.RTon2F)-var.cao)/(exp((x[0]-15.0)/var.RTon2F)-1.0);

}

// Na-Ca Exchanger NCX
void comp_inaca (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ncx.hca=ncx.Thca[iV]*d2 + ncx.Thca[iV+1]*d1;
	ncx.hna=ncx.Thna[iV]*d2 + ncx.Thna[iV+1]*d1;

	ncx.h1 = x[17]*x[17]*x[17]*var.cao;
	ncx.h2 = var.nao*var.nao*var.nao*x[14]*ncx.a;
	ncx.h3 = ncx.kmna*ncx.kmna*ncx.kmna+var.nao*var.nao*var.nao;
	ncx.h4 = ncx.kmca+var.cao;

	ncx.j = ncx.max*(ncx.hca*ncx.h1-ncx.hna*ncx.h2)/(ncx.h3*ncx.h4*(1.0+ncx.ksat*ncx.hna));

}

// Na-K Pump
void comp_inak (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	inak.knai = inak.Tknai[iV]*d2 + inak.Tknai[iV+1]*d1;
	inak.knao = inak.Tknao[iV]*d2 + inak.Tknao[iV+1]*d1;

	inak.fnak = 1.0/(1.0 + inak.knai + inak.knao);

	inak.inak = inak.max*inak.fnak*(var.ko*x[17]/((var.ko + inak.kmk)*(x[17]+inak.kmna)));

}

// Sarcolemmal Ca Pump 
void comp_ipca (double x[])
{
	ipca.ca = ipca.G*x[14]/(ipca.km + x[14]);
}

// Plateau K Current
void comp_ikp (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ikp.u = ikp.Tu[iV]*d2 + ikp.Tu[iV+1]*d1;
	ikp.k = ikp.G*ikp.u*(x[0] - var.Ek);
}

// Ca Background Current 
void comp_icab (double x[])
{
	icab.ca = icab.G*(x[0] - var.Eca);
}

// Na Background Current 
void comp_inab (double x[])
{
	inab.na = inab.G*(x[0] - var.Ena);
}

void comp_jleak (double x[])
{
	// Ca2+ leak current from the SR
	jleak.ca = jleak.vleak*(x[15]-x[14]);

}

void comp_jup (double x[])
{
	// Ca2+ uptake in the SR
	jup.ca = jup.vmaxup/(1.0+(jup.kup/x[14])*(jup.kup/x[14]));
}

void comp_jxfer (double x[])
{
	jxfer.ca = jxfer.vxfer*(x[16]-x[14]);
}

void comp_jrel (double x[])
{
	
	jrel.kcasr = jrel.maxsr - ((jrel.maxsr - jrel.minsr)/(1.0+(jrel.EC/x[15])*(jrel.EC/x[15])));
	jrel.k1    = jrel.k1_prime/jrel.kcasr;
	jrel.k2    = jrel.k2_prime*jrel.kcasr;

	jrel.o = (jrel.k1*x[16]*x[16]*x[13])/(jrel.k3+jrel.k1*x[16]*x[16]);
	jrel.ca = jrel.vrel*jrel.o*(x[15]-x[16]);

}

void comp_concentration (double x[])
{
	buf.cai =  1.0/(1.0+buf.bufc*buf.kcai/((buf.kcai+x[14])*(buf.kcai+x[14])));
	buf.casr = 1.0/(1.0+buf.bufsr*buf.kcasr/((buf.kcasr+x[15])*(buf.kcasr+x[15])));
	buf.cass = 1.0/(1.0+buf.bufss*buf.kcass/((buf.kcass+x[16])*(buf.kcass+x[16])));

}


// Reversal potentials */

void comp_reversal_potential(double x[])
{
	var.Ena = var.RTonF*log(var.nao/x[17]);
	var.Ek = var.RTonF*log(var.ko/x[18]);
	var.Eks = var.RTonF*log((var.ko+var.prnak*var.nao)/(x[18]+var.prnak*x[17]));
	var.Eca = var.RTon2F*log(var.cao/x[14]);
	
	//printf("Ena=%lf, Ek=%lf, Eks=%lf\n",var.Ena,var.Ek,var.Eks);
}

