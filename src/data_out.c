#include "syspara.h"

void data_out(FILE *fp2, double t, double u[])
{

	fprintf(fp2,"%e %e %10.9e %e %e %e %e %e %e\n",
		t,u[0],u[14],u[15],u[16],u[17],u[18],var.Itotal,var.dvdt);

}

void current(FILE *fp4, FILE *fp5, FILE *fp6, FILE *fp7, FILE *fp8, FILE *fp9, double t, double u[])
{

	out_ikr(fp4,t,u);
	out_iks(fp5,t,u);
	out_ical(fp6,t,u);
	out_inaca(fp7,t,u);
	out_inak(fp8,t,u);
	out_cicr(fp9,t,u);
	//printf("t=%lf\n",t);

}

void out_ina(double p[])
{
	//ina.fast = ina.Gna*p[1]*p[1]*p[1]*p[2]*p[2]*(p[0] - var.RTF*log(ina.u1/(ina.u2 + p[15])));

}

// Ito Transient Outward Current
void out_ito (double p[])
{

	//ito.ik = ito.Gto*p[5]*p[6]*(p[0] - var.RTF*log(ito.k1/(var.ki + ito.pnato*p[15])));

}


// L-type calcium current
void out_ical(FILE *fp6, double time, double p[])
{
	fprintf(fp6,"%e %e %e %e %e %e\n",time,ical.ica,p[9],p[10],p[11],p[12]);

}

// Rapidly Activating Potassium Current 
void out_ikr (FILE *fp4, double time, double p[])
{
	//ikr.ik = -(ikr.Gkr*ikr.rate*p[7]*(var.Ek - p[0]))/(var.uni + exp((ikr.k1 + p[0])/ikr.k2));
	fprintf(fp4,"%lf %lf\n",time,ikr.ik);

}

// Slowly Activating Potassium Current 
void out_iks (FILE *fp5, double time, double p[])
{
	
	//iks.ik = iks.rate*iks.Gks*p[8]*p[8]*(p[0] - var.RTF*log(iks.k1/(var.ki + var.prnak*p[15])));
	fprintf(fp5,"%lf %lf\n",time,iks.ik);

}

// Inward rectifier potassium current (Ik1)
void out_ik1 (double p[])
{

	//ik1.ik = -(ik1.Gk1*ik1.c5*(var.Ek - p[0])*(var.uni + exp(ik1.c6*(var.Ek - p[0]))))/((var.uni + exp(-ik1.c1*(var.Ek + ik1.c2 - p[0])))*(exp(-ik1.c5*(var.Ek + ik1.c7 - p[0])) + ik1.c8*exp(ik1.c3*(ik1.c4 - var.Ek + p[0]))) + ik1.c5*(var.uni + exp(ik1.c6*(var.Ek - p[0]))));

}

// Sodium-Calcium Exchanger V-S

void out_inaca (FILE *fp7, double time, double p[])
{
	//ncx.j = -(ncx.c1*(ncx.c2*p[12]*exp(var.FRT*ncx.c3*p[0]) - var.cao*p[15]*p[15]*p[15]*exp(var.FRT*ncx.rncx*p[0])))/(var.uni + ncx.ksat*exp(var.FRT*ncx.c3*p[0]));
	fprintf(fp7,"%lf %lf\n",time,ncx.j);

}

// Sodium-Potassium Pump

void out_inak (FILE *fp8, double time, double p[])
{
	//inak.inak= (inak.G*inak.c1)/((var.uni + (inak.c4*sqrt(inak.c4/p[15]))/p[15])*(var.uni + inak.c3*exp(-var.FRT*p[0]) + inak.c2*exp(-var.FRT*inak.c5*p[0])));
	fprintf(fp8,"%lf %lf\n",time,inak.inak);

}

// Sarcolemmal Ca Pump 

void out_ipca (double p[])
{

	//ipca.ca = (ipca.G*p[12])/(ipca.km + p[12]);

}

// Ca Background Current 

void out_icab (double p[])
{
	
	//icab.ca = icab.pcab*(p[0] - (var.RTF*log(var.cao/p[12]))/var.zca);

}

// Na Background Current 

void out_inab (double p[])
{

	//inab.na = inab.pnab*(p[0] - var.RTF*log(var.nao/p[15]));

}

void out_cicr (FILE *fp9, double time, double p[])
{

	fprintf(fp9,"%e %e %e %e %e %e %e\n",time,jrel.ca,jrel.o,p[13],jrel.kcasr,jrel.k1,jrel.k2);

}

