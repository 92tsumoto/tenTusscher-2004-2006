//#ifndef __SYSPARA_H_INCLUDE 
//#define __SYSPARA_H_INCLUDE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "mkl.h"
#include "/home/tsumoto/lib/xhplot.h"

#define NN 19
#define BUF 100
#define NUM 60

//#define R 8314.472		// J/mmol/K
//#define F 96485.33771638995	// C/mol
//#define T 310.0		// K
#define R 8314.4	// J/mol/K
#define F 96485.3	// C/mmol
#define T 310.0		// K

#define VNMAX 400*5+1
#define dvm 5

struct varstruct {

    int datas;
    int line_wid[NUM];
	
	int n;
    double Istim;
	double dIstim;

	// An invariant constant
	double RTonF,RTon2F,RTon4F;

	// Cell tupe
	// 0: Endo, 1: Mid, 2: Epi
	int celltype;

	// Cell Geometry
	double length,a,vcell,acap;
	double vmyo,vsr,vss;
	double vr1,vr2,vr3,vr4,vr5;

	// Ion Valences 
	double zna,zk,zca;

	// Q10
	double K_Q10;	
	// Reversal potential
	double Ena,Ek,Eks,Eca;
	double prnak;
			
	// Total Ion currents 
	double Ina_total;
	double Ik_total;
	double Ica_total;
	double Itotal;

	// leak and xfer
	double vleak, vxfer;

	// Extracellular ion concentrations
	double nao,ko,cao;

	// Base Currnt Stimulus
	double Istim_base;

	// test variable
	double dt,dvdt;
	// Sttimulus parameters
	double BCL;  // Base cycle length = stimulus period
	int beat; // Number of stimulus

    int m;
    int l;

    double x0[NUM][NN];
    double tsign[NUM];
    double tend[NUM];

    int pflag;
    int write, graph;
    int write0;
    int half;
	int deb;

} var;

// Fast sodium currnets
struct inastruct {

	double Gna,fast;
	double *Tmss,*Ttaum,*Thss,*Ttauh,*Tjss,*Ttauj;
	double mss,taum,hss,tauh,jss,tauj;

} ina;

// Transient Outward Current (Ito)
struct itostruct {

	double ik,Gto;
	double rss,taur,sss,taus;
	double *Trss,*Ttaur,*Tsss,*Ttaus;

} ito;

// L-type Calcium channel current (IcaL)
struct icalstruct {

	double dss,taud,fss,tauf,f2ss,tauf2,fcass,taufcass;
	double *Tdss,*Ttaud,*Tfss,*Ttauf,*Tf2ss,*Ttauf2;
	double ica,gca;

} ical;

// Rapid activating potassium current (Ikr)
struct ikrstruct {

	double ik,Gkr,rategkr;
	double xr1,xr1ss,tauxr1,xr2,xr2ss,tauxr2;
	double *Txr1ss,*Ttauxr1,*Txr2ss,*Ttauxr2;

} ikr;

// Slowlactivating potassium current (Iks)
struct iksstruct {

	double ik,Gks;
	double xsss,tauxs;
	double *Txsss,*Ttauxs;
		
} iks;

// Inward rectifier potassium current (Ik1)
struct ik1struct {

	double ik,Gk1,rategk1;
	double k1ss,ak1,bk1;

} ik1;

// Na-K Pump
struct inakstruct {

	double kmna,kmk;
	double knai,knao,*Tknai,*Tknao;
	double fnak,max,inak;

} inak;

// Na-Ca exchanger
struct ncxstruct {

	double hca,hna;
	double *Thca,*Thna;
	double h1,h2,h3,h4;
	double ksat,gamma;
	double kmna,kmca;
	double j,max,a;

} ncx;

// Sarcolemmal Ca Pump
struct ipcastruct {

	double G,km,ca;

} ipca;

// plateau K Current
struct ikpstruct {
	double G,k;
	double u,*Tu;

} ikp;

// Na Background Current
struct inabstruct {
	double G,na;
} inab;

// Ca Background Current
struct icabstruct {
	double G,ca;
} icab;

// Calcium leak via SERCA pump
struct jleakstruct {

	double vleak,ca;

} jleak;

// Calcium uptake via SERCA pump
struct jupstruct {

	double vmaxup,kup,ca;

} jup;

// SR calcium release flux, via RyR (Jrel)
struct jrelstruct {

	double vrel,o,k1_prime,k2_prime,k1,k2,k3,k4,kcasr,EC,maxsr,minsr,ca;

} jrel;

// Ixfer 
struct jxferstruct {

	double vxfer,ca;

} jxfer;

// Ca buffer
struct bufstruct {

	double bufc,kcai,bufsr,kcasr,bufss,kcass;
	double cai,casr,cass;

} buf;

// Translocation of Ca Ions from NSR to JSR
struct jtrstruct {

	double tau,ca;

} jtr;

void val_consts(FILE *);
void make_ExPTable();

void eular(int n,double h,double x[],double t);
void function(double x[],double f[],double t);
void input_para(FILE *);
//void initial_mem(int tMAX);
void initial_mem();
void closed_mem();

void eventloop(FILE *, int *mode, int *P, double m[]);
void orbit(int *mode, double m[], double x2);
void draw_p(int *mode, int P, double x[], double x2);
void mouse(int *mode, double x[], double x2);

void data_out(FILE *, double t, double u[]);
void out_ikr (FILE *, double time, double p[]);
void out_iks (FILE *, double time, double p[]);
void out_ical(FILE *, double time, double p[]);
void out_inaca (FILE *f, double time, double p[]);
void out_inak (FILE *f, double time, double p[]);
void out_cicr (FILE *f, double time, double p[]);

void comp_reversal_potential(double x[]);
void comp_ina(double x[]);
void comp_ito(double x[]);
void comp_ical(double x[]);
void comp_ikr(double x[]);
void comp_iks(double x[]);
void comp_ik1(double x[]);
void comp_inaca(double x[]);
void comp_inak(double x[]);
void comp_ipca(double x[]);
void comp_ikp(double x[]);
void comp_icab(double x[]);
void comp_inab(double x[]);
void comp_jleak(double x[]);
void comp_jup(double x[]);
void comp_jxfer(double x[]);
void comp_jrel(double x[]);
void comp_concentration (double x[]);

main(int argc,char **argv);
