/* produced by Tsumoto. K 2008.10.27 */

#include <string.h>
#include <stdlib.h>
#include "syspara.h"

FILE *fopen(), *fpin, *fp0, *fp1, *fp2, *fp3;
FILE *fp4, *fp5, *fp6, *fp7, *fp8, *fp9;
int mode = 1;
int P = 2;
int beats = 10;

typedef double Number;

main(argc,argv)
int argc;
char **argv;
{
	int i,w;
	int ii=0;
	double x[NN];
	double t = 0.0;
	double time=0.0;
	double h;
	double v_old,dvdt,dvdt_new;
	double t_stok;
	char *tmpname;
	char cmd[BUFSIZ];
	double tend;

/* Action Potential Duration and Max. Info */
	double *vmax ; // Max. Voltage (mV)
	double *dvdtmax ; // Max. dv/dt (mV/ms)
	double *apd; // Action Potential Duration
	double *toneapd; // Time of dv/dt Max.
	double *ttwoapd; // Time of 90% Repolarization
	double *rmbp; // Resting Membrane Potential
	double *nair; // Intracellular Na At Rest
	double *cair; // Intracellular Ca At Rest
	double *kir ; // Intracellular K At Rest
	double caimax [beats] ; // Peak Intracellular Ca

	vmax=(Number *)calloc(beats,sizeof(Number));
	dvdtmax=(Number *)calloc(beats,sizeof(Number));
	apd=(Number *)calloc(beats,sizeof(Number));
	toneapd=(Number *)calloc(beats,sizeof(Number));
	ttwoapd=(Number *)calloc(beats,sizeof(Number));
	rmbp=(Number *)calloc(beats,sizeof(Number));
	nair=(Number *)calloc(beats,sizeof(Number));
	cair=(Number *)calloc(beats,sizeof(Number));
	kir=(Number *)calloc(beats,sizeof(Number));
	if(vmax==NULL || dvdtmax==NULL || apd==NULL || toneapd==NULL || ttwoapd==NULL 
		|| rmbp==NULL || nair==NULL || cair==NULL || kir==NULL
		) exit(1);

//int i; // Stimulation Counter

	tmpname = "temp";

	sprintf(cmd, "/usr/bin/cpp -P %s > %s", argv[1],tmpname);
	if(system(cmd) == -1){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if((fpin=fopen(tmpname,"r"))==NULL){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if ((fp1 = fopen("para.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp2 = fopen("data.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp3 = fopen("ndata.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}

// parameter inputs
	input_para(fpin);

	if (var.write){
		if ((fp0 = fopen(argv[2],"w"))==NULL){
			fprintf(stderr, "%s cannot open.\n",argv[2]);
			exit(-1);
		}
	}
	if (var.write0){
		if ((fp4 = fopen("ikr.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp5 = fopen("iks.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp6 = fopen("ical.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp7 = fopen("incx.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp8 = fopen("inak.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp9 = fopen("jrel.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
	}

	xhplot(WINDOW, 700.0, 700.0, WHITE);
	xhplot(DIRECT, 0.0, 0.0, WHITE);

	for (ii = 0; ii < var.datas; ii++){
		long j;
		time = 0.0;
		tend = var.tend[ii];
		for (i = 0; i < NN; i++){ 
			x[i] = var.x0[ii][i];
		}
		if (var.half){
			h = 1.0 / var.m;
		}
		else {
			h = 1.0 / var.m;
		}
		h *= var.tsign[ii];
		xddp.line_wid = var.line_wid[ii];
		xhplot(LINEATT,0,0,WHITE);

		
		// initial values input.
		val_consts(fp1);
		printf("exit consts\n");
		printf("cell volume=%e\n",M_PI*var.a*var.a*var.length);
		printf("Vcell=%e\n",var.vcell);
		printf("Cm=%e\n",var.acap);
				
	// initial values input.
		initial_mem();
		printf("exit memory initialization\n");

		printf("Istim=%lf\n",var.Istim_base);

	// Tablize exp functions.	
	printf("start tablization\n");
	make_ExpTable();
	printf("finished tablization\n");

	// Initialization time
		time -= h;
		var.dt = h;
		var.beat = 0;

		if(var.deb==1){	
			printf("%d %lf %lf %lf %lf\n",var.beat,apd[var.beat],toneapd[var.beat],ttwoapd[var.beat],rmbp[var.beat]);
			printf("%lf %lf %f %f %e\n", x[0],x[1],x[2],x[3],x[4]);
			printf("%lf %lf %f %f %e\n", x[5],x[6],x[7],x[8],x[9]);
			printf("%lf %e %lf %lf %lf\n", x[10],x[11],x[12],x[13],x[14]);
			printf("%lf %e %lf %lf\n", x[15],x[16],x[17],x[18]);
			printf("time=%lf,Istim=%lf\n",time,var.Istim);
			printf("dvdtmax[%d]=%lf\n",var.beat,dvdt);
			printf("ENa=%lf, EK=%lf ECa=%lf\n", var.Ena,var.Ek,var.Eca);
			printf("vr1=%lf, vr2=%lf vr3=%lf\n", var.vr1,var.vr2,var.vr3);
		}

		while (var.beat < beats){
			eventloop(fp1,&mode,&P,x);

			for (j = 0; j< (var.m * var.BCL ); j++){
				t = h*j;
				time += h;

				if(var.deb==1){
					if( time-(var.BCL*var.beat+50.0) >= 0.0 && time-(var.BCL*var.beat+50.0) < h ){
						printf("%d %lf %lf %lf %lf\n",var.beat,apd[var.beat],toneapd[var.beat],ttwoapd[var.beat],rmbp[var.beat]);
						printf("%lf %lf %f %f %e\n", x[0],x[1],x[2],x[3],x[4]);
						printf("%lf %lf %f %f %e\n", x[5],x[6],x[7],x[8],x[9]);
						printf("%lf %e %lf %lf %lf\n", x[10],x[11],x[12],x[13],x[14]);
						printf("%lf %e %lf %lf\n", x[15],x[16],x[17],x[18]);
						printf("time=%lf,Istim=%lf\n",time,var.Istim);
						printf("dvdtmax[%d]=%lf\n",var.beat,dvdt);
						printf("ENa=%lf, EK=%lf ECa=%lf\n", var.Ena,var.Ek,var.Eca);
					}
				}

				if( time-(var.BCL*var.beat+50.0) >= 0.0 && time-(var.BCL*var.beat+50.0) < h ){
					apd[var.beat] =0; toneapd[var.beat] =0; ttwoapd[var.beat] =0;
					rmbp[var.beat] =x[0]; nair[var.beat] = x[17]; kir[var.beat] = x[18]; cair[var.beat] = x[14];
					caimax[var.beat] = x[14]; vmax[var.beat] = -90.0; dvdtmax[var.beat] = 0.0;
				}

				if (time-(var.BCL*var.beat+50.0) >= 0.0 && time-(var.BCL*var.beat+50.0) < 1.0){
					var.Istim = var.Istim_base;
				} else {
					var.Istim = 0;
				}

				if (fabs(time) > tend &&  tend != 0.0) break;

				v_old = x[0];

				eular(NN,h,x,t);
				
				dvdt_new = (x[0]-v_old)/h;
				//printf("dvdt_new=%lf\n",dvdt_new);

				if(var.beat>=0){
					if (x[0] > vmax[var.beat] )
						vmax[var.beat] = x[0];
					if (x[14] > caimax[var.beat] )
						caimax[var.beat] = x[14];
					if (dvdt_new > dvdtmax[var.beat] ){
						dvdtmax[var.beat] = dvdt_new;
						toneapd[var.beat] = time;
					}
					if (dvdt_new < 0 && x[0] >= (vmax[var.beat] -0.9*(vmax[var.beat]-rmbp[var.beat]) ) )
						ttwoapd[var.beat] = time;
				}

				if (var.pflag) orbit(&mode,x,dvdt_new);

				if (time>= (beats-5)*var.BCL && time < beats*var.BCL){
					data_out(fp2,time,x);
					//fprintf(fp2,"%lf %lf\n",time,x[0]);
					if(var.write0){
						current(fp4,fp5,fp6,fp7,fp8,fp9,time,x);
					}
				}

				dvdt = dvdt_new;
				var.dvdt=dvdt;

			} // end j-loop
			var.beat++;

			printf("%e ",time);
			for (i=0; i<NN; i++) printf("%e ",x[i]);
			printf("\n");

			draw_p(&mode,P,x,dvdt);
			mouse(&mode,x,dvdt);
			if (fabs(time) > tend &&  tend != 0.0) break;

		} // end while loop

		for(w=0;w<beats;w++){
			apd[w] = ttwoapd [w] -toneapd [w] ;
			fprintf(fp3,"%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%e\t%e\t%g\n",w,
				vmax[w],dvdtmax[w],apd[w],toneapd[w],ttwoapd[w],nair[w],kir[w],cair[w],caimax[w],rmbp[w]);
			printf("%d %lf %lf %lf %lf %lf %lf %lf %e %e %lf\n",w,
				vmax[w],dvdtmax[w],apd[w],toneapd[w],ttwoapd[w],nair[w],kir[w],cair[w],caimax[w],rmbp[w]);
		}
		fclose(fp1);
		fclose(fp2);
		fclose(fp3);
		if(var.write0){
			fclose(fp4);fclose(fp5);fclose(fp6);fclose(fp7);fclose(fp8);fclose(fp9);
		}
		free(vmax);free(dvdtmax);free(apd);free(toneapd);free(ttwoapd);
		free(rmbp);free(nair);free(cair);free(kir);
		closed_mem();
	} // end i-loop
}

