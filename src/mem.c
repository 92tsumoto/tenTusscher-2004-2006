
#include "syspara.h"

typedef double Number;

void initial_mem()
{

	// ina_fast
	ina.Tmss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttaum=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Thss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttauh=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Tjss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttauj=(Number *)calloc(VNMAX,sizeof(Number));
	if( ina.Tmss==NULL || ina.Ttaum==NULL || ina.Thss==NULL || ina.Ttauh==NULL || ina.Tjss==NULL || ina.Ttauj==NULL) exit(1);
	
	// ito
	ito.Trss=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Ttaur=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Tsss=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Ttaus=(Number *)calloc(VNMAX,sizeof(Number));
	if( ito.Trss==NULL || ito.Ttaur==NULL || ito.Tsss==NULL || ito.Ttaus==NULL) exit(1);
	
	// ikr
	ikr.Txr1ss=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.Ttauxr1=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.Txr2ss=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.Ttauxr2=(Number *)calloc(VNMAX,sizeof(Number));
	if( ikr.Txr1ss==NULL || ikr.Ttauxr1==NULL || ikr.Txr2ss==NULL || ikr.Ttauxr2==NULL ) exit(1);

	// iks
	iks.Txsss=(Number *)calloc(VNMAX,sizeof(Number));
	iks.Ttauxs=(Number *)calloc(VNMAX,sizeof(Number));
	if( iks.Txsss==NULL || iks.Ttauxs==NULL ) exit(1);

	// ical
	ical.Tdss=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttaud=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Tfss=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttauf=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Tf2ss=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttauf2=(Number *)calloc(VNMAX,sizeof(Number));
	if( ical.Tdss==NULL || ical.Ttaud==NULL || ical.Tfss==NULL || ical.Ttauf==NULL || ical.Tf2ss==NULL || ical.Ttauf2==NULL) exit(1);

	// inak
	inak.Tknai=(Number *)calloc(VNMAX,sizeof(Number));
	inak.Tknao=(Number *)calloc(VNMAX,sizeof(Number));
	if( inak.Tknai==NULL || inak.Tknao==NULL ) exit(1);
	
	// inaca
	ncx.Thca=(Number *)calloc(VNMAX,sizeof(Number));
	ncx.Thna=(Number *)calloc(VNMAX,sizeof(Number));
	if( ncx.Thca==NULL || ncx.Thna==NULL ) exit(1);
	
	// iks
	ikp.Tu=(Number *)calloc(VNMAX,sizeof(Number));
	if( ikp.Tu==NULL) exit(1);

}


void closed_mem()
{

	free(ina.Tmss); free(ina.Ttaum); free(ina.Thss); free(ina.Ttauh); free(ina.Tjss); free(ina.Ttauj);
	free(ito.Trss); free(ito.Ttaur); free(ito.Tsss); free(ito.Ttaus);
	free(ikr.Txr1ss); free(ikr.Ttauxr1); free(ikr.Txr2ss); free(ikr.Ttauxr2); 
	free(iks.Txsss); free(iks.Ttauxs);
	free(ical.Tdss); free(ical.Ttaud); free(ical.Tfss); free(ical.Ttauf); free(ical.Tf2ss); free(ical.Ttauf2);
	free(inak.Tknai); free(inak.Tknao);
	free(ncx.Thna); free(ncx.Thca);
	free(ikp.Tu);

}

