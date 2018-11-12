#include "syspara.h"

void make_ExpTable()
{

	int vindex,kiindex;
	double v,ki;
	double am,bm,ah,bh,aj,bj;
	double ao,bo,ai,bi;
	double ua,ub,ia,ib;
	double axr,bxr;
	double axs,bxs;
    
	for(vindex=0;vindex<VNMAX;vindex++){

        v = (double)vindex/dvm-200.0;
		
		/** m-gate for ina **/
		ina.Tmss[vindex] = 1.0/((1.0+exp((-v-56.86)/9.03))*(1.0+exp((-v-56.86)/9.03)));
		ina.Ttaum[vindex] = (1.0/(1.0+exp((-60.0-v)/5)))*(0.1/(1.0+exp((v+35.0)/5.0))+0.1/(1.0+exp((v-50)/200)));

		/** h-gate for ina **/
		if(v<-40.0){
			ah = 0.057*exp(-(v+80.0)/6.8);
			bh = 2.7*exp(0.079*v)+3.1E+5*exp(0.3485*v);
		} else {
			ah = 0.0;
			bh = 0.77/(0.13*(1.0+exp(-(v+10.66)/11.1)));
		}
		ina.Thss[vindex] = 1.0/((1.0+exp((v+71.55)/7.43))*(1.0+exp((v+71.55)/7.43)));
		ina.Ttauh[vindex] = 1.0/(ah+bh);

		/** j-gate for ina **/
		if(v<-40.0){
			aj = (v+37.78)*(-2.5428E+4*exp(0.2444*v)-6.948E-6*exp(-0.04391*v))/(1.0+exp(0.311*(v+79.23)));
			bj = 0.02424*exp(-0.01052*v)/(1.0+exp(-0.1378*(v+40.14)));
		} else {
			aj = 0.0;
			bj = 0.6*exp(0.057*v)/(1.0+exp(-0.1*(v+32.0)));
		}
		ina.Tjss[vindex] = 1.0/((1.0+exp((v+71.55)/7.43))*(1.0+exp((v+71.55)/7.43)));
		ina.Ttauj[vindex] = 1.0/(aj+bj);

		// ito
		ito.Trss[vindex] = 1.0/(1.0+exp((v+20.0)/6.0));
		ito.Ttaur[vindex] = 9.5*exp(-(v+40.0)*(v+40.0)/1800.0)+0.8;
		
		if(var.celltype==0){								// for Endo
			ito.Tsss[vindex] = 1.0/(1.0+exp((v+28.0)/5.0));
			ito.Ttaus[vindex] = 1.000*exp(-(v+67.0)*(v+67.0)/1000)+8.0;
		} else if(var.celltype==1 || var.celltype==2){		// for Mid and Epi
			ito.Tsss[vindex] = 1.0/(1.0+exp((v+20.0)/5.0));
			ito.Ttaus[vindex] = 85.0*exp(-(v+45.0)*(v+45.0)/320)+5.0/(1.0+exp((v-20.0)/5.0))+3.0;
		}

		// for ikr 
		ikr.Txr1ss[vindex] = 1.0/(1.0+exp((-26.0-v)/7.0));
		ikr.Ttauxr1[vindex] = (450/(1.0+exp((-45.0-v)/10.0)))*(6.0/(1.0+exp((v+30.0)/11.5)));

		ikr.Txr2ss[vindex] = 1.0/(1.0+exp((v+88.0)/24.0));
		ikr.Ttauxr2[vindex] = (3.0/(1.0+exp((-60-v)/20.0)))*(1.12/(1.0+exp((v-60.0)/20.0)));

		// for iks 
		iks.Txsss[vindex] = 1.0/(1.0+exp((-5.0-v)/14.0));
		iks.Ttauxs[vindex] = (1400.0/sqrt(1.0+exp((5.0-v)/6.0)))*(1.0/(1.0+exp((v-35.0)/15.0)))+80.0;

		// for ical
		ical.Tdss[vindex] = 1.0/(1.0+exp((-8.0-v)/7.5));
		ical.Ttaud[vindex] = (0.25+1.4/(1.0+exp((-35.0-v)/13.0)))*(1.4/(1.0+exp((v+5.0)/5.0)))+(1.0/(1.0+exp((50.0-v)/20.0)));

		ical.Tfss[vindex] = 1.0/(1.0+exp((v+20.0)/7.0));
		ical.Ttauf[vindex] = 1102.5*exp(-((v+27.0)*(v+27.0)/225.0))+200.0/(1.0+exp((13.0-v)/10.0))+180.0/(1.0+exp((v+30.0)/10.0))+20.0;

		ical.Tf2ss[vindex] = 0.67/(1.0+exp((v+35.0)/7.0))+0.33;
		ical.Ttauf2[vindex] = 600.0*exp(-((v+25.0)*(v+25.0)/170.0))+31.0/(1.0+exp((25.0-v)/10.0))+16.0/(1.0+exp((v+30.0)/10.0));

		// inak 
		inak.Tknai[vindex] = 0.1245*exp((-0.1*v)/var.RTonF);
		inak.Tknao[vindex] = 0.0353*exp((-1.0*v)/var.RTonF);

		// incx
		ncx.Thca[vindex] = exp(ncx.gamma*v/var.RTonF);
		ncx.Thna[vindex] = exp((ncx.gamma-1.0)*v/var.RTonF);
		
		// ikp
		ikp.Tu[vindex] = 1.0/(1.0+exp((25.0-v)/5.98));

	} //for i loop end

}
