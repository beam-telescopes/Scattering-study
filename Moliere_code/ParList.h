

// define input array of variables (old version)
/* 
	* par[0]:   Expected beam energy at u,v=0;
	* par[1]:   Beam particle charge
	* par[2]:   Beam particle mass
	* par[3]:   Target material density
	* par[4]:   Target material atomic number
	* par[5]:   Target material atomic weight
	* par[6]:   Thickness
	* par[7]:   Expected angle reconstruction error
	* par[8]:   reco error calibration factor
	* par[9]:   Normalization
	* par[10]:  u coordinate
	* par[11]:  v coordinate
	* par[12]:  u BE gradient
	* par[13]:  v BE gradient
*/
double par[7] = {1, 1, 0.000512, 2.7, 13, 26.98, 0.025};
//Al sample

	struct dataSet {
	 	int E;
		int width;
		int run;
		int run_bkg;
		};
dataSet al[30];
dataSet* initSampleAl(void){
	
	for(int i = 0; i < 5; i++){
		al[i].width = 25;
		al[i+5].width = 50;
		al[i+10].width = 100;
		al[i+15].width = 200;
		al[i+20].width = 1000;
		al[i+25].width = 10000;
	}
// 25mkm

al[0].E = 5;
al[0].run = 7;
al[0].run_bkg = 5;

al[1].E = 4;
al[1].run = 9;
al[1].run_bkg = 4;

al[2].E = 3;
al[2].run = 11;
al[2].run_bkg = 3;

al[3].E = 2;
al[3].run = 12;
al[3].run_bkg = 1; 

al[4].E = 1;
al[4].run = 13;
al[4].run_bkg = 2;
// 50mkm

al[5].E = 1;
al[5].run = 15;
al[5].run_bkg = 2;

al[6].E = 2;
al[6].run = 17;
al[6].run_bkg = 1;

al[7].E = 3;
al[7].run = 18;
al[7].run_bkg = 3;

al[8].E = 4;
al[8].run = 19;
al[8].run_bkg = 4; 

al[9].E = 5;
al[9].run = 20;
al[9].run_bkg = 5;

// 100mkm

al[10].E = 5;
al[10].run = 22;
al[10].run_bkg = 5;

al[11].E = 4;
al[11].run = 23;
al[11].run_bkg = 4;

al[12].E = 3;
al[12].run = 24;
al[12].run_bkg = 3;

al[13].E = 2;
al[13].run = 26;
al[13].run_bkg = 1; 

al[14].E = 1;
al[14].run = 28;
al[14].run_bkg = 2;

// 200mkm

al[15].E = 1;
al[15].run = 30;
al[15].run_bkg = 2;

al[16].E = 2;
al[16].run = 32;
al[16].run_bkg = 1;

al[17].E = 3;
al[17].run = 33;
al[17].run_bkg = 3;

al[18].E = 4;
al[18].run = 35;
al[18].run_bkg = 4; 

al[19].E = 5;
al[19].run = 36;
al[19].run_bkg = 5;

//1 mm 

al[20].E = 5;
al[20].run = 38;
al[20].run_bkg = 5;

al[21].E = 4;
al[21].run = 39;
al[21].run_bkg = 4;

al[22].E = 3;
al[22].run = 41;
al[22].run_bkg = 3;

al[23].E = 2;
al[23].run = 42;
al[23].run_bkg = 1; 

al[24].E = 1;
al[24].run = 44;
al[24].run_bkg = 2;

//10 mm 

al[25].E = 1;
al[25].run = 46;
al[25].run_bkg = 2;

al[26].E = 2;
al[26].run = 47;
al[26].run_bkg = 1;

al[27].E = 3;
al[27].run = 48;
al[27].run_bkg = 3;

al[28].E = 4;
al[28].run = 49;
al[28].run_bkg = 4; 

al[29].E = 5;
al[29].run = 50;
al[29].run_bkg = 5;
return al;
}

