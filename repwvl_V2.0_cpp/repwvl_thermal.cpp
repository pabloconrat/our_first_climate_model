
#include <iostream>
#include <vector>
#include <algorithm> 
#include <netcdf>
#include <string>
#include "repwvl_thermal.h"


using namespace std;
using namespace netCDF;

template <typename T> int sgn(T val){
    return (T(0) < val) - (val < T(0));
};



size_t LowerPos(vector<double> &tempsOnLayer, double &currT){
    const size_t P_REFSIZE = tempsOnLayer.size();
    double delta;
    size_t lowerPos=0;
    int sign;

    sign = sgn<double>(tempsOnLayer[0] - currT);

    for(size_t p_ref_pos = 1; p_ref_pos != P_REFSIZE; ++ p_ref_pos ) {

        delta = tempsOnLayer[p_ref_pos] - currT;

        if (sign != sgn(delta)) {
            lowerPos = p_ref_pos - 1;
            return lowerPos;
        }
        else {
            if (p_ref_pos == P_REFSIZE - 1) {
                lowerPos = p_ref_pos - 1;
                return lowerPos;
            }
        }
        sign = sgn<double>(delta);

    }
    return lowerPos;
}



void read_tau (const char *reducedLkpPath, int nLev, vector<double> &plevel, vector<double> &Tvector, double *H20_VMR,
	       double *CO2_VMR, double *O3_VMR, double *N2O_VMR, double *CO_VMR, double *CH4_VMR,
	       double *O2_VMR, double *HNO3_VMR, double *N2_VMR,
	       double ***tau, double **wvl, double **weight, int *nWvl,
	       int T_at_Lev) {
    //Optionally we could use default values for some trace gases and reduce input.
    const int numOfSpecies = 9;
    int nTemp=0;
    
    vector<vector<double>> tauMat;

    if (T_at_Lev==0)  // temperature for layers
      nTemp = nLev - 1;
    else            // temperature for levels
      nTemp = nLev;
    
    double* p = plevel.data();
    double* T = Tvector.data();
    
    double* plevPa = new double[nLev];
    
    for (int ilev=0; ilev<nLev; ilev++) {
      /* convert pressure from hPa to Pa */
      plevPa[ilev] = p[ilev] * 100.0;
    }

    //Convert data in arrays to vectors first.
    vector<double> currentTemp(T, T + nTemp);
    vector<double> currentPress(plevPa, plevPa + nLev);
    vector<double> H20_VMR_VEC(H20_VMR, H20_VMR + nLev);
    vector<double> CO2_VMR_VEC(CO2_VMR, CO2_VMR + nLev);
    vector<double> O3_VMR_VEC(O3_VMR, O3_VMR + nLev);
    vector<double> N2O_VMR_VEC(N2O_VMR, N2O_VMR + nLev);
    vector<double> CO_VMR_VEC(CO_VMR, CO_VMR + nLev);
    vector<double> CH4_VMR_VEC(CH4_VMR, CH4_VMR + nLev);
    vector<double> O2_VMR_VEC(O2_VMR, O2_VMR + nLev);
    vector<double> HNO3_VMR_VEC(HNO3_VMR, HNO3_VMR + nLev);
    vector<double> N2_VMR_VEC(N2_VMR, N2_VMR + nLev);

    // need to reverse order of vectors
    reverse(currentPress.begin(), currentPress.end());
    reverse (currentTemp.begin(),  currentTemp.end());
    reverse (H20_VMR_VEC.begin(),  H20_VMR_VEC.end());
    reverse (CO2_VMR_VEC.begin(),  CO2_VMR_VEC.end());
    reverse (O3_VMR_VEC.begin(),   O3_VMR_VEC.end());
    reverse (N2O_VMR_VEC.begin(),  N2O_VMR_VEC.end());
    reverse (CO_VMR_VEC.begin(),   CO_VMR_VEC.end());
    reverse (CH4_VMR_VEC.begin(),  CH4_VMR_VEC.end());
    reverse (O2_VMR_VEC.begin(),   O2_VMR_VEC.end());
    reverse (HNO3_VMR_VEC.begin(), HNO3_VMR_VEC.end());
    reverse (N2_VMR_VEC.begin(),   N2_VMR_VEC.end());
    
    vector<vector<double>> VMRS;
    VMRS.emplace_back(H20_VMR_VEC);
    VMRS.emplace_back(CO2_VMR_VEC);
    VMRS.emplace_back(O3_VMR_VEC);
    VMRS.emplace_back(N2O_VMR_VEC);
    VMRS.emplace_back(CO_VMR_VEC);
    VMRS.emplace_back(CH4_VMR_VEC);
    VMRS.emplace_back(O2_VMR_VEC);
    VMRS.emplace_back(HNO3_VMR_VEC);
    VMRS.emplace_back(N2_VMR_VEC);
    
    //Variables for getting Data from Lookuptable
    NcVar  xsecVar, VMRs, tempVar, pressVar, tPertVar, currTempVar, currPressVar, wvlVar, weightVar;
    NcFile dataFile(reducedLkpPath,NcFile::read);
    
    const int xsec_nbooks = dataFile.getDim("xsec_nbooks").getSize();
    const int xsec_npages = dataFile.getDim("xsec_npages").getSize();
    const int xsec_nrows = dataFile.getDim("xsec_nrows").getSize();
    const int xsec_ncols = dataFile.getDim("xsec_ncols").getSize();
    const int vmrs_ref_nrows = dataFile.getDim("vmrs_ref_nrows").getSize();
    const int vmrs_ref_ncols = dataFile.getDim("vmrs_ref_ncols").getSize();
    const int t_pert_nelem = dataFile.getDim("t_pert_nelem").getSize();
    const int num_of_nodes = dataFile.getDim("numOfNodes").getSize();
    const double avog = 6.02214076e23; // mol^⁻1
    const double molMassAir = 0.0289647; // (kg/mol)
    const double earthAccel = 9.80665; // (m/s²)

    const int currPnElem = currentPress.size();


    vector<size_t> startp,countp;

    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(0);

    countp.push_back(xsec_nbooks);
    countp.push_back(xsec_npages);
    countp.push_back(1);
    countp.push_back(xsec_ncols);

    xsecVar = dataFile.getVar("xsec");
    VMRs = dataFile.getVar("vmrs_ref");
    tempVar = dataFile.getVar("t_ref");
    pressVar = dataFile.getVar("p_grid");
    tPertVar = dataFile.getVar("t_pert");
    wvlVar = dataFile.getVar("ChosenWvls");
    weightVar = dataFile.getVar("ChosenWeights");


    
    double xsec[xsec_nbooks][xsec_npages][xsec_ncols];
    double vmrs_ref[vmrs_ref_nrows][vmrs_ref_ncols];
    double t_pert[t_pert_nelem];
    double press[vmrs_ref_ncols];
    double t_ref[vmrs_ref_ncols];
    
    double numDens[currPnElem];

    tauMat.resize(currPnElem - 1, vector<double>(xsec_nrows));

    *nWvl = num_of_nodes;
    *wvl    = (double *)  calloc (*nWvl, sizeof(double));
    *weight = (double *)  calloc (*nWvl, sizeof(double));
    *tau    = (double **) calloc (*nWvl, sizeof(double *));
    for (int iwvl=0; iwvl<*nWvl; iwvl++)
      (*tau)[iwvl] = (double *) calloc (nLev-1, sizeof(double));

      
    VMRs.getVar(vmrs_ref);
    tPertVar.getVar(t_pert);
    tempVar.getVar(t_ref);
    pressVar.getVar(press);
    weightVar.getVar(*weight);
    wvlVar.getVar(*wvl);

    vector<double> p_ref;

    for(int i = 0; i != xsec_ncols; ++i){
        p_ref.emplace_back(press[i]);
    }

    vector<double> tempsOnLayer;
    tempsOnLayer.resize(t_pert_nelem);

    double tempXsec;

    size_t lowPosT;
    size_t lowPosP;
    double c0, cP, cT, cPT;
    double delP, delT;
    double midT;
    double midP;
    double midVMR;

    for (int i = 0; i != currPnElem - 1; ++i)
        for (int j = 0; j!= xsec_nrows; ++j)
            tauMat[i][j] = 0;

    for(int lyr = 0; lyr != currPnElem - 1; ++lyr)
        numDens[lyr] = (currentPress[lyr] - currentPress[lyr + 1]) * avog / molMassAir / earthAccel;

    //Calculation of optical thickness
    for (int wvl = 0; wvl != xsec_nrows; ++wvl) {
        startp[2] = xsec_nrows - wvl - 1 ;

        xsecVar.getVar(startp, countp, xsec);

        for(int lyr = 0; lyr < currPnElem - 1; ++lyr){

            for(int spec = 0; spec != numOfSpecies; ++spec){

                midP = (currentPress[lyr + 1] + currentPress[lyr]) / 2;
		if (T_at_Lev == 0)  // temperature for layers
		  midT = currentTemp[lyr];
		else              // temperature for layers
		  midT = (currentTemp[lyr + 1] + currentTemp[lyr]) / 2;

		  
                midVMR = (VMRS[spec][lyr] + VMRS[spec][lyr + 1]) / 2;

                lowPosP = LowerPos(p_ref,midP);


                for(int pertNum = 0; pertNum != t_pert_nelem; ++pertNum){
                    tempsOnLayer[pertNum] = t_ref[lowPosP] + t_pert[pertNum];
                }
                lowPosT = LowerPos(tempsOnLayer,midT);

                // 2D Interpolation in p and T
                c0 = xsec[lowPosT][spec][lowPosP];
                cT = xsec[lowPosT + 1][spec][lowPosP] - c0;
                cP = xsec[lowPosT][spec][lowPosP + 1] - c0;
                cPT = xsec[lowPosT + 1][spec][lowPosP + 1] -cP - cT - c0;
                delT = (midT- tempsOnLayer[lowPosT]) / (tempsOnLayer[lowPosT + 1] - tempsOnLayer[lowPosT]);
                delP = (midP - p_ref[lowPosP]) / (p_ref[lowPosP + 1] - p_ref[lowPosP]);
                tempXsec = c0 + cT * delT + cP * delP + cPT * delT * delP;

		/* tau stored backwards in wavelength dimension! */
                tauMat[currPnElem - 2 - lyr][xsec_nrows - 1 - wvl] += tempXsec * midVMR;
            }
            tauMat[currPnElem - 2 - lyr][xsec_nrows - 1 - wvl] *= numDens[lyr] ;
        }
    }

    // convert Matrix to double **
    if (tauMat.size() != nLev-1)
      cerr << "Error, number of levels inconsistent" << endl;

    if (tauMat[0].size() != *nWvl)
      cerr << "Error, number of wavelengths inconsistent" << endl;

    for (int iwvl=0; iwvl<*nWvl; iwvl++)
      for (int ilyr=0; ilyr<nLev-1; ilyr++)
	(*tau)[iwvl][ilyr] = tauMat[ilyr][iwvl];
    
    delete[] plevPa;
}
