#include <iostream>
#include <string>

#include "fitdipole.h"
#include "harmotools.h"
#include "maptools.h"

#include "STClibrary.h"

// Events Map
THealpixMap gEvtMap;

// Coverage Map
THealpixMap gCovMap;

// Size of my Map
unsigned int gMapSize;

// Factor of form due to the integration of the coverage map through the beam and its normalization 
// (healpixmap.cc -> l.189).
double gCorrection = 1.;

// Galactic coordinates of the pixels of my map
vector<double> gLpix, gBpix;

TFitDipole::TFitDipole(const THealpixMap & evtMap, const THealpixMap & covMap, const double ampDipole, const double RAdipole, const vector<double> & lobe, const vector<double> & thetaLobe)
  : fLobe(lobe), fThetaLobe(thetaLobe), fAmpDipole(ampDipole), fRAdipole(RAdipole)
{
  gEvtMap = evtMap;
  gCovMap = covMap/covMap.Max();
  gMapSize = evtMap.NPix();
  fMinuit = 0x0;
  Init();
}



TFitDipole::~TFitDipole()
{
  if(fMinuit) delete fMinuit;
}



void TFitDipole::Init()
{
  /* Initialization of gLpix and gBpix */
  vector<long> iPix(gMapSize);
  for(unsigned int i=0; i<gMapSize; i++) iPix[i] = i;
  gEvtMap.GiveLB(iPix, gLpix, gBpix);
  iPix.clear();

  /* Initialization of gCorrection */
  if(fLobe.size() != 0)
    {
      // Correct for the normalization of the map
      double max = *max_element(fLobe.begin(),fLobe.end());
      vector<double> lobeNorm;
      lobeNorm.resize(fLobe.size());
      for(unsigned int i = 0; i < fLobe.size(); i++) lobeNorm[i] = fLobe[i]*1./max*sin(fThetaLobe[i]*DTOR);
      double area = integrate_nc5(fThetaLobe,lobeNorm)*2*M_PI/DTOR;
      double normCorrect = sqrt((double)(gMapSize*area/(4.*M_PI/(DTOR*DTOR))));
      
      // Correct for the integration through the beam
      vector<double> lobe2Theta(fLobe.size()), lobeTheta(fLobe.size());
      for(unsigned int i = 0; i < fLobe.size(); i++) lobe2Theta[i] = fLobe[i]*fLobe[i]*sin(fThetaLobe[i]*DTOR);
      for(unsigned int i=0; i<fLobe.size(); i++) lobeTheta[i] = fLobe[i]*sin(fThetaLobe[i]*DTOR);
      double mWeight = integrate_nc5(fThetaLobe, lobeTheta);
      double varWeight = integrate_nc5(fThetaLobe, lobe2Theta);
      double lobeCorrection = sqrt(varWeight/mWeight);
      
      gCorrection = normCorrect*lobeCorrection;
    }
}



void TFitDipole::FitProcedure()
{
  /* Minimization with ROOT using TMinuit
     Tminuit expects the function to minimize to be always a static external function
     The function pointer member in TMinuit is :
     void (*)(int & npar, double * gin, double & f, double * par, int flag) fFCN */

  // Instantiate Minuit for 3 parameters
  fMinuit = new TMinuit(3);
  
  // Error flag
  int isOK = 0; 

  // Set the address of the minimization function
  fMinuit->SetFCN(FCN);

  // Define parameters
  string parName[3] = {"Amp (%)","RA", "DEC"};
  double parBaselineValue[3] = {fAmpDipole,fRAdipole,0.};
  double parMinValue[3] = {0.,0.,-90.};
  double parMaxValue[3] = {1.,360.,90.};
  double step[3] = {0.001,0.001,0.001};
  fMinuit->mnparm(0,parName[0],parBaselineValue[0],step[0],parMinValue[0],parMaxValue[0],isOK);
  if( isOK != 0 ) {cout << "In FitProcedure() : unable to implement " << parName[0] << endl; exit(0);}
  fMinuit->mnparm(1,parName[1],parBaselineValue[1],step[1],parMinValue[1],parMaxValue[1],isOK);
  if( isOK != 0 ) {cout << "In FitProcedure() : unable to implement " << parName[1] << endl; exit(0);}
  fMinuit->mnparm(2,parName[2],parBaselineValue[2],step[2],parMinValue[2],parMaxValue[2],isOK);
  if( isOK != 0 ) {cout << "In FitProcedure() : unable to implement " << parName[1] << endl; exit(0);}

  // Set Output : 0 for minimum output. 1 for normal output.
  fMinuit->mncomd("SET PRI 0",isOK);
  if( isOK != 0 ) {cout << "In FitProcedure() : unable to execute command SET PRIntout" << endl; exit(0);}

  // No Warnings
  fMinuit->mncomd("SET NOW",isOK);
  if( isOK != 0 ) {cout << "In FitProcedure() : unable to execute command SET NOWarnings" << endl; exit(0);}

  // Set strategy : 2 to improve minimization
  fMinuit->mncomd("SET STR 2",isOK);
  if( isOK != 0 ) {cout << "In FitProcedure() : unable to execute command SET STRategy" << endl; exit(0);}

  // Set Error definition : 1 for Chi square minimization
  fMinuit->mncomd("SET ERR 1", isOK);
  if( isOK != 0 ) {cout << "In FitProcedure() : unable to execute command SET ERRordef" << endl; exit(0);}

  // Minimization itself. Using MIGRAD with 500 iterations. Default tolerance is 0.1.
  fMinuit->mncomd("MIGRAD 500", isOK);
  if( isOK != 0 ) {cout << "In FitProcedure() : unable to execute command MIGRAD" << endl; exit(0);}

  // Get the information on the parameters
  fNParameters = fMinuit->GetNumPars();
  fParameters.resize(fNParameters);
  fParametersErrors.resize(fNParameters);
  for(int i = 0; i < fNParameters; i++) fMinuit->GetParameter(i,fParameters[i],fParametersErrors[i]);
}



void FCN(int& npar, double *gin, double& f, double *par, int iflag) 
{
  double l, b;
  radec2gal(par[1]/15., par[2], &l, &b);
  vector<double> uvDip = ll2uv(l,b);
  THealpixMap modelMap(gEvtMap.NSide());
  for(unsigned int i = 0; i < gMapSize; i++) 
    {
      if(gCovMap[i] != 0)
        {
          vector<double> uvPix = ll2uv(gLpix[i], gBpix[i]);
          modelMap[i] = gCovMap[i]*(1.+par[0]*(uvPix[0]*uvDip[0]+uvPix[1]*uvDip[1]+uvPix[2]*uvDip[2]));
        }
    }
  modelMap = modelMap*(gEvtMap.Total()*1./modelMap.Total());

  THealpixMap chiMap(gEvtMap.NSide());
  double  error;
  f = 0;
  for(unsigned int i = 0; i < gMapSize; i++) 
    {
      if(gCovMap[i] != 0)
        {
          error = sqrt(gCovMap[i])*gCorrection;
          chiMap[i] = (gEvtMap[i]-modelMap[i])/error;
          f += chiMap[i]*chiMap[i];
        }
    }
}
