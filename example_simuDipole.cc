#include <iostream>
#include <fstream>
#include <vector>
#include <sys/stat.h>

// Files
#include "events.h"
#include "coverage.h"
#include "angdist.h"
#include "maptools.h"
#include "projmap.h"
#include "simuevents.h"
#include "common.h"
#include "Cl.h"
#include "rayleigh.h"
#include "fitdipole.h"


// ROOT
#include "TRint.h"
#include "TROOT.h"
#include "TStyle.h"

#ifdef gcc323
char* operator+( std::streampos&, char* );
#endif

using namespace std;

void Usage(string myname)
{
  cout << endl;
  cout << " Synopsis : " << endl;
  cout << myname << " -h, --help to obtain this message" << endl;
  cout << myname << endl << endl;

  cout << " Description :" << endl;
  cout << "Events are simulated following a dipolar modulation (amplitude and orientation are fixed in " 
       << myname << ") in addition to a constant acceptance. The coverage map assuming a flat acceptance is "
       << "computed. Using the events and coverage maps, the power spectrum is drawn. A Rayleigh analysis to "
       << "the Right Ascension of the events and a direct fit are also performed." 
       << endl;

  exit(0);
}

int main(int argc, char* argv[])
{
  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                        To start (initialization)                       //
  //                                                                        //
  //////////////////////////////////////////////////////////////////////////// 
  
  // Command line
  if(argc != 1) Usage(argv[0]);

  int fargc = 1;
  string extension;
  TRint *rint = new TRint("simuDipole", &fargc, argv);
  gROOT->SetStyle("Plain");
  gStyle->SetTitleFont(30,"TITLE");
  extension = ".png";

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                Simulation of the dipole & the events                   //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  // Simulation of my dipole
  unsigned int nSide = 64;
  int sizeX = 800, sizeY = 400;
  double decLimit = 25.;
  double longStep = 60., latStep = 30.;
  THealpixMap dipoleMap(nSide);
  
  vector<long> iPix(dipoleMap.NPix());
  vector<double> lPix; // Deg
  vector<double> bPix; // Deg
 
  for(unsigned int i = 0; i < dipoleMap.NPix(); i++) iPix[i] = i;
  dipoleMap.GiveLB(iPix,lPix,bPix);
 
  double raDipole, decDipole, lDipole, bDipole, ampDipole;
  ampDipole = 0.1;
  raDipole = 0.;
  decDipole = -20.;
  radec2gal(raDipole/15., decDipole, &lDipole, &bDipole);
  vector<double> uvDipole = ll2uv(lDipole, bDipole);

  for(unsigned int i = 0; i < dipoleMap.NPix(); i++)
    {
      vector<double> uvPix = ll2uv(lPix[i], bPix[i]);
      dipoleMap[i] = ampDipole*(uvPix[0]*uvDipole[0]+uvPix[1]*uvDipole[1]+uvPix[2]*uvDipole[2]);
    }
  dipoleMap = 1.+dipoleMap;
  TProjMap dipoleMapProj(dipoleMap, sizeX, sizeY, 90.);
  dipoleMapProj.SkyMap("Full Sky Dipole");
  dipoleMapProj.SetPalette(PaletteRGB, 255); 
  dipoleMapProj.ShowGrid(longStep,latStep);
  dipoleMapProj.Save("fullSkyDipole"+extension);
  
  // Simulation of the events : no acceptance law
  double thetaMax = 60.;
  unsigned int nVal = 5000;
  vector<double> thVal(nVal);
  vector<double> pthVal(nVal);
  for(unsigned int i=0; i<nVal; i++)
    {
      thVal[i] = i*180./(nVal-1);
      pthVal[i] = sin(thVal[i]*M_PI/180)*cos(thVal[i]*M_PI/180);
      if (thVal[i] > thetaMax) pthVal[i] = 0;
    }
  double latSite = kConstantsTK::AugerSouthLatitude;
  double lonSite = kConstantsTK::AugerSouthLongitude;

  unsigned int nb = 300000;
  vector<TEvent> simData = SimulateEvents(dipoleMap, nb, thVal, pthVal, latSite, lonSite);
  cout << "Number of events simulated : " << simData.size() << endl;
  ShowLocalCoord(simData, thVal, pthVal); ShowEqCoord(simData), ShowArrivalTimesCoord(simData);

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                       Compute the coverage map                         //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  THealpixMap covMap = GetAnalyticalCoverage(nSide, thetaMax, latSite);
  covMap *= simData.size()*1./covMap.Total();

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                         Compute the events map                         //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  DECLARE_VECTOR(double, lData, simData, fL);
  DECLARE_VECTOR(double, bData, simData, fB);

  THealpixMap evtMap = map_events(nSide, lData, bData);

  TProjMap evtMapProj(evtMap, sizeX, sizeY, decLimit);
  evtMapProj.SkyMap("Events Map");
  evtMapProj.SetPalette(PaletteRGB, 255); 
  evtMapProj.ShowGrid(longStep,latStep);
  evtMapProj.ShowFOV(kConstantsTK::AugerSouthLatitude+60.); 
  evtMapProj.Save("evtMap"+extension);
  
  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                        Plot the Power Spectrum                         //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  // Binning of the Power Spectrum in l.
  const unsigned int lmax = 20;
  vector<unsigned int> lbins;
  for(unsigned int i = 0; i < lmax+2; i++) lbins.push_back(i);
  vector<vector<double> > lvalue = lvalues(lbins);
  vector<double> ErrorCl;
  
  vector<double> Cl = compute_Cl(simData, covMap, evtMap, lmax, lbins, ErrorCl);

  double* lgraph = new double[lvalue[0].size()];
  double* Clgraph = new double[lvalue[0].size()];
  double* lerrorgraph = new double[lvalue[0].size()];
  double* Clerrorgraph = new double[lvalue[0].size()];
  for(unsigned int i = 0; i < lvalue[0].size(); i++)
    {
      lgraph[i] = lvalue[0][i];
      lerrorgraph[i] = lvalue[1][i];
      Clgraph[i] = Cl[i];
      Clerrorgraph[i] = ErrorCl[i];
      if(lgraph[i] == 0) {lgraph[i] = 0.; lerrorgraph[i] = 0.; Clgraph[i] = 0.; Clerrorgraph[i] = 0.;}
      cout.precision(4);
      cout << " l = " << lgraph[i] << " and Cl = " << Clgraph[i] << " +/- " << Clerrorgraph[i] << endl;
    }

  // Plot the Power Spectrum
  string name = "Power Spectrum";
  string Xaxis = "l";
  string Yaxis = "C_{l}";
  string save = "Cl"+extension;
  TCanvas* cPS = new TCanvas("cPS", name.c_str(), 700, 700);
  TGraphErrors* PS = new TGraphErrors(lvalue[0].size(), lgraph, Clgraph, lerrorgraph, Clerrorgraph);
  PlotXY(cPS, PS, 0., lmax+1, name, Xaxis, Yaxis);
  PS->Draw("AP");
  cPS->Update();
  cPS->SaveAs(save.c_str());

  delete [] lgraph;
  delete [] Clgraph;
  delete [] lerrorgraph;
  delete [] Clerrorgraph;

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                            Rayleigh Analysis                           //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  unsigned int nHarmonic = 1;
  TRayleigh* rayleigh = new TRayleigh(simData, covMap, nHarmonic);

  rayleigh->ComputeRAweight();
  rayleigh->ComputeDist();
  rayleigh->DrawEvtsDist();
  rayleigh->DrawCovDist();
  rayleigh->DrawTrue();
  rayleigh->ComputeAmplitude();
  rayleigh->ComputePhase();
  rayleigh->ComputeSignificance();
  rayleigh->ComputeChanceProbability();
  
  double amplitude, phase, significance, probability;
  amplitude = rayleigh->GetAmplitude();
  phase = rayleigh->GetPhase();
  significance = rayleigh->GetSignificance();
  probability = rayleigh->GetChanceProbability();
  cout << "Rayleigh" << endl;
  cout << "--------" << endl;
  cout << " amplitude = " << amplitude << endl;
  cout << " phase = " << phase << endl;
  cout << " probability = " << probability << endl;
  cout << " significance = " << significance << endl << endl;

  //////////////////////////////////////////////////////////////////
  //                                                              //
  //                  Fit my dipole on the sky                    //
  //                                                              //
  //////////////////////////////////////////////////////////////////
  
  TFitDipole* fitDipole = new TFitDipole(evtMap, covMap, amplitude, phase);
  fitDipole->FitProcedure();
  vector<double> parameter = fitDipole->GetParameters();
  vector<double> errorparameter = fitDipole->GetErrorParameters();
  
  double ampFit = parameter[0];
  double raFit = parameter[1];
  double decFit = parameter[2];
  double ampFitError = errorparameter[0];
  double raFitError = errorparameter[1];
  double decFitError = errorparameter[2];

  cout << "Direct Fit" << endl;
  cout << "----------" << endl;  
  cout << "RA : " << raFit << " +/- " << raFitError << " "
       << "dec : " << decFit << " +/- " << decFitError << " "
       << "A : " << ampFit << " +/- " << ampFitError << endl;

  cout << "Program Finished Normally" << endl;
  rint->Run(kTRUE);
}
