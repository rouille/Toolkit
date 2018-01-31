#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
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
  cout << myname << " <lobe file>" << endl << endl;

  cout << " Description :" << endl;
  cout << "The events are simulated following an acceptance that vary with time. We compute the coverage map "
       << "using the UTC and the JD distributions of the simulated events in order to absorb the acceptance "
       << "modulation. Finally, the Li & Ma map and the power spectrum are computed. The <lobe file> must "
       << "contain two colums: theta (in deg.) and the lobe value normalized to one at maximum. You can easily "
       << "produce one of this file using the compute_lobe executable." << endl;

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
  if(argc != 2) Usage(argv[0]);
  string lobeFile = argv[1];
  if( !CheckFile(lobeFile) ) {cerr << "File: " << lobeFile << " not found." << endl; exit(0);}

  int fargc = 1;
  string extension;
  TRint *rint = NULL;
  rint = new TRint("SimuAcceptance", &fargc, argv);
  extension = ".png";
  gROOT->SetStyle("Plain");
  gStyle->SetTitleFont(30,"TITLE");


  // Read the lobe file
  ifstream ifsLobe(lobeFile.c_str());
  vector<double> thetaLobe, lobe;
  while( !ifsLobe.eof() )
    {
      float th, lo;
      ifsLobe >> th >> lo;
      thetaLobe.push_back(th);
      lobe.push_back(lo);
    }

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //             Simulation of my acceptance law & the events               //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  unsigned int nSide = 64;
  THealpixMap map(nSide);
  // Isotropic full sky
  map = map+1;
 
  // Simulation of the events with an acceptance law
  double thetaMin = 0.;
  double thetaMax = 60.;
  unsigned int nVal = 10000;
  vector<double> thVal(nVal);
  vector<double> pthVal(nVal);
  for(unsigned int i=0; i<nVal; i++)
    {
      thVal[i] = i*180./(nVal-1);
      pthVal[i] = sin(thVal[i]*M_PI/180)*cos(thVal[i]*M_PI/180);
      if (thVal[i]>thetaMax) pthVal[i] = 0;
    }
    
  // JD law to be taken by the simulator
  double minJD, maxJD;
  date2jd(2005, 1, 1, 0., &minJD);
  date2jd(2005, 3, 1, 0., &maxJD);
  unsigned int nbJD = (unsigned int)((maxJD-minJD)*10);
  double xJD, yJD;
  string lawFileJD = "jdLaw.txt";
  ofstream fileJD(lawFileJD.c_str());
  for(unsigned int i=0; i<nbJD; i++)
    {
      xJD = minJD+(maxJD-minJD)*i*1./(nbJD-1);
      yJD = (1.+0.5*sin(xJD/365*2*M_PI))*(1.+(xJD-minJD)/(maxJD-minJD)*2);
      fileJD << setprecision(12) << xJD << " " << yJD << endl;
    }
  fileJD.close();

  // UTC law to be taken by the simulator
  unsigned int nbUTC = 1000;
  double xUTC, yUTC;
  string lawFileUTC = "utcLaw.txt";
  ofstream fileUTC(lawFileUTC.c_str());
  for(unsigned int i=0; i<nbUTC; i++)
    {
      xUTC = 24.*i*1./(nbUTC-1);
      yUTC = 1.+0.6*sin(xUTC/24*2*M_PI);
      fileUTC << xUTC << " " << yUTC << endl;
    }
  fileUTC.close();
  
  double latSite = kConstantsTK::AugerSouthLatitude;
  double lonSite = kConstantsTK::AugerSouthLongitude;
  unsigned int nb = 100000;

  vector<TEvent> simData;
  simData = SimulateEvents(map, nb, thVal, pthVal, latSite, lonSite, lawFileUTC, lawFileJD);
  cout << "Number of events simulated : " << simData.size() << endl;
  ShowLocalCoord(simData, thVal, pthVal); ShowEqCoord(simData), ShowArrivalTimesCoord(simData);

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                       Compute the coverage map                         //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  // Compute the Coverage Map
  unsigned int nBins = 30;
  double chi2Theta = 1.5;

  TCoverage coverage(nSide);
  coverage.SetCoordSystem('G');
  coverage.SetLatitude(latSite);
  coverage.SetLongitude(lonSite);

  // Zenith angle modulation
  coverage.fThetaDist.SetExtension(extension);
  coverage.fThetaDist.SetAngleName("theta");
  DECLARE_VECTOR(double, thetaData, simData, fTheta);
  coverage.fThetaDist.SetData(thetaData);
  coverage.fThetaDist.SetOptions(thetaMin, thetaMax, nBins, chi2Theta);

  // Fitting with a simple geometric function
  coverage.fThetaDist.fAngFitFunc = geosplFunction;
  coverage.fThetaDist.SetNbParameters(4);
  // param 0 : dummy (about the Fermi-Dirac)
  // param 1 : dummy (about the Fermi-Dirac)
  // param 2 : thetamin
  // param 3 : thetamax
  coverage.fThetaDist.fDegreeMax = 4;
  coverage.fThetaDist.fDataMin = thetaMin;
  coverage.fThetaDist.fDataMax = thetaMax;
  coverage.fThetaDist.fParameters[0].SetParameter(eFixed,0);
  coverage.fThetaDist.fParameters[1].SetParameter(eFixed,0);
  coverage.fThetaDist.fParameters[2].SetParameter(eFixed,thetaMin);
  coverage.fThetaDist.fParameters[3].SetParameter(eFixed,thetaMax);
  bool fitok = coverage.fThetaDist.Run();

  if( !fitok )
    {
      cout << "Impossible to fit theta distribution. Try another fitting function. Exiting." << endl;
      exit(0);
    }

  // Take into account the zenith angle distribution to compute the 
  // coverage map
  coverage.CorrectForAngularModulation("theta");

  // Time modulation
  string model = "UTC+JD";
  coverage.fTimeMod.SetAccTimeModel(model);
  DECLARE_VECTOR(double,UTCh,simData,fUTCh);
  DECLARE_VECTOR(double,UTCs,simData,fUTCs);
  double chi2TimeLimit = 5.;
  double timeStep = 7200;
  coverage.fTimeMod.ComputeAccTime(UTCh,UTCs,timeStep,chi2TimeLimit);
  
  // Take into account the time distribution to compute the coverage map
  coverage.CorrectForTimeModulation(model);
   
  
  // VERY IMPORTANT (initializes many constants to save time)
  // To be called once latSite, thetaMax and thetaMin are set
  coverage.ComputeDeclinationLimits();
  cout << "declination limits : " << coverage.GetDecMin() << " " << coverage.GetDecMax() << endl;
  
  // Compute the coverage map
  coverage.ComputeCoverage();
  
  // Put the coverage map into a Healpix map
  THealpixMap covMap = coverage.GetMap();
  // Normalize it to the total number of events
  covMap *= thetaData.size()*1./covMap.Total();
  // Integrate the coverage map in the beam
  THealpixMap intCovMap = covMap.IntBeam(thetaLobe, lobe);


  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                         Compute the events map                         //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  DECLARE_VECTOR(double, lData, simData, fL);
  DECLARE_VECTOR(double, bData, simData, fB);

  THealpixMap evtMap = map_events(nSide, lData, bData);
  THealpixMap intEvtMap = evtMap.IntBeam(thetaLobe, lobe);

  // Assign 0 to pixels that were at zero in the initial coverage map
  long nPix = covMap.size();
  for(long i = 0; i < nPix; i++) if(covMap[i] == 0.)
    {
      intCovMap[i] = 0.;
      intEvtMap[i] = 0.;
    }

  int sizeX = 800, sizeY = 500;
  double decLimit = 25.;
  double longStep = 60., latStep = 30.;

  TProjMap intCovMapProj(intCovMap, sizeX, sizeY, decLimit);
  intCovMapProj.SkyMap("Coverage Map");
  intCovMapProj.SetPalette(PaletteOrange, 255); 
  intCovMapProj.ShowGrid(longStep,latStep);
  intCovMapProj.ShowFOV(kConstantsTK::AugerSouthLatitude+60.); 
  intCovMapProj.Save("covMap"+extension);

  TProjMap intEvtMapProj(intEvtMap, sizeX, sizeY, decLimit);
  intEvtMapProj.SkyMap("Events Map"); 
  intEvtMapProj.SetPalette(PaletteOrange, 255);
  intEvtMapProj.ShowGrid(longStep,latStep);
  intEvtMapProj.ShowFOV(kConstantsTK::AugerSouthLatitude+60.); 
  intEvtMapProj.Save("evtMap"+extension);

  THealpixMap diffMap = intEvtMap-intCovMap;
  TProjMap diffMapProj(diffMap, sizeX, sizeY, decLimit);
  diffMapProj.SkyMap("Difference Map"); 
  diffMapProj.SetPalette(PaletteOrange, 255);
  diffMapProj.ShowGrid(longStep,latStep);
  diffMapProj.ShowFOV(kConstantsTK::AugerSouthLatitude+60.); 
  diffMapProj.Save("diffMap"+extension);

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                        Plot the Power Spectrum                         //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  // Binning of the Power Spectrum in l.
  const unsigned int lmax = 20;
  vector<unsigned int> lbins;
  for(unsigned int i=0; i<lmax+2; i++) lbins.push_back(i);
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

  cout << "Program Finished Normally" << endl;
  rint->Run(kTRUE);
}
