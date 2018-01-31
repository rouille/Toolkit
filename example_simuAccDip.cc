#include <iostream>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <iomanip>

// Files
#include "events.h"
#include "coverage.h"
#include "angdist.h"
#include "maptools.h"
#include "projmap.h"
#include "simuevents.h"
#include "common.h"
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
  cout << myname << " <lobe file>" << endl << endl;
  
  cout << " Description :" << endl;
  cout << "Events are simulated following a dipolar modulation in addition to an acceptance that vary with time. "
       << "We compute the coverage map using the UTC and JD distributions of the simulated events in order to "
       << "absorb the acceptance modulation. A Rayleigh analysis to the Right Ascension of the events and a "
       << "direct fit are then performed. The <lobe file> must contain two colums: theta (in deg.) and the lobe "
       << "value normalized to one at maximum. You can easily produce one of this file using the compute_lobe "
       << "executable." << endl;
  
  cout << endl;  
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
  if(strncmp(argv[1],"-",1) == 0) Usage(argv[0]);
  
  // Checking that we can find the lobe file.
  string lobeFile = argv[1];
  struct stat fileStat;
  if(stat(lobeFile.c_str(),&fileStat) == -1) 
    {
      cerr << "The file : " << lobeFile << " has not been found." << endl;
      exit(0);
    }

  int fargc = 1;
  string extension;
  TRint *rint = NULL;
  rint = new TRint("SimuAccDip", &fargc, argv);
  extension = ".png";
  gROOT->SetStyle("Plain");
  gStyle->SetTitleFont(30,"TITLE");

  
  // Read the lobe file
  ifstream ifsLobe(lobeFile.c_str());
  vector<double> thetaLobe, lobe;
  while(!ifsLobe.eof())
    {
      float thetaLobeTmp, lobeTmp;
      ifsLobe >> thetaLobeTmp >> lobeTmp;
      thetaLobe.push_back(thetaLobeTmp);
      lobe.push_back(lobeTmp);
    }

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //           Simulation of the events according to the dipolar            //
  //                  modulation and the acceptance effect                  //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  // Simulation of my dipole
  unsigned int nSide = 64;
  int sizeX = 800, sizeY = 500;
  double decLimit = 25.;
  double longStep = 60., latStep = 30.;

  THealpixMap dipoleMap(nSide);
  vector<long> iPix(dipoleMap.NPix());
  vector<double> lPix, bPix; // Deg
  for(unsigned int i = 0; i < dipoleMap.NPix(); i++) iPix[i] = i;
  dipoleMap.GiveLB(iPix, lPix, bPix);
 
  double raDipole, decDipole, lDipole, bDipole, ampDipole;
  ampDipole = 0.1;
  raDipole = 120.;
  decDipole = -30.;
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
  dipoleMapProj.SetPalette(PaletteOrange, 255); 
  dipoleMapProj.ShowGrid(longStep,latStep);
  dipoleMapProj.Save("fullSkyDipole"+extension);
  
  // Simulation of the acceptance effect
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
 
  // UTC law to be taken by the simulator
  unsigned int nbUTC = 1000;
  double xUTC, yUTC;
  string lawFileUTC = "utcLaw.txt";
  ofstream fileUTC(lawFileUTC.c_str());
  for(unsigned int i = 0; i < nbUTC; i++)
    {
      xUTC = 24.*i*1./(nbUTC-1);
      yUTC = 1+0.3*sin(xUTC/24*2*M_PI);
      fileUTC << xUTC << " " << yUTC << endl;
    }
  fileUTC.close();

  // JD law to be taken by the simulator
  double minJD, maxJD;
  date2jd(2004, 1, 1, 0., &minJD);
  date2jd(2004, 7, 1, 0., &maxJD);
  unsigned int nbJD = (unsigned int)((maxJD-minJD)*10);
  double xJD, yJD;
  string lawFileJD = "jdLaw.txt";
  ofstream fileJD(lawFileJD.c_str());
  for(unsigned int i = 0; i < nbJD; i++)
    {
      xJD = minJD+(maxJD-minJD)*i*1./(nbJD-1);
      yJD = (1.+0.4*sin(xJD/365*2*M_PI))*(1.+(xJD-minJD)/(maxJD-minJD)*5);
      fileJD << setprecision(12) << xJD << " " << yJD << endl;
    }
  fileJD.close();

  double latSite = kConstantsTK::AugerSouthLatitude;
  double lonSite = kConstantsTK::AugerSouthLongitude;
  unsigned int nb = 500000;
  
  vector<TEvent> simData;
  simData = SimulateEvents(dipoleMap, nb, thVal, pthVal, latSite, lonSite, lawFileUTC, lawFileJD);
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
  
  // First try : fitting with a simple geometric function
  coverage.fThetaDist.fAngFitFunc = geosplFunction;
  coverage.fThetaDist.SetNbParameters(4);
  // param 0 : dummy (about the Fermi-Dirac)
  // param 1 : dummy (about the Fermi-Dirac)
  // param 2 : thetamin
  // param 3 : thetamax
  coverage.fThetaDist.fDegreeMax = 3;
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

  // Take into account the zenith angle distribution to compute the 
  // coverage map
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
  for(long i = 0; i < nPix; i++) if(covMap[i] == 0.) {intCovMap[i] = 0.; intEvtMap[i] = 0.;}

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
  
  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                         Map of the difference                          //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////
  
  THealpixMap diffMap = intEvtMap-intCovMap;

  TProjMap diffMapProj(diffMap, sizeX, sizeY, decLimit);
  diffMapProj.SkyMap("Difference Map");
  diffMapProj.SetPalette(PaletteOrange, 255); 
  diffMapProj.ShowGrid(longStep,latStep);
  diffMapProj.ShowFOV(kConstantsTK::AugerSouthLatitude+60.); 
  diffMapProj.Save("diffMap"+extension);
 
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
  cout << "amplitude = " << amplitude << endl;
  cout << " phase = " << phase << endl;
  cout << " probability = " << probability << endl;
  cout << " significance = " << significance << endl;

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
