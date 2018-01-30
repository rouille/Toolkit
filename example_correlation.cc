#include <iostream>

#include "events.h"
#include "coverage.h"
#include "angdist.h"
#include "maptools.h"
#include "correlation.h"
#include "simuevents.h"
#include "common.h"

// ROOT
#include "TCanvas.h"
#include "TRint.h"
#include "TROOT.h"
#include "TStyle.h"


#ifdef gcc323
char* operator+( std::streampos&, char* );
#endif

using namespace std;



void Usage(string myName)
{
  cout << endl;
  cout << " Synopsis : " << endl;
  cout << myName << " <events file> <simu flag>" << endl << endl;
  
  cout << " Description :" << endl;  
  cout << myName << " computes the 2-points correlation function. The <simu flag> must be set to 0 if one wants "
       << "to compute the 2-points correlation function of the events stored in <events file> and set to 1 if "
       << "one wants to use simulated events. In the latter case <events file> will be ignored. The "
       << "<events file> must contain the following fields :" << endl;
  
  DumpFields();
  
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
  if(argc != 3) Usage(argv[0]);
  string eventFile = argv[1];
  int isSimu = atoi(argv[2]);
  if( !CheckFile(eventFile) ) {cerr << "File: " << eventFile << " not found" << endl; exit(0);}

  // ROOT
  int fargc = 1;
  string extension;
  TRint* rint = new TRint("2pt Correlation", &fargc, argv);
  extension = ".png";
  gROOT->SetStyle("Plain");
  gStyle->SetTitleFont(30,"TITLE");

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                      Events (file or simulation)                       //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////
 
  unsigned int nSide = 64;
  double latSite = kConstantsTK::AugerSouthLatitude;
  double lonSite = kConstantsTK::AugerSouthLongitude;
  double thetaMax = 60.;

  unsigned int nVal = 5000;
  vector<double> thVal(nVal), pthVal(nVal);

  for(unsigned int i = 0; i < nVal; i++)
    {
      thVal[i] = i*180./(nVal-1);
      pthVal[i] = sin(thVal[i]*DTOR)*cos(thVal[i]*DTOR);
      if (thVal[i] > thetaMax) pthVal[i] = 0;
    }

  vector<TEvent> events;
  
  if(isSimu == 1)
    {
      THealpixMap mapSimu(nSide);
      mapSimu = 1; // Isotropic case
      unsigned int nb = 69;
      events = SimulateEvents(mapSimu, nb, thVal, pthVal, latSite, lonSite);
      cout << events.size() << " events simulated" << endl;
    }
  else
    {
      cout << "Reading events file " << eventFile << endl;
      events = GetEvents(eventFile);
      if( !events.size() ) {cout << "Program Failed: No events read. Exiting." << endl; exit(0);}
    }

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                                 Analysis                               //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  /* Three way to compute the coverage correlation function :
     + using coverage map (events or analytic) and then the relation existing 
     between C(l) and C(theta)
     + using simulation to generate events
     + one can just scrambled the arrival times of the events to get new 
     samples of events 
  */

#define simulation
  /*
    #define cMapAnalytic
    #define scrambling
    #define simulation
  */
  unsigned int dAlpha = 1;
  double alphaCut = 60.;
  double dispersion = 0.9;

#ifdef simulation
  // Using Simulation
  TCorrelation correlation(events, dAlpha, alphaCut);
  correlation.SetExtension(extension);
  cout << "Events Correlation Function Computation" << endl;
  correlation.ComputeEventsCorrelation();
  
  cout << "Coverage Correlation Function Computation using simulation" << endl;
  unsigned int nSimu = 250;
  THealpixMap map(nSide);
  map = 1;
  correlation.ComputeCoverageCorrelationSimu(lonSite, latSite, map, nSimu, dispersion, thVal, pthVal);
  vector<double> alpha, Np, errLow, errHigh;
  correlation.GetCoverageCorrelation(alpha,Np,errLow,errHigh);
  for(unsigned int i = 0; i < Np.size(); i++)
    cout << "Angle: " << alpha[i] << " Expected Np: " << Np[i] << " Error: " << errLow[i] << "/" << errHigh[i] << endl;
#endif


#ifdef scrambling  
  // Using Monte Carlo
  TCorrelation correlation(events, dAlpha, alphaCut);
  correlation.SetExtension(extension);
  cout << "Events Correlation Function Computation" << endl;
  correlation.ComputeEventsCorrelation();
  
  cout << "Coverage Correlation Function Computation using scrambled data Set" << endl;
  unsigned int nMC = 10;
  string binningType = "EVENTS";
  string scramblingType = "UTC+JD";
  correlation.ComputeCoverageCorrelationMC(lonSite, latSite, nMC, dispersion, binningType, scramblingType, thetaMax);
  vector<double> alpha, Np, errLow, errHigh;
  correlation.GetCoverageCorrelation(alpha,Np,errLow,errHigh);
  for(unsigned int i = 0; i < Np.size(); i++)
    cout << "Angle: " << alpha[i] << " Expected Np: " << Np[i] << " Error: " << errLow[i] << "/" << errHigh[i] << endl;
#endif


#ifdef cMapAnalytic
  // Using Coverage Map (Analytic)
  TCoverage coverage(nSide);
  coverage.SetCoordSystem('G');
  coverage.fMap = GetAnalyticalCoverage(nSide, thetaMax, kConstantsTK::AugerSouthLatitude);
  THealpixMap covMap = coverage.GetMap();
  covMap *= 1./covMap.Max(); // Normalization: maximum equal to 1
 
  // Events and coverage correlation function 
  TCorrelation correlation(events, dAlpha, alphaCut);
  correlation.SetExtension(extension);
  cout << "Events Correlation Function Computation" << endl;
  correlation.ComputeEventsCorrelation();
  
  cout << "Coverage Correlation Function Computation using the coverage map" << endl;
  correlation.ComputeCoverageCorrelationMap(covMap);
  vector<double> alpha, Np, errLow, errHigh;
  correlation.GetCoverageCorrelation(alpha,Np,errLow,errHigh); 
  for(unsigned int i = 0; i < Np.size(); i++)
    cout << "Angle: " << alpha[i] << " Expected Np: " << Np[i] << " Error: " << errLow[i] << "/" << errHigh[i] << endl;
#endif
  

  // 2 pts correlation function
  correlation.DrawEvents();
  correlation.DrawCov();
  correlation.DrawBoth();

  cout << "Program Finished Normally" << endl;
  rint->Run(kTRUE);
}
