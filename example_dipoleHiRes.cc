#include <iostream>
#include <cmath>

#include "dipoleHiRes.h"
#include "simuevents.h"
#include "common.h"
#include "maptools.h"
#include "projmap.h"

#include "TROOT.h"
#include "TRint.h"
#include "TCanvas.h"
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
  cout << myName << " is based upon comparison between data and a large quantity of events generated "
       << "isotropically by Monte Carlo simulation to determine the presence of a dipole. A binning technisque "
       << "that considers the event counts for the full range of opening angles from the center of each proposed "
       << "dipole distribution is used. The <simu flag> must be set to 0 if one wants to use the events stored "
       << "in <events file> and set to 1 if one wants to use simulated events according to a dipole pattern. In "
       << "the latter case <events file> will be ignored. The <events file> must contain the following fields :" 
       << endl;
  
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
  TRint* rint = new TRint("Fit Dipole HiRes", &fargc, argv);
  extension = ".png";
  gROOT->SetStyle("Plain");
  gStyle->SetTitleFont(30,"TITLE");

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //            dipole (simulation) & events (file or simulation)           //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////
  
  unsigned int nSide = 128;
  
  double latSite = kConstantsTK::AugerSouthLatitude;
  double lonSite = kConstantsTK::AugerSouthLongitude;

  double thetaMax = 60.;
  unsigned int nVal = 10000;
  vector<double> thVal(nVal);
  vector<double> pthVal(nVal);
  for(unsigned int i=0; i<nVal; i++)
    {
      thVal[i] = i*180./(nVal-1);
      pthVal[i] = sin(thVal[i]*M_PI/180)*cos(thVal[i]*M_PI/180);
      if (thVal[i] > thetaMax) pthVal[i] = 0;
    }

  double lDipole, bDipole, ampDipole;
  ampDipole = 0.05;
  lDipole = 0.;
  bDipole = 45.;
  vector<TEvent> events;

  if(isSimu == 1)
    {
      // Simulation of my dipole
      unsigned int nSide = 128;
      THealpixMap mapSimu(nSide);
      
      vector<long> iPix(mapSimu.NPix());
      vector<double> lPix; // Deg
      vector<double> bPix; // Deg
      
      for(unsigned int i = 0; i < mapSimu.NPix(); i++) iPix[i] = i;
      mapSimu.GiveLB(iPix,lPix,bPix);
      
      vector<double> uvDipole = ll2uv(lDipole, bDipole);
      vector<double> dipoleMap(mapSimu.NPix());
      
      for(unsigned int i=0; i<mapSimu.NPix(); i++)
        {
          vector<double> uvPix = ll2uv(lPix[i], bPix[i]);
          dipoleMap[i] = ampDipole*(uvPix[0]*uvDipole[0]+uvPix[1]*uvDipole[1]+uvPix[2]*uvDipole[2]);
        }
      
      mapSimu = 1.+dipoleMap;
      unsigned int nb = 50000;
      events = SimulateEvents(mapSimu, nb, thVal, pthVal, latSite, lonSite);
      cout << "Number of simulated events: " << events.size() << endl;
    }
  else
    {
      cout << "Reading events file " << eventFile << endl;
      events = GetEvents(eventFile);
      long nEvents = events.size();
      if(nEvents == 0) {cout << "Program Failed: No events read. Exiting." << endl; exit(0);}
    }


  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                                 Analysis                               //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  double dCosTheta = 0.05;
  TDipoleHiRes dipole(events, lDipole, bDipole, dCosTheta);

  // Events Opening Angle
  dipole.ComputeEventsOpeningAngle();

  // Background Opening Angle
  unsigned int nSimu = 100;
  THealpixMap map(nSide);
  map = map+1;
  dipole.ComputeCoverageOpeningAngle(lonSite, latSite, map, nSimu, thVal, pthVal);
  
  dipole.ComputeTrueOpeningAngle();
  dipole.fitDipole();

  dipole.DrawEvents();
  dipole.DrawCov();
  dipole.DrawBoth();
  dipole.DrawTrueFit();

  vector<double> fitParameters = dipole.GetFitParameters();
  double yIntercept = fitParameters[0];
  double slope = fitParameters[1];
  vector<double> fitParametersErrors = dipole.GetFitParametersErrors();
  double yInterceptError = fitParametersErrors[0];
  double slopeError = fitParametersErrors[1];
  double amplitude, varAmplitude, dAmplitude;
  amplitude = slope/yIntercept;
  varAmplitude = (slope*yInterceptError/pow(yIntercept,2))*(slope*yInterceptError/pow(yIntercept,2))+(slopeError/yIntercept)*(slopeError/yIntercept);
  dAmplitude = sqrt(varAmplitude);

  cout << "amplitude = " << amplitude << " +/- " << dAmplitude << endl;

  cout << "Program Finished Normally" << endl;
  rint->Run(kTRUE);
}
