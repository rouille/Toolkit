#include <iostream>
#include <fstream>
#include <sys/stat.h>


#include "events.h"
#include "coverage.h"
#include "angdist.h"
#include "maptools.h"
#include "Cl.h"
#include "rayleigh.h"
#include "common.h"

// ROOT
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"

#ifdef gcc323
char* operator+(std::streampos&, char*);
#endif

using namespace std;


void Usage(string myName)
{
  cout << endl;
  cout << " Synopsis : " << endl;
  cout << myName << " <events file> <lobe file>" << endl << endl;
  
  cout << " Description :" << endl;
  cout << myName << " extracts data from <events file> and computes the related coverage map. The events and "
       << "coverage maps are then integrated in the lobe specified by <lobe file>. The Angular power spectrum "
       << "(Cl) is then computed debiasing from the partial sky. The <lobe file> must contain two colums: theta "
       << "(in deg.) and the lobe value normalized to one at maximum. You can easily produce one of this file "
       << "using the compute_lobe executable. The <events file> must contain this fields :" << endl;

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
  string eventFile = argv[1], lobeFile = argv[2];
  if( !CheckFile(eventFile) || !CheckFile(lobeFile) )
    {
      cerr << "File(s): " << lobeFile << " or/and " << eventFile << " not found" << endl;
      exit(0);
    }


  // ROOT
  int fargc = 1;
  string extension;
  TRint* rint = new TRint("AngularPowerSpectrum", &fargc, argv);
  extension = ".png";
  gROOT->SetStyle("Plain");
  gStyle->SetTitleFont(30,"TITLE");
  
  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                          Reading events file                           //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////
 
  cout << "Reading events file " << eventFile << endl;

  vector<TEvent> events = GetEvents(eventFile);
  long nEvents = events.size();
  if( nEvents == 0) {cout << "Program Failed : No events read. Exiting." << endl; exit(0);}

  // Min and Max dates
  DECLARE_VECTOR(double,utcsTime,events,fUTCs);
  int minYear, minMonth, minDay, maxYear, maxMonh, maxDay;
  double minUtcs, maxUtcs, minhourUtc, maxhourUtc;
  minUtcs = *min_element(utcsTime.begin(),utcsTime.end());
  maxUtcs = *max_element(utcsTime.begin(),utcsTime.end());
  utcs2date(minUtcs, &minYear, &minMonth, &minDay, &minhourUtc);
  utcs2date(maxUtcs, &maxYear, &maxMonh, &maxDay, &maxhourUtc);
  cout << "First event : " << minDay << "/" << minMonth << "/" << minYear << " UTC = " << minhourUtc << endl;
  cout << "Last event : " << maxDay << "/" << maxMonh << "/" << maxYear << " UTC = " << maxhourUtc << endl;
  
  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                                 Analysis                               //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  unsigned int nSide = 64;

#define cMapWithData

#ifdef cMapWithData
  double thetaMin = 0.;
  double thetaMax = 60.;
  unsigned int nBins = 20;
  double chi2Theta = 1.5;
  string accTimeModel = "TIME_FLAT";

  // Compute the Coverage Map
  TCoverage coverage(nSide);
  coverage.SetCoordSystem('G');
  coverage.SetLatitude(kConstantsTK::AugerSouthLatitude);
  coverage.SetLongitude(kConstantsTK::AugerSouthLongitude);

  // Time modulation
  coverage.fTimeMod.SetAccTimeModel(accTimeModel);

  // Zenith angle modulation
  coverage.fThetaDist.SetExtension(extension);
  coverage.fThetaDist.SetAngleName("theta");
  DECLARE_VECTOR(double,vTheta,events,fTheta);
  coverage.fThetaDist.SetData(vTheta);
  coverage.fThetaDist.SetOptions(thetaMin,thetaMax,nBins, chi2Theta);

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

  // Second try : fitting with a Fermi-Dirac + geometric function
  if( !fitok )
    {
      cout << "Add Fermi-Dirac function in order to correctly fit zenith "
	   << "angle distribution" << endl;
      coverage.fThetaDist.fAngFitFunc = fdsplFunction;
      coverage.fThetaDist.SetNbParameters(4);
      // param 0 : FD angle cutoff (degrees)
      // param 1 : FD width
      // param 2 : thetamin
      // param 3 : thetamax
      coverage.fThetaDist.fDegreeMax = 10;
      coverage.fThetaDist.fDataMin = thetaMin;
      coverage.fThetaDist.fDataMax = thetaMax;
      coverage.fThetaDist.fParameters[0].SetParameter(eFree,50.);
      coverage.fThetaDist.fParameters[1].SetParameter(eFree,6.);
      coverage.fThetaDist.fParameters[2].SetParameter(eFixed,thetaMin);
      coverage.fThetaDist.fParameters[3].SetParameter(eFixed,thetaMax);
      fitok = coverage.fThetaDist.Run();
    }

  if( !fitok )
    {
      cout << "Impossible to fit theta distribution. Try another fitting function. Exiting." << endl;
      exit(0);
    }

  // Take into account the zenith angle distribution to compute the 
  // coverage map
  coverage.CorrectForAngularModulation("theta");


  // Take into account time modulation (UTC and/or JD)
  if(accTimeModel != "TIME_FLAT")
    {
      DECLARE_VECTOR(double, UTCh, events, fUTCh);
      DECLARE_VECTOR(double, UTCs, events, fUTCs);
      double chi2Lim = 5;
      double timeStep = 7200.;
      coverage.fTimeMod.ComputeAccTime(UTCh, UTCs, timeStep, chi2Lim);
      coverage.CorrectForTimeModulation(accTimeModel);
    }
  
  // VERY IMPORTANT (initializes many constants to save time)
  // To be called once latsite, thetaMax and thetaMin are set
  coverage.ComputeDeclinationLimits();
  cout << "declination limits: " << coverage.GetDecMin() << " " << coverage.GetDecMax() << endl;
  
  // Compute the coverage map
  coverage.ComputeCoverage();
#endif

#ifdef cMapAnalytical
  TCoverage coverage(nSide);
  double thetaMax = 60.;
  coverage.fMap = GetAnalyticalCoverage(nSide, thetaMax, kConstantsTK::AugerSouthLatitude);
#endif


  // Convert the coverage map into an Healpix map
  THealpixMap covMap = coverage.GetMap();
  covMap *= (nEvents*1./covMap.Total());

  // Making events map
  DECLARE_VECTOR(double,ll,events,fL);
  DECLARE_VECTOR(double,bb,events,fB);
  THealpixMap evtMap = map_events(nSide, ll, bb);

  // Binning of the Power Spectrum in l.
  const unsigned int lmax = 20;
  vector<unsigned int> lbins;
  for(unsigned int i = 0; i < lmax+2; i++) lbins.push_back(i);
  vector<vector<double> > lvalue = lvalues(lbins);
  vector<double> ErrorCl;
  vector<double> Cl = compute_Cl(events, covMap, evtMap, lmax, lbins, ErrorCl);

  // Plot the power spectrum
  double* lgraph = new double[lvalue[0].size()];
  double* Clgraph = new double[lvalue[0].size()];
  double* lerrorgraph = new double[lvalue[0].size()];
  double* Clerrorgraph = new double[lvalue[0].size()];
  for(unsigned int i=0; i<lvalue[0].size(); i++)
    {
      lgraph[i] = lvalue[0][i];
      lerrorgraph[i] = lvalue[1][i];
      Clgraph[i] = Cl[i];
      Clerrorgraph[i] = ErrorCl[i];
      if(lgraph[i] == 0)
	{
	  lgraph[i] = 0.;
	  lerrorgraph[i] = 0.;
	  Clgraph[i] = 0.;
	  Clerrorgraph[i] = 0.;
	}
      cout.precision(4);
      cout << " l = " << lgraph[i] << " and Cl = " << Clgraph[i] << " +/- " << Clerrorgraph[i] << endl;
    }

  // Plot the Power Spectrum
  string name = "Power Spectrum"; string Xaxis = "l"; string Yaxis = "C_{l}"; string save = "Cl"+extension;
  TCanvas* cPS = new TCanvas("cPS", name.c_str(), 700, 700);
  TGraphErrors* PS = new TGraphErrors(lvalue[0].size(), lgraph, Clgraph, lerrorgraph, Clerrorgraph);
  PlotXY(cPS, PS, 0., lmax+1, name, Xaxis, Yaxis);
  PS->SetMarkerStyle(20); PS->SetMarkerSize(0.8);
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
