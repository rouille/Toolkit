#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sys/stat.h>

// Files
#include "healpixmap.h"
#include "events.h"
#include "coverage.h"
#include "angdist.h"
#include "maptools.h"
#include "projmap.h"
#include "simuevents.h"
#include "common.h"
#include "lima.h"

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
  cout << "Two sources are simulated on a full and uniform sky. Their coordinates are those of the Galactic "
       << "Center and Centaurus A. Events are then simulated following this anisotropic pattern in addition to "
       << "a constant acceptance. The coverage map assuming a flat acceptance is then computed. Finally, the "
       << "Li & Ma map and the histogram of the significance of the pixels are computed. We also show a zoom "
       << "on both candidates. The <lobe file> must contain two colums: theta (in deg.) and the lobe value "
       << "normalized to one at maximum. You can easily produce one of this file using the compute_lobe "
       << "executable." << endl;

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
  if( !CheckFile(lobeFile) ) {cerr << "File: " << lobeFile << " not found" << endl; exit(0);}

  int fargc = 1;
  string extension;
  TRint * rint = NULL;
  rint = new TRint("SimuSources", &fargc, argv);
  gROOT->SetStyle("Plain");
  gStyle->SetTitleFont(30,"TITLE");
  extension = ".png";

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                Simulation of the sources & the events                  //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  // Simulation of the sources
  unsigned int nSide = 128;
  int sizeX = 800, sizeY= 500;
  double decLimit = 25.;
  double longStep = 60., latStep = 30.;

  THealpixMap sourceMap(nSide), sourceMapTmp(nSide);

  // Pixels
  vector<long> iPix(sourceMap.NPix());
  vector<double> lPix, bPix; // Deg
  for(unsigned int i = 0; i < sourceMap.NPix(); i++) iPix[i] = i;
  sourceMap.GiveLB(iPix, lPix, bPix);
  iPix.clear();
 
  // Sources
  double lSource[2] = {0.,309.5}, bSource[2] = {0,19.5}, ampSource[2] = {0.6,0.9}, sigmaSource[2] = {4.,2.};
  
  for(unsigned int i = 0; i < 2; i++)
    {
      for(unsigned int j = 0; j < sourceMap.NPix(); j++)
        {
          vector<double> uvPix = ll2uv(lPix[j], bPix[j]);
          vector<double> uvSource = ll2uv(lSource[i],bSource[i]);
          double cosAngle = uvPix[0]*uvSource[0]+uvPix[1]*uvSource[1]+uvPix[2]*uvSource[2];
          sourceMapTmp[j] = exp( (cosAngle-1)/(2.*(1-cos(sigmaSource[i]*DTOR))) ); 
        }
      sourceMapTmp *= ampSource[i]/sourceMapTmp.Max();
      sourceMap = sourceMap + sourceMapTmp;
    }
  
  sourceMap = 1.+sourceMap;
  TProjMap mapSourcesProj(sourceMap, sizeX, sizeY, 90.);
  mapSourcesProj.SkyMap("Simulated Sources");
  mapSourcesProj.SetPalette(PaletteOrange, 255); 
  mapSourcesProj.ShowGrid(longStep,latStep);
  mapSourcesProj.Save("fullSkySources"+extension);
  
  // Simulation of the events : no acceptance law
  double thetaMax = 60.;
  unsigned int nVal = 5000;
  vector<double> thVal(nVal);
  vector<double> pthVal(nVal);
  for(unsigned int i=0; i<nVal; i++)
    {
      thVal[i] = i*180./(nVal-1);
      pthVal[i] = sin(thVal[i]*M_PI/180)*cos(thVal[i]*M_PI/180);
      if (thVal[i]>thetaMax) pthVal[i] = 0;
    }
  double latSite = kConstantsTK::AugerSouthLatitude;
  double lonSite = kConstantsTK::AugerSouthLongitude;
  unsigned int nb = 100000;

  vector<TEvent> simData = SimulateEvents(sourceMap, nb, thVal, pthVal, latSite, lonSite);
  cout << "Number of events simulated : " << simData.size() << endl;
  ShowLocalCoord(simData, thVal, pthVal); ShowEqCoord(simData), ShowArrivalTimesCoord(simData);

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                                 Analysis                               //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  TCoverage coverage(nSide);
  coverage.SetCoordSystem('G');
  coverage.fMap = GetAnalyticalCoverage(nSide, thetaMax, kConstantsTK::AugerSouthLatitude);

  // Li & Ma analysis
  double threshold = 5;
  TLiMa LiMa(coverage, simData, lobeFile, threshold);
  LiMa.ComputeLiMaMap();
  LiMa.ComputeMaxima();
  int nbinslima = 31;
  double minlima = -5;
  double maxlima = 5;
  TH1F* hlima = LiMa.GetLiMaHistogram(nbinslima,minlima,maxlima);
  LiMa.DrawLiMaHistogram(hlima);
  LiMa.PrintResults(false);

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                                  Maps                                  //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  // Coverage map
  THealpixMap covMap = coverage.fMap;
  covMap *= 1./covMap.Max();

  /* Coverage Map : Band of equal exposure
  covMap *= 1./covMap.Total();
  THealpixMap covMapIntExposure(nSide);
  double raPix[covMap.NPix()], decPix[covMap.NPix()];
  for(unsigned int i = 0; i < covMap.NPix(); i++)
    {
      gal2radec(lPix[i],bPix[i],&raPix[i],&decPix[i]);
      raPix[i] = raPix[i]*15;
    }
  int * index = new int[covMap.NPix()];
  TMath::Sort(covMap.NPix(),decPix,index,false);
  double exposure = 0;
  double band = 1.;
  for(unsigned int i = 0; i < covMap.NPix(); i++)
    {
      exposure += covMap[index[i]];
      if( exposure > band/6. && band < 6 ) band++; 
      covMapIntExposure[index[i]] = (7-band)/6;
    }
  cout << exposure << " " << band << endl;
  
  for(unsigned int i = 0; i < covMap.NPix(); i++) if( decPix[i] > 25. ) covMapIntExposure[i] = 0.;
  */

  TProjMap covMapProj(covMap, sizeX, sizeY, decLimit);
  covMapProj.SkyMap("Coverage Map");
  covMapProj.SetPalette(PaletteRGB, 255); 
  covMapProj.ShowGrid(longStep,latStep);
  covMapProj.ShowFOV(kConstantsTK::AugerSouthLatitude+60., 1); 
  covMapProj.Save("covMap"+extension);

  // Events map
  THealpixMap eventsMap = LiMa.GetEventsMap();
  TProjMap eventsMapProj(eventsMap, sizeX, sizeY, decLimit);
  eventsMapProj.SkyMap("Events Map");
  eventsMapProj.SetPalette(PaletteRGB, 255); 
  eventsMapProj.ShowGrid(longStep,latStep);
  eventsMapProj.ShowFOV(kConstantsTK::AugerSouthLatitude+60., 1); 
  eventsMapProj.Save("eventsMap"+extension);

  // LiMa map
  THealpixMap LiMaMap = LiMa.GetLiMaMap();
  vector<long> ipAboveThreshold;
  vector<double> valueAbovethreshold;
  LiMa.GetMaxima(ipAboveThreshold,valueAbovethreshold);

  TProjMap LiMaMapProj(LiMaMap, sizeX, sizeY, decLimit);
  LiMaMapProj.SkyMap("Blind Search Candidates",-10,10);
  LiMaMapProj.SetPalette(PaletteBlueAndRed, 255); 
  LiMaMapProj.ShowGrid(longStep,latStep);
  LiMaMapProj.ShowFOV(kConstantsTK::AugerSouthLatitude+60., 1); 
  LiMaMapProj.ShowSGP(1, kRed, 2);
  LiMaMapProj.PlotMaxima(ipAboveThreshold);
  LiMaMapProj.Save("LiMaMapMollweide"+extension);


  cout << "Program Finished Normally" << endl;
  rint->Run(kTRUE);
}
