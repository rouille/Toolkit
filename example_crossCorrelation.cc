#include <iostream>
#include <fstream>

#include "events.h"
#include "coverage.h"
#include "maptools.h"
#include "crosscorrelation.h"
#include "simuevents.h"
#include "common.h"
#include "projmap.h"

// ROOT
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
  cout << myName << " <events file> <catalog file>" << endl << endl;
  
  cout << " Description :" << endl;  
  cout << myName << " computes the cross correlation between the events in <events file> and astrophysical "
       << "objects in <catalog file>. The <catalog file> must have four columns : RA (right ascension), "
       << "DEC (declination), Z (redshift) and F (flux of the sources). The <events file> must contain the "
       << "following fields :" << endl;
  
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
  char * catalogueFile = argv[2];
  if( !CheckFile(eventFile) ) {cerr << "File: " << eventFile << " not found" << endl; exit(0);}
  if( !CheckFile(catalogueFile) ) {cerr << "File: " << catalogueFile << " not found" << endl; exit(0);}
 
  // ROOT
  int fargc = 1;
  string extension = ".png";
  TRint* rint = new TRint("Cross Correlation", &fargc, argv);
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

  // We want to remove events and galaxies found to be at less than 18° of Centaurus A
  vector<double> uvCenA = ll2uv(309.512,19.4179);
  double sepAng = 0.;

  // Events
  vector<TEvent> events;
  cout << "Reading events file " << eventFile << endl;
  events = GetEvents(eventFile);
  if( !events.size() ) {cout << "Program Failed: No events read. Exiting." << endl; exit(0);}
  vector<TEvent> eventsWithoutCenA;
  for(unsigned int i = 0; i < events.size(); i++) 
    {
      vector<double> uvEvents = ll2uv(events[i].fL,events[i].fB);
      double ang = acos(uvEvents[0]*uvCenA[0]+uvEvents[1]*uvCenA[1]+uvEvents[2]*uvCenA[2]);
      if( ang > sepAng*DTOR ) eventsWithoutCenA.push_back(events[i]);
    }
  cout << eventsWithoutCenA.size() << " events excluding CenA" << endl;

  // Catalogue
  cout << "Reading catalog file " << catalogueFile << endl;
  ifstream catalogue(catalogueFile);
  double raSourceTmp, decSourceTmp, lSourceTmp, bSourceTmp, redshiftSourceTmp, fluxSourceTmp;
  vector<double> lSources, bSources, fluxSources;
  vector<string> nameSources;
  catalogue.ignore(1000,'\n');
  while( !catalogue.eof() )
    {
      catalogue >> raSourceTmp >> decSourceTmp >> redshiftSourceTmp >> fluxSourceTmp;
      radec2gal(raSourceTmp/15.,decSourceTmp, &lSourceTmp, &bSourceTmp);
      vector<double> uvCatalogue = ll2uv(lSourceTmp,bSourceTmp);
      double ang = acos(uvCatalogue[0]*uvCenA[0]+uvCatalogue[1]*uvCenA[1]+uvCatalogue[2]*uvCenA[2]);
      if( ang > sepAng*DTOR )
        {
          lSources.push_back(lSourceTmp);
          bSources.push_back(bSourceTmp);
          fluxSources.push_back(fluxSourceTmp);
          nameSources.push_back("");
        }
    }
  catalogue.close();
  cout << fluxSources.size() << " astrophysical objects found in " << catalogueFile << endl;

  // Map
  double decLimit = 90.;
  int sizeX = 800, sizeY = 500;
  double longStep = 60., latStep = 30.;
  char mapTitle[500]; sprintf(mapTitle,"SWIFT and events above 55 EeV");

  THealpixMap mapEvents2MASS(nSide);
  
  TProjMap mapProjEvents2MASS(mapEvents2MASS, sizeX, sizeY, decLimit);
  mapProjEvents2MASS.SkyMap(mapTitle);
  mapProjEvents2MASS.SetPalette(PaletteOrange, 255); 
  mapProjEvents2MASS.ShowGrid(longStep,latStep);
  mapProjEvents2MASS.ShowSGP(1, kBlue, 2);
  mapProjEvents2MASS.ShowFOV(25.);
  mapProjEvents2MASS.PlotSources(lSources, bSources,nameSources,false,kCircle);
  mapProjEvents2MASS.PlotEvents(eventsWithoutCenA,kFullCircle,kRed);
  mapProjEvents2MASS.Save("crossCorrelationMap"+extension);

  // Coverage map
  THealpixMap covMap = GetAnalyticalCoverage(nSide, thetaMax, kConstantsTK::AugerSouthLatitude);
  covMap *= 1./covMap.Max();
  char mapName[500]; sprintf(mapName,"CoverageMap.fits");

  TProjMap covMapProj(covMap, sizeX, sizeY, decLimit);
  covMapProj.SkyMap("Coverage Map");
  covMapProj.SetPalette(PaletteOrange, 255); 
  covMapProj.ShowGrid(longStep,latStep);
  covMapProj.ShowFOV(kConstantsTK::AugerSouthLatitude+60.);
  covMapProj.Save("covMap"+extension);
  
  for(unsigned int i = 0; i < fluxSources.size(); i++) fluxSources[i] = 1.;//covMap.Value(lSources[i],bSources[i]);

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                                 Analysis                               //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  double dAlpha = 1.;
  double alphaCut = 30.;
  double dispersion = 0.9;

  TCrossCorrelation crossCorrelation(eventsWithoutCenA, lSources, bSources, fluxSources, dAlpha, alphaCut);
  cout << "Events Correlation Function Computation" << endl;
  crossCorrelation.ComputeEventsCrossCorrelation();
  
  cout << "Coverage Correlation Function Computation using simulation" << endl;
  unsigned int nSimu = 1000;
  THealpixMap map(nSide);
  map = 1;
  crossCorrelation.ComputeCoverageCrossCorrelation(lonSite, latSite, map, nSimu, dispersion, thVal, pthVal);

  crossCorrelation.DrawEvents();
  crossCorrelation.DrawCov();
  crossCorrelation.DrawBoth();
  crossCorrelation.DrawBothRelativeExcess();
 
  cout << "Program Finished Normally" << endl;
  rint->Run(kTRUE);
}
