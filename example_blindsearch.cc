#include <iostream>
#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstring>
#include <iomanip>

#include "STClibrary.h"
#include "common.h"
#include "events.h"
#include "coverage.h"
#include "angdist.h"
#include "maptools.h"
#include "projmap.h"
#include "lima.h"

// ROOT
#include "TRint.h"
#include "TROOT.h"
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
	<< "coverage maps are then integrated in the lobe specified by <lobe_file> and a Li & Ma significance "
	<< "map is computed. The <lobe file> must contain two colums: theta (in deg.) and the lobe value "
	<< "normalized to one at maximum. You can easily produce one of this file using the compute_lobe "
	<< "executable. The <events file> must contain the following fields :" << endl;
	
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
	TRint* rint = new TRint("BlindSearch", &fargc, argv);
	extension = ".png";
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFont(30,"TITLE");
	
	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                          Reading events file                           //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////
	
	cout << "Reading events file " <<  eventFile << endl;
	vector<TEvent> events = GetEvents(eventFile);
	long nEvents = events.size();
	if( nEvents == 0 ) {cout << "Program Failed: No events read. Exiting." << endl; exit(0);}


	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                                 Analysis                               //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////
	
	unsigned int nSide = 128;
	double threshold = 3.5;
	double lCenter = kSTC::GalacticSouthPoleL;
	double bCenter = kSTC::GalacticSouthPoleB;
	
#define cMapWithData
	
#ifdef cMapWithData
	double thetaMin = 0.;
	double thetaMax = 60.;
	unsigned int nBins = 60;
	double chi2Theta = 1.5;
	string accTimeModel = "TIME_FLAT";
	TCoverage coverage(nSide);
	coverage.SetExtension(extension);
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
	coverage.fThetaDist.SetOptions(thetaMin,thetaMax,nBins,chi2Theta);
	
	// First try : fitting with a simple geometric function
	coverage.fThetaDist.fAngFitFunc = geopolyFunction;
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
		  cout << "Add Fermi-Dirac function in order to correctly fit zenith angle distribution" << endl;
		  coverage.fThetaDist.fAngFitFunc = fdpolyFunction;
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
	// Take into account the zenith angle distribution to compute the coverage map
	coverage.CorrectForAngularModulation("theta");
	
	// Take into account time modulation (UTC and/or JD)
	if(accTimeModel != "TIME_FLAT")
    {
		  DECLARE_VECTOR(double,UTCh,events,fUTCh);
		  DECLARE_VECTOR(double,UTCs,events,fUTCs);
		  double chi2Lim = 5;
		  double timeStep = 3600.;
		  coverage.fTimeMod.ComputeAccTime(UTCh,UTCs,timeStep,chi2Lim);
		  coverage.CorrectForTimeModulation(accTimeModel);
    }
	
	// Initializes many constants to save time. To be called once latitude, thetaMax and thetaMin are set
	coverage.ComputeDeclinationLimits();
	cout << "declination limits: " << coverage.GetDecMin() << " " << coverage.GetDecMax() << endl;
	
	// COMPUTE THE COVERAGE MAP
	coverage.ComputeCoverage();
#endif
	
#ifdef cMapAnalytical
	TCoverage coverage(nSide);
	double thetaMax = 60.;
	coverage.fMap = GetAnalyticalCoverage(nSide, thetaMax, kConstantsTK::AugerSouthLatitude);
#endif
	
	// Blind search
	TLiMa LiMa(coverage, events, lobeFile, threshold);
	LiMa.SetExtension(extension);
	LiMa.ComputeLiMaMap();
	LiMa.ComputeMaxima();
	int nbinslima = 30;
	double minlima = -6;
	double maxlima = 6;
	TH1F* hlima = LiMa.GetLiMaHistogram(nbinslima,minlima,maxlima);
	LiMa.DrawLiMaHistogram(hlima);
	
	bool detailed = true;  
	LiMa.PrintResults(detailed);
	
	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                                    Plot                                //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////
	
	double decLimit = 25.;
	int sizeX = 800, sizeY = 400;
	double longStep = 60., latStep = 30.;
	
	// Coverage map
	THealpixMap covMap = LiMa.GetCovMap();
	TProjMap covMapProj(covMap, sizeX, sizeY, decLimit);
	covMapProj.SkyMap("Coverage Map");
	covMapProj.SetPalette(PaletteRGB, 255); 
	covMapProj.ShowGrid(longStep,latStep);
	covMapProj.Save("covMap"+extension);
	
	// Events map
  THealpixMap eventsMap = LiMa.GetEventsMap();
	TProjMap eventsMapProj(eventsMap, sizeX, sizeY, decLimit);
	eventsMapProj.SkyMap("Events Map"); 
	eventsMapProj.SetPalette(PaletteRGB, 255); 
	eventsMapProj.ShowGrid(longStep,latStep);
	eventsMapProj.Save("eventsMap"+extension);

	// LiMa map
	THealpixMap LiMaMap = LiMa.GetLiMaMap();
	vector<long> ipAboveThreshold;
	vector<double> valueAbovethreshold;
	LiMa.GetMaxima(ipAboveThreshold,valueAbovethreshold);
	
	TProjMap LiMaMapProjMollweide(LiMaMap, sizeX, sizeY, decLimit);
	LiMaMapProjMollweide.SkyMap("Blind Search Candidates"); 
	LiMaMapProjMollweide.SetPalette(PaletteBlueAndRed, 255); 
	LiMaMapProjMollweide.ShowGrid(longStep,latStep);
	LiMaMapProjMollweide.ShowSGP(1, kBlack, 2);
	
	ifstream catalogue("Catalogue/hess.dat");
	double lCat, bCat;
	string nameCat;
	vector<double> lSources, bSources;
	vector<string> nameSources;
	while(catalogue >> lCat)
    {
		  catalogue >> bCat >> nameCat;
		  lSources.push_back(lCat);
		  bSources.push_back(bCat);
		  nameSources.push_back((string)nameCat);
    }
	catalogue.close();      
	
	LiMaMapProjMollweide.PlotSources(lSources, bSources, nameSources);
	LiMaMapProjMollweide.Save("LiMaMapMollweide"+extension);

	TProjMap LiMaMapProjLambert(LiMaMap, 600, 600, lCenter, bCenter, 100., decLimit);
	LiMaMapProjLambert.SkyMap("Blind Search Candidates");
	LiMaMapProjLambert.SetPalette(PaletteBlueAndRed, 255);
	LiMaMapProjLambert.ShowGrid(30.,30.);
	LiMaMapProjLambert.ShowFOV(kConstantsTK::AugerSouthLatitude+60.); 
	LiMaMapProjLambert.ShowSGP(1, kBlack, 2);
	LiMaMapProjLambert.PlotSources(lSources, bSources, nameSources);
	LiMaMapProjLambert.Save("LiMaMapLambert"+extension);

#ifdef File
	// File for Paolo
	THealpixMap covMapOriginal = coverage.GetMap();
	char covMapFile[64]; sprintf(covMapFile,"covMap.fits");
	covMapOriginal.WriteFits(covMapFile);
	covMapOriginal *= events.size()/covMapOriginal.Total();
	DECLARE_VECTOR(double,lEvents,events,fRa);
	DECLARE_VECTOR(double,bEvents,events,fB);
	THealpixMap evtMapOriginal = map_events(nSide, lEvents, bEvents);

	ofstream Paolo("HotSpot.txt");
	double lPix, bPix, raPix, decPix;
	for(unsigned int i = 0; i < covMap.NPix(); i++)
	  {
		  covMap.GiveLB(i,lPix,bPix);
		  vector<double> uvPix = ll2uv(lPix,bPix);
		  gal2radec(lPix,bPix,&raPix,&decPix);
		  Paolo << i << " " << LiMaMap[i] << " " << uvPix[0] << " " << uvPix[1] << " " << uvPix[2] << " " 
		        << lPix << " " << bPix << " " << raPix << " " << decPix << " " << evtMapOriginal[i] << " " 
		        << eventsMap[i] << " " << covMapOriginal[i] << " " << covMap[i] << endl; 
    }
	Paolo.close();
#endif
	
	cout << "Program Finished Normally" << endl;
	rint->Run(kTRUE);
}
