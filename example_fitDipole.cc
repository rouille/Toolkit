#include <iostream>
#include <fstream>
#include <vector>
#include <sys/stat.h>

// Files
#include "events.h"
#include "coverage.h"
#include "angdist.h"
#include "maptools.h"
#include "simuevents.h"
#include "common.h"
#include "fitdipole.h"
#include "rayleigh.h"


// ROOT
#include "TRint.h"
#include "TROOT.h"
#include "TH1F.h"
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
  cout << myname << " <raDipole> <decDipole> <ampDipole> <nEvt> <nMC>" << endl << endl;

  cout << " Description :" << endl;
  cout << "<nEvt> events are simulated following a dipolar modulation with amplitude <ampDipole> and orientation "
       << "(<raDipole>,<decDipole>) in addition to a constant acceptance. The coverage map assuming a flat "
       << "acceptance is computed. We directly fit the product of the coverage map with a general dipole (three "
       << "degrees of freedom). The operation is repeated <nMC> times." << endl;
  
  exit(0);
}


int main(int argc, char* argv[])
{
  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                        To start (initialization)                       //
  //                                                                        //
  //////////////////////////////////////////////////////////////////////////// 
    
  // ROOT in batch
  gROOT->SetBatch(true);

  // Style
  gROOT->SetStyle("Plain");
  gStyle->SetTitleFont(30,"TITLE");

  // just to have color table
  TCanvas* cTmp = new TCanvas("cTmp");
  delete cTmp;
  string extension = ".png";
  
  // Command line
  if(argc != 6) Usage(argv[0]);
  if(strncmp(argv[1],"-",1) == 0) Usage(argv[0]);

  double raDipole = atof(argv[1]);
  double decDipole = atof(argv[2]);
  double ampDipole = atof(argv[3]);
  unsigned int nb = atoi(argv[4]);
  unsigned int nMC = atoi(argv[5]);

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                Simulation of the dipole & the events                   //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  unsigned int nSide = 8;

  // Simulation of my dipole
  unsigned int nSideSimu = 64;
  THealpixMap map(nSideSimu);
  
  vector<long> iPix(map.NPix());
  vector<double> lPix; // Deg
  vector<double> bPix; // Deg
 
  for(unsigned int i=0; i<map.NPix(); i++) iPix[i] = i;
  map.GiveLB(iPix,lPix,bPix);
  
  double lDipole, bDipole;
  radec2gal(raDipole/15., decDipole, &lDipole, &bDipole);
  vector<double> uvDipole = ll2uv(lDipole, bDipole);
  
  for(unsigned int i=0; i<map.NPix(); i++)
    {
      vector<double> uvPix = ll2uv(lPix[i], bPix[i]);
      map[i] = 1.+ampDipole*(uvPix[0]*uvDipole[0]+uvPix[1]*uvDipole[1]+uvPix[2]*uvDipole[2]);
    }
  
  // Simulation of the events : no acceptance law
  double thetaMin = 0.;
  double thetaMax = 60.;
  unsigned int nVal = 10000;
  vector<double> thVal(nVal);
  vector<double> pthVal(nVal);
  for(unsigned int i=0; i<nVal; i++)
    {
      thVal[i] = i*180./(nVal-1);
      pthVal[i] = sin(thVal[i]*DTOR)*cos(thVal[i]*DTOR);
      if (thVal[i] > thetaMax) pthVal[i] = 0;
    }
  double latSite = kConstantsTK::AugerSouthLatitude;
  double lonSite = kConstantsTK::AugerSouthLongitude;

  double* raFit = new double[nMC];
  double* decFit = new double[nMC];
  double* ampFit = new double[nMC];
  double* raFitError = new double[nMC];
  double* decFitError = new double[nMC];
  double* ampFitError = new double[nMC];
  char fileName[1024];
  sprintf(fileName, "dipole_RA%i_dec%i_amp%5.3f_N%i.txt", (int)raDipole, (int)decDipole, ampDipole, (int)nb);
  ofstream fitDipoleFile(fileName);
  for(unsigned int i = 0; i < nMC; i++)
    {
      vector<TEvent> simData = SimulateEvents(map, nb, thVal, pthVal, latSite, lonSite);

      //////////////////////////////////////////////////////////////////
      //                                                              //
      //                   Compute the  coverage map                  //
      //                                                              //
      //////////////////////////////////////////////////////////////////

      // Compute the Coverage Map
      unsigned int nBinsTheta = 30;
      double chi2Theta = 1.5;
      
      TCoverage coverage(nSide);
      coverage.SetCoordSystem('G');
      coverage.SetLatitude(latSite);
      coverage.SetLongitude(lonSite);

      // Time modulation
      coverage.fThetaDist.SetExtension(extension);
      coverage.fThetaDist.SetAngleName("theta");
      DECLARE_VECTOR(double, thetaData, simData, fTheta);
      coverage.fThetaDist.SetData(thetaData);
      coverage.fThetaDist.SetOptions(thetaMin,thetaMax,nBinsTheta,chi2Theta);

      // First try : fitting with a simple geometric function
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
      
      // To be called once latitude, thetaMax and thetaMin are set
      coverage.ComputeDeclinationLimits();
      cout << "declination limits : " << coverage.GetDecMin() << " " << coverage.GetDecMax() << endl;
      
      // Compute the coverage map
      coverage.ComputeCoverage();
      
      // Put the coverage map into a Healpix map
      THealpixMap covMap = coverage.GetMap();
      // Normalize it to the total number of events
      covMap *= thetaData.size()*1./covMap.Total();
       
      //////////////////////////////////////////////////////////////////
      //                                                              //
      //                   Compute the  events map                    //
      //                                                              //
      //////////////////////////////////////////////////////////////////
      
      DECLARE_VECTOR(double, lData, simData, fL);
      DECLARE_VECTOR(double, bData, simData, fB);
      
      THealpixMap evtMap = map_events(nSide, lData, bData);
      
      //////////////////////////////////////////////////////////////////
      //                                                              //
      //                  Fit my dipole on the sky                    //
      //                                                              //
      //////////////////////////////////////////////////////////////////

      TRayleigh* rayleigh = new TRayleigh(simData, covMap);
      rayleigh->ComputeAmplitude();
      rayleigh->ComputePhase();
      double amplitude = rayleigh->GetAmplitude();
      double phase = rayleigh->GetPhase();

      TFitDipole* fitDipole = new TFitDipole(evtMap, covMap, amplitude, phase);
      fitDipole->FitProcedure();
      vector<double> parameter = fitDipole->GetParameters();
      vector<double> errorparameter = fitDipole->GetErrorParameters();

      ampFit[i] = parameter[0];
      raFit[i] = parameter[1];
      decFit[i] = parameter[2];
      ampFitError[i] = errorparameter[0];
      raFitError[i] = errorparameter[1];
      decFitError[i] = errorparameter[2];

      fitDipoleFile << raFit[i] << " " <<raFitError[i] << " "
		                << decFit[i] << " " << decFitError[i] << " "
		                << ampFit[i] << " " << ampFitError[i] << endl;
      delete fitDipole;
    }
  fitDipoleFile.close();  

  for(unsigned int i=0; i<nMC; i++)
    cout << "MC : " << i << "   " << "RA : " << raFit[i] << " +/- " << raFitError[i] << " "
	       << "dec : " << decFit[i] << " +/- " << decFitError[i] << " "
         << "A : " << ampFit[i] << " +/- " << ampFitError[i] << endl;
  
  double raMean = 0, raMean2 = 0, raRMS = 0, decMean = 0, decMean2 = 0, decRMS = 0, ampMean = 0, ampMean2 = 0, ampRMS = 0;
  for(unsigned int i = 0; i < nMC; i++)
    {
      raMean += raFit[i]/nMC;
      raMean2 += raFit[i]*raFit[i]/nMC;
      decMean += decFit[i]/nMC;
      decMean2 += decFit[i]*decFit[i]/nMC;
      ampMean += ampFit[i]/nMC;
      ampMean2 += ampFit[i]*ampFit[i]/nMC;
    }
  raRMS = sqrt(-raMean*raMean+raMean2);
  decRMS = sqrt(-decMean*decMean+decMean2);
  ampRMS = sqrt(-ampMean*ampMean+ampMean2);
  
  cout << "Mean and RMS " << endl;
  cout << "------------" << endl;
  cout << "RA : " << raMean << " +/- " << raRMS << endl;
  cout << "dec : " << decMean << " +/- " << decRMS << endl;
  cout << "A : " << ampMean << " +/- " << ampRMS << endl;

  unsigned int nBins = 60;
  double minHisto = -3;
  double maxHisto = 3;
  string Xaxis = "Significance";
  string Yaxis = "";

  TCanvas* cRA = new TCanvas("cRA", "RA", 700, 700);
  TH1F* raHisto = new TH1F("raHisto", "RA", nBins, minHisto, maxHisto);
  char raSave[1024];
  sprintf(raSave, "hRA_RA%i_dec%i_amp%4.2f_N%i.png", (int)raDipole, (int)decDipole, ampDipole, (int)nb);
  for(unsigned int i = 0; i < nMC; i++) raHisto->Fill((raFit[i]-raMean)/raRMS);
  DrawHisto(cRA, raHisto, Xaxis, Yaxis, raSave);
  delete [] raFit;
  delete [] raFitError;
  delete cRA;
  delete raHisto;

  TCanvas* cDec = new TCanvas("cDec", "Dec", 700, 700);
  TH1F* decHisto = new TH1F("decHisto", "Dec", nBins, minHisto, maxHisto);
  char decSave[1024];
  sprintf(decSave, "hDec_RA%i_dec%i_amp%4.2f_N%i.png", (int)raDipole, (int)decDipole, ampDipole, (int)nb);
  for(unsigned int i = 0; i < nMC; i++) decHisto->Fill((decFit[i]-decMean)/decRMS);
  DrawHisto(cDec, decHisto, Xaxis, Yaxis, decSave);
  delete [] decFit;
  delete [] decFitError;
  delete cDec;
  delete decHisto;

  TCanvas* cAmplitude = new TCanvas("cAmplitude", "Amplitude", 700, 700);
  TH1F* aHisto = new TH1F("aHisto", "Amplitude", nBins, minHisto, maxHisto);
  char ampSave[1024];
  sprintf(ampSave, "hAmp_RA%i_dec%i_amp%4.2f_N%i.png", (int)raDipole, (int)decDipole, ampDipole, (int)nb);
  for(unsigned int i = 0; i < nMC; i++) aHisto->Fill((ampFit[i]-ampMean)/ampRMS);
  DrawHisto(cAmplitude, aHisto, Xaxis, Yaxis, ampSave);
  delete [] ampFit;
  delete [] ampFitError;
  delete cAmplitude;
  delete aHisto;

  cout << "Program Finished Normally" << endl;
}
