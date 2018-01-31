#include <iostream>
#include <algorithm>
#include <cmath>

// Toolkit
#include "STClibrary.h"
#include "rayleigh.h"
#include "simuevents.h"
#include "common.h"
#include "maptools.h"
#include "coverage.h"

// ROOT
#include "TROOT.h"
#include "TRint.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TStyle.h"

#ifdef gcc323
char* operator+( std::streampos&, char* );
#endif

static char gObjName[1024];
static int gObjNumber = 0;

using namespace std;



void Usage(string myName)
{
  cout << endl;
  cout << " Synopsis : " << endl;
  cout << myName << " <events file>" << endl << endl;
  
  cout << " Description :" << endl;
  cout << myName << " performs a Rayleigh analysis to the Right Ascension of the events. The <events file> must "
       << "contain the following fields :" << endl;
  
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
  if(argc != 2) Usage(argv[0]);
  string eventFile = argv[1];
  if( !CheckFile(eventFile) ) {cerr << "File: " << eventFile << " not found" << endl; exit(0);}


  // ROOT
  int fargc = 1;
  string extension;
  TRint* rint = new TRint("Rayleigh", &fargc, argv);
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
  if(events.size() == 0) {cout << "Program Failed : No events read. Exiting." << endl; exit(0);}
  
  // Min and Max dates
  DECLARE_VECTOR(double, utcstime, events, fUTCs);
  int minyear, minmonth, minday, maxyear, maxmonth, maxday;
  double minutcs, maxutcs, minhourutc, maxhourutc;
  minutcs = *min_element(utcstime.begin(),utcstime.end());
  maxutcs = *max_element(utcstime.begin(),utcstime.end());
  utcs2date(minutcs, &minyear, &minmonth, &minday, &minhourutc);
  utcs2date(maxutcs, &maxyear, &maxmonth, &maxday, &maxhourutc);
  cout << "First event : " << minday << "/" << minmonth << "/" << minyear << " UTC = " << minhourutc << endl;
  cout << "Last event : " << maxday << "/" << maxmonth << "/" << maxyear << " UTC = " << maxhourutc << endl;
  
  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                          Energy binning                                //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  unsigned int nBinsEnergy = 5;
  double * eMin = new double[nBinsEnergy];
  double eMax;
  vector<TEvent> eventsBin;
  vector<TEvent> eventsTmp1;
  vector<vector<TEvent> > allEvents;
  allEvents.resize(nBinsEnergy);
  
  DECLARE_VECTOR(double, energy, events,fEnergy);
  sort(energy.begin(),energy.end());
  eMin[0] = energy[0];
  eMax = 5.;
  for(unsigned int i = 1; i < nBinsEnergy; i++) eMin[i] = eMin[0]+i*(eMax-eMin[0])/(nBinsEnergy-1);

  for(unsigned int i = 0; i < nBinsEnergy; i++)
    {
      for(long j = 0; j < nEvents; j++) if(events[j].fEnergy >= eMin[i]) eventsTmp1.push_back(events[j]);
      allEvents[i].resize(eventsTmp1.size());
      for(long k = 0; k < (long)eventsTmp1.size(); k++) allEvents[i][k] = eventsTmp1[k];
      eventsTmp1.clear();
    }

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //               Rayleigh analysis for each energy bin                    //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////

  unsigned int nSide = 64;
  double thetaMin = 0.;
  double thetaMax = 60.;
  string accTimeModel = "TIME_FLAT";
 
  vector<TEvent> eventsTmp;
  double * amplitude = new double[nBinsEnergy];
  double * phase = new double[nBinsEnergy];
  double * significance = new double[nBinsEnergy];
  double * probability = new double[nBinsEnergy];
  double * ampLimit = new double[nBinsEnergy]; // Amp. corresponding to a chance probability of 5%

  for(unsigned int i = 0; i < nBinsEnergy; i++)
    {
      unsigned int nEvents = allEvents[i].size();
      eventsTmp.resize(nEvents);
      for(unsigned int j = 0; j < nEvents; j++) eventsTmp[j] = allEvents[i][j];

      // Compute the Coverage Map
      unsigned int nBins = 20;
      double chi2Theta = 1.5;
      double chi2Phi = 19;

      TCoverage coverage(nSide);
      coverage.SetCoordSystem('G');
      coverage.SetLatitude(kConstantsTK::AugerSouthLatitude);
      coverage.SetLongitude(kConstantsTK::AugerSouthLongitude);

      // Time modulation
      coverage.fTimeMod.SetAccTimeModel(accTimeModel);

      // Zenith angle modulation
      coverage.fThetaDist.SetExtension(extension);
      coverage.fThetaDist.SetAngleName("theta");
      DECLARE_VECTOR(double, thetaData, eventsTmp, fTheta);
      coverage.fThetaDist.SetData(thetaData);
      coverage.fThetaDist.SetOptions(thetaMin,thetaMax,nBins, chi2Theta); 

      // First : fitting with a simple geometric function
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
          cout << "Add Fermi-Dirac function in order to correctly fit zenith angle distribution" << endl;
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
          coverage.fThetaDist.fParameters[2].SetParameter(eFixed, thetaMin);
          coverage.fThetaDist.fParameters[3].SetParameter(eFixed, thetaMax);
          fitok = coverage.fThetaDist.Run();
        }
      if( !fitok )
        {
          cout << "Impossible to fit theta distribution. Try another fitting function. Exiting." << endl;
          exit(0);
        }
      
      // Take into account the zenith angle acceptance to compute the coverage map
      coverage.CorrectForAngularModulation("theta");
      
      // Azimut angle modulation
      double phimin = -180;
      double phimax = 180;
      unsigned int nThetaBins = 10;
      coverage.fPhiMod.SetAngleName("phi");
      coverage.fPhiMod.SetOptions(phimin, phimax, nBins, chi2Phi);
      coverage.fPhiMod.SetExtension(extension);
      // Set the function used to fit the azimut angle distribution in each theta bin
      // warning : we assume that the modulation follows the same law in each theta bin
      coverage.fPhiMod.fAngFitFunc = AugerPhiFunction;
      coverage.fPhiMod.SetNbParameters(2);
      // a(1+(b/a)cos(6*phi))
      // param 0 : a
      // param 1 : b/a
      coverage.fPhiMod.fDataMin = -180;
      coverage.fPhiMod.fDataMax = 180;
      coverage.fPhiMod.fParameters[0].SetParameter(eFree,1.);
      coverage.fPhiMod.fParameters[1].SetParameter(eFree,1.);
      coverage.fPhiMod.fDegreeMax = 0;
      
      // Set the function used to fit the parameters of the modulation as a function of theta
      TF1 * philaw = coverage.fPhiMod.SetPhiLaw(ModPhiThetaLaw, thetaMin, thetaMax,2);
      philaw->SetParameter(1,(thetaMin+thetaMax)/2);
      coverage.fPhiMod.ComputePhiModulation(events,nThetaBins);
      
      // Take into account the azimut angle distribution to compute the coverage map
      coverage.CorrectForAngularModulation("phi");
      
      // Take into account time modulation (UTC and/or JD)
      if(accTimeModel != "TIME_FLAT")
        {
          DECLARE_VECTOR(double,UTCh,eventsTmp,fUTCh);
          DECLARE_VECTOR(double,UTCs,eventsTmp,fUTCs);
          double chi2lim = 5;
          double timestep = 7200;
          coverage.fTimeMod.ComputeAccTime(UTCh,UTCs,timestep,chi2lim);
          coverage.CorrectForTimeModulation(accTimeModel);
        }

      // VERY IMPORTANT (initializes many constants to save time)
      // To be called once latsite, thetaMax and thetaMin are set
      coverage.ComputeDeclinationLimits();
      cout << "declination limits : " << coverage.GetDecMin() << " " << coverage.GetDecMax() << endl;
      
      // Compute the coverage map
      coverage.ComputeCoverage();
      
      // Put the coverage map into a Healpix map
      THealpixMap covMap = coverage.GetMap();
      // Normalize it to the total number of events
      covMap *= thetaData.size()*1./covMap.Total();
      
      //-------------------- Rayleigh --------------------//
      unsigned int nHarmonic = 1;
      TRayleigh* rayleigh = new TRayleigh(eventsTmp, covMap, nHarmonic);

      rayleigh->ComputeRAweight();
      rayleigh->ComputeAmplitude();
      rayleigh->ComputePhase();
      rayleigh->ComputeSignificance();
      rayleigh->ComputeChanceProbability();

      amplitude[i] = rayleigh->GetAmplitude();
      phase[i] = rayleigh->GetPhase();
      significance[i] = rayleigh->GetSignificance();
      probability[i] = rayleigh->GetChanceProbability();
      double exponent = -4.*log(0.05)*1./nEvents;
      double value = sqrt(exponent);
      ampLimit[i] = value;
      cout << "eMin = " << eMin[i] << endl;
      cout << "amplitude = " << amplitude[i] << endl; 
      cout  << "phase = " << phase[i] << endl;
      cout << "probability = " << probability[i] << endl;
      cout << "significance = " << significance[i] << endl;
      
      delete rayleigh;
      thetaData.clear();
      eventsTmp.clear();
    }

  ////////////////////////////////////////////////////////////////////////////
  //                                                                        //
  //                                  Plot                                  //
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////
 
  string Xaxis = "Energy [EeV]";

  // Amplitude
  string YaxisAmplitude = "Amplitude";
  string nameAmplitude = "Rayleigh Amplitude";
  string saveAmplitude = "rayleighAmplitude"+extension;
  sprintf(gObjName, "cTmpExRayleigh%d", gObjNumber++);
  TCanvas* cAmplitude = new TCanvas(gObjName, "Rayleigh Amplitude", 700, 700);
  cAmplitude->SetLogy();

  sprintf(gObjName, "cTmpExRayleigh%d", gObjNumber++);
  TH1F* hAmplitude = cAmplitude->DrawFrame(eMin[0], 0.001, eMax, 1.);
  hAmplitude->SetNameTitle(gObjName, nameAmplitude.c_str());
  DrawHisto(cAmplitude, hAmplitude, Xaxis, YaxisAmplitude, "");

  TGraph* plotAmplitude = new TGraph(nBinsEnergy, eMin, amplitude);
  plotAmplitude->SetMarkerStyle(20);
  plotAmplitude->SetMarkerColor(kBlack);
  plotAmplitude->SetMarkerSize(0.6);
  plotAmplitude->Draw("P");

  TGraph* plotAmpLimit = new TGraph(nBinsEnergy, eMin, ampLimit);
  plotAmpLimit->SetLineStyle(2);
  plotAmpLimit->SetLineWidth(1);
  plotAmpLimit->SetLineColor(kRed);
  plotAmpLimit->Draw("L");

  cAmplitude->Update();
  cAmplitude->SaveAs(saveAmplitude.c_str());

  // Phase
  string YaxisPhase = "Phase";
  string namePhase = "Rayleigh Phase";
  string savePhase = "rayleighPhase"+extension;
  sprintf(gObjName, "cTmpExRayleigh%d", gObjNumber++);
  TCanvas *cPhase = new TCanvas(gObjName, namePhase.c_str(), 700, 700);
  cPhase->SetGrid();
  TGraphErrors* plotPhase = new TGraphErrors(nBinsEnergy, eMin, phase);

  PlotXY(cPhase, plotPhase, eMin[0], eMax, namePhase, Xaxis, YaxisPhase);
  plotPhase->Draw("AP");
  cPhase->Update();
  cPhase->SaveAs(savePhase.c_str());

  // Probability
  string YaxisProbability = "Probability";
  string nameProbability = "Rayleigh Probability";
  string saveProbability = "rayleighProbability"+extension;
  sprintf(gObjName, "cTmpExRayleigh%d", gObjNumber++);
  TCanvas *cProbability = new TCanvas(gObjName, nameProbability.c_str(), 700, 700);
  cProbability->SetLogy();
  cProbability->SetGrid();
  TGraphErrors* plotProbability = new TGraphErrors(nBinsEnergy, eMin, 
						   probability);

  PlotXY(cProbability, plotProbability, eMin[0], eMax, nameProbability, Xaxis, YaxisProbability);
  plotProbability->Draw("AP");
  cProbability->Update();
  cProbability->SaveAs(saveProbability.c_str());

  // Significance
  string YaxisSignificance = "Significance";
  string nameSignificance = "Rayleigh Significance";
  string saveSignificance = "rayleighSignificance"+extension;
  sprintf(gObjName, "cTmpExRayleigh%d", gObjNumber++);
  TCanvas *cSignificance = new TCanvas(gObjName, nameSignificance.c_str(), 700, 700);
  cSignificance->SetGrid();
  TGraphErrors* plotSignificance = new TGraphErrors(nBinsEnergy, eMin, significance);

  PlotXY(cSignificance, plotSignificance, eMin[0], eMax, nameSignificance, Xaxis, YaxisSignificance);
  plotSignificance->Draw("AP");
  cSignificance->Update();
  cSignificance->SaveAs(saveSignificance.c_str());

  delete [] eMin;
  delete [] amplitude;
  delete [] phase;
  delete [] probability;
  delete [] significance;
  delete [] ampLimit;

  cout << "Program Finished Normally" << endl;
  rint->Run(kTRUE);
}
