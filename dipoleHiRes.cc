#include <iostream>

#include "dipoleHiRes.h"
#include "common.h"
#include "TCanvas.h"
#include "simuevents.h"
#include "maptools.h"
#include "harmotools.h"

#include "TMath.h"

static TProgressBar gProgB;
static char gObjName[1024];
static int gObjNumber = 0;
static char * GetObjName() {sprintf(gObjName, "DipoleHiRes%d", gObjNumber++); return gObjName;}



TDipoleHiRes::TDipoleHiRes(const vector<TEvent> &events, double lDipole, double bDipole, double dCosTheta)
  : fEvents(events), fdCosTheta(dCosTheta), fLdipole(lDipole), fBdipole(bDipole)
{
  fNevents = fEvents.size();
  fExtension = ".png";
  InitHisto();
  fFitFcn = 0x0;
}



TDipoleHiRes::~TDipoleHiRes()
{
  if(fFitFcn     ) delete fFitFcn;
  if(fHistoEvents) delete fHistoEvents;
  if(fHistoCov   ) delete fHistoCov;
  if(fHistoTrue  ) delete fHistoTrue;
  if(fCosTheta   ) delete [] fCosTheta;
}



void TDipoleHiRes::InitHisto()
{
  fNbins = (unsigned int)floor(2./fdCosTheta);
  fHistoEvents = new TH1F(GetObjName(), "Events Opening Angle", fNbins, -1, 1);
  fHistoCov = new TH1F(GetObjName(), "Coverage Opening Angle", fNbins, -1, 1);
  fHistoTrue = new TH1F(GetObjName(), "True Opening Angle", fNbins, -1, 1);
  fCosThetaBins.resize(fNbins);
  for(unsigned int i = 0; i < fNbins; i++) fCosThetaBins[i] = fHistoEvents->GetBinCenter(i+1);
}



void TDipoleHiRes::ComputeOpeningAngle(const vector<TEvent> &events)
{
  vector<double> uvDipole = ll2uv(fLdipole, fBdipole);
  fCosTheta = new double[fNevents];
  for(unsigned int i = 0; i < fNevents; i++) 
    {
      vector<double> uvEvt = ll2uv(events[i].fL, events[i].fB);
      fCosTheta[i] = uvEvt[0]*uvDipole[0]+uvEvt[1]*uvDipole[1]+uvEvt[2]*uvDipole[2];
    }
}



void TDipoleHiRes::ComputeEventsOpeningAngle()
{
  ComputeOpeningAngle(fEvents);
  for(unsigned int i = 0; i < fNevents; i++) fHistoEvents->Fill(fCosTheta[i]);
}



void TDipoleHiRes::ComputeCoverageOpeningAngle(double longitude, double latitude, const THealpixMap &map, unsigned int nSimu, const vector<double> &thVal, const vector<double> &pthVal, string utcFile, string jdFile, string globalFile)
{
  vector<vector<unsigned int> > binContent(fNbins);
  gProgB.Zero();
  gProgB.fBegin = 0;
  gProgB.fEnd = nSimu;
  gProgB.InitPercent();
  for(unsigned int i = 0; i < nSimu; i++)
    {
      vector<TEvent> simEvents = SimulateEvents(map, fNevents, thVal, pthVal, latitude, longitude, utcFile, jdFile, globalFile);
   
      ComputeOpeningAngle(simEvents);
      TH1F* hCoverage = new TH1F(GetObjName(), "", fNbins, -1, 1);
      for(unsigned int j = 0; j < simEvents.size(); j++) hCoverage->Fill(fCosTheta[j]);

      for(unsigned int j = 0; j < fNbins; j++)
      {
          if(i == 0) binContent[j].resize(nSimu);
          binContent[j][i] = hCoverage->GetBinContent(j+1);
      }
      delete hCoverage;
      gProgB.PrintPercent(i);
    }
  gProgB.EndPercent();


  long * binContentTmp = new long[nSimu];
  double avgSimu, stdSimu;
  for(unsigned int i = 0; i < fNbins; i++)
    {
      avgSimu = 0., stdSimu = 0.;
      for(unsigned int j = 0; j < nSimu; j++) {avgSimu += binContent[i][j]; binContentTmp[j] = binContent[i][j];}
      avgSimu /= nSimu;
      fHistoCov->SetBinContent(i+1,avgSimu);
      stdSimu = TMath::RMS(nSimu,binContentTmp);
      fHistoCov->SetBinError(i+1,stdSimu);
    }
  delete [] binContentTmp;
}




void TDipoleHiRes::ComputeTrueOpeningAngle()
{
  for(unsigned int i = 0; i < fNbins; i++) fHistoTrue->SetBinContent(i+1,(*fHistoEvents)[i+1]);  
  fHistoTrue->Divide(fHistoCov);
}



void TDipoleHiRes::fitDipole()
{
  if( !fHistoTrue )
  {
      cout << "The fit cannot be performed. The opening angle histogram have to be filled first" << endl;
      exit(0);
  }
  unsigned int nPar = 2;
  fFitFcn = new TF1(GetObjName(), linearFunction, -1, 1, nPar);

  // Set Fit default parameters
  fFitFcn->SetParameter(0,0); // y intercept
  fFitFcn->SetParameter(1,1); // slope
  // Set parameter names
  string parNames[2];
  parNames[0] = "y intercept"; parNames[1] = "slope";

  fHistoTrue->Fit(fFitFcn,"QN");
  double chi2perNDF = (fFitFcn->GetChisquare())/(fFitFcn->GetNDF());
  cout << "Chi2/NDF = " << chi2perNDF << endl;

  // Fitted parameters
  double getParameter[nPar];
  fFitFcn->GetParameters(getParameter);
  double error[nPar];
  for(unsigned int i = 0; i < nPar; i++) error[i] = fFitFcn->GetParError(i);

  fFitParameters.resize(nPar); fFitParametersErrors.resize(nPar);
  for(unsigned int i = 0; i < nPar; i++)
    {
      fFitParameters[i] = getParameter[i];
      fFitParametersErrors[i] = error[i];
      cout << parNames[i].c_str() << " : " << fFitParameters[i] << " error : " << fFitParametersErrors[i] << endl;
    }
}


void TDipoleHiRes::DrawEvents() const
{
  string Xaxis = "Cos #alpha"; string Yaxis = "Count"; string save = "";
  TCanvas* cEvents = new TCanvas(GetObjName(), "Events Opening Angle", 700, 700);

  DrawHisto(cEvents, fHistoEvents, Xaxis.c_str(), Yaxis.c_str(), save.c_str());
  fHistoEvents->SetMarkerStyle(20); fHistoEvents->SetMarkerSize(0.8); fHistoEvents->SetMarkerColor(kRed);
  fHistoEvents->Draw("P");
  string saveEvents = "dipoleHiResEvents"+fExtension;
  cEvents->SaveAs(saveEvents.c_str());
}



void TDipoleHiRes::DrawCov() const
{
  string Xaxis = "Cos #alpha"; string Yaxis = "Count"; string save = "";

  TCanvas* cCoverage = new TCanvas(GetObjName(), "Coverage Opening Angle", 700, 700);

  DrawHisto(cCoverage, fHistoCov, Xaxis.c_str(), Yaxis.c_str(), save.c_str());
  fHistoCov->SetMarkerStyle(20); fHistoCov->SetMarkerSize(0.8);
  fHistoCov->Draw("e1p");
  string saveCov = "dipoleHiResCov"+fExtension;
  cCoverage->SaveAs(saveCov.c_str());
}



void TDipoleHiRes::DrawBoth() const
{
  string Xaxis = "Cos #alpha"; string Yaxis = "Count"; string save = "";
  TCanvas* cBoth = new TCanvas(GetObjName(), "Both", 700, 700);
  
  fHistoEvents->SetMarkerStyle(20); fHistoEvents->SetMarkerSize(0.8); fHistoEvents->SetMarkerColor(kRed);
  DrawHisto(cBoth, fHistoEvents, Xaxis.c_str(), Yaxis.c_str(), save.c_str());
  fHistoCov->SetMarkerStyle(20); fHistoCov->SetMarkerSize(0.8);
  fHistoCov->Draw("e1p"); fHistoEvents->Draw("same P");
  string saveBoth = "dipoleHiResBoth"+fExtension;
  cBoth->SaveAs(saveBoth.c_str());
}



void TDipoleHiRes::DrawTrue() const
{
  string Xaxis = "Cos #alpha"; string Yaxis = "Normalized Event Count"; string save = "dipoleHiResTrue"+fExtension;

  TCanvas * cTrue = new TCanvas(GetObjName(),"True Opening Angle", 700, 700);
  fHistoTrue->SetMarkerStyle(20); fHistoTrue->SetMarkerSize(0.8);
  DrawHisto(cTrue, fHistoTrue, Xaxis.c_str(), Yaxis.c_str(), save.c_str());
}



void TDipoleHiRes::DrawTrueFit() const
{
  string Xaxis = "Cosinus Opening Angle";
  string Yaxis = "Normalized Event Count";
  string save = "";
  
  TCanvas * cTrueFit = new TCanvas(GetObjName(),"True Opening Angle", 700, 700);

  DrawHisto(cTrueFit, fHistoTrue, Xaxis.c_str(), Yaxis.c_str(), save.c_str());
  fHistoTrue->SetMarkerStyle(20); fHistoTrue->SetMarkerSize(0.8);
  fHistoTrue->Draw("e1p");
  fFitFcn->SetLineWidth(1); fFitFcn->SetLineColor(kRed);
  fFitFcn->Draw("same");
  cTrueFit->Update();
  string saveTrueFit = "dipoleHiResTrueFit"+fExtension;
  cTrueFit->SaveAs(saveTrueFit.c_str());
}



vector<double> TDipoleHiRes::GetEventsCounts() const
{
  vector<double> tmp(fNbins);
  for(unsigned int i = 0; i < fNbins; i++) tmp[i] = fHistoEvents->GetBinContent(i+1);
  return tmp;
}



vector<double> TDipoleHiRes::GetCoverageCounts() const
{
  vector<double> tmp(fNbins);
  for(unsigned int i = 0; i < fNbins; i++) tmp[i] = fHistoCov->GetBinContent(i+1);
  return tmp;
}



vector<double> TDipoleHiRes::GetTrueValues() const
{
  vector<double> tmp(fNbins);
  for(unsigned int i = 0; i < fNbins; i++) tmp[i] = fHistoTrue->GetBinContent(i+1);
  return tmp;
}



double linearFunction(double *x, double *par)
{
  return par[0]+par[1]*x[0];
}
