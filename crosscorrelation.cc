#include <iostream>
#include <fstream>
#include <iomanip>

#include "crosscorrelation.h"
#include "TCanvas.h"
#include "TMath.h"
#include "simuevents.h"
#include "common.h"


static TProgressBar gProgB;
static char gObjName[1024];
static int gObjNumber = 0 ;
static char * GetObjName() {sprintf(gObjName,"CrossCorrelation%d",gObjNumber++); return gObjName;}


TCrossCorrelation::TCrossCorrelation(const vector<TEvent>& events, const vector<double>& lSources, const vector<double>& bSources, const vector<double>& wSources, double alphaStep, double alphaMax)
  : fEvents(events), fLsources(lSources), fBsources(bSources), fWsources(wSources), fAlphaStep(alphaStep), fAlphaMax(alphaMax)
{
  fNevents = fEvents.size();
  fNsources = fLsources.size();
  fExtension = ".pdf";
  InitGraph();
  InitAnglesSources();
}



TCrossCorrelation::~TCrossCorrelation()
{
  if( fEventsCrossCorrelation   ) delete fEventsCrossCorrelation;
  if( fCoverageCrossCorrelation ) delete fCoverageCrossCorrelation;
  if( fThetaEvents              ) delete [] fThetaEvents;
  if( fPhiEvents                ) delete [] fPhiEvents;
  if( fCosThetaEvents           ) delete [] fCosThetaEvents;
  if( fSinThetaEvents           ) delete [] fSinThetaEvents;
}



void TCrossCorrelation::InitGraph()
{
  fNstep = (int)floor(fAlphaMax/fAlphaStep);
  fEventsCrossCorrelation = new TGraph(fNstep);
  fCoverageCrossCorrelation = new TGraphAsymmErrors(fNstep);
}



void TCrossCorrelation::InitAnglesEvents(const vector<TEvent>& events)
{
  DECLARE_POINTER(double, thetaEvents, events, fB);
  DECLARE_POINTER(double, phiEvents, events, fL);

  fThetaEvents = new double[fNevents]; fPhiEvents = new double[fNevents];
  for(unsigned int i = 0; i < fNevents; i++) {fThetaEvents[i] = 90.-thetaEvents[i]; fPhiEvents[i] = phiEvents[i];}

  fCosThetaEvents = new double[fNevents]; fSinThetaEvents = new double[fNevents];
  for(unsigned int i = 0; i < fNevents; i++) 
    {
      fCosThetaEvents[i] = cos(fThetaEvents[i]*DTOR);
      fSinThetaEvents[i] = sin(fThetaEvents[i]*DTOR);
    }

  delete [] thetaEvents;
  delete [] phiEvents;
}



void TCrossCorrelation::InitAnglesSources()
{
  fThetaSources = new double[fNsources]; fPhiSources = new double[fNsources];
  for(unsigned int i = 0; i < fNsources; i++) {fThetaSources[i] = 90.-fBsources[i]; fPhiSources[i] = fLsources[i];}

  fCosThetaSources = new double[fNsources]; fSinThetaSources = new double[fNsources];
  for(unsigned int i = 0; i < fNsources; i++)
    {
      fCosThetaSources[i] = cos(fThetaSources[i]*DTOR);
      fSinThetaSources[i] = sin(fThetaSources[i]*DTOR);
    }
}



void TCrossCorrelation::ComputeEventsCrossCorrelation()
{
  InitAnglesEvents(fEvents);
  double dt;
  TH1F * histoEvents = new TH1F(GetObjName(), "", fNstep, 0, fAlphaMax);
  for(unsigned int i = 0; i < fNevents; i++)
    {
      for(unsigned int j = 0; j < fNsources; j++)
        {
          dt = acos(fSinThetaEvents[i]*fSinThetaSources[j]*cos((fPhiEvents[i]-fPhiSources[j])*DTOR)+fCosThetaEvents[i]*fCosThetaSources[j]);
          if(dt/DTOR <= fAlphaMax) histoEvents->Fill(dt/DTOR,fWsources[j]);
        }	
    }
  for(unsigned int i = 0; i < fNstep; i++)
    fEventsCrossCorrelation->SetPoint(i,histoEvents->GetBinCenter(i+1),histoEvents->Integral(1,i+1));

  delete histoEvents;
}



void TCrossCorrelation::ComputeCoverageCrossCorrelation(double longitude, double latitude, const THealpixMap& map, unsigned int nSimu, double dispersion, const vector<double>& thVal, const vector<double>& pthVal, string utcFile, string jdFile, string globalFile)
{
  vector<vector<double> > Np(fNstep);
  vector<double> BinCenter(fNstep);
  gProgB.Zero();
  gProgB.fBegin = 0;
  gProgB.fEnd = nSimu;
  gProgB.InitPercent();
  for(unsigned int i = 0; i < nSimu; i++)
    {
      vector<TEvent> simEvents = SimulateEvents(map, fNevents, thVal, pthVal, latitude, longitude, utcFile, jdFile, globalFile);
      
      InitAnglesEvents(simEvents);
      TH1F * histoCoverage = new TH1F(GetObjName(), "", fNstep, 0, fAlphaMax);
      double dt;
      for(unsigned int j = 0; j < fNevents; j++)
        {
          for(unsigned int k = 0; k < fNsources; k++)
            {
              dt = acos(fSinThetaEvents[j]*fSinThetaSources[k]*cos((fPhiEvents[j]-fPhiSources[k])*DTOR)+fCosThetaEvents[j]*fCosThetaSources[k]);
              if(dt/DTOR <= fAlphaMax) histoCoverage->Fill(dt/DTOR,fWsources[k]);
            }
        }
      
      for(unsigned int j = 0; j < fNstep; j++)
        {
          if(i == 0) {Np[j].resize(nSimu); BinCenter[j] = histoCoverage->GetBinCenter(j+1);}
          if(j == 0) Np[j][i] = histoCoverage->GetBinContent(j+1);
          else Np[j][i] = Np[j-1][i]+histoCoverage->GetBinContent(j+1);
        }
      delete histoCoverage;
      gProgB.PrintPercent(i);
    }
  gProgB.EndPercent();

  vector<double> NpTmp(nSimu);
  double mean, dispersionHigh, dispersionLow;
  for(unsigned int i = 0; i < fNstep; i++)
    {
      mean = 0., dispersionHigh = 0., dispersionLow = 0.;
      for(unsigned int j = 0; j < nSimu; j++) {mean += Np[i][j]; NpTmp[j] = Np[i][j];}
      mean /= nSimu;
      sort(NpTmp.begin(),NpTmp.end());
      unsigned int indexMean = 0;
      for(unsigned int j = 0; j < nSimu; j++) if(NpTmp[j] >= mean) {indexMean = j; break;}
      dispersionHigh = NpTmp[indexMean+(unsigned int)ceil((nSimu-indexMean)*dispersion)];
      dispersionLow = NpTmp[indexMean-(unsigned int)ceil(indexMean*dispersion)];
      fCoverageCrossCorrelation->SetPoint(i,BinCenter[i],mean);
      fCoverageCrossCorrelation->SetPointError(i,0.,0.,fabs(mean-dispersionLow),fabs(mean-dispersionHigh));
    }
}



void TCrossCorrelation::DrawEvents() const
{
  string Xaxis = "#alpha"; string Yaxis = "n_{p}"; string save = "";
  TCanvas * cEvents = new TCanvas(GetObjName(), "Events", 800, 800);
  double minX = 0., maxX, minY = 0, maxY;
  fEventsCrossCorrelation->GetPoint(fNstep-1,maxX,maxY);
  TH1F * hEvents = cEvents->DrawFrame(minX,minY,maxX+fAlphaStep,maxY+10);
  DrawHisto(cEvents, hEvents, Xaxis.c_str(), Yaxis.c_str(), save.c_str());
   
  fEventsCrossCorrelation->SetMarkerStyle(kFullCircle); fEventsCrossCorrelation->SetMarkerSize(1);
  fEventsCrossCorrelation->Draw("P");

  string saveEvents = "crossCorrelationEvents"+fExtension;
  cEvents->SaveAs(saveEvents.c_str());

}



void TCrossCorrelation::DrawCov() const
{
  string Xaxis = "#alpha"; string Yaxis = "n_{p}"; string save = "";
  TCanvas* cCoverage = new TCanvas(GetObjName(), "Coverage", 800, 800);
  double minX = 0., maxX, minY = 0, maxY;
  fCoverageCrossCorrelation->GetPoint(fNstep-1,maxX,maxY);
  TH1F * hCoverage = cCoverage->DrawFrame(minX,minY,maxX+fAlphaStep,maxY+fCoverageCrossCorrelation->GetErrorYhigh(fNstep-1)+10);
  DrawHisto(cCoverage, hCoverage, Xaxis.c_str(), Yaxis.c_str(), save.c_str());

  fCoverageCrossCorrelation->SetMarkerStyle(kMultiply);
  fCoverageCrossCorrelation->Draw("P");

  string saveCov = "crossCorrelationCov"+fExtension;
  cCoverage->SaveAs(saveCov.c_str());

}



void TCrossCorrelation::DrawBoth() const
{
  string Xaxis = "#alpha"; string Yaxis = "n_{p}"; string save = "";
  TCanvas* cBoth = new TCanvas(GetObjName(),"Both", 800, 800);
  double minX = 0., maxX, minY = 0, maxY;
  fCoverageCrossCorrelation->GetPoint(fNstep-1,maxX,maxY);
  TH1F* hBoth = cBoth->DrawFrame(minX,minY,maxX+fAlphaStep,maxY+fCoverageCrossCorrelation->GetErrorYhigh(fNstep-1)+10);
  DrawHisto(cBoth, hBoth, Xaxis.c_str(), Yaxis.c_str(), save.c_str());
  
  fCoverageCrossCorrelation->SetMarkerStyle(kMultiply); fCoverageCrossCorrelation->Draw("P");
  fEventsCrossCorrelation->SetMarkerStyle(kFullCircle); fEventsCrossCorrelation->SetMarkerSize(1); fEventsCrossCorrelation->SetMarkerColor(kRed);
  fEventsCrossCorrelation->Draw("same P");

  string saveBoth = "crossCorrelationBoth"+fExtension;
  cBoth->SaveAs(saveBoth.c_str());
}



void TCrossCorrelation::DrawBothRelativeExcess() const
{
  double alpha, NpEvents, NpCoverage, errLow, errHigh;
  TGraph * eventsRelativeExcess = new TGraph(fNstep);
  TGraphAsymmErrors * coverageRelativeExcess = new TGraphAsymmErrors(fNstep);
  for(unsigned int i = 0; i < fNstep; i++)
    {
      fEventsCrossCorrelation->GetPoint(i,alpha,NpEvents);
      fCoverageCrossCorrelation->GetPoint(i,alpha,NpCoverage);
      errLow = fCoverageCrossCorrelation->GetErrorYlow(i);
      errHigh = fCoverageCrossCorrelation->GetErrorYhigh(i);
      
      eventsRelativeExcess->SetPoint(i,alpha,(NpEvents-NpCoverage)/NpCoverage);
      coverageRelativeExcess->SetPoint(i,alpha,0.);
      coverageRelativeExcess->SetPointError(i,0.,0.,errLow/NpCoverage,errHigh/NpCoverage);
    }

  TCanvas * cBothRelativeExcess = new TCanvas(GetObjName(),"BothRelativeExcess", 800, 800);
  TH1F* hBoth = cBothRelativeExcess->DrawFrame(0.,-2,fAlphaMax+fAlphaStep,2);
  
  string Xaxis = "#alpha"; string Yaxis = "Relative Excess of pairs"; string save = "";
  DrawHisto(cBothRelativeExcess, hBoth, Xaxis.c_str(), Yaxis.c_str(), save.c_str());
  eventsRelativeExcess->SetMarkerColor(kRed); eventsRelativeExcess->SetMarkerStyle(kFullCircle); eventsRelativeExcess->SetMarkerSize(1);
  coverageRelativeExcess->SetFillColor(18);
  coverageRelativeExcess->Draw("3"); eventsRelativeExcess->Draw("same P");
  string saveBoth = "crossCorrelationBothRelativeExcess"+fExtension;
  cBothRelativeExcess->SaveAs(saveBoth.c_str());
}



void TCrossCorrelation::GetEventsCrossCorrelation(vector<double>& alpha, vector<double>& Np) const
{
  Np.resize(fNstep); alpha.resize(fNstep);
  for(unsigned int i = 0; i < fNstep;i++) fEventsCrossCorrelation->GetPoint(i,alpha[i],Np[i]);
}



void TCrossCorrelation::GetCoverageCrossCorrelation(vector<double>& alpha,vector<double>& Np,vector<double>& errLow,vector<double>& errHigh) const
{
  Np.resize(fNstep); alpha.resize(fNstep); errLow.resize(fNstep); errHigh.resize(fNstep);
  for(unsigned int i = 0; i < fNstep; i++)
    {
      fCoverageCrossCorrelation->GetPoint(i,alpha[i],Np[i]);
      errLow[i] = fCoverageCrossCorrelation->GetErrorYlow(i);
      errHigh[i] = fCoverageCrossCorrelation->GetErrorYhigh(i);
    
    }
}

