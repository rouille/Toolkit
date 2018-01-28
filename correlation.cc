#include <iostream>
#include <fstream>
#include <iomanip>

#include "TCanvas.h"
#include "TMath.h"

#include "harmotools.h"
#include "Cl.h"
#include "simuevents.h"
#include "maptools.h"
#include "common.h"
#include "correlation.h"


static TProgressBar gProgB;
static char gObjName[1024];
static int gObjNumber = 0 ;
static char * GetObjName() {sprintf(gObjName,"Correlation%d",gObjNumber++); return gObjName;}


TCorrelation::TCorrelation(const vector<TEvent>& events, unsigned int alphaStep, double alphaMax)
  : fEvents(events), fAlphaStep(alphaStep), fAlphaMax(alphaMax)
{
  fNevents = fEvents.size();
  fExtension = ".pdf";
  InitGraph();
}



TCorrelation::~TCorrelation()
{
  if( fEventsCorrelation   ) delete fEventsCorrelation;
  if( fCoverageCorrelation ) delete fCoverageCorrelation;
  if( fTheta               ) delete [] fTheta;
  if( fPhi                 ) delete [] fPhi;
  if( fCosTheta            ) delete [] fCosTheta;
  if( fSinTheta            ) delete [] fSinTheta;
}



void TCorrelation::InitGraph()
{
  fNstep = (int)floor(fAlphaMax/fAlphaStep);
  fEventsCorrelation = new TGraph(fNstep); 
  fCoverageCorrelation = new TGraphAsymmErrors(fNstep);
}



void TCorrelation::InitAngles(const vector<TEvent> &events)
{
  DECLARE_POINTER(double, thetaEvents, events, fB);
  DECLARE_POINTER(double, phiEvents, events, fL);

  fTheta = new double[fNevents]; fPhi = new double[fNevents];
  for(unsigned int i = 0; i < fNevents; i++) {fTheta[i] = 90.-thetaEvents[i]; fPhi[i] = phiEvents[i];}

  fCosTheta = new double[fNevents]; fSinTheta = new double[fNevents];
  for(unsigned int i = 0; i < fNevents; i++) {fCosTheta[i] = cos(fTheta[i]*DTOR); fSinTheta[i] = sin(fTheta[i]*DTOR);}

  delete [] thetaEvents;
  delete [] phiEvents;
}



void TCorrelation::ComputeEventsCorrelation()
{
  InitAngles(fEvents);
  double dt;
  TH1F * histoEvents = new TH1F(GetObjName(), "", fNstep, 0, fAlphaMax);
  for(unsigned int i = 0; i < fNevents; i++)
    {
      for(unsigned int j = i+1; j < fNevents; j++)
        {
          dt = acos(fSinTheta[i]*fSinTheta[j]*cos((fPhi[i]-fPhi[j])*DTOR)+fCosTheta[i]*fCosTheta[j]);
          if(dt/DTOR <= fAlphaMax) histoEvents->Fill(dt/DTOR);
        }
    }
  
  for(unsigned int i = 0; i < fNstep; i++)
    fEventsCorrelation->SetPoint(i,histoEvents->GetBinCenter(i+1),histoEvents->Integral(1,i+1));
  
  delete histoEvents;
}



void TCorrelation::ComputeCoverageCorrelationMap(const THealpixMap& covMap)
{
  int nSide = covMap.NSide();
  int lMax = 3*nSide-1;
  vector<double> Cl = covMap.Map2Cl(lMax);
  vector<double> pixWindow = GetPixWindow(nSide,lMax);
  for(unsigned int i = 0; i < Cl.size(); i++) Cl[i] /= (pixWindow[i]*pixWindow[i]);

  TH1F * histoCoverage = new TH1F(GetObjName(), "", fNstep, 0, fAlphaMax); 
  ComputeCtheta(lMax, &histoCoverage, Cl);
  
  double alpha, normEvents;
  fEventsCorrelation->GetPoint(fNstep-1,alpha,normEvents);
  double normCoverage = histoCoverage->Integral(1,fNstep);
  for(unsigned int i = 0; i < fNstep; i++)
    { 
      fCoverageCorrelation->SetPoint(i,histoCoverage->GetBinCenter(i+1),(normEvents/normCoverage)*histoCoverage->Integral(1,i+1));
      fCoverageCorrelation->SetPointError(i,0.,0.,0.,0.);
    }
 delete histoCoverage; 
}



void TCorrelation::ComputeCoverageCorrelationSimu(double longitude, double latitude, const THealpixMap &map, unsigned int nSimu, double dispersion, const vector<double> &thVal, const vector<double> &pthVal, string utcFile, string jdFile, string globalFile)
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
      
      InitAngles(simEvents);
      TH1F * histoCoverage = new TH1F(GetObjName(), "", fNstep, 0, fAlphaMax);
      double dt;
      for(unsigned int j = 0; j < fNevents; j++)
        {
          for(unsigned int k = j+1; k < fNevents; k++)
            {
              dt = acos(fSinTheta[j]*fSinTheta[k]*cos((fPhi[j]-fPhi[k])*DTOR)+fCosTheta[j]*fCosTheta[k]);
              if(dt/DTOR <= fAlphaMax) histoCoverage->Fill(dt/DTOR);
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
      fCoverageCorrelation->SetPoint(i,BinCenter[i],mean);
      fCoverageCorrelation->SetPointError(i,0.,0.,fabs(mean-dispersionLow),fabs(mean-dispersionHigh));
    }
}



void TCorrelation::ComputeCoverageCorrelationMC(double longitude, double latitude, unsigned int nMC, double dispersion, string binningType, string scramblingType, double thetaMax)
{
  vector<vector<double> > Np(fNstep);
  vector<double> BinCenter(fNstep);
  gProgB.Zero();
  gProgB.fBegin = 0;
  gProgB.fEnd = nMC;
  gProgB.InitPercent();
  for(unsigned int i = 0; i < nMC; i++)
    {
      unsigned int nBins = 6;
      vector<TEvent> mcEvents = ScrambleData(fEvents, nBins, binningType, latitude, longitude, thetaMax, scramblingType);
      InitAngles(mcEvents);
      double dt;
      TH1F * histoCoverage = new TH1F(GetObjName(), "", fNstep, 0, fAlphaMax);
      for(unsigned int j = 0 ; j < fNevents; j++)
        {
          for(unsigned int k = j+1; k < fNevents; k++)
            {
              dt = acos(fSinTheta[j]*fSinTheta[k]*cos((fPhi[j]-fPhi[k])*DTOR)+fCosTheta[j]*fCosTheta[k]);
              if(dt/DTOR <= fAlphaMax) histoCoverage->Fill(dt/DTOR);
            }
        }
        
      for(unsigned int j = 0; j < fNstep; j++)
        {
          if(i == 0) {Np[j].resize(nMC); BinCenter[j] = histoCoverage->GetBinCenter(j+1);}
          if(j == 0) Np[j][i] = histoCoverage->GetBinContent(j+1);
          else Np[j][i] = Np[j-1][i]+histoCoverage->GetBinContent(j+1);
        }
      delete histoCoverage;
      gProgB.PrintPercent(i);
    }
  gProgB.EndPercent();
  
  vector<double> NpTmp(nMC);
  double mean, dispersionHigh, dispersionLow;
  for(unsigned int i = 0; i < fNstep; i++)
    {
      mean = 0., dispersionHigh = 0., dispersionLow = 0.;
      for(unsigned int j = 0; j < nMC; j++) {mean += Np[i][j]; NpTmp[j] = Np[i][j];}
      mean /= nMC;
      sort(NpTmp.begin(),NpTmp.end());
      unsigned int indexMean = 0;
      for(unsigned int j = 0; j < nMC; j++) if(NpTmp[j] >= mean) {indexMean = j; break;}
      dispersionHigh = NpTmp[indexMean+(unsigned int)ceil((nMC-indexMean)*dispersion)];
      dispersionLow = NpTmp[indexMean-(unsigned int)ceil(indexMean*dispersion)];
      fCoverageCorrelation->SetPoint(i,BinCenter[i],mean);
      fCoverageCorrelation->SetPointError(i,0.,0.,fabs(mean-dispersionLow),fabs(mean-dispersionHigh));
    }
}



void TCorrelation::DrawEvents() const
{
  string Xaxis = "Separation Angle"; string Yaxis = "Number of Pairs"; string save = "";
  TCanvas * cEvents = new TCanvas(GetObjName(), "Events", 600, 600);
  double minX = 0., maxX, minY = 0, maxY;
  fEventsCorrelation->GetPoint(fNstep-1,maxX,maxY);
  TH1F * hEvents = cEvents->DrawFrame(minX,minY,maxX+fAlphaStep,maxY+10);
  DrawHisto(cEvents, hEvents, Xaxis.c_str(), Yaxis.c_str(), save.c_str());
  
  fEventsCorrelation->SetMarkerStyle(kFullCircle); fEventsCorrelation->SetMarkerSize(1);
  fEventsCorrelation->Draw("P");

  string saveEvents = "correlationEvents"+fExtension;
  cEvents->SaveAs(saveEvents.c_str());
}



void TCorrelation::DrawCov() const
{
  string Xaxis = "Separation Angle"; string Yaxis = "Number of Pairs"; string save = "";
  TCanvas* cCoverage = new TCanvas(GetObjName(), "Coverage", 600, 600);
  double minX = 0., maxX, minY = 0, maxY;
  fCoverageCorrelation->GetPoint(fNstep-1,maxX,maxY);
  TH1F * hCoverage = cCoverage->DrawFrame(minX,minY,maxX+fAlphaStep,maxY+fCoverageCorrelation->GetErrorYhigh(fNstep-1)+10);
  DrawHisto(cCoverage, hCoverage, Xaxis.c_str(), Yaxis.c_str(), save.c_str());
  
  fCoverageCorrelation->SetMarkerStyle(kMultiply);
  fCoverageCorrelation->Draw("P");

  string saveCov = "correlationCov"+fExtension;
  cCoverage->SaveAs(saveCov.c_str());
}



void TCorrelation::DrawBoth() const
{
  string Xaxis = "Separation Angle"; string Yaxis = "Number of Pairs"; string save = "";
  TCanvas* cBoth = new TCanvas(GetObjName(),"Both", 600, 600);
  double minX = 0., maxX, minY = 1., maxY;
  fCoverageCorrelation->GetPoint(fNstep-1,maxX,maxY);
  TH1F* hBoth = cBoth->DrawFrame(minX,minY,maxX+fAlphaStep,maxY+fCoverageCorrelation->GetErrorYhigh(fNstep-1)+10);
  DrawHisto(cBoth, hBoth, Xaxis.c_str(), Yaxis.c_str(), save.c_str());
  
  fCoverageCorrelation->SetMarkerStyle(kMultiply); fCoverageCorrelation->Draw("P");
  fEventsCorrelation->SetMarkerStyle(kFullCircle); fEventsCorrelation->SetMarkerSize(1); fEventsCorrelation->SetMarkerColor(kRed);
  cBoth->SetLogy();
  fEventsCorrelation->Draw("same P");

  string saveBoth = "correlationBoth"+fExtension;
  cBoth->SaveAs(saveBoth.c_str());
}



void TCorrelation::GetEventsCorrelation(vector<double>& alpha, vector<double>& Np) const
{
  Np.resize(fNstep); alpha.resize(fNstep);
  for(unsigned int i = 0; i < fNstep;i++) fEventsCorrelation->GetPoint(i,alpha[i],Np[i]);
}



void TCorrelation::GetCoverageCorrelation(vector<double>& alpha,vector<double>& Np,vector<double>& errLow,vector<double>& errHigh) const
{
  Np.resize(fNstep); alpha.resize(fNstep); errLow.resize(fNstep); errHigh.resize(fNstep);
  for(unsigned int i = 0; i < fNstep; i++)
    {
      fCoverageCorrelation->GetPoint(i,alpha[i],Np[i]);
      errLow[i] = fCoverageCorrelation->GetErrorYlow(i);
      errHigh[i] = fCoverageCorrelation->GetErrorYhigh(i);
    }
}



vector<double> Cl2Ctheta(int lMax, const vector<double>& cl, const vector<double>& thetaInRad)
{
  unsigned int size = thetaInRad.size();
  vector< vector<double> > poln;
  poln.resize(lMax+1);
  for(long i = 0; i < lMax+1; i++) poln[i].resize(size);
  
  for(unsigned int i = 0; i < size; i++) poln[0][i] = 1.;
  
  if(lMax >= 1) for(unsigned int i = 0; i < size; i++) poln[1][i] = cos(thetaInRad[i]);
  if(lMax >= 2)
    {
      for(long k = 2; k < lMax+1; k++)
        for(unsigned int i = 0; i < size; i++) 
          poln[k][i] = (1./k)*((2.*k-1.)*poln[k-1][i]*cos(thetaInRad[i])-(k-1)*poln[k-2][i]);
    }
  vector<double> ct(size);
  for(unsigned int i = 0; i < size; i++)
    {
      ct[i] = 0;
      for(int l = 0; l <= lMax; l++)
        {
          double norm = (2*l+1)/(4.*M_PI);
          ct[i] += norm*cl[l]*poln[l][i];
        }
    }
  return ct;
}



double GetCtheta(int lMax, const vector<double>& cl, double thetaInRad)
{
  vector<double> poln;
  poln.resize(lMax+1);
  
  if(lMax >= 1) poln[1] = cos(thetaInRad);
  if(lMax >= 2) for(long k = 2; k < lMax+1; k++) poln[k] = (1./k)*((2.*k-1.)*poln[k-1]*cos(thetaInRad)-(k-1)*poln[k-2]);
  double ct = 0;
  double norm;
  for( int l = 0 ;l <= lMax; l++ ) {norm = (2*l+1)/(4.*M_PI); ct += norm*cl[l]*poln[l];}
  return ct;  
}



void ComputeCtheta(int lMax, TH1F** h, const vector<double>& Cl)
{
  int nBins = (*h)->GetNbinsX();
  vector<double> theta(nBins);
  for(int i = 0; i < nBins; i++) theta[i] = (*h)->GetBinCenter(i+1)*DTOR;
  vector<double> Ctheta = Cl2Ctheta(lMax,Cl,theta);
  for(int i = 0; i < nBins; i++) (*h)->SetBinContent(i+1,Ctheta[i]*sin(theta[i]));
}
