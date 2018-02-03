#include <cmath>
#include <iostream>

#include "STClibrary.h"

#include "rayleigh.h"
#include "common.h"
#include "angdist.h"
#include "harmotools.h"
#include "userfcn.h"

#include "TCanvas.h"

static char gObjName[1024];
static int gObjNumber = 0;
static char * GetObjName() {sprintf(gObjName,"Rayleigh%d",gObjNumber++);return gObjName;}


TRayleigh::TRayleigh(const vector<TEvent>& events, const THealpixMap& covMap, const unsigned int nHarmonic, double dRA)
  : fEvts(events), fCovMap(covMap), fNharmonic(nHarmonic), fdRA(dRA)
{
  fNbins = (unsigned int)floor(360./fdRA);
  fNevts = events.size();
  fHistoEvts = 0x0;
  fPlotCov = 0x0;
  fHistoTrue = 0x0;
  fFitFcnTrue = 0x0;
  fComponent = false;
  Init();
}



TRayleigh::~TRayleigh()
{
  if( fHistoEvts  ) delete fHistoEvts;
  if( fPlotCov    ) delete fPlotCov;
  if( fFitFcnTrue ) delete fFitFcnTrue;
  if( fHistoTrue  ) delete fHistoTrue;
}



void TRayleigh::Init()
{
  // fRA
  DECLARE_POINTER(double,raEvts,fEvts,fRa);
  fRA.resize(fNevts);
  for(unsigned int i = 0; i < fNevts; i++) fRA[i] = raEvts[i];
  delete [] raEvts;

  // fHistoEvts
  fHistoEvts = new TH1F(GetObjName(), "RA Events", fNbins, 0, 360);

  // fBinCenter
  fBinCenter.resize(fNbins);
  for(unsigned int i = 0; i < fNbins; i++) fBinCenter[i] = fHistoEvts->GetBinCenter(i+1);
  
  // fHistoTrue
  fHistoTrue = new TH1F(GetObjName(), "RA True", fNbins, 0, 360);

  // fRAweight - TIME_FLAT case
  fRAweight.resize(fNevts);
  for(unsigned int i = 0; i < fNevts; i++) fRAweight[i] = 1;

  // fRAweightBins - TIME_FLAT case
  fRAweightBins.resize(fNbins);
  for(unsigned int i = 0; i < fNbins; i++) fRAweightBins[i] = 1;
}



void TRayleigh::ComputeRAweight()
{
  fRAweight.clear();
  fRAweightBins.clear();
  // We integer all the declinations corresponding to one RA.
  fCovMap = fCovMap.Map2Map(512);
  unsigned int nPix = fCovMap.NPix();

  vector<long> iPix;
  iPix.resize(nPix);
  vector<double> lPix, bPix, raPix, decPix, pixVal;
  raPix.resize(nPix);
  decPix.resize(nPix);
  pixVal.resize(nPix);

  for(unsigned int i = 0; i < nPix; i++) iPix[i] = i;
  fCovMap.GiveLB(iPix, lPix, bPix);
  pixVal = fCovMap.Values(lPix, bPix); 
  // Warning : gal2radec converts galactic coordinates in equatorial coordinates(WITH RIGHT ASCENSION IN HOURS)
  for(unsigned int i = 0; i < nPix; i++) gal2radec(lPix[i], bPix[i], &raPix[i], &decPix[i]);

  // We want to know the bin corresponding to the raPix
  unsigned int index = 0;
  vector<double> raw;
  fRAweightBins.resize(fNbins);
  raw.resize(fNbins);
  for(unsigned int i = 0; i < nPix; i++)
    {
      index = (unsigned int)floor(fNbins*raPix[i]*15./360);
      raw[index] += 1;
      fRAweightBins[index] += pixVal[i]*cos(decPix[i]*DTOR);
    }
  for(unsigned int i = 0; i < fNbins; i++) fRAweightBins[i] = fRAweightBins[i]/raw[i];

  // Interpolation at the RA of the events given by fRA.
  fRAweight.resize(fNevts);
  fRAweight = linear_interp(fBinCenter, fRAweightBins, fRA);

  // We need to normalize fRAweight to the mean of fRAweight to have fRAweight equal to 1 when the coverage 
  // does not depend on RA
  double raWeightMean = 0.;
  for(unsigned int i = 0; i < fNevts; i++) raWeightMean += fRAweight[i]/fNevts;
  for(unsigned int i = 0; i < fNevts; i++) fRAweight[i] = fRAweight[i]/raWeightMean;
}



void TRayleigh::ComputeDist()
{
  // Evts
  for(unsigned int i = 0; i < fNevts; i++) fHistoEvts->Fill(fRA[i]);
  vector<double> tmp(fNbins);
  for(unsigned int i = 0; i < fNbins; i++) 
    {
      tmp[i] = fHistoEvts->GetBinContent(i+1);
      fHistoEvts->SetBinError(i+1,sqrt(fHistoEvts->GetBinContent(i+1)));
    }

  // Normalization Evts - Cov
  fNormEvts = integrate_nc5(fBinCenter,tmp);
  fNormCov = integrate_nc5(fBinCenter,fRAweightBins);
  
  // Cov
  for(unsigned int i = 0; i < fNbins; i++) fRAweightBins[i] = fRAweightBins[i]*fNormEvts/fNormCov;
  
  // True
  fRAtrue.resize(fNbins);
  for(unsigned int i = 0; i < fNbins; i++) fRAtrue[i] = tmp[i]/fRAweightBins[i];
  double errorTrue;
  for(unsigned int i = 0; i < fNbins; i++)
    {
      fHistoTrue->SetBinContent(i+1,fRAtrue[i]);
      errorTrue = fHistoEvts->GetBinError(i+1)/fRAweightBins[i];
      fHistoTrue->SetBinError(i+1,errorTrue);
    }
}



void TRayleigh::FitDistTrue()
{
  unsigned int nPar = 2;
  double chi2Fit;
  fFitFcnTrue = new TF1(GetObjName(), fitFunction, 0., 360., nPar);
  fFitFcnTrue->SetParNames("Amplitude","Phase");
  fHistoTrue->Fit(fFitFcnTrue, "R");
  chi2Fit = (fFitFcnTrue->GetChisquare())/(fFitFcnTrue->GetNDF());
  cout << "Chi2/NDF = " << chi2Fit << endl;	  

  fFitParameter.resize(nPar);
  double* getParameter = new double[nPar];
  fFitFcnTrue->GetParameters(getParameter);
  fFitParameter[0] = abs(getParameter[0]);
  fFitParameter[1] = getParameter[1]*RTOD;
  fFitParameter[1] = mod(720.+fFitParameter[1],360);
}



void TRayleigh::DrawEvtsDist() const
{
  string Xaxis = "Right Ascension";
  string Yaxis = "Raw Event Count";
  string nameEvts = "RA Evts";
  string save = "rayleighRAevt.png";

  TCanvas *cRAevt =  new TCanvas(GetObjName(), nameEvts.c_str(), 700, 700);

  DrawHisto(cRAevt, fHistoEvts, Xaxis, Yaxis, save);
}



void TRayleigh::DrawCovDist()
{
  double* binCenterTmp  = new double[fNbins];
  double* raWeightBinsTmp = new double[fNbins];
  for(unsigned int i = 0; i < fNbins; i++)
    {
      binCenterTmp[i] = fBinCenter[i];
      raWeightBinsTmp[i] = fRAweightBins[i];
    }
  string Xaxis = "Right Ascension";
  string Yaxis = "";
  string nameCov = "RA Coverage";

  TCanvas *cRAcov = new TCanvas(GetObjName(), nameCov.c_str(), 700, 700);
  fPlotCov = new TGraphErrors(fNbins, binCenterTmp, raWeightBinsTmp);

  PlotXY(cRAcov, fPlotCov, 0., 360., nameCov, Xaxis, Yaxis);
  fPlotCov->Draw("AP");
  cRAcov->Update();
  cRAcov->SaveAs("rayleighRAcov.png");

  delete [] binCenterTmp;
  delete [] raWeightBinsTmp;
}



void TRayleigh::DrawTrue()
{
  string Xaxis = "Right Ascension";
  string Yaxis = "Normalized Event Count";
  string nameTrue = "RA True";
  string save = "";

  TCanvas *cRAtrue = new TCanvas(GetObjName(), nameTrue.c_str(), 700, 700);

  DrawHisto(cRAtrue, fHistoTrue, Xaxis, Yaxis, save);
  FitDistTrue(); // Fit
  fHistoTrue->Draw("e1p");
  fFitFcnTrue->SetLineWidth(2);
  fFitFcnTrue->SetLineColor(kRed);
  fFitFcnTrue->Draw("same");
  cRAtrue->Update();
  cRAtrue->SaveAs("rayleighRAtrue.png");
}



void TRayleigh::ComputeComponent()
{
  // fCompFirst & fCompSec
  double compFirstTmp = 0;
  double compSecTmp = 0;
  for(unsigned int i=0; i<fNevts; i++)
    {
      compFirstTmp += (1./fRAweight[i])*cos(fNharmonic*fRA[i]*DTOR);
      compSecTmp += (1./fRAweight[i])*sin(fNharmonic*fRA[i]*DTOR);
    }
  fCompFirst = (2./fNevts)*compFirstTmp;
  fCompSec = (2./fNevts)*compSecTmp;
}



void TRayleigh::ComputeAmplitude()
{
  if(! fComponent) ComputeComponent();
  fComponent = true;
  fAmplitude = sqrt(fCompFirst*fCompFirst+fCompSec*fCompSec);
}



void TRayleigh::ComputePhase()
{
  if(! fComponent) ComputeComponent();
  fComponent = true;
  fPhase = atan2(fCompSec,fCompFirst)*RTOD; // [-180,180]
  fPhase = mod(fPhase+720,360); // [0,360]
}



void TRayleigh::ComputeSignificance()
{
  fSignificance = (fNevts*fAmplitude*fAmplitude)/4.;
}



void TRayleigh::ComputeChanceProbability()
{
  double significance = (fNevts*fAmplitude*fAmplitude)/4.;
  fChanceProbability = exp(-1.*significance);
}



double fitFunction(double* x, double* par)
{
  static const double factor = 2.*M_PI/360.;
  double t = x[0];
  double F = 1+par[0]*cos(t*factor+par[1]);
  return F;
}
