#include <cmath>
#include <iomanip>
#include <fstream>

#include "angdist.h"
#include "harmotools.h"
#include "common.h"
#include "userfcn.h"

// ROOT
#include "TCanvas.h"
#include "TGraph.h"
#include "TProfile2D.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"
#include "TLatex.h"

static char gObjName[1024];
static int gObjNumber = 0;
static char * GetObjName() {sprintf(gObjName,"Phi%d",gObjNumber++); return gObjName;}


using namespace std;

// About the TAngularDistribution class
TAngularDistribution::TAngularDistribution() : TFitFunction()
{
  fExtension = ""; // no figure is the default
  fPhiModulation = false;
}


void TAngularDistribution::SetOptions(double xMin, double xMax, unsigned int nBins, double chi2)
{
  fDataMin = xMin;
  fDataMax = xMax;
  fNBins = nBins;
  fChi2Limit = chi2;
}



double TAngularDistribution::GetAccAngle(double* angle) const
{
  double val;
  if( *angle <= fDataMin || *angle >= fDataMax ) val = 0;
  else val = fFitFcn->EvalPar(angle,fPFitParameters);
  return val;
}



// About the TPhiModulation class
TPhiModulation::TPhiModulation() : TAngularDistribution(),
				   fPhiLawFunc(0x0),
				   fModPhi(0x0),
				   fPhiLaw(0x0),
				   fModParams(0x0)
{
}


TPhiModulation::~TPhiModulation()
{
  if( fModPhi ) delete fModPhi;
  if( fPhiLaw ) delete fPhiLaw;
  if( fModParams ) delete [] fModParams;
}



TF1 * TPhiModulation::SetPhiLaw(double (*f)(double *t, double *par), double thetamin, double thetamax, int npars)
{
  fPhiLawFunc = f;
  fThetaMin = thetamin;
  fThetaMax = thetamax;
  fNParsPhiLaw = npars;
  fPhiLaw = new TF1(GetObjName(),fPhiLawFunc,fThetaMin,fThetaMax,fNParsPhiLaw);
  return fPhiLaw;
}



void TPhiModulation::ComputePhiModulation(const vector<TEvent>& events, unsigned int nthetabins)
{
  if( fModPhi != 0x0 ) delete fModPhi;
  if( fPhiLaw == 0x0 )
    {
      cout << "Set the law to be used to fit the phi modulation by calling TPhiModulation::SetPhiLaw." << endl;
      cout << "Returning." << endl;
      return;
    }
  fPhiModulation = true;  
  TCanvas * cPhi = new TCanvas(GetObjName(), "Phi Modulation", 1);
  fModPhi = new TH1F(GetObjName(), "Phi modulation", nthetabins, fThetaMin, fThetaMax);
  DrawHisto(cPhi, fModPhi, "Zenith Angle", "#phi Modulation Coefficient", "");  
  vector<double> thetavals(nthetabins+1);
  unsigned int nbinsphi = 9;
  SetNBins(nbinsphi);
  vector< vector<double> > pars;
  vector< vector<double> > parserr;
  vector<double> vnull(fNPars,0);
  for(unsigned int i = 0; i <= nthetabins; i++) thetavals[i] = fThetaMin+(fThetaMax-fThetaMin)*i/nthetabins;
  string extSave = fExtension;
  
  //  fExtension = ""; // do not draw all distributions
  for(unsigned int i = 0; i < nthetabins; i++)
    {
      vector<double> vphi;
      // fill data
      for(unsigned int j = 0; j < events.size(); j++)
        {
          if( events[j].fTheta < thetavals[i+1] && events[j].fTheta >= thetavals[i] ) vphi.push_back(events[j].fPhi);
        }
      if( vphi.size() > 0 )
        {
          SetData(vphi);
          // fit phi distribution in this zenith angle bin
          Run();
          pars.push_back(GetFitParameters());
          parserr.push_back(GetFitParametersErrors());
          fModPhi->SetBinContent(i+1,pars[i][1]);
          fModPhi->SetBinError(i+1,parserr[i][1]);
        }
      else
        {
          pars.push_back(vnull);
          parserr.push_back(vnull);
          fModPhi->SetBinContent(i+1,pars[i][1]);
          fModPhi->SetBinError(i+1,parserr[i][1]);
        }
    }
  fExtension = extSave;
  fModPhi->Fit(fPhiLaw,"Q");
  if( fModParams ) delete [] fModParams;
  else fModParams = new double[fNParsPhiLaw];
  fPhiLaw->GetParameters(fModParams);
  for(unsigned int i = 0; i < 2; i++) cout << "phi law : " << fModParams[i] << endl;
  if( fExtension != "" )
    {
      fModPhi->Draw("E");
      fPhiLaw->SetLineWidth(2);
      fPhiLaw->SetLineColor(kRed);
      fPhiLaw->Draw("same");
      cPhi->Update();
      string fileName = "phiLaw"+fExtension;
      cPhi->SaveAs(fileName.c_str());
    }
}



double TPhiModulation::GetAccAngle(double* angle, double * pars)
{
  if( !fPhiLaw || !fFitFcn )
    {
      cout << "Cannot compute phi modulation, previous fit must have failed" << endl;
      cout << "PhiLaw = " << fPhiLaw << " " << " FitFcn = " << fFitFcn << endl;
      cout << "Exiting." << endl;
      exit(0);
    }
  if( *angle < fDataMin || *angle > fDataMax ) return 0;
  if( pars[0] < fThetaMin || pars[0] > fThetaMax ) return 0;
  // fPhiLaw describes fFitFcn parameters law as a function of theta
  double amplitude = fPhiLaw->EvalPar(pars,fModParams);
  double thepars[2];
  thepars[0] = 1;
  thepars[1] = amplitude;
  double modulation = fFitFcn->EvalPar(angle,thepars);
  return modulation;
}
