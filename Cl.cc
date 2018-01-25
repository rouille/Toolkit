#include <iostream>
#include <cmath>

#include "Cl.h"
#include "harmotools.h"

// ROOT
#include "TMath.h"

using namespace std;

vector<double> compute_Cl(const vector<TEvent>& events, const THealpixMap& incovmap, const THealpixMap& inevtmap, unsigned int lmax, const vector<unsigned int>& lbins, vector<double>& ErrorCl)
{
  THealpixMap covmap = incovmap;
  THealpixMap evtmap = inevtmap;
  // WARNING
  if(lmax >= 3*covmap.NSide())
    {
      cout << "Program Failed : lmax should not be greater than 3*nside-1" << endl;
      exit(0);
    }
  for(unsigned int i = 0; i < lbins.size(); i++)
    {
      if(lbins[i] > lmax+1)
        {
          cout << "Program Failed : lbins should not be greater than lmax+1" << endl;
          exit(0);
        }
    }
  
  // Normalization of the coverage map to the total number of events
  covmap = covmap*(events.size()*1./covmap.Total());

  // First order moment of the window
  double f1 = covmap.Total()/covmap.NPix();
  
  // Upgrade the coverage map. Window field that measures the relative exposure in the direction n in the sky.
  THealpixMap window = covmap.Map2Map(8*covmap.NSide());
  
  // We expand the window field on the spherical harmonics
  vector<unsigned int> wl;
  vector<double> wCl;
  for(unsigned int i=0; i<lmax*2+1; i++) wl.push_back(i);
  wCl = window.Map2Cl(lmax*2);
  
  /*
    Experimentally, we only have access to the power spectrum of the product of:
    Window field * stochastic field (measuring the distance to isotropy in the case of a uniform and full sky coverage)
    Since the harmonic transform of the product of two fields is the convolution of the harmonic transforms of these 
    fields, the mode-coupled power spectrum is related to the true underlying one through a convolution kernel.
    Computation of the coupling matrix.
  */
  double detM;
  TMatrixD Mll(lmax+1,lmax+1);
  Mll = compute_Mll(wl, wCl, lmax);

  TMatrixD invertMll = Mll;
  invertMll.Invert(&detM);
  cout << "Determinant of the Mll matrix = " << detM << endl;

  // Get the pixel window function of the coverage map.
  vector<double> PixWindow = GetPixWindow(covmap.NSide(),lmax);

  // Binning of the Power Spectrum in l.
  TMatrixD matP;
  TMatrixD matQ;
  makePQ(lbins, lmax, matP, matQ);

  // Pixellisation, filtering and binning -> Matrix K
  TMatrixD mat(lmax+1,lmax+1);
  unsigned int lbs = lbins.size()-1;
  TMatrixD Kll(lbs,lbs);
  TMatrixD Temp(lmax+1,lbs);
  for(unsigned int i = 0; i < lmax+1; i++)
    for(unsigned int j = 0; j < lmax+1; j++) mat(j,i) = Mll(j,i)*(PixWindow[j]*PixWindow[j]);
  
  Temp.Mult(mat,matP);
  Kll.Mult(matQ,Temp);

  // Inversion of the Kll matrix.
  TMatrixD invertKll = Kll;
  double detK;
  invertKll.Invert(&detK);
  cout << "Determinant of the Kll matrix = " << detK << endl;
  
  // Variance of the angular power spectrum estimate.
  TMatrixD CovClEstimatedInit(lmax+1,lmax+1);
  TMatrixD CovClEstimated(lbs,lbs);
  TMatrixD tmatP(lbs,lmax+1);
  TMatrixD tinvertKll(lbs,lbs);
  for(unsigned int i = 0; i < lmax+1; i++) for(unsigned int j = 0; j < lbs; j++) tmatP[j][i] = matP[i][j];
  for(unsigned int i = 0; i < lbs; i++) for(unsigned int j = 0; j < lbs; j++) tinvertKll[j][i] = invertKll[i][j];

  CovClEstimatedInit = covariance_Cltilde(wl, wCl, lmax, f1, events.size());
  Temp.Mult(CovClEstimatedInit,matP);
  CovClEstimated.Mult(tmatP,Temp);

  // Variance on the experimental power spectrum.
  TMatrixD CovCl(lbs,lbs);
  Temp.ResizeTo(lbs,lbs);
  Temp.Mult(CovClEstimated,invertKll);
  CovCl.Mult(tinvertKll,Temp);

  // RMS of the experimental power spectrum.
  ErrorCl.resize(lbs);
  for(unsigned int i = 0; i < lbs; i++) ErrorCl[i] = sqrt(CovCl[i][i]);

  // Bias induced by the finite number of arrival directions available.
  TMatrixD Bias = noise_bias(matP, f1, events.size(), lmax, lbs);

  // Expected coverage.
  THealpixMap expcovmap(covmap.NSide());
  expcovmap = covmap/covmap.Total()*events.size()*covmap.NPix()*1./(4.*M_PI);

  // Normalization of the events map.
  evtmap = evtmap*covmap.NPix()*1./(4.*M_PI);

  // Excess map.
  THealpixMap excessmap(covmap.NSide());
  excessmap = (evtmap-expcovmap)/(events.size()/(4.*M_PI*f1));

  // We expand the excessmap on the spherical harmonics.
  vector<unsigned int> expl(lmax+1);
  vector<double> expCltemp;
  Temp.ResizeTo(1,lmax+1);
  TMatrixD expCl(1,lbs);

  for(unsigned int i = 0; i < lmax+1; i++) expl[i] = i;
  expCltemp = excessmap.Map2Cl(lmax);
  for(unsigned int i = 0; i < lmax+1; i++) Temp[0][i] = expCltemp[i];
  expCl.Mult(Temp,matP);
  
  // True Cl.
  TMatrixD Cl(1,lbs);
  vector<double> Cltrue(lbs);
  Temp.ResizeTo(1,lbs);
  for(unsigned int i = 0; i < lbs; i++) Temp[0][i] = expCl[0][i]-Bias[0][i];
  Cl.Mult(Temp,invertKll);
  for(unsigned int i = 0; i < lbs; i++) Cltrue[i]=Cl[0][i];

  return Cltrue;
}



TMatrixD compute_Mll(vector<unsigned int>& wl, vector<double>& wCl, unsigned int lmax)
{
  vector<double> c(wl.size());
  for(unsigned int i = 0; i < wl.size(); i++) c[i] = (2.*wl[i]+1.)*wCl[i];
  double b;
  
  vector<double> lngammaW;
  lngammaW = wigner_init(wl.size());
  vector<double> Wsquare;
  
  TMatrixD Mll(lmax+1,lmax+1);
  for(unsigned int l1 = 0; l1 < lmax+1; l1++)
    {
      for(unsigned int l2 = 0; l2 < lmax+1; l2++)
        {
          Wsquare = wigner_3j(l1, l2, wl, lngammaW);
          b = c[0]*Wsquare[0];
          for(unsigned int i = 1; i < Wsquare.size(); i++) b += c[i]*Wsquare[i];
          // line - column.
          Mll(l2,l1) = (2.*l2+1)*b/(4.*M_PI);
        }
    }
  return Mll;
}



vector<double> wigner_init(unsigned int lmax)
{
  vector<double> lngammaW(2+lmax*3);
  for(unsigned int i = 0; i < lngammaW.size(); i++) lngammaW[i] = TMath::LnGamma(2.*i+1.)-2.*TMath::LnGamma(i+1.);
  return lngammaW;
}



vector<double> wigner_3j(unsigned int l1, unsigned int l2, const vector<unsigned int>& l3, const vector<double>& lngammaW)
{
  unsigned int size = l3.size();
  vector<unsigned int> L(size);
  vector<unsigned int> L_2(size);
  for(unsigned int i = 0; i < size; i++) {L[i] = l1+l2+l3[i];L_2[i] = L[i]/2;}

  // Triangular inequality.
  unsigned int lmin, lmax;
  lmin = abs((int)(l1-l2));
  lmax = l1+l2;
  
  double lngammaW1, lngammaW2, lngammaW3, lngammaWl, lngammaWsquare;
  vector<double> Wsquare(size,0.);
  for(unsigned int i = 0; i < size; i++)
    {
      // Parity and triangular inequality conditions.
      if(L_2[i]*2-L[i] == 0 && l3[i] >= lmin && l3[i] <= lmax)
        {
          lngammaW1 = lngammaW[L_2[i]-l1];
          lngammaW2 = lngammaW[L_2[i]-l2];
          lngammaW3 = lngammaW[L_2[i]-l3[i]];
          lngammaWl = lngammaW[L_2[i]];
          lngammaWsquare = -log(L[i]+1.)-lngammaWl+lngammaW1+lngammaW2+lngammaW3;
          Wsquare[i] = exp(lngammaWsquare);
        }
    }
  return Wsquare;
}



vector<double> GetPixWindow(unsigned int nside, unsigned int lmax)
{
  // Get the pixel window function of the coverage map.
  string healpixdir(getenv("HEALPIX_DIR"));
  char PixWindowFunction[100];
  if(nside > 64) sprintf(PixWindowFunction, "%s/data/pixel_window_n0%d.fits", healpixdir.c_str(), nside);
  else sprintf(PixWindowFunction, "%s/data/pixel_window_n00%d.fits", healpixdir.c_str(), nside);
 
  vector<double> PixWindow;

  char comment[81];
  fitsfile *fptr;
  int status = 0, hdutyp, anynul;
  long _Nele;
  double * _donnees;
 
  ffopen(&fptr,PixWindowFunction, 0, &status);
  ffmahd(fptr, 2, &hdutyp, &status);
  ffgky(fptr, TLONG, (char*)"NAXIS2", &_Nele, comment, &status);
  _donnees = new double [_Nele];
  // Read first column for Temperature window function
  ffgcvd(fptr, 1, 1 , 1, _Nele,-9999, _donnees, &anynul, &status);  
  ffclos(fptr, &status);
  PixWindow.resize(lmax+1);
  for (unsigned int j = 0;j < lmax+1; j++) PixWindow[j]=_donnees[j];

  delete [] _donnees;
  return(PixWindow);
}



void makePQ(const vector<unsigned int>& lbins, unsigned int lmax, TMatrixD& matP, TMatrixD& matQ)
{
  unsigned int nbins = lbins.size()-1;
  matP.ResizeTo(lmax+1,nbins);
  matQ.ResizeTo(nbins,lmax+1);
  
  for(unsigned int i = 0; i < nbins; i++)
    {
      for(unsigned int j = max(0,(int)lbins[i]); j < lbins[i+1]; j++)
        {
          matP(j,i) = 1./(double)(lbins[i+1]-lbins[i]);
          matQ(i,j) = 1.;
        }
    }
}



vector<vector<double> > lvalues(const vector<unsigned int>& lbins)
{
  vector<vector<double> > lvalues(2);
  lvalues[0].resize(lbins.size()-1);
  lvalues[1].resize(lbins.size()-1);
  vector<unsigned int> dec;
  dec.resize(lbins.size());
   
  for(unsigned int i = 0; i < lbins.size(); i++)
    {
      if(i == 0) dec[i] = lbins[i]-lbins[lbins.size()-1]-1;
      else dec[i] = lbins[i]-lbins[i-1]-1;
    }

  for(unsigned int j = 0; j < lvalues[0].size(); j++)
    {
      lvalues[1][j] = (dec[j+1]+1)*1./2;
      lvalues[0][j] = lbins[j]-0.5+lvalues[1][j];
    }
  return lvalues;
}



TMatrixD covariance_Cltilde(const vector<unsigned int>& wl, const vector<double>& wCl, unsigned int lmax, double f1, unsigned int nevt)
{
  unsigned int size=wl.size();
  vector<double> c(size);
  for(unsigned int i = 0; i < size; i++) c[i] = (2.*wl[i]+1.)*wCl[i];
  double b;
  
  vector<double> lngammaW = wigner_init(wl.size());
  vector<double> Wsquare;
  
  TMatrixD Mll(lmax+1,lmax+1);
  for(unsigned int l1 = 0; l1 < lmax+1; l1++)
    {
      for(unsigned int l2 = 0; l2 < lmax+1; l2++)
        {
          Wsquare = wigner_3j(l1, l2, wl, lngammaW);
          b = c[0]*Wsquare[0];
          for(unsigned int i = 1; i < Wsquare.size(); i++) b += c[i]*Wsquare[i];
          // line - column
          Mll(l2,l1) = pow(4.*M_PI*f1/nevt,2)*(1./(2.*M_PI))*b;
        }
    }
  return Mll;
}



TMatrixD noise_bias(TMatrixD matP, double f1, unsigned int nevt, unsigned int lmax, unsigned int nbin)
{
  TMatrixD Bias(1,nbin);
  TMatrixD Temp(1,lmax+1);
  
  for(unsigned int i = 0; i < lmax+1; i++) Temp[0][i] = (4.*M_PI/nevt)*pow(f1,2);
  Bias.Mult(Temp,matP);

  return Bias;
}


