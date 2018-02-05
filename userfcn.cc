#include "userfcn.h"
#include <fstream>

static char gObjName[1024];
static int gObjNumber = 0;
static char * GetObjName() {sprintf(gObjName,"UserFcn%d",gObjNumber++); return gObjName;}


// About the TParameter class
TParameter::TParameter() {fStatus = eFree; fValue = 0;}



// About the TFitFunction class
TFitFunction::TFitFunction()
{
	fNBins = 10;
	fDegreeMax = 10;
	fFitFcn = 0x0;
	fPFitParameters = 0x0;
}



TFitFunction::~TFitFunction()
{
	if( fFitFcn ) delete fFitFcn;
	if( fPFitParameters ) delete [] fPFitParameters;
}



void TFitFunction::ComputeFitFunction()
{
	fDataBins.clear();
	fDataFit.clear();
	if( fFitFcn == 0x0 )
    {
      cout << "TFitFunction::ComputeFitFunction: Impossible to compute fit function. Pointer not allocated." << endl;
      return;
    }
	unsigned int mybins=1000;
	fDataBins.resize(mybins,0.);
	fDataFit.resize(mybins,0.);
	double step = (fDataMax-fDataMin)/(mybins-1);
	for(unsigned int i=0; i<mybins; i++)
    {
      fDataBins[i] = fDataMin+i*step;
      fDataFit[i] = fFitFcn->EvalPar(&fDataBins[i],fPFitParameters);
    }
	string fn = fDataName+"distribution.txt";
	ofstream outangdist(fn.c_str());
	for( unsigned int i=0;i<mybins;i++ ) outangdist << fDataBins[i] << " " << fDataFit[i] << endl;
	outangdist.close();
}



void TFitFunction::PlotDistribution() const
{
	if( fExtension == "" ) return;
	unsigned int size = fDataBins.size();
	double* th  = new double[size];
	double* pth = new double[size];
	for(unsigned int i=0; i<size; i++)
    {
      th[i] = fDataBins[i];
      pth[i] = fDataFit[i];
    }
	string Xaxis = "#"+fDataName;
	string Yaxis = "";
	string save = fDataName+"Fit"+fExtension;
	TCanvas* cAngle = new TCanvas(GetObjName(), fDataName.c_str(),700,700);
	TGraphErrors* angle = new TGraphErrors(size, th, pth);
	PlotXY(cAngle, angle, fDataMin, fDataMax, "", Xaxis, Yaxis);
	angle->SetLineWidth(2);
	angle->Draw("AL");
	cAngle->Update();
	cAngle->SaveAs(save.c_str());
	delete [] th;
	delete [] pth;
}



bool TFitFunction::Run()
{
	string filename= "";
	string xAxis = "#"+fDataName;
	string name = "";
	TH1F* Distribution = new TH1F(GetObjName(), name.c_str(), fNBins, fDataMin, fDataMax);
	TCanvas * cDistribution = 0x0;  
	GetHisto(Distribution,fData);
	if( fExtension != "" )
    {
      cDistribution = new TCanvas(GetObjName(), "Fit", 700, 700);
      filename = fDataName+"RawDistribution"+fExtension;
      DrawHisto(cDistribution,Distribution,xAxis,"",filename);		
    }
	// Now let's fit the histogram !
	unsigned int degree = 0;
	unsigned int nPar = 0;
	fChi2perNDF = fChi2Limit+1; // be sure that chi2perNDF > fChi2Limit
	bool status = false;
	cout << "Chi2/NDF max : " << fChi2Limit << endl;
	int addpar = (fDegreeMax == 0) ? 0 : 1 ;
	while( fChi2perNDF >= fChi2Limit )
    {
      nPar = degree+fNPars+addpar;
      fFitFcn = new TF1(GetObjName(),fAngFitFunc,fDataMin,fDataMax,nPar);
      for( unsigned int i = 0; i < fNPars;i++ ) // restore values
        {
          if( fParameters[i].fStatus == eFixed ) fFitFcn->FixParameter(i,fParameters[i].fValue);
          else fFitFcn->SetParameter(i,fParameters[i].fValue);
        }
      fFitFcn->FixParameter(fNPars,degree);
      {
        double getParameter[nPar];
        fFitFcn->GetParameters(getParameter);
      }
      if( fFitFcn->GetNumberFreeParameters() == 0 )
        {
          cout << "No free parameters ! Try increasing degree..." << endl;
          delete fFitFcn;
          fFitFcn = 0x0;
          degree++;
          if( ++degree > fDegreeMax )
            {
              status = false;
              break;
            }
          continue;
        }
      string opts = "Q";
      if( fExtension == "" ) opts += "N";
      Distribution->Fit(fFitFcn,opts.c_str());
      fChi2 = fFitFcn->GetChisquare();
      fChi2perNDF = (fFitFcn->GetChisquare())/(fFitFcn->GetNDF());
      cout << "Chi2/NDF = " << fChi2perNDF << endl;
      if( fChi2perNDF >= fChi2Limit )
        {
          delete fFitFcn;
          fFitFcn = 0x0;
          if( ++degree > fDegreeMax )
            {
              status = false;
              break;
            }
        }
      else status = true;
    }
	if( !status )
    {      
      cout << "Program Failed : Too much splines/poly used." << endl;
      cout << "Fit failed. Returning." << endl;
      return status;
    }
	if( fExtension != "" )
    {
      //Plot Histogram + Fit
      Distribution->Draw("e1p");
      fFitFcn->SetLineWidth(2);
      fFitFcn->SetLineColor(kRed);
      fFitFcn->Draw("same");
      cDistribution->Update();
      string filename = fDataName+"DistributionFitted"+fExtension;
      cDistribution->SaveAs(filename.c_str());
    }
	// Fitted parameters (to be returned)
	double getParameter[nPar];
	fFitFcn->GetParameters(getParameter);
	double error[nPar];
	for( unsigned int i=0;i<nPar;i++ ) error[i] = fFitFcn->GetParError(i);
	fFitParameters.resize(nPar);
	fPFitParameters = new double[nPar];
	fFitParametersErrors.resize(nPar);
	for( unsigned int i=0;i<nPar;i++ )
    {
      fFitParameters[i] = getParameter[i];
      fPFitParameters[i] = getParameter[i];
      fFitParametersErrors[i] = error[i];
      cout << "parameter " << i << " : " << fFitParameters[i] << " error : " << fFitParametersErrors[i] << endl;
    }
	ComputeFitFunction();
	if( fExtension != "" ) PlotDistribution();
	
	return status;
}



vector<double> spline_init(const vector<double>& x, const vector<double>& y)
{
	unsigned int size = x.size();
	vector<double> u(size-1), y2(size);
	double sig, p;
	// The so-called natural cubic spline has zero second derivative 
	// on one or both of its boudaries
	y2[0] = u[0] = 0.;
	// This is the decomposition loop of the tridiagonal algorithm
	for( unsigned int i=1;i<size-1;i++ )
    {
      sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p = sig*y2[i-1]+2.;
      y2[i] = (sig-1.)/p;
      u[i] = (y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i] = (6.*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
	y2[size-1] = 0.;
	// Back substitution loop of the tridiagonal algorithm
	for( unsigned int i=size-2;i>0;i-- ) y2[i] = y2[i]*y2[i+1]+u[i];
	
	return y2;
}



vector<double> spline_interp(const vector<double>& xa, const vector<double>& ya, const vector<double>& y2a, const vector<double>& x)
{
	vector<double> y(x.size(),0.);
	unsigned int k, klow, khigh;
	float h, b, a;
	klow = 0;
	khigh = xa.size()-1;
	for( unsigned int i=0;i<x.size();i++ )
    {
      while( khigh-klow>1 )
        {
          k = (khigh+klow) >> 1;
          if(xa[k]>x[i]) khigh = k;
          else klow = k;
        }
      h = xa[khigh]-xa[klow];
      if( h == 0. )
        {
          cout << "Bad xa input to routine spline_interp" << endl;
          return y;
        }
      a = (xa[khigh]-x[i])/h;
      b = (x[i]-xa[klow])/h;
      y[i] = a*ya[klow]+b*ya[khigh]+((a*a*a-a)*y2a[klow]+(b*b*b-b)*y2a[khigh])*(h*h)/6.0;
    }
	
	return y;
}



vector<double> linear_interp(const vector<double>& x, const vector<double>& y, const vector<double>& u)
{
	unsigned int size = u.size();
	unsigned int xsize = x.size();
	vector<double> v(size);
	unsigned int k, klow, khigh;
	for( unsigned int i=0;i<size;i++ )
    {
      klow = 0;
      khigh = xsize-1;
      while(khigh-klow>1)
        {
          k = (khigh+klow)/2.;
          if( x[k]>u[i] ) khigh = k;
          else klow = k;
        }
      v[i] = y[klow]+((y[klow]-y[khigh])/(x[klow]-x[khigh]))*(u[i]-x[klow]);
    }
	
	return v;
}



double linear_interp(const vector<double>& x, const vector<double>& y, double u)
{
	vector<double> v(1);
	vector<double> uu(1);
	uu[0] = u;
	v = linear_interp(x,y,uu);
	
	return v[0];
}



double splineFunction(double* t, double* par)
{
	// par[0] : tMin
	// par[1] : tMax
	// par[2] : number of spline
	// par[3+i] : coeff. deg i
	static double xmin = par[0];
	static double xmax = par[1];
	static unsigned int nbspl = (unsigned int)par[2]+100; // be sure that at first call, nbspl != newnbspl
	double nxmin = par[0];
	double nxmax = par[1];
	unsigned int newnbspl = (unsigned int)par[2];
	if( newnbspl == 1 ) return 1;
	static vector<double> xspl(newnbspl);
	bool splinit = false;
	if( newnbspl != nbspl || nxmin != xmin || nxmax != xmax )
    {
      xmin = nxmin;
      xmax = nxmax;
      nbspl = newnbspl;
      splinit = true;
      xspl.resize(nbspl);
      for( unsigned int i = 0;i < nbspl;i++ ) xspl[i] = xmin+i*(xmax-xmin)/(nbspl-1.);
    }
	// We create a base of spline function
	vector<double> fcsplines(nbspl);
	static vector< vector<double> > vy2;
	static vector< vector<double> > vyspl;
	if( splinit )
    {
      vy2.resize(nbspl);
      vyspl.resize(nbspl);
      for( unsigned int i = 0;i < nbspl;i++ )
        {
          vector<double> yspl(nbspl,0);
          yspl[i] = 1.;
          vyspl[i] = yspl;
          vy2[i] = spline_init(xspl,yspl);
        }
    }
	vector<double> tt(1,t[0]);
	for( unsigned int i = 0;i < nbspl;i++ )
    {
      vector<double> tmp = spline_interp(xspl,vyspl[i],vy2[i],tt);
      fcsplines[i] = tmp[0];
    }
	// Linear combination of the spline function of the base
	double value = 0;
	for( unsigned int i = 0;i<nbspl;i++ ) value += par[3+i]*fcsplines[i];
	
	return value;
}



double polyFunction(double* t, double* par)
{
	// t[0] : theta
	// par[0] : tMin
	// par[1] : tMax
	// par[2] : polynom degree
	// par[3+i] : coeff. deg i
	unsigned int degree = (unsigned int)par[2];
	double value = 0.;
	for(unsigned int i = 0; i < degree; i++) value = value+par[3+i]*pow(t[0],i*1.);
	
	return value;
}



double geoFunction(double* t, double *par)
{
	// t[0] = theta
	return cos(t[0]*DTOR)*sin(t[0]*DTOR);
}



double fdFunction(double* t, double* par)
{
  // t[0] : theta
  // par[0] : FD cutoff
  // par[1] : FD width 
  double arg = (t[0]-par[0])/par[1];
  
  return 1./(1.+exp(arg));
}



double fdsplFunction(double* t, double *par)
{
	double val;
	if(*t < par[2] || *t > par[3])
    {
      TF1::RejectPoint();
      val = 0;
    }
	else val = fdFunction(t,par)*geoFunction(t)*splineFunction(t,&par[2]);
	
	return val;
}



double fdpolyFunction(double* t, double* par)
{
	double val;
	if(*t < par[2] || *t > par[3]) val = 0;  
	else val = fdFunction(t,par)*geoFunction(t)*polyFunction(t,&par[2]);
	
	return val;
}



double geosplFunction(double* t, double* par)
{
	double val;
	if(*t < par[2] || *t > par[3]) val = 0;  
	else val = geoFunction(t)*splineFunction(t,&par[2]);
	
	return val;
}



double geopolyFunction(double* t, double* par)
{
	double val;
	if(*t < par[2] || *t > par[3]) val = 0;  
	else val = geoFunction(t)*polyFunction(t,&par[2]);
	return val;
}



double scintillatorFunction(double *t, double *par)
{
	// polysfd works very well
	// par[0] : FD cutoff
	// par[1] : FD width
	// par[2] : thetamin
	// par[3] : thetamax
	double val;
	if(*t < par[2] || *t > par[3]) val = 0;
	else val = fdsplFunction(t,par)*cos(*t*DTOR);
	
	return val;
}



double GaussFunction(double *t, double *par)
{
	double arg = (t[0]-par[2])*DTOR/par[3];
	
	return par[0]+par[1]*exp(-arg*arg/2.);
}



double AugerPhiFunction(double *t, double *par)
{
	return par[0]*(1+par[1]*cos(t[0]*DTOR*1));
}



double ModPhiThetaLaw(double *t, double *par)
{
	// force small zenith angle to 0 correction
	return par[0]*t[0]*exp(t[0]/par[1]);
}



double PoissonFluctuation(double mean, double value)
{
	double probability = 0.;
	double valueTmp = 0;
	while(valueTmp <= value)
    {
      probability += TMath::PoissonI(valueTmp, mean);
      valueTmp++;
    }
	
	return(1.0 - probability);
}



double Binomial(int NTot, int k, double p)
{
	double CkN = TMath::Binomial(NTot, k);
	
	return ( CkN * pow( p, 1.*k ) * pow( 1.-p, 1.*(NTot-k) ) );
}



double bigP(int NTot, int k, double p)
{
	double P = 0.;
	for(int j = k; j <= NTot; j++) P += Binomial(NTot, j, p);
	
	return P;
}



double OneSigmaOpeningAngle(double theta, double dTheta, double dPhi)
{
	return ( sqrt( dTheta*dTheta + sin(theta*DTOR)*sin(theta*DTOR)*dPhi*dPhi) * 1/sqrt(2.0) );
}


double integrate_nc5(const vector<double> & x, const vector<double> & y)
{
  // This is a five points Newton-Cotes (Bode's formula) integrator
  // See NumRec for details
  // We assume that the data is regularly gridded
  unsigned int npts = x.size();
  double h = x[1]-x[0];
  vector<unsigned int> ii;
  unsigned int nbii = (unsigned int)floor((npts-1.)/4);
  unsigned int rest = (npts-1)-nbii*4;  
  unsigned int nbii2;
  if(rest == 1 || rest == 2) nbii2 = nbii-1;
  else nbii2 = nbii;
  ii.resize(nbii2);
  for(unsigned int i=0; i<nbii2; i++) ii[i] = (i+1)*4;
                    
  double integral = 0;
  for(unsigned int i=0; i<nbii2; i++)
    {
      integral += 2.*h*(7.*(y[ii[i]-4]+y[ii[i]])+32.*(y[ii[i]-3]+y[ii[i]-1])+12.*y[ii[i]-2])/45.;
    }

  if( rest+1 == 2 )
    {
      // decoupage 4-3
      unsigned int shift = npts-2;
      // 4 points
      integral += 3*h*(y[shift-4]+3*y[shift-3]+3*y[shift-2]+y[shift-1])/8.; 
      // 3 points
      integral += h*(y[npts-3]+4*y[npts-2]+y[npts-1])/3.;
    }
  else if( rest+1 == 3 )
    {
      // decoupage 4-4
      unsigned int shift = npts-3;
      // 4 points
      integral += 3*h*(y[shift-4]+3*y[shift-3]+3*y[shift-2]+y[shift-1])/8.;
      // 3 points
      shift = npts;
      integral += 3*h*(y[shift-4]+3*y[shift-3]+3*y[shift-2]+y[shift-1])/8.;
    }
  else if(rest+1 == 4)
    {
      integral += 3*h*(y[npts-4]+3*y[npts-3]+3*y[npts-2]+y[npts-1])/8.;
    }
  return integral;
}

