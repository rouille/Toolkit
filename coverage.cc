#include "coverage.h"
#include "userfcn.h"

#include <algorithm>

// ROOT
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

#ifdef gcc323
char* operator+( std::streampos&, char* );
#endif

static char gObjName[1024];
static int gObjNumber = 0;
static char * GetObjName() {sprintf(gObjName,"Coverage%d",gObjNumber++);return gObjName;}

static TProgressBar gProgB;
static vector<double> gLpix;
static vector<double> gBpix;
static vector<double> gRApix;
static vector<double> gDECpix;
static vector<unsigned int> gPixelIndexesInFOV;
static unsigned int gNPixInFOV;

// Pointer on the function used for the zenithal distribution
static double (*gFitFcn)(double*,double*);

using namespace std;


double cov_ra_function(double *t, double *par)
{
	// initialization at first call
	static double thelatitude = par[1];
	static double clat = cos(thelatitude*DTOR);
	static double slat = sin(thelatitude*DTOR);
	if( par[1] != thelatitude )// new latitude (should not happen in normal conditions of use)
    {
      thelatitude = par[1];
      clat = cos(thelatitude*DTOR);
      slat = sin(thelatitude*DTOR);
    }
	double dec = par[0]*DTOR;
	double* parameter = par+3;
	double thetamin = parameter[2];
	double thetamax = parameter[3];
	double costheta = sin(dec)*slat+cos(dec)*cos((*t)*DTOR)*clat;
	double thetheta = acos(costheta)*RTOD;
	double val = 0;
	if( thetheta < thetamax && thetheta > thetamin ) val = (*gFitFcn)(&thetheta,parameter)/sin(thetheta*DTOR);
	return val;
}


TCoverage::TCoverage(unsigned int nside)
{
	fMap.SetNSide(nside);
	fCoordSystem = 'G';
	fThetaModulation = false;
	fPhiModulation = false;
	fTimeModulation = false;
	fExtension = "";
}



void TCoverage::ComputeConstants() const
{
	// pertinent only if coordsys = 'G'
	gPixelIndexesInFOV.clear();
	unsigned int npix = fMap.NPix();  
	gRApix.resize(npix);
	gDECpix.resize(npix);
	fMap.GiveLB(fMap.fIpix,gLpix,gBpix);
	for(unsigned int i = 0; i < npix; i++)
    {
      gal2radec(gLpix[i],gBpix[i],&gRApix[i],&gDECpix[i]);
      // Warning : gal2radec converts galactic coordinates in equatorial coordinates with Right Ascension in hours
      gRApix[i] = gRApix[i]*15.;
      if(IsInFOV(gDECpix[i])) gPixelIndexesInFOV.push_back(i);
    }
	gNPixInFOV = gPixelIndexesInFOV.size();
}



vector<double> TCoverage::GetZenithalRA(unsigned int nptsra) const
{
	vector<double> zen;
	if( !fTimeModulation )// no time modulation
    {
      // Get array of zenith R.A. as a function of time
      zen.resize(nptsra);
      for(unsigned int i = 0; i < nptsra; i++) zen[i] = 360.*i/nptsra; // fake utctimes in hours
    }
	else // time modulation
    {
      // Get array of zenith R.A. as a function of time     
      vector<double> utcTime = fTimeMod.GetUTCs();
      unsigned int size = utcTime.size();
      if( !size )
        {
          cout << "TCoverage::GetZenithalRA : You must call before TCoverage::ComputeAcceptance. Exiting." << endl;
          exit(0);
        }
      zen.resize(size);
      double zenithAz = 180.; // Whatever
      double zenithEl = 90.;// zenith elevation in Deg
      double lJD, lRA, lDec;
      for(unsigned int j = 0; j < size; j++)
        {
          utcs2jd(utcTime[j],&lJD);
          azel2radec(fLatitude,fLongitude,zenithAz,zenithEl,lJD,&lRA,&lDec);
          // Warning : azel2radec converts elevation and azimuth (WITH AZIMUTH MEASURED WESTWARDS FROM SOUTH) in equatorial coordinates (WITH RA IN HOURS)
          lRA *= 360/24.; // conversion hour -> Deg
          zen[j] = lRA;
        }
    }  
	return zen;
}



void TCoverage::ComputeDeclinationLimits()
{
	double thedec = fLatitude+fThetaDist.GetMaxAngle();
	fDecMax = (thedec > 90 ? 90 : thedec);
	thedec = fLatitude-fThetaDist.GetMaxAngle();
	fDecMin = (thedec < -90 ? -90 : thedec);
	ComputeConstants();
}



void TCoverage::ZenithToDeclination(unsigned int nb)
{
	fDeclination.clear();
	fDecDistribution.clear();
	fDeclination.resize(nb);
	fDecDistribution.resize(nb);
	vector<double> theproba;
	for(unsigned int i = 0; i < nb; i++) fDeclination[i] = -90.+180.*i/(nb-1.);// Deg
	int degree = (int)fThetaDist.fPFitParameters[4];
	// Fixed : FD cutoff + FD width + thetamin + thetamax + degree
	// Free : degree (free parameters related to spline or poly function)
	// New : longitude + latitude + declination
	int nbpar = degree+5+3;
	double * theparam = new double[nbpar];
	theparam[1] = fLatitude;
	theparam[2] = fLongitude;
	for(int k = 0; k < degree+5; k++) theparam[k+3] = fThetaDist.fPFitParameters[k];
	TF1 * cov_func = new TF1("cov_func",cov_ra_function,0.,360.,nbpar);
	gFitFcn = fThetaDist.fAngFitFunc;
	double precision = 1e-8;
	double cosdec, sindec, cosomega, omegamin(0), omegamax(0), otmp;
	double cthetamax = cos(fThetaDist.GetMaxAngle()*DTOR);
	cout << "ZenithToDeclination : Integrating Coverage Map : " << endl;
	gProgB.Zero();
	gProgB.fBegin = 0;
	gProgB.fEnd = nb;
	gProgB.InitPercent();
	for(unsigned int i = 0; i < nb; i++)
    {
      // Compute min and max value of the hour angle
      cosdec = cos(fDeclination[i]*DTOR);
      sindec = sin(fDeclination[i]*DTOR);
      cosomega = (cthetamax-sindec*fsinLatitude)/(cosdec*fcosLatitude);
      omegamin = 0;
      omegamax = 360.;
      if(fabs(cosomega) <= 1)
        {
          omegamax = acos(cosomega)*RTOD;
          omegamin = -omegamax;
          if(omegamin > omegamax)
            {
              otmp = omegamax;
              omegamax = omegamin;
              omegamin = otmp;
            }
        }
      theparam[0] = fDeclination[i];   
      cov_func->SetParameters(theparam);
      fDecDistribution[i] = cov_func->Integral(omegamin,omegamax,precision)*cosdec;
      gProgB.PrintPercent(i);
    }
	gProgB.EndPercent();
	cout << "Coverage Map done" << endl;
	delete [] theparam;
	delete cov_func;
}



void TCoverage::ComputeCoverage(unsigned int npts)
{
	cout << "Coverage Map Computation" << endl;
	
	vector<double> accTime = fTimeMod.GetAccTime();
	
	if( fTimeMod.GetAccTimeModel() == "TIME_FLAT" )
    {
      if( !fThetaModulation )
        {
          cout << "Take into account the ZENITHAL modulation to continue. Exiting" << endl;
          exit(0);
        }
      if( fThetaModulation && !fPhiModulation )
        {
          cout << "Acceptance modulation(s) used to compute the coverage map:" << endl;
          cout << "ZENITHAL : YES" << endl;
          cout << "AZIMUTAL : NO" << endl;
          cout << "TIME : NO" << endl;
          ZenithToDeclination(npts);
          // Linear interpolation
          double factor = integrate_nc5(fDeclination,fDecDistribution);  
          double* y = new double[npts];
          double* x = new double[npts];
          for( unsigned int i = 0; i < npts;i++ )
            {
              y[i] = fDecDistribution[i]/factor;
              x[i] = fDeclination[i];
            }
          if( fExtension != "" )
            {
              string name = "";
              string Xaxis = "#delta";
              string Yaxis = "";
              string save = "declination.png";
              TCanvas* cDeclination = new TCanvas(GetObjName(),"declination", 700, 700);
              TGraphErrors* declination = new TGraphErrors(npts, x, y);
              double xMin = *min_element(fDeclination.begin(),fDeclination.end());
              double xMax = *max_element(fDeclination.begin(),fDeclination.end());
              PlotXY(cDeclination, declination, xMin, xMax, name, Xaxis, Yaxis); 
              declination->Draw("AL");
              cDeclination->Update();
              cDeclination->SaveAs(save.c_str());
            }
			
          // Healpixmap of the exposure
          fMap = linear_interp(fDeclination,fDecDistribution,gDECpix);
          unsigned int idx;
          gProgB.Zero();
          gProgB.fBegin = 0;
          gProgB.fEnd = gNPixInFOV;
          gProgB.InitPercent();
          for(unsigned int i = 0; i < gNPixInFOV; i++)
            {
              idx = gPixelIndexesInFOV[i];
              fMap[idx] = fMap[idx]/cos(gDECpix[idx]*DTOR)/factor;
              gProgB.PrintPercent(i);
            }
          gProgB.EndPercent();
          cout << endl;
          delete [] x;
          delete [] y;
          return; // nothing more to do
        }
		
      if( fThetaModulation && fPhiModulation )
        {
          cout << "Acceptance mobulation(s) used to compute the coverage map:" << endl;
          cout << "ZENITHAL : YES" << endl;
          cout << "AZIMUTAL : YES" << endl;
          cout << "TIME : NO" << endl;
        }
    }
	else if( fTimeModulation )
    {
      if( accTime.size() == 0 )
        {
          cout << "You must call before TCoverage::ComputeAcceptance. Exiting." << endl;
          exit(0);
        }
      if( !fThetaModulation )
        {
          cout << "Take into account the ZENITHAL modulation to continue. Exiting" << endl;
          exit(0);
        }
      if( accTime.size() != 0 && fThetaModulation && !fPhiModulation )
        {
          cout << "Acceptance mobulation(s) used to compute the coverage map:" << endl;
          cout << "ZENITHAL : YES" << endl;
          cout << "AZIMUTAL : NO" << endl;
          cout << "TIME : YES" << endl;
        }
      if( accTime.size() != 0 && fThetaModulation && fPhiModulation )
        {
          cout << "Acceptance mobulation(s) used to compute the coverage map:"  << endl;
          cout << "ZENITHAL : YES" << endl;
          cout << "AZIMUTAL : YES" << endl;
          cout << "TIME : YES" << endl;
        }
    }
  vector<double> lZenithRA = GetZenithalRA();
  unsigned int nptsra = lZenithRA.size();
  cout << "TCoverage::ComputeCoverage : Loop over pixels." << endl;
  double acc_theta(1), acc_phi(1), acc_time(1);
  double thetheta, costheta, sintheta;
	double thephi, drad, cdrad, sdrad, dalpha;
	vector<double> acceptance(nptsra);
	vector<double> timevec;
	if( !fTimeModulation ) timevec = lZenithRA;
	else timevec = fTimeMod.GetUTCs();
	unsigned int idx;
	gProgB.Zero();
	gProgB.fBegin = 0;
	gProgB.fEnd = gNPixInFOV;
	gProgB.InitPercent();
	for(unsigned int i = 0; i < gNPixInFOV; i++)
    {
      idx = gPixelIndexesInFOV[i];
      drad = gDECpix[idx]*DTOR;
      cdrad = cos(drad);
      sdrad = sin(drad);
      for(unsigned int j = 0; j < nptsra; j++)
        {
          dalpha = (gRApix[idx]-lZenithRA[j])*DTOR;
          costheta = sdrad*fsinLatitude+cdrad*cos(dalpha)*fcosLatitude;
          thetheta = acos(costheta);
          sintheta = sin(thetheta);
          thephi = atan2(fcosLatitude*sdrad-fsinLatitude*cdrad*cos(dalpha),cdrad*sin(dalpha));// in [-pi,pi]
          thephi /= DTOR;
          thetheta /= DTOR;
          if( fThetaModulation ) acc_theta = fThetaDist.GetAccAngle(&thetheta)/sintheta;
          if( fPhiModulation ) acc_phi = fPhiMod.GetAccAngle(&thephi,&thetheta);
          if( fTimeModulation ) acc_time  = accTime[j];
			
          // probability products
          acceptance[j] = acc_time*acc_theta*acc_phi;
        }
      fMap[idx] = integrate_nc5(timevec,acceptance);
      gProgB.PrintPercent(i);
    }
	gProgB.EndPercent();
	cout << endl;
}



void TCoverage::CorrectForAngularModulation(string angleName)
{
	if(angleName == "theta") fThetaModulation = true;
	else if(angleName == "phi") fPhiModulation = true;
	else cout << "No angle with name : " << angleName << endl;
}



void TCoverage::CorrectForTimeModulation(string timeModel)
{
	if(timeModel != "TIME_FLAT") fTimeModulation = true;
	else
    {
      cout << "TCoverage::CorrectForTimeModulation : TIME_FLAT means that no time variation is expected" << endl;
      cout << "-> The coverage map is now computed following this assumption" << endl;
    } 
}   



THealpixMap GetAnalyticalCoverage(int nside, double thetamax, double latsite)
{
	// theta law in sin*cos
	THealpixMap anacoverage(nside,'G');
	double clat = cos(latsite*DTOR);
	double slat = sin(latsite*DTOR);
	double th, ph, ll, bb;
	double ctm = cos(thetamax*DTOR);
	double ksi, ra, dec, am, cdec, sdec;
	for(unsigned int i = 0; i < anacoverage.NPix(); i++)
    {
      pix2ang_ring(nside, i, &th, &ph);
      bb = 90.-th/DTOR;
      ll = ph/DTOR;
      gal2radec(ll, bb, &ra, &dec);
      ra = ra*DTOR;
      dec = dec*DTOR;
      cdec = cos(dec);
      sdec = sin(dec);
      ksi = (ctm-slat*sdec)/(clat*cdec);
      if( ksi > 1 ) am = 0;
      else if(ksi < -1) am = M_PI;
      else am = acos(ksi);
      anacoverage[i] = clat*cdec*sin(am)+am*slat*sdec;
    }
	return anacoverage;
}


THealpixMap ComputeCoverageScrambling(const vector<TEvent> &eventsInput, unsigned int nSide, double thetaMax, double latitude, double longitude, unsigned int nScramble, unsigned int nBins, string binning, string variable)
{
	cout << "You're using the scrambling method to compute the coverage map" << endl;
	vector< vector<TEvent> > eventsBinned;
	eventsBinned = BinEvents(eventsInput, nBins, binning, thetaMax);
	cout << variable << " arrival times variable is scrambled" << endl;
	
	vector<TEvent> eventsScrambled;
	THealpixMap coverageScrambling(nSide,'G');
	gProgB.Zero();
	gProgB.fBegin = 0;
	gProgB.fEnd = nScramble;
	gProgB.InitPercent();
	for(unsigned int i = 0; i < nScramble; i++)
    {
      eventsScrambled = DoTheScrambling(eventsBinned, variable, latitude, longitude);
      DECLARE_VECTOR(double,lMap,eventsScrambled,fL);
      DECLARE_VECTOR(double,bMap,eventsScrambled,fB);
      
      if( i == 0 ) coverageScrambling = map_events(coverageScrambling.NSide(), lMap, bMap);
      else coverageScrambling = coverageScrambling+map_events(coverageScrambling.NSide(), lMap, bMap);
      
      eventsScrambled.clear();
      lMap.clear();
      bMap.clear();
      gProgB.PrintPercent(i);
    }
	gProgB.EndPercent();
	coverageScrambling = coverageScrambling/nScramble;
	
	return coverageScrambling;
}
