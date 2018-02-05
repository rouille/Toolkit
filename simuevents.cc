#include <iostream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>

#include "simuevents.h"
#include "angdist.h"
#include "common.h"
#include "userfcn.h"

#include "TCanvas.h"
#include "TH1F.h"
#include "STClibrary.h"
#include "TGraph.h"
#include "harmotools.h"

#ifdef gcc323
char* operator+( std::streampos&, char* );
#endif


vector<TEvent> SimulateEvents(const THealpixMap &map, unsigned int nEvents, const vector<double> &thetaLaw, const vector<double> & pThetaLaw, double latSite, double lonSite, string utcFile, string jdFile, string globalFile)
{
  // Initialize random seed to the machine clock
  gRandom->SetSeed(0);
  unsigned int seed;
  struct timeval myTimeVal;
  struct timezone myTimeZone;
  gettimeofday(&myTimeVal, &myTimeZone);
  seed = (unsigned int) (myTimeVal.tv_usec+(myTimeVal.tv_sec % 1000)*1000000);
  TRandom trandom(seed);
  
  // Pixel
  unsigned int nPix = map.NPix();
  double szPix = map.GetPixSize()*DTOR;
  vector<double> pixNumber(nPix);
  for(unsigned int i = 0; i < nPix; i++) pixNumber[i] = i;

  // Simulated events
  vector<TEvent> simEvents;
  
  // Compute thetaMax
  double thetaMax = 0.;
  for(unsigned int i = 0; i < pThetaLaw.size(); i++) if(pThetaLaw[i] != 0 && thetaLaw[i] > thetaMax) thetaMax = thetaLaw[i];

  // Temporal variables
  vector<double> utc, pUtc;
  vector<double> jd, pJd;
  if( globalFile == "" ) {UtcGenerator(utcFile,utc,pUtc); JdGenerator(jdFile,jd,pJd);}
  else GlobalTimeGenerator(globalFile,jd,pJd);
  double jdMin = *min_element(jd.begin(),jd.end());
  double jdMax = *max_element(jd.begin(),jd.end());

  unsigned int nDraw = 10000;
  simEvents.reserve(nEvents);
  while(simEvents.size() < nEvents)
    {
      // Drawing pixel numbers according to map
      vector<long> pixNumberDraw(nDraw);
      vector<double> pixNumberDrawTmp = distri(pixNumber,map,nDraw);
      for(unsigned int i = 0; i < nDraw; i++) pixNumberDraw[i] = ceil(pixNumberDrawTmp[i]);

      // Get Pixel coordinates
      vector<double> lPix, bPix;
      map.GiveLB(pixNumberDraw, lPix, bPix);
      
      // Celestial coordinates of the simulated events
      vector<double> lSimEvents(nDraw), bSimEvents(nDraw), raSimEvents(nDraw), decSimEvents(nDraw);
      // All events are located in the center of the pixels. We move them randomly inside their pixel
      double lSimEventsTmp, bSimEventsTmp, raSimEventsTmp, decSimEventsTmp;
      double phiSimEventsTmp, thetaSimEventsTmp;
      long pixNumberMove;
      for(unsigned int i = 0; i < nDraw; i++)
        {
          double thetaPix = PiOver2-bPix[i]*DTOR;
          double phiPix = lPix[i]*DTOR;
          double szPixProj = 2*sin(thetaPix)*szPix;
          bool status = false;
          while( !status )
            {
              thetaSimEventsTmp = acos(cos(thetaPix)+(trandom.Rndm()-0.5)*szPixProj);
              phiSimEventsTmp = phiPix+(trandom.Rndm()-0.5)*szPix;
              ang2pix_ring(map.NSide(),thetaSimEventsTmp,phiSimEventsTmp,&pixNumberMove);
              if(pixNumberMove == pixNumberDraw[i]) status = true;
            }
          lSimEventsTmp = phiSimEventsTmp*RTOD;
          bSimEventsTmp = 90.-thetaSimEventsTmp*RTOD;
          gal2radec(lSimEventsTmp,bSimEventsTmp,&raSimEventsTmp,&decSimEventsTmp);
          lSimEvents[i] = lSimEventsTmp;
          bSimEvents[i] = bSimEventsTmp;
          raSimEvents[i] = raSimEventsTmp*15; // we want RA in degrees
          decSimEvents[i] = decSimEventsTmp;
        }
      // Note that the number of elements in lSimEvents, bSimEvents, raSimEvents and decSimEvents is not nEvents now
      // Only some satisfying events have been kept            
      // Generation of the times of the events according to the input files
      vector<double> jdSimEvents(nDraw);
      vector<double> utcSimEvents(nDraw);
      if(globalFile == "")
        {
          // Now we generate times
          // First UTC
          utcSimEvents = distri(utc,pUtc,nDraw);
          // Integer julian days (in fact it is not integer but .5)
          vector<double> jdSimEventsTmp = distri(jd,pJd,nDraw);
          for(unsigned int i = 0; i < nDraw; i++)
            {
              if(jdSimEventsTmp[i] > (jdMax-1)) jdSimEventsTmp[i] = jdMax-1;
              if(jdSimEventsTmp[i] < jdMin) jdSimEventsTmp[i] = jdMin;
              jdSimEventsTmp[i] = (long)(jdSimEventsTmp[i]+0.5)-0.5;
            }
          // Now define the true double julian day
          for(unsigned int i = 0; i < nDraw; i++) jdSimEvents[i] = jdSimEventsTmp[i]+(utcSimEvents[i]/24);
        }
      else
        {
          jdSimEvents = distri(jd, pJd, nDraw);
          for(unsigned int i = 0; i < nDraw; i++)
            {
              int yearSimEventsTmp, monthSimEventsTmp;
              double daySimEventsTmp;
              jd2date(jdSimEvents[i], &yearSimEventsTmp, &monthSimEventsTmp, &daySimEventsTmp);
              utcSimEvents[i] = (daySimEventsTmp-(int)daySimEventsTmp)*24.;
            }
        }
      vector<int> yearSimEvents(nDraw), monthSimEvents(nDraw), daySimEvents(nDraw);
      for(unsigned int i = 0; i < nDraw; i++)
        {
          double daySimEventsTmp;
          jd2date(jdSimEvents[i], &yearSimEvents[i], &monthSimEvents[i], &daySimEventsTmp);
          daySimEvents[i] = (long)daySimEventsTmp;
        }      
      // Convert to Azimuth and Elevation
      vector<double> azimutSimEvents(nDraw), elevationSimEvents(nDraw), thetaSimEvents(nDraw), phiSimEvents(nDraw);
      for(unsigned int i = 0; i < nDraw; i++)
        {
          double azimutSimEventsTmp, elevationSimEventsTmp;
          // RA are considered to be in hours. Azimut is given westwards from South
          radec2azel(latSite,lonSite,raSimEvents[i]/15.,decSimEvents[i],jdSimEvents[i],&azimutSimEventsTmp,&elevationSimEventsTmp);
          azimutSimEvents[i] = mod((180.+azimutSimEventsTmp),360);
          elevationSimEvents[i] = elevationSimEventsTmp;
          thetaSimEvents[i] = 90.-elevationSimEventsTmp;
          phiSimEvents[i] = mod((90.-azimutSimEvents[i]+720.),360);
          if(phiSimEvents[i] > 180.) phiSimEvents[i] = phiSimEvents[i]-360;
        }	
      
      /* 
         The simulated events have a distribution in theta that follows sin(theta) modulated by the anisotropy. 
         But we know that the true zenith angle distribution is different. It is given by the input vectors 
         thetaLaw and pThetaLaw which come from the coverage map estimation or from a prior knowledge of the 
         acceptance. We therefore have to perform an acceptance rejection to keep only events matching this 
         realistic distribution. 
         Set the law on which acceptance rejection is performed (one at max)
      */
      vector<double> law(thetaLaw.size());
      double fact = 0.;
      double maxVal = 0.;
      for(unsigned int i = 0; i < thetaLaw.size(); i++)
        {
          law[i] = 0;
          fact = sin(thetaLaw[i]*DTOR);
          if(fact != 0) law[i] = pThetaLaw[i]/fact;
          if(law[i] > maxVal) maxVal = law[i];
        }
      for(unsigned int i = 0; i < thetaLaw.size(); i++) law[i] = law[i]/maxVal;
      
      // Acceptance-rejection
      double jdRef;
      double lutcSimEvents;
      for(unsigned int i = 0; i < nDraw; i++)
        {
          double lawValue;
          lawValue = linear_interp(thetaLaw,law,thetaSimEvents[i]);
          if( lawValue >= gRandom->Rndm() && simEvents.size() < nEvents && thetaSimEvents[i] < thetaMax )
            {
              // this event is kept and added to the events vector
              TEvent event;
              event.fTheta = thetaSimEvents[i];
              event.fElevation = elevationSimEvents[i];
              event.fPhi = phiSimEvents[i];
              event.fAzimuth = azimutSimEvents[i];
              event.fL = lSimEvents[i];
              event.fB = bSimEvents[i];
              event.fRa = raSimEvents[i];
              event.fDec = decSimEvents[i];
              date2jd(1970, 1, 1, 0., &jdRef);
              lutcSimEvents = (jdSimEvents[i]-jdRef)*3600.*24;
              event.fUTCs = lutcSimEvents;
              event.fYear = yearSimEvents[i];
              event.fDay = daySimEvents[i];
              event.fMonth = monthSimEvents[i];
              event.fUTCh = utcSimEvents[i];
              simEvents.push_back(event);
            }
        }
    }
  return simEvents;
}



vector<TEvent> SimulateIsoEvents(unsigned int nEvents, const vector<double> & thetaLaw, const vector<double> & pThetaLaw, double latSite, double lonSite)
{
	// Initialize random seed to the machine clock
	unsigned int seed;
	struct timeval myTimeVal;
	struct timezone myTimeZone;
	gettimeofday(&myTimeVal, &myTimeZone);
	seed = (unsigned int) (myTimeVal.tv_usec+(myTimeVal.tv_sec % 1000)*1000000);

	gRandom->SetSeed(seed);
	
	// Simulated events
	vector<TEvent> simEvents(nEvents);

	vector<double> theta;
	long utcs_start = 1072915200; // 2004-01-01 00:00:00
	long utcs_stop  = 1104537600; // 2005-01-01 00:00:00
	vector<double> mtheta(thetaLaw.size());
	theta = distri(thetaLaw,pThetaLaw,nEvents);
	double refjd;
	utcs2jd(utcs_start,&refjd);
	double jd;
	for(unsigned int i = 0; i < nEvents; i++)
	{
		simEvents[i].fUTCs = double(utcs_start+long((utcs_stop-utcs_start)*gRandom->Rndm()));
		utcs2jd(simEvents[i].fUTCs,&jd);
		simEvents[i].fTheta = theta[i];
		simEvents[i].fElevation = 90.-simEvents[i].fTheta;
		simEvents[i].fAzimuth = gRandom->Rndm()*2.*M_PI*RTOD;
		simEvents[i].fPhi = mod((90.-simEvents[i].fAzimuth+720.),360);
		utcs2date(simEvents[i].fUTCs,&simEvents[i].fYear,&simEvents[i].fMonth,&simEvents[i].fDay,&simEvents[i].fUTCh);
		// Warning : azel2gal converts elevation and azimuth (WITH AZIMUTH MEASURED WESTWARDS FROM SOUTH) in galactic coordinates
		azel2gal(latSite,lonSite,simEvents[i].fAzimuth+180,simEvents[i].fElevation,jd,&simEvents[i].fL,&simEvents[i].fB,refjd);
		gal2radec(simEvents[i].fL,simEvents[i].fB,&simEvents[i].fRa,&simEvents[i].fDec);
	}
	return simEvents;
}



void UtcGenerator(string utcFile, vector<double> & utc, vector<double> & pUtc)
{
  // Constant law
  if( utcFile == "" )
    {
      unsigned int nUtc = 1000;
      utc.resize(nUtc); pUtc.resize(nUtc);
      for(unsigned int i = 0; i < nUtc; i++) {utc[i] = 24.*i*1./nUtc; pUtc[i] = 1.;}
    }
  // From utcFile
  else
    {
      struct stat fileStat;
      if(stat(utcFile.c_str(),&fileStat) == -1) 
        {
          cerr << "The file : " << utcFile << " has not been found." << endl;
          exit(0);
        }
      cout << "Reading UTC model file: " << utcFile << endl;
      ifstream utcModel(utcFile.c_str());
      double utcTmp, pUtcTmp;
      while(utcModel >> utcTmp)
        {
          utcModel >> pUtcTmp;
          utc.push_back(utcTmp);
          pUtc.push_back(pUtcTmp);
        }
      cout << "Read " << utc.size() << " elements" << endl;
      utcModel.close();
    }
}


void JdGenerator(string jdFile, vector<double> & jd, vector<double> & pJd)
{
  // Constant Law
  if( jdFile == "" )
    {
      double minJd, maxJd;
      date2jd(2004, 1, 1, 0., &minJd);
      date2jd(2005, 1, 1, 0., &maxJd);
      unsigned int nJd = (unsigned int)((maxJd-minJd)*100);
      jd.resize(nJd);
      pJd.resize(nJd);
      for(unsigned int i = 0; i < nJd; i++) {jd[i] = minJd+(maxJd-minJd)*i*1./(nJd-1); pJd[i] = 1.;}
    }
  // From jdFile
  else
    {
      struct stat fileStat;
      if(stat(jdFile.c_str(),&fileStat) == -1) 
        {
          cerr << "The file : " << jdFile << " has not been found." << endl;
          exit(0);
        }
      cout << "Reading JD model fil: " << jdFile << endl;
      ifstream jdModel(jdFile.c_str());
      double jdTmp, pJdTmp;
      while(jdModel >> jdTmp)
        {
          jdModel >> pJdTmp;
          jd.push_back(jdTmp);
          pJd.push_back(pJdTmp);
        }
      cout << "Read " << jd.size() << " elements" << endl;
      jdModel.close();
    }
}


void GlobalTimeGenerator(string globalFile, vector<double> & jd, vector<double> & pJd)
{
  struct stat fileStat;
  if(stat(globalFile.c_str(),&fileStat) == -1) 
    {
      cerr << "The file : " << globalFile << " has not been found." << endl;
      exit(0);
    }
  cout << "Reading Global model file :" << globalFile << endl;
  ifstream globalModel(globalFile.c_str());
  double utcTime, putcTime, jdTmp;
  while(globalModel >> utcTime)
    {
      globalModel >> putcTime;
      utcs2jd(utcTime, &jdTmp);
      jd.push_back(jdTmp);
      pJd.push_back(putcTime);
    }
  cout << "Read " << jd.size() << " elements" << endl;
  globalModel.close();
}


vector<double> distri(const vector<double> &x, const vector<double> &y, unsigned int nb)
{
  // This is an invert transform method

  // Get the seed from microsecond precision time. The root option for 
  // setting the seed (gRandom->SetSeed(0)) is only at the second precision...
  // Thanks guys ...
  unsigned int seed;
  struct timeval mytimeval;
  struct timezone mytimezone;
  gettimeofday(&mytimeval,&mytimezone);
  seed = (unsigned int)(mytimeval.tv_usec+(mytimeval.tv_sec % 1000)*1000000);
  // Set the seed to the current machine clock
  gRandom->SetSeed(seed);

  // vector size
  unsigned int sz = x.size();

  // characteristic function normalized to one at max
  vector<double> integ(sz);
  integ[0] = y[0];
  for(unsigned int i = 1; i < sz; i++) integ[i] = integ[i-1]+y[i];
  double maxi = *max_element(integ.begin(),integ.end());
  for(unsigned int i = 1; i < sz; i++) integ[i] = integ[i]/maxi;

  // draw uniform values
  vector<double> distinit;
  distinit.resize(nb);

  // get the random values
  for (unsigned int i = 0; i < nb; i++) distinit[i] = gRandom->Rndm();

  // interpolate to get the random numbers according to the distribution
  vector<double> xcorresp = linear_interp(integ, x, distinit);

  // return random numbers
  return xcorresp;
}



double distri(const vector<double> &x, const vector<double> &y)
{
  // This is an invert transform method

  // Get the seed from microsecond precision time. The root option for 
  // setting the seed (gRandom->SetSeed(0)) is only at the second precision...
  // Thanks guys ...
  unsigned int seed;
  struct timeval mytimeval;
  struct timezone mytimezone;
  gettimeofday(&mytimeval,&mytimezone);
  seed = (unsigned int)(mytimeval.tv_usec+(mytimeval.tv_sec % 1000)*1000000);
  // Set the seed to the current machine clock
  gRandom->SetSeed(seed);

  // vector size
  unsigned int sz = x.size();

  // characteristic function normalized to one at max
  vector<double> integ(sz);
  integ[0] = y[0];
  for(unsigned int i = 1; i < sz; i++) integ[i] = integ[i-1]+y[i];
  double maxi = *max_element(integ.begin(),integ.end());
  for(unsigned int i = 1; i < sz; i++) integ[i] = integ[i]/maxi;

  // draw a uniform value
  double distinit = gRandom->Rndm();
  
  // interpolate to get the random numbers according to the distribution
  double xcorresp = linear_interp(integ, x, distinit);

  // return random numbers
  return xcorresp;
}



vector<long> distriDiscrete(const vector<long> &x, const vector<double> &y, unsigned int nb)
{
  // This is an invert transform method

  // Get the seed from microsecond precision time. The root option for 
  // setting the seed (gRandom->SetSeed(0)) is only at the second precision...
  // Thanks guys ...
  unsigned int seed;
  struct timeval mytimeval;
  struct timezone mytimezone;
  gettimeofday(&mytimeval,&mytimezone);
  seed = (unsigned int)(mytimeval.tv_usec+(mytimeval.tv_sec % 1000)*1000000);
  // Set the seed to the current machine clock
  gRandom->SetSeed(seed);

  // vector size
  unsigned int sz = x.size();

  // characteristic function normalized to one at max
  vector<double> integ(sz);
  integ[0] = y[0];
  for(unsigned int i=1; i<sz; i++) integ[i] = integ[i-1]+y[i];
  double maxi = *max_element(integ.begin(),integ.end());
  for(unsigned int i=1;i<sz;i++) integ[i] = integ[i]/maxi;

  // draw uniform values
  vector<double> distinit;
  distinit.resize(nb);

  // get the random values
  for(unsigned int i = 0; i < nb; i++) distinit[i] = gRandom->Rndm();

  // interpolate to get the random numbers according to the distribution
  vector<long> xcorresp(nb);
  for(unsigned int i = 0; i < nb; i++)
      for(unsigned int j = 0; j < sz; j++) if(distinit[i] <= integ[j]) {xcorresp[i] = x[j]; break;}
  
  // return random numbers
  return xcorresp;
}



long distriDiscrete(const vector<long> &x, const vector<double> &y)
{
  // This is an invert transform method

  // Get the seed from microsecond precision time. The root option for 
  // setting the seed (gRandom->SetSeed(0)) is only at the second precision...
  // Thanks guys ...
  unsigned int seed;
  struct timeval mytimeval;
  struct timezone mytimezone;
  gettimeofday(&mytimeval,&mytimezone);
  seed = (unsigned int)(mytimeval.tv_usec+(mytimeval.tv_sec % 1000)*1000000);
  // Set the seed to the current machine clock
  gRandom->SetSeed(seed);

  // vector size
  unsigned int sz = x.size();

  // characteristic function normalized to one at max
  vector<double> integ(sz);
  integ[0] = y[0];
  for(unsigned int i = 1; i < sz; i++) integ[i] = integ[i-1]+y[i];
  double maxi = *max_element(integ.begin(),integ.end());
  for(unsigned int i = 1; i < sz; i++) integ[i] = integ[i]/maxi;

  // draw a uniform value
  double distinit = gRandom->Rndm();
  
  // interpolate to get the random numbers according to the distribution
  long xcorresp = 0;
  for(unsigned int i = 0; i < sz; i++) if(distinit <= integ[i]) {xcorresp = x[i]; break;}

  // return random numbers
  return xcorresp;
}
