#include <iostream>
#include <fstream>
#include <algorithm>

#include "STClibrary.h"
#include "TMath.h"
#include "angdist.h"
#include "events.h"
#include "harmotools.h"
#include <sys/time.h>
#include "TRandom2.h"
#include "userfcn.h"

static char gObjName[1024];
static int gObjNumber=0;
static char * GetObjName() {sprintf(gObjName,"Events%d",gObjNumber++);return gObjName;}


using namespace std;


vector<TEvent> GetEvents(string fileNameInit, string cutFile)
{
  // Check that the file exists
  FILE * evtFile = fopen(fileNameInit.c_str(),"r");  
  if( !evtFile )
    {
      cerr << "Error: Cannot access to the events file : " << fileNameInit << endl;
      vector<TEvent> bad;
      return bad;
    }
  int nlines = 0;
  int size = 2048;
  fseek(evtFile,0,SEEK_END);// go to EOF
  const long lastpos = ftell(evtFile);
  fseek(evtFile,0,SEEK_SET);// rewind
  while( ftell(evtFile) < lastpos )
    {
      char line[size];
      fgets(line,size,evtFile);
      nlines++;
    }
  fseek(evtFile,0,SEEK_SET);// rewind
  TEvent event;
  vector<double> utcTime(nlines);
  vector<TEvent> tmp(nlines);
  char ctmp[100];
  char *pline, *ptmp;
  int count = 0;
  while( ftell(evtFile) < lastpos )
    {
      char line[size];
      fgets(line,size,evtFile);
      pline = line;
      
      ptmp = ctmp;
      while( *pline != ' ' ) *ptmp++ = *pline++;
      *ptmp = '\0';
      event.fId = strtoull(ctmp,0x0,10);

      pline = SkipSpaces(pline);

      ptmp = ctmp;
      pline++;
      while( *pline != ' ' ) *ptmp++ = *pline++;
      *ptmp = '\0';
      event.fTheta  = atof(ctmp);

      pline = SkipSpaces(pline);

      ptmp = ctmp;
      pline++;
      while( *pline != ' ' ) *ptmp++ = *pline++;
      *ptmp = '\0';
      event.fdTheta  = atof(ctmp);

      pline = SkipSpaces(pline);

      ptmp = ctmp;
      pline++;
      while( *pline != ' ' ) *ptmp++ = *pline++;
      *ptmp = '\0';
      event.fPhi  = atof(ctmp);

      pline = SkipSpaces(pline);

      ptmp = ctmp;
      pline++;
      while( *pline != ' ' ) *ptmp++ = *pline++;
      *ptmp = '\0';
      event.fdPhi  = atof(ctmp);

      pline = SkipSpaces(pline);

      ptmp = ctmp;
      pline++;
      while( *pline != ' ' ) *ptmp++ = *pline++;
      *ptmp = '\0';
      event.fL  = atof(ctmp);

      pline = SkipSpaces(pline);

      ptmp = ctmp;
      pline++;
      while( *pline != ' ' ) *ptmp++ = *pline++;
      *ptmp = '\0';
      event.fB  = atof(ctmp);

      pline = SkipSpaces(pline);

      ptmp = ctmp;
      pline++;
      while( *pline != ' ' ) *ptmp++ = *pline++;
      *ptmp = '\0';
      event.fRa  = atof(ctmp);

      pline = SkipSpaces(pline);

      ptmp = ctmp;
      pline++;
      while( *pline != ' ' ) *ptmp++ = *pline++;
      *ptmp = '\0';
      event.fDec  = atof(ctmp);

      pline = SkipSpaces(pline);

      ptmp = ctmp;
      pline++;
      while( *pline != ' ' ) *ptmp++ = *pline++;
      *ptmp = '\0';
      event.fUTCs  = atof(ctmp);

      pline = SkipSpaces(pline);

      ptmp = ctmp;
      pline++;
      while( *pline != ' ' ) *ptmp++ = *pline++;
      *ptmp = '\0';
      event.fTcore  = ctmp;

      pline = SkipSpaces(pline);

      ptmp = ctmp;
      pline++;
      while( *pline != '\n' ) *ptmp++ = *pline++;
      *ptmp = '\0';
      event.fEnergy  = atof(ctmp);

      
      event.fElevation = 90-event.fTheta;
      event.fdElevation = event.fdTheta;
      event.fAzimuth = mod(90-event.fPhi+720, 360.);
      event.fdAzimuth = event.fdPhi;
      utcs2date(event.fUTCs,&event.fYear,&event.fMonth,&event.fDay,&event.fUTCh);
      utcTime[count] = event.fUTCs;
      tmp[count] = event;
      count++;
    }
  fclose(evtFile);
  cout << tmp.size() << " Events were read in file : " << fileNameInit << endl << endl;

  if( cutFile == "" ) return tmp;
  else
    {
      vector<TEvent> tmp2;
      tmp2 = KillBadEvents(tmp,cutFile);
      return tmp2;
    }
}


vector<TEvent> SelectPeriod(const vector<TEvent> &input, int yearMin, int monthMin, int dayMin, double UTChMin, int yearMax, int monthMax, int dayMax, double UTChMax)
{
	cout << "Selecting Time Period from " << dayMin << "/"  << monthMin << "/" << yearMin << " UTC = " << UTChMin 
	<< " To " 
	<< dayMax << "/" << monthMax << "/" << yearMax << " UTC = " << UTChMax << endl;
	cout << "Initially : " << input.size() << " Events" << endl;
	double UTCsMin, UTCsMax;
	double JDref, JDmax, JDmin;
	date2jd(1970, 1, 1, 0., &JDref);
	date2jd(yearMin, monthMin, dayMin, UTChMin, &JDmin);
	date2jd(yearMax, monthMax, dayMax, UTChMax, &JDmax);
	UTCsMin = (JDmin-JDref)*3600.*24.;
	UTCsMax = (JDmax-JDref)*3600.*24.;
	cout << "Min UTC : " << UTCsMin << endl;
	cout << "Max UTC : " << UTCsMax << endl;
	vector<TEvent> output;
	for(unsigned int i = 0; i < input.size(); i++) if(input[i].fUTCs >= UTCsMin && input[i].fUTCs < UTCsMax)	output.push_back(input[i]);
	cout << output.size() << " events selected" << endl << endl;  
	return output;
}


vector<TEvent> KillBadEvents(const vector<TEvent> & input, string cutFileName)
{
	vector<double> UTCsCut;
	vector<double> cut;
	double UTCsCutTmp;
	double cutTmp;
	ifstream cutFile( cutFileName.c_str() );
	if( !cutFile.is_open() )
    {
		  cerr << "Program Failed: Error: Cannot access to cut file " << cutFileName << endl;
		  exit(-1);
    }
	while( cutFile >> UTCsCutTmp )
    {
		  cutFile >> cutTmp;
		  UTCsCut.push_back(UTCsCutTmp);
		  cut.push_back(cutTmp);
    }
	cout << "Cut file loaded" << endl ;
	cout << cut.size() << " line were read" << endl;
	
	double UTCsCutMin,UTCsCutMax;
	UTCsCutMin = *min_element(UTCsCut.begin(),UTCsCut.end());
	UTCsCutMax = *max_element(UTCsCut.begin(),UTCsCut.end());

  DECLARE_VECTOR(double,eventsUTCs,input,fUTCs);  
	
	vector<double> eventsUTCsCut;
	eventsUTCsCut = linear_interp(UTCsCut,cut,eventsUTCs);
	
	cout << "Performing cuts based on file : " << cutFileName << endl;
	cout << input.size() << " events on input" << endl;
	vector<TEvent> output;
	for(unsigned int i = 0; i < input.size(); i++) if(eventsUTCsCut[i] != 0 && input[i].fUTCs >= UTCsCutMin && input[i].fUTCs <= UTCsCutMax) output.push_back(input[i]);
	cout << output.size() << " events were kept" << endl << endl;
	
	return output;
} 


vector<TEvent> ScrambleData(const vector<TEvent> & input, unsigned int nBins, string binning, double latitude, double longitude, double thetaMax, string variable)
{
	vector< vector<TEvent> > eventsBinned;
	eventsBinned = BinEvents(input, nBins, binning, thetaMax);
	vector<TEvent> output;
	output = DoTheScrambling(eventsBinned, variable, latitude, longitude);
	
	return output;
}


vector<vector<TEvent> > BinEvents(const vector<TEvent> & input, unsigned int nBins, string binning, double thetaMax)
{
	double * thetaMinBin = new double[nBins];
	double * thetaMaxBin = new double[nBins];
	vector<TEvent> eventsInBin;
	vector<vector<TEvent> > output;
	output.resize(nBins);
	if( binning == "THETA" )
    {
		  thetaMinBin[0] = 0.;
		  thetaMaxBin[nBins-1] = thetaMax;
		  for(unsigned int i = 1; i < nBins; i++)
        {
			    thetaMinBin[i] = i*thetaMax/nBins; 
			    thetaMaxBin[i-1] = i*thetaMax/nBins;
        }
    }
	else if( binning == "EVENTS" )
    {
		  vector<double> theta(input.size());
		  double nInBin = floor((double)input.size()/nBins);
		  for(unsigned int i = 0; i < input.size(); i++) theta[i] = input[i].fTheta;
		  sort(theta.begin(), theta.end());
		  thetaMinBin[0] = 0.;
		  thetaMaxBin[nBins-1] = theta[input.size()-1];
		  for(unsigned int i = 1; i < nBins; i++)
        {
			    thetaMinBin[i] = theta[i*(int)nInBin];
			    thetaMaxBin[i-1] = theta[i*(int)nInBin];
        }
    }
	else
    {
		  cout << "Program Failed : " << binning << ": binning method not adapted in TEvent::BinEvents" << endl;
		  exit(0);
    }
	
	for(unsigned int i = 0; i < nBins; i++)
    {
		  for(long j = 0; j < (long)input.size(); j++)
			  if(input[j].fTheta <= thetaMaxBin[i] && input[j].fTheta > thetaMinBin[i]) eventsInBin.push_back(input[j]);
		  output[i].resize(eventsInBin.size());
		  for(long k = 0; k < (long)eventsInBin.size(); k++) output[i][k] = eventsInBin[k];
		  eventsInBin.clear();
    }
	delete [] thetaMinBin;
	delete [] thetaMaxBin;
	
	return output;
}



vector<TEvent> DoTheScrambling(const vector<vector<TEvent> >& eventsBinned, string variable, double latitude, double longitude)
{
	unsigned int nBins = eventsBinned.size();
	unsigned int nEvents = 0;
	for(unsigned int i = 0; i < nBins; i++) nEvents += eventsBinned[i].size();
	
	// Initialize random seed to the machine clock
	gRandom->SetSeed(0);
	unsigned int seed;
	struct timeval myTimeVal;
	struct timezone myTimeZone;
	gettimeofday(&myTimeVal, &myTimeZone);
	seed = (unsigned int) (myTimeVal.tv_usec+(myTimeVal.tv_sec % 1000)*1000000);
	
	TRandom2 * randomDay = new TRandom2(seed);
	double * randomDayArray;
	int * dayIndex;
	
	TRandom2 * randomUTCh = new TRandom2(randomDay->Integer(nEvents));
	double * randomUTChArray;
	int * UTCindex;
	
	vector<TEvent> eventsInBin;
	vector<TEvent> eventsInBinSave;
	vector<TEvent> eventsScrambled;
	for(unsigned int j = 0; j < nBins; j++)
    {
      int nEventsInBin = eventsBinned[j].size();
      eventsInBin.resize(nEventsInBin); eventsInBinSave.resize(nEventsInBin);
      for(int k = 0; k < nEventsInBin; k++) {eventsInBin[k] = eventsBinned[j][k]; eventsInBinSave[k] = eventsBinned[j][k];}
		
      randomDayArray = new double[nEventsInBin];
      dayIndex = new int[nEventsInBin];
      randomDay->RndmArray(nEventsInBin, randomDayArray);
			// carefull with the prototype of Sort, either Sort(long...) or Sort(int...)
      TMath::Sort(nEventsInBin, randomDayArray, dayIndex, kFALSE);
		
      randomUTChArray = new double[nEventsInBin];
      UTCindex = new int[nEventsInBin];
      randomUTCh->RndmArray(nEventsInBin, randomUTChArray);
      TMath::Sort(nEventsInBin, randomUTChArray, UTCindex, kFALSE);
		
      delete [] randomDayArray;
      delete [] randomUTChArray;
      double jd;
      for(int k = 0; k < nEventsInBin; k++)
        {
          if( variable == "UTC+JD" ) 
            {
              eventsInBin[k].fYear  = eventsInBinSave[dayIndex[k]].fYear;
              eventsInBin[k].fMonth = eventsInBinSave[dayIndex[k]].fMonth;
              eventsInBin[k].fDay   = eventsInBinSave[dayIndex[k]].fDay;
              eventsInBin[k].fUTCh   = eventsInBinSave[UTCindex[k]].fUTCh;
            }
          else if( variable == "UTC" ) eventsInBin[k].fUTCh = eventsInBinSave[UTCindex[k]].fUTCh;
          else
            {
				      cout << "Program Failed : " << variable << " is not adapted in TEvent::DoTheScrambling" << endl;
				      exit(0);
            }
			    date2jd(eventsInBin[k].fYear, eventsInBin[k].fMonth, eventsInBin[k].fDay+eventsInBin[k].fUTCh/24., &jd);
			    // Warning : azel2gal converts elevation and azimuth (WITH AZIMUTH MEASURED WESTWARDS FROM SOUTH) in galactic coordinates
			    azel2gal(latitude,longitude,eventsInBin[k].fAzimuth+180.,eventsInBin[k].fElevation,jd,&eventsInBin[k].fL,&eventsInBin[k].fB,0,2000);
			    gal2radec(eventsInBin[k].fL,eventsInBin[k].fB,&eventsInBin[k].fRa,&eventsInBin[k].fDec);
			    // Warning : gal2radec converts galactic coordinates in equatorial coordinates (WITH RIGHT ASCENSION IN HOURS)
			    eventsInBin[k].fRa = eventsInBin[k].fRa*15.;
			    eventsScrambled.push_back(eventsInBin[k]);
        }
		  delete [] dayIndex;
		  delete [] UTCindex;
		  eventsInBin.clear();
    }
	return eventsScrambled;
}



char * SkipSpaces(char* input)
{
	bool ok = false;
	while( *input == ' ' ) {input++; ok = true;}
	if( ok ) return input-1;
	else return input;
}



void ShowLocalCoord(const vector<TEvent> & events, const vector<double> & theta, const vector<double> & pTheta)
{
	TCanvas * cLocalCoord = new TCanvas(GetObjName(),"Local Coordinates",1100, 700);
	cLocalCoord->Divide(2,1);
	
  /* Zenith Angle */
	double thetaMax = 0.;
	for(unsigned int i = 0; i < theta.size(); i++) if(pTheta[i] != 0 && theta[i] > thetaMax) thetaMax = theta[i];
	
	TH1F * hTheta = new TH1F(GetObjName(), "Zenith Angle Distribution", 60, 0, round(thetaMax));
	hTheta->SetMinimum(0.);
	for(unsigned int i = 0; i < events.size(); i++) hTheta->Fill(events[i].fTheta);
	vector<double> binCenterTheta(60), binContentTheta(60);
	for(unsigned int i = 0; i < 60; i++)
    {
		  binCenterTheta[i] = hTheta->GetBinCenter(i+1);
		  binContentTheta[i] = hTheta->GetBinContent(i+1);
    }
	double normEvents = integrate_nc5(binCenterTheta,binContentTheta);
	
	// Law given by the arguments of the function
	double normLaw = integrate_nc5(theta,pTheta);
	double * thetaPlot = new double[60]; 
	double * pThetaPlot = new double[60];
	for(unsigned int i = 0; i < 60; i++)
    {
		  thetaPlot[i] = binCenterTheta[i];
		  pThetaPlot[i] = linear_interp(theta,pTheta,binCenterTheta[i])*normEvents/normLaw;
    }
	
	TGraph * lawPlot = new TGraph(60, thetaPlot, pThetaPlot);
	lawPlot->SetLineColor(kRed); lawPlot->SetLineWidth(2);
	cLocalCoord->cd(1);
	DrawHisto(cLocalCoord,hTheta,"#theta","","");
	hTheta->Draw("e1p");
	lawPlot->Draw("same");
	cLocalCoord->Update();
	
	delete [] thetaPlot;
	delete [] pThetaPlot;
	
	/* Azimuth Angle */
	TH1F * hAzimuth = new TH1F(GetObjName(),"Azimuth Angle Distribution",90,0,360.);
	hAzimuth->SetMinimum(0.);
	for(unsigned int i = 0; i < events.size(); i++) hAzimuth->Fill(events[i].fAzimuth);
	cLocalCoord->cd(2);
	DrawHisto(cLocalCoord,hAzimuth,"#phi is measured eastwards from North","", "");
	hAzimuth->Draw("e1p");
	cLocalCoord->Update();
	cLocalCoord->SaveAs("LocalCoord.png");
}



void ShowLocalCoord(const vector<TEvent> & events)
{
	TCanvas * cLocalCoord = new TCanvas(GetObjName(),"Local Coordinates",1100, 700);
	cLocalCoord->Divide(2,1);
	
  /* Zenith Angle */	
	TH1F * hTheta = new TH1F(GetObjName(), "Zenith Angle Distribution", 60, 0, 60);
	hTheta->SetMinimum(0.);
	for(unsigned int i = 0; i < events.size(); i++) hTheta->Fill(events[i].fTheta);
	cLocalCoord->cd(1);
	DrawHisto(cLocalCoord,hTheta,"#theta","", "");
	hTheta->Draw("e1p");
	cLocalCoord->Update();
	
	
	/* Azimuth Angle */
	TH1F * hAzimuth = new TH1F(GetObjName(),"Azimuth Angle Distribution",90,0,360.);
	hAzimuth->SetMinimum(0.);
	for(unsigned int i = 0; i < events.size(); i++) hAzimuth->Fill(events[i].fAzimuth);
	cLocalCoord->cd(2);
	DrawHisto(cLocalCoord,hAzimuth,"#phi is measured eastwards from North","", "");
	hAzimuth->Draw("e1p");
	cLocalCoord->Update();
	
	cLocalCoord->SaveAs("LocalCoord.png");
}



void ShowEqCoord(const vector<TEvent> & events)
{
	TCanvas * cEqCoord = new TCanvas(GetObjName(),"Equatorial Coordinates",1100,700);
	cEqCoord->Divide(2,1);
	
  /* Right Ascension */	
	TH1F * hRa = new TH1F(GetObjName(),"Right Ascension Distribution",90, 0, 360);
	hRa->SetMinimum(0.);
	for(unsigned int i = 0; i < events.size(); i++) hRa->Fill(events[i].fRa);
	cEqCoord->cd(1);
	DrawHisto(cEqCoord,hRa,"RA","","");
	hRa->Draw("e1p");
	cEqCoord->Update();
	
	/* Declination */
	DECLARE_VECTOR(double,dec,events,fDec);
	double decMin = *min_element(dec.begin(),dec.end());
	double decMax = *max_element(dec.begin(),dec.end());
	TH1F * hDec = new TH1F(GetObjName(),"Declination Distribution", 90, decMin, decMax);
	hDec->SetMinimum(0.);
	for(unsigned int i = 0; i < events.size(); i++) hDec->Fill(events[i].fDec);
	cEqCoord->cd(2);
	DrawHisto(cEqCoord,hDec,"Dec","","");
	hDec->Draw("e1p");
	cEqCoord->Update();
	cEqCoord->SaveAs("EqCoord.png");
}



void ShowArrivalTimesCoord(const vector<TEvent> & events)
{
	TCanvas * cArrivalTime = new TCanvas(GetObjName(),"Arrival Time Coordinates",1100,700);
	cArrivalTime->Divide(2,1);
	
  /* UTCh */
	TH1F * hUTCh = new TH1F(GetObjName(),"UTC Distribution",48,0,24.);
	hUTCh->SetMinimum(0.);
	for(unsigned int i = 0; i < events.size(); i++) hUTCh->Fill(events[i].fUTCh);
	cArrivalTime->cd(1);
	DrawHisto(cArrivalTime,hUTCh,"UTC Hour","","");
	hUTCh->Draw("e1p");
	cArrivalTime->Update();
	
	/* Julian Days */
	vector<double> JD(events.size());
	for(unsigned int i = 0; i < JD.size(); i++)
    {
	  	double JDtmp;
		  date2jd(events[i].fYear,events[i].fMonth,events[i].fDay,events[i].fUTCh,&JDtmp);
		  JD[i] = JDtmp;
    }
	double JDmin = *min_element(JD.begin(),JD.end());
	double JDmax = *max_element(JD.begin(),JD.end());
	TH1F * hJD = new TH1F(GetObjName(),"Julian Days Distribution",100, JDmin, JDmax);
	hJD->SetMinimum(0.);
	for(unsigned int i = 0; i < JD.size(); i++) hJD->Fill(JD[i]);
	cArrivalTime->cd(2);
	DrawHisto(cArrivalTime,hJD,"JD","","");
	hJD->Draw("e1p");
	cArrivalTime->Update();
	cArrivalTime->SaveAs("ArrivalTime.png");
}
