#include <iostream>
#include <fstream>
#include <cmath>

#include "lima.h"
#include "common.h"
#include "harmotools.h"
#include "maptools.h"
#include "projmap.h"

#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"

static TProgressBar gProgB;
static char gObjName[1024];
static int gObjNumber = 0 ;
static char * GetObjName() {sprintf(gObjName,"LiMa%d",gObjNumber++);return gObjName;}


using namespace std;



TLiMa::TLiMa(const TCoverage & coverage, const vector<TEvent>& events, string lobeFile, double threshold)
{
	fEvents = events;
	fNevents = events.size();
	fThreshold = threshold;
	fExtension = ".png";
  InitLobe(lobeFile);
	InitMap(coverage);
	fWindowRadius = integrate_nc5(fThetaLobe,fLobe);
	ComputeLobeCorrection();	
}



TLiMa::TLiMa(const TCoverage & coverage, const vector<TEvent> & events, const vector<double>& thetalobe, const vector<double>& lobe, double threshold)
{
	fEvents = events;
	fNevents = events.size();
	fThreshold = threshold;
	fExtension = ".png";
	fThetaLobe = thetalobe;
	fLobe = lobe;
	InitMap(coverage);
	fWindowRadius = integrate_nc5(fThetaLobe,fLobe);
	ComputeLobeCorrection();
}



void TLiMa::InitLobe(string lobe_file)
{
	ifstream lobeFile(lobe_file.c_str());
	while( !lobeFile.eof() )
    {
      float thetaLobe, lobe;
      lobeFile >> thetaLobe >> lobe;
      fThetaLobe.push_back(thetaLobe);
      fLobe.push_back(lobe);
    }
	fThetaLobe.pop_back();
	fLobe.pop_back();
	lobeFile.close();
}



void TLiMa::InitMap(const TCoverage & coverage)
{
  // Coverage map
  fRawCovMap = coverage.GetMap();
  fRawCovMap *= fNevents/fRawCovMap.Total();
  fCovMap = fRawCovMap.IntBeam(fThetaLobe, fLobe);
  
  // Prepare LiMa map
  fLiMaMap.SetNSide(fRawCovMap.NSide());
  
  // Events Map
  DECLARE_VECTOR(double,lEvents,fEvents,fL);
  DECLARE_VECTOR(double,bEvents,fEvents,fB);
  fEventsMap = map_events(coverage.NSide(),lEvents, bEvents);

  fEventsMap = fEventsMap.IntBeam(fThetaLobe, fLobe);

  // Assign 0 to pixels that were at zero in the initial coverage map
  for(unsigned int i = 0; i < fCovMap.NPix(); i++) if(fRawCovMap[i] == 0.) { fCovMap[i] = 0.; fEventsMap[i] = 0; }

  // Difference map
  fDiffMap = fEventsMap-fCovMap;
}



void TLiMa::ComputeLobeCorrection()
{
	vector<double> weight2(fLobe.size());
	vector<double> weight1(fLobe.size());
	weight2.resize(fLobe.size()), weight1.resize(fLobe.size());
	
	for( unsigned int i = 0; i < fLobe.size(); i++ )
    {
      weight2[i] = fLobe[i]*fLobe[i]*sin(fThetaLobe[i]*DTOR);
      weight1[i] = fLobe[i]*sin(fThetaLobe[i]*DTOR);
    }
	
	double meanWeight, varWeight;
	meanWeight = 2.*M_PI*integrate_nc5(fThetaLobe,weight1);
	varWeight = 2.*M_PI*integrate_nc5(fThetaLobe,weight2);
	fLobeCorrection = sqrt(varWeight/meanWeight);
}



double TLiMa::LiMa(double nOn, double nOff, double alpha)
{
	double excess = nOn-nOff*alpha;
	double n1, n2, sign, a;
	if(alpha == 0) return sqrt(1.0*nOn);
	if(excess > 0)
    {
      n1 = nOn;
      n2 = nOff;
      sign = 1.0;
      a = alpha;
    }
	else
    {
      n2 = nOn;
      n1 = nOff;
      sign = -1.0;
      a = 1.0/alpha;
    }
	if(n2 == 0) return sign*sqrt(2*n1*log((1+a)/a));
	if(n1 == 0) return sign*sqrt(2*n2*log(1+a));
	double nTot = n1+n2;
	double t1 = n1*log(((1+a)/a)*(n1/nTot));
	double t2 = n2*log((1+a)*(n2/nTot));
	double significance = sqrt(2.)*sqrt(t1+t2);  
	
	return sign*significance;
}



void TLiMa::ComputeLiMaMap()
{
	cout << "Computing Li & Ma map" << endl;
		
	double nOff = fNevents*1.;
  for(unsigned int i = 0; i < fLiMaMap.NPix(); i++)
    {
		if(fCovMap[i] != 0)
      {
        double alpha = max(fCovMap[i]*1./nOff,0.);
        double nOn = max(fEventsMap[i],0.);
        double significance = LiMa(nOn,nOff,alpha);
        if(isfinite(significance) == 0) fLiMaMap[i] = 0;
        else fLiMaMap[i] = significance/fLobeCorrection;
      }
		else fLiMaMap[i] = 0.;
    }
}



TH1F * TLiMa::GetLiMaHistogram(int nbins, double minLiMa, double maxLiMa)
{
	TH1F * histoLiMa = new TH1F(GetObjName(),"Li & Ma", nbins, minLiMa, maxLiMa);
	if( fLiMaMap.NPix() == 0 )
    {
      cout << "The map is empty. Call TLiMa::ComputeLiMaMap first. Returning." << endl;
      return histoLiMa;
    }	
	for(unsigned int i = 0; i < fLiMaMap.size(); i++) if(fLiMaMap[i] != 0) histoLiMa->Fill(fLiMaMap[i]);		
	histoLiMa->SetMinimum(0.1);
	histoLiMa->Sumw2();
	return histoLiMa;
}


void TLiMa::DrawLiMaHistogram(TH1F * histoLiMa)
{
	if( fExtension == "" )
	{
		cout << "DrawLiMaHistogram:: LiMa extension filename not defined. Returning." << endl;	
		return;
	}
	string name = "Li & Ma Distribution";
	TCanvas * cLiMa = new TCanvas(GetObjName(), name.c_str(), 700, 700);  
	string Xaxis = "#sigma";
	string Yaxis = "";
	string saveTmp = "";
	DrawHisto(cLiMa,histoLiMa,Xaxis.c_str(),Yaxis.c_str(),saveTmp.c_str());
	cLiMa->SetLogy();
	double norm = histoLiMa->Integral("width");
	histoLiMa->Draw("e1p");
	TF1 * normalLaw = new TF1("NormalLaw", "([0]/sqrt(2.*3.1415926))*exp(-0.5*x*x)",histoLiMa->GetBinCenter(1),histoLiMa->GetBinCenter(histoLiMa->GetNbinsX()));
	normalLaw->SetParameter(0,norm);
	normalLaw->SetLineWidth(1);
	normalLaw->SetLineStyle(2);
	normalLaw->SetLineColor(kRed);
	normalLaw->Draw("same");
	
	cLiMa->Update();
	string save = "histoLima"+fExtension;
	cLiMa->SaveAs(save.c_str());
}


void TLiMa::ComputeMaxima()
{
	if( fLiMaMap.NPix() == 0 )
    {
      cout << "The map does not contain pixel. Call TLiMa::ComputeLiMaMap size first. Returning." << endl;
      return;
    }
	vector<double> maxima; 
	vector<long> ipMax;
	maxima = fLiMaMap.FindMaxima(ipMax);
	for(unsigned int i = 0; i < ipMax.size(); i++)
    {
		if(maxima[i] >= fThreshold)
      {
        fIpAboveThreshold.push_back(ipMax[i]);
        fValueAboveThreshold.push_back(maxima[i]);
      }
    }
}



void TLiMa::GetMaxima(vector<long> & ipAboveThreshold, vector<double> & valueAboveThreshold)
{
	ipAboveThreshold = fIpAboveThreshold;
	valueAboveThreshold = fValueAboveThreshold;
}



void TLiMa::PrintResults(bool detailed)
{
	if( fIpAboveThreshold.size() )
    {
      vector<double> lMax, bMax, raMax, decMax;
      fLiMaMap.GiveLB(fIpAboveThreshold, lMax, bMax);
      raMax.resize(lMax.size());
      decMax.resize(lMax.size());
      for(unsigned int i = 0; i < lMax.size(); i++)
        {
          gal2radec(lMax[i], bMax[i], &raMax[i], &decMax[i]);
          // Warning : gal2radec converts G coordinates in Q coordinates (WITH RIGHT ASCENSION IN HOURS)
          raMax[i] = raMax[i]*15.;
        }
		
      // Get events dataset informations
      DECLARE_VECTOR(double, utcTime, fEvents, fUTCs);
      DECLARE_VECTOR(double, theta, fEvents, fTheta);
      DECLARE_VECTOR(double, energy, fEvents, fEnergy);
      int yearMin, yearMax, monthMin, monthMax, dayMin, dayMax;
      double utcsMin, utcsMax, utchourMin, utchourMax, thetaMin, thetaMax, eMin, eMax;
      utcsMin = *min_element(utcTime.begin(),utcTime.end());
      utcsMax = *max_element(utcTime.begin(),utcTime.end());
      utcs2date(utcsMin, & yearMin, & monthMin, & dayMin, & utchourMin);
      utcs2date(utcsMax, & yearMax, & monthMax, & dayMax, & utchourMax);
      thetaMin = *min_element(theta.begin(),theta.end());
      thetaMax = *max_element(theta.begin(),theta.end());
      eMin = *min_element(energy.begin(), energy.end());
      eMax = *max_element(energy.begin(), energy.end());
      
      cout << "---------------------------------------------------" << endl;
      cout << "Blind Hotspots Search :" << endl;
      cout << "---------------------------------------------------" << endl;
      cout << "Significance threshold used : " << fThreshold << endl;
      cout << "Number of events in dataset : " << fNevents << endl;
      cout << "First event : " << dayMin << "/" << monthMin << "/" << yearMin << " UTC = " << utchourMin << endl;
      cout << "Last event : " << dayMax << "/" << monthMax << "/" << yearMax << " UTC = " << utchourMax << endl;
      cout << "Min Energy : " << eMin << " Max Energy : " << eMax << endl;
      cout << "Min Theta : " << thetaMin << " Max Theta : " << thetaMax << endl << endl;
      if( detailed == true )
        {
          for(unsigned int i = 0; i < fIpAboveThreshold.size(); i++)
            {
              vector<TEvent> eventAboveThreshold;
              eventAboveThreshold = GetEventsFromLB(fEvents, lMax[i], bMax[i], fWindowRadius);
              printf("Maxima %2d : l = %5.2f b = %5.2f significance = %5.2f \n", i+1,lMax[i],bMax[i],fLiMaMap[fIpAboveThreshold[i]]);
              printf("nobs = %8.2f, nexp = %8.2f\n",fEventsMap[fIpAboveThreshold[i]],fCovMap[fIpAboveThreshold[i]]);
              printf("list of events in a %5.3f deg. circular window \n",fWindowRadius);
              printf("%12s | %4s | %5s | %3s | %5s | %6s | %5s | %6s | %6s | %5s \n", "EventsID", "year", "month", "day", "UTC", "Energy", "Theta", "Phi", "l", "b");
              for(unsigned int j = 0; j < eventAboveThreshold.size(); j++)
                {
                  int year, month, day;
                  double utc;
                  utcs2date(eventAboveThreshold[j].fUTCs, &year, &month, &day, &utc);
                  printf("%12llu | %4d | %5d | %3d | %5.2f | %6.2f | %5.1f | %6.1f | %6.1f | %6.1f \n", eventAboveThreshold[j].fId, year, month, day, utc, eventAboveThreshold[j].fEnergy, eventAboveThreshold[j].fTheta, eventAboveThreshold[j].fPhi, eventAboveThreshold[j].fL, eventAboveThreshold[j].fB);
                }
              cout << endl;
            }
        }
		else
      {
        cout << "Num |    Ip  |    l   |     b  |   ra   |  dec   | events |   Bg   |  LiMa " << endl;
			  for(unsigned int i = 0; i < fIpAboveThreshold.size(); i++) printf("%3d |%7d |%7.2f |%7.2f |%7.2f |%7.2f |%7.1f |%7.1f |%7.2f\n", i+1, (int)fIpAboveThreshold[i], lMax[i], bMax[i], raMax[i], decMax[i], fEventsMap[fIpAboveThreshold[i]], fCovMap[fIpAboveThreshold[i]], fValueAboveThreshold[i]);
      }
    }
  else cout << "No candidates found." << endl;
	cout << "---------------------------------------------------" << endl;
}




double LiMa(double nOn, double nOff, double alpha)
{
	double excess = nOn-nOff*alpha;
	double n1, n2, sign, a;
	if(alpha == 0) return sqrt(1.0*nOn);
	if(excess > 0)
    {
      n1 = nOn;
      n2 = nOff;
      sign = 1.0;
      a = alpha;
    }
	else
    {
      n2 = nOn;
      n1 = nOff;
      sign = -1.0;
      a = 1.0/alpha;
    }
	if(n2 == 0) return sign*sqrt(2*n1*log((1+a)/a));
	if(n1 == 0) return sign*sqrt(2*n2*log(1+a));
	double nTot = n1+n2;
	double t1 = n1*log(((1+a)/a)*(n1/nTot));
	double t2 = n2*log((1+a)*(n2/nTot));
	double significance = sqrt(2.)*sqrt(t1+t2);  
	
	return sign*significance;
}

