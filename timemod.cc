#include "timemod.h"
#include "userfcn.h"

// ROOT
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

#ifdef gcc323
char* operator+( std::streampos&, char* );
#endif

static char gObjName[1024];
static int gObjNumber=0;
static char * GetObjName() {sprintf(gObjName,"Time%d",gObjNumber++); return gObjName;}



TTimeModulation::TTimeModulation() {fAccTimeModel = "TIME_FLAT";}



void TTimeModulation::ComputeAccTime(const vector<double> & UTCh, const vector<double> & UTCs, double timeStep, double chi2lim)
{
	/* Read an acceptance file */
	if(fAccTimeModel != "TIME_FLAT" && fAccTimeModel != "UTC" && fAccTimeModel != "UTC+JD")
    {
		  cout << "Acceptance file : " << fAccTimeModel << endl;
		  fUTCs.clear();
		  fAccTime.clear();
		  ifstream AccFile(fAccTimeModel.c_str());
		  double UTCsTmp, AccTimeTmp;
		  if( !AccFile.is_open() )
		    {
			    cerr << "Program Failed : Error : Cannot access to acceptance file" << endl;
			    exit(-1);
		    }
		  while( AccFile >> UTCsTmp )
		    {
			    AccFile >> AccTimeTmp;
			    fUTCs.push_back(UTCsTmp);
			    fAccTime.push_back(AccTimeTmp);
		    }
		  cout << "Acceptance file loaded" << endl ;
		  cout << fUTCs.size() << " line were read" << endl;
      AccFile.close();
		  return;
    }
	
  /* From the distribution */
	// UTCs time
	double UTCsMin = *min_element(UTCs.begin(),UTCs.end());
	double UTCsMax = *max_element(UTCs.begin(),UTCs.end());
	unsigned int nElement = (int)floor((UTCsMax-UTCsMin)/timeStep);
	fUTCs.clear();
	fUTCs.resize(nElement);
	for(unsigned int i = 0; i < nElement; i++) fUTCs[i] = UTCsMin+(UTCsMax-UTCsMin)*i*1./(nElement-1);
	
	// Julian Days
	vector<double> JD(UTCs.size());
	double JDmin, JDmax;
	for(unsigned int j = 0; j < UTCs.size(); j++) utcs2jd(UTCs[j], &JD[j]);
	JDmax = *max_element(JD.begin(), JD.end());  
	JDmin = *min_element(JD.begin(), JD.end());
	
	
	 /**********************************\
	 | Plot the UTCh & JD distributions |
	 \**********************************/
	
	// UTCh
	TCanvas * cUTCh = new TCanvas(GetObjName(), "UTC Hour", 1100, 700);
	unsigned int nBinsUTCh = 48; 
	TH1F * hUTCh = new TH1F(GetObjName(),"UTC Hour Distribution", nBinsUTCh, 0, 24.);
	GetHisto(hUTCh,UTCh);
	hUTCh->SetMinimum(0.);
	DrawHisto(cUTCh,hUTCh,"UTC Hour","","UTCh.pdf");
	
  // Julian Days
  TCanvas * cJD = new TCanvas(GetObjName(), "JD", 1100, 700);
  unsigned int nBinsJD = (unsigned int)(JDmax-JDmin);
  TH1F * hJD = new TH1F(GetObjName(),"Julian Days Distribution", nBinsJD, JDmin, JDmax);
  hJD->SetMinimum(0.);
  GetHisto(hJD,JD);
  DrawHisto(cJD,hJD,"Julian Days","","JulianDays.pdf");
	
	// Interpolation of the UTCh histogram
	vector<double> accUTChBin(nBinsUTCh);
	vector<double> UTChBin(nBinsUTCh);
	for(unsigned int i = 0; i < nBinsUTCh; i++)
	  {
		  UTChBin[i] = hUTCh->GetBinCenter(i+1);
		  accUTChBin[i] = hUTCh->GetBinContent(i+1);
	  }	  
	
  // Fill the acceptance vector with UTCh variation
	vector<double> accUTCh(fUTCs.size());
	double UTChTmp;
	int yearTmp, monthTmp, dayTmp;
	for(unsigned int k = 0; k < fUTCs.size(); k++)
	  {
		  utcs2date(fUTCs[k], &yearTmp, &monthTmp, &dayTmp, &UTChTmp);
		  accUTCh[k] = linear_interp(UTChBin, accUTChBin, UTChTmp);
	  }
	
	fAccTime.resize(fUTCs.size());
	if( fAccTimeModel == "UTC" )
	  {
		  cout << "Acceptance computed from the UTC distribution" << endl;
		  for(unsigned int j = 0; j < fUTCs.size(); j++) fAccTime[j] = accUTCh[j];
	  }
	else if(fAccTimeModel == "UTC+JD")
	  {
		  cout << "Acceptance computed from the JD & UTC distributions" << endl;
		  vector<double> accJDbin(nBinsJD);
		  vector<double> JDbin(nBinsJD);
		  for(unsigned int i = 0; i < nBinsJD; i++)
		    {
			    JDbin[i] = hJD->GetBinCenter(i+1);
			    accJDbin[i] = hJD->GetBinContent(i+1);
		    }
		  // Fill the acceptance vector with JD+UTC variation
		  double accJDtmp, JDtmp;
		  for(unsigned int k = 0; k < fUTCs.size(); k++)
		    {
			    accJDtmp = 0.;
			    utcs2jd(fUTCs[k], &JDtmp);
			    accJDtmp = linear_interp(JDbin, accJDbin, JDtmp);
			    fAccTime[k] = accJDtmp*accUTCh[k];
		    }
	  }
}

