#ifndef _LIMA_H_
#define _LIMA_H_

#include <vector>
#include <string>

#include "coverage.h"
#include "healpixmap.h"
#include "events.h"
#include "TH1.h"

/*!
  This class allows to perform a blind search of excess in any direction of the sky. For each pixel of the map, we 
  compute the number of expected events integrating the coverage map in a given window (top hat, gaussian, ...) 
  centered on this pixel. Then the signal is obtained by applying the same filtering on the event map. The 
  significance of the signal with respect to the expected background is computed following : 
  T. Li & Y. Ma, The Astrophysical Journal, \b 272:317-324, 1983. 
*/  


//! Li & Ma Analysis.
class TLiMa
{
 public:
  //! Constructor.
  TLiMa(const TCoverage & coverage, const vector<TEvent> & events, string lobeFile, double threshold);

  //! Constructor  
	TLiMa(const TCoverage & coverage, const vector<TEvent> & events, const vector<double>& thetalobe, const vector<double>& lobe, double threshold);
	
  //! Set #fExtension
  void SetExtension(string extension) {fExtension = extension;}

  //! Set norm of #fRawCovMap
  void SetRawCovMapNorm(unsigned int norm)
  {
    fNevents = norm;
    fRawCovMap *= fNevents/fRawCovMap.Total(); fCovMap = fRawCovMap.IntBeam(fThetaLobe,fLobe);
    for(unsigned int i = 0; i < fCovMap.NPix(); i++ ) if( fRawCovMap[i] == 0.) fCovMap[i] = 0.;
  }

  //! Computes the significance map.
  void ComputeLiMaMap();

  //! Computes #fIpAboveThreshold and #fValueAboveThreshold.
  void ComputeMaxima();

  //! Returns #fLiMaMap.
  THealpixMap GetLiMaMap() const {return fLiMaMap;}

  //! Returns #fCovMap.
  THealpixMap GetCovMap() const {return fCovMap;}

  //! Returns #fEventsMap.
  THealpixMap GetEventsMap() const {return fEventsMap;}

  //! Returns #fDiffMap.
  THealpixMap GetDiffMap() const {return fDiffMap;}

  //! Returns #fThetaLobe
  vector<double> GetThetaLobe() const {return fThetaLobe;}

  //! Returns #fLobe
  vector<double> GetLobe() const {return fLobe;}
  
  //! Returns #fLobeCorrection.
  double GetLobeCorrection() const {return fLobeCorrection;}
  
  //! Compute the Li & Ma histogram.
	TH1F * GetLiMaHistogram(int nbins = 40, double minLiMa = -5, double maxLiMa = 5);
	
  //! Draw the Li & Ma histogram
  void DrawLiMaHistogram(TH1F*);
	
  //! Returns #fIpAboveThreshold and #fValueAboveThreshold.
  void GetMaxima(vector<long> & ipAboveThreshold, vector<double> & valueAboveThreshold);

  //! Prints results of the search.
  void PrintResults(bool detailed = true);

 private:

  //! Figure extension.
  string fExtension;

  //! Initialize #fThetaLobe and #fLobe.
  void InitLobe(string lobe_file);

  //! Computes #fCovMap, #fEventsMap and #fDiffMap.
  void InitMap(const TCoverage & coverage);

  //! Computes #fLobeCorrection.
  void ComputeLobeCorrection();

  //! This function returns the significance.
  double LiMa(double nOn, double nOff, double alpha);

  //! Events.
  vector<TEvent> fEvents;

  //! Number of events.
  unsigned int fNevents;

  //! Significance threshold.
  double fThreshold;

  //! Events map.
  THealpixMap fEventsMap;

  //! Initial coverage map (not smoothed)
  THealpixMap fRawCovMap;

  //! Coverage map.
  THealpixMap fCovMap;

  //! Difference map.
  THealpixMap fDiffMap;

  //! Significance map.
  THealpixMap fLiMaMap;

  //! Opening angle at a given direction.
  vector<double> fThetaLobe;

  //! Values of the lobe.
  vector<double> fLobe;

  //! Effective radius of the window.
  double fWindowRadius;

  //! Correction to apply due to the smoothing.
  double fLobeCorrection;

  //! Pixels of #fLiMaMap with significance above #fThreshold.
  vector<long> fIpAboveThreshold;

  //! Values corresponding to #fIpAboveThreshold.
  vector<double> fValueAboveThreshold;
};

//! The Li & Ma significance
double LiMa(double nOn, double nOff, double alpha);

#endif
