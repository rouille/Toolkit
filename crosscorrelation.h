#ifndef _CROSSCORR_H_
#define _CROSSCORR_H_

#include "events.h"
#include "TGraph.h"
#include "TH1F.h"
 #include <TGraphAsymmErrors.h>
#include "healpixmap.h"

#define WITHSINE


//! Cross correlation between events and astrophysical objects.
class TCrossCorrelation
{
 public:
  //! Constructor.
  TCrossCorrelation(const vector<TEvent>& events, const vector<double>& lSources, const vector<double>& bSources, const vector<double>& wSources, double alphaStep = 1, double alphaMax = 180);

  //! Destructor.
  ~TCrossCorrelation();

  //! Set fExtension. 
  void SetExtension(string ext) {fExtension = ext;}

  /*! 
    Computes for each events the separation angle with the astrophysical objects. For each event direction, 
    \f$ \vec{u_i} \f$ (0 \f$ \leq i \leq \f$ #fNevents), and each astrophysical sources \f$ \vec{v_j} \f$
    (0 \f$ \leq j \leq \f$ #fNsources), the separation angle \f$ \alpha = acos(\vec{u_i}.\vec{v_j})\f$ is 
    computed. #fEventsCrossCorrelation is filled according to #fWsources.
 */ 
  void ComputeEventsCrossCorrelation();

  /*! 
    Computes for each simulated data set of #fNevents the separation angle with the astrophysical objects. The 
    operation is repeated nSimu times and the average separation angle is computed. The error bars correspond 
    to the dispersion of disp% of the simulations. 
 */ 
  void ComputeCoverageCrossCorrelation(double longitude, double latitude, const THealpixMap& map, unsigned int nSimu, double dispersion, const vector<double>& thVal, const vector<double>& pthVal, string utcFile = "", string jdFile = "", string globalFile = "");
  
  //! Draw #fEventsCrossCorrelation.
  void DrawEvents() const;

  //! Draw #fCoverageCrossCorrelation.
  void DrawCov() const;

  //! Draw #fEventsCrossCorrelation and #fCoverageCrossCorrelation in the same canvas.
  void DrawBoth() const;

  //! Draw #fHistoEvents and #fHistoCov in Relative excess in the same canvas.
  void DrawBothRelativeExcess() const;

  //! Get the number of event pairs with angular separation smaller than\f$ \alpha \f$.
  void GetEventsCrossCorrelation(vector<double>& alpha, vector<double>& Np) const;

  //! Get the mean number of pairs expected and their dispersion with angular separation smaller than\f$ \alpha \f$. 
  void GetCoverageCrossCorrelation(vector<double>& alpha,vector<double>& Np,vector<double>& errLow,vector<double>& errHigh) const;

private:
  //! Initializes #fThetaEvents and #fPhiEvents.
  void InitAnglesEvents(const vector<TEvent> & events);

  //! Initializes #fThetaSources and #fPhiSources.
  void InitAnglesSources();
  
  //! Initializes #fNstep, #fEventsCrossCorrelation and #CoverageCrossCorrelation.
  void InitGraph();

  //! Computes the relative excess of cumulative number of pairs in each separation angle bin with respect to 
  //! #fCoverageCrossCorrelation.
  void RelativeExcess();

  //! Events (actual or simulated).
  const vector<TEvent> & fEvents;

  //! Galactic longitude of the astrophysical objects.
  vector<double> fLsources;

  //! Galactic latitude of the astrophysical objects.
  vector<double> fBsources;

  //! Statistical weight of the astrophysical objects.
  vector<double> fWsources;

  //! Number of bins of #fHistoEvents and #fHistoCov.
  unsigned int fNstep;

  //! Size of the bins of #fHistoEvents and #fHistoCov.
  double fAlphaStep;

  //! Upper limit for which we are looking for the cross correlation.
  double fAlphaMax;

  //! Events cross correlation angle distribution.
  TGraph * fEventsCrossCorrelation;

  //! Background cross correlation angle distribution.
  TGraphAsymmErrors * fCoverageCrossCorrelation;

  //! \f$ \theta \f$ spherical angle (on the sky) of fEvents.
  double * fThetaEvents;

  //! \f$ \phi \f$ spherical angle (on the sky) of fEvents.
  double * fPhiEvents;

  //! \f$ \cos \theta \f$. Cosinus of the spherical angle of fEvents.
  double * fCosThetaEvents;

  //! \f$ \sin \theta \f$. Cosinus of the spherical angle of fEvents.
  double * fSinThetaEvents;
  
  //! \f$ \theta \f$ spherical angle (on the sky) of the sources.
  double * fThetaSources;

  //! \f$ \phi \f$ spherical angle (on the sky) of the sources.
  double * fPhiSources;

  //! \f$ \cos \theta \f$. Cosinus of the spherical angle of the sources.
  double * fCosThetaSources;

  //! \f$ \sin \theta\f$. Cosinus of the spherical angle of the sources.
  double * fSinThetaSources;

  //! Number of events.
  unsigned int fNevents;

  //! Number of astrophysical objects
  unsigned int fNsources;

  //! Figure extension.
  string fExtension;
};

#endif
