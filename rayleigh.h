#ifndef _RAYLEIGH_H_
#define _RAYLEIGH_H_

#include <vector>

#include "events.h"
#include "healpixmap.h"
#include "TGraphErrors.h"

#include "TH1F.h"
#include "TF1.h"

/*!
  Harmonic analysis to the Right Ascension distribution of the events. To determine the amplitude and the phase of 
  the anisotropy, we fit the distribution to a sine wave with period \f$ \frac{2\pi}{m} (m^{th} harmonic)\f$.
*/

//! Harmonic analysis to the right ascension of the events
class TRayleigh
{
 public:
  //! Constructor.
  TRayleigh(const vector<TEvent> &events, const THealpixMap &covMap, const unsigned int nHarmonic = 1, double dRA = 3.);

  //! Destructor.
  ~TRayleigh();

  //! Computes #fRAweight and #fRAweightBins.
  void ComputeRAweight();

  //! Computes #fRAtrue and fill #fHistoTrue. 
  void ComputeDist();

  //! Draw #fHistoEvts.
  void DrawEvtsDist() const;

  //! Draw #fPlotCov
  void DrawCovDist();

  //! Draw #fHistoTrue.
  void DrawTrue();

  //! Returns #fFitParameter.
  vector<double> GetFitParameter() const {return fFitParameter;}

  //! Computes #fAmplitude. 
  void ComputeAmplitude();

  //! Computes #fPhase.
  void ComputePhase();

  //! Computes #fSignificance.
  void ComputeSignificance();

  //! Computes #fChanceProbability.
  void ComputeChanceProbability();

  //! Returns #fAmplitude.
  double GetAmplitude() const {return fAmplitude;}
  
  //! Returns #fPhase.
  double GetPhase() const {return fPhase;}
  
  //! Returns #fSignificance.
  double GetSignificance() const {return fSignificance;}

  //! Returns #fChanceProbability.
  double GetChanceProbability() const {return fChanceProbability;}
  
 private:
  //! Initializes #fRA, #fRAweight, #fRAweightBins, #fBinCenter, #fHistoEvts and #fHistoTrue.
  void Init();

  //! Computes #fCompFirst and #fCompSec.
  void ComputeComponent();

  //! Fit #fHistoTrue.
  void FitDistTrue();

  //! To know if the #fCompFirst and #fCompSec have been already computed.
  bool fComponent;

  //! Events
  const vector<TEvent> & fEvts;

  //! Coverage Map.
  THealpixMap fCovMap;

  //! Number of harmonic.
  unsigned int fNharmonic;

  //! Size of the bins of #fHistoEvts and #fHistoTrue.
  double fdRA;

  //! Number of bins of #fHistoEvts and #fHistoTrue.
  unsigned int fNbins;

  //! Number of events.
  unsigned int fNevts;

  //! Right Ascension of the events.
  vector<double> fRA;

  //! Right Ascension distribution of the events.
  TH1F * fHistoEvts;

  //! Plot #fRAweightBins as a function of #fBinCenter.
  TGraphErrors * fPlotCov;

  //! Plot #fRAtrue. 
  TH1F * fHistoTrue;

  //! Function used to fit #fHistoTrue.
  TF1 * fFitFcnTrue;

  //! Parameter of the fit.
  vector<double> fFitParameter;

  //! Equal to the integral of #fHistoEvts.
  double fNormEvts;

  //! Equal to the integral of (#fBinCenter,#fRAweightBins).
  double fNormCov;

  //! Bins center of #fHistoEvts and #fHistoTrue.
  vector<double> fBinCenter;

  //! Right Ascension values of the coverage map interpolated at the Right Ascension of the events.
  vector<double> fRAweight;

  //! Right Ascension values of the coverage map in each bins of size #fdRA.
  vector<double> fRAweightBins;

  /*!
    True values of the Right Ascension in each bins of size #fdRA. It's computed in each bin of size #fdRA 
    as the ratio of the bin content of #fHistoEvts by #fRAweightBins once normalized to the same integral.
  */
  vector<double> fRAtrue;

  //! First component. 
  double fCompFirst;

  //! Second component.
  double fCompSec;

  //! Amplitude.
  double fAmplitude;

  //! Phase.
  double fPhase;

  //! Significance.
  double fSignificance;

  //! Chance probability.
  double fChanceProbability;
};

//! Harmonic modulation.
double fitFunction(double* t, double* par);

#endif
