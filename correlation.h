#ifndef _CORR_H_
#define _CORR_H_

#include "healpixmap.h"
#include "events.h"

#include "TGraph.h"
#include <TGraphAsymmErrors.h>

#define WITHSINE



//! 2 points angular correlation function.
class TCorrelation
{
 public:
  /*! 
    Constructor.
    Default : number of pairs as a function of the maximum separation angle \f$ \alpha \f$  will be computed from 
    \f$0^\circ up to 180^\circ\f$ by step of \f$1^\circ\f$.
   */
  TCorrelation(const vector<TEvent>& events, unsigned int alphaStep = 1, double alphaMax = 180);

  //! Destructor.
  ~TCorrelation();

  //! Set #fExtension.
  void SetExtension(string ext) {fExtension = ext;}

  /*! 
    Number of event pairs with angular separation smaller than a given value. For each event \f$ \vec{u_i} \f$, the angle 
    \f$ \alpha = acos(\vec{u_i}.\vec{u_j}) \f$ with \f$ 0 \leq i \leq \f$ #fNevents, \f$ i+1 \leq j \leq \f$ #fNevents is 
    computed. The #fEventsCorrelation graph is filled.
   */ 
  void ComputeEventsCorrelation();

  /*! 
    Computes the mean number of pairs expected. The \f$ C_\ell \f$ of the coverage map are computed and the \f$ C(\theta) \f$ 
    are then derived computed. Only the mean number of pairs is computed. The #fCoverageCorrelation graph is filled.
    Important: #ComputeEventsCorrelation must be called first. 
   */
  void ComputeCoverageCorrelationMap(const THealpixMap& covMap);

  /*!
    Computes the mean number of pairs expected in simulated isotropic sets of #fNevents, distributed in proportion to 
    the exposure of the Observatory. The errors bars correspond to the dispersion of disp% of the simulations. 
   */
  void ComputeCoverageCorrelationSimu(double longitude, double latitude, const THealpixMap & map, unsigned int nSimu, double dispersion, const vector<double> & thVal, const vector<double> & pthVal, string utcFile = "", string jdFile = "", string globalFile = "");

  /*!
    Computes the mean number of pairs expected by scrambling the original data set. The errors bars correspond to 
    the dispersion of disp% of the simulations.
   */
  void ComputeCoverageCorrelationMC(double longitude, double latitude, unsigned int nMC, double dispersion, string binningType, string scramblingType, double thetaMax);

  //! Draw #fEventsCorrelation.
  void DrawEvents() const;

  //! Draw #fCoverageCorrelation.
  void DrawCov() const;

  //! Draw #fEventsCorrelation and #fCoverageCorrelation in the same canvas.
  void DrawBoth() const;

  //! Get the number of event pairs with angular separation smaller than\f$ \alpha \f$.
  void GetEventsCorrelation(vector<double>& alpha, vector<double>& Np) const;

  //! Get the mean number of pairs expected and their dispersion with angular separation smaller than\f$ \alpha \f$.
  void GetCoverageCorrelation(vector<double>& alpha,vector<double>& Np,vector<double>& errLow,vector<double>& errHigh) const;

private:
  //! Initializes #fTheta and #fPhi.
  void InitAngles(const vector<TEvent> & events);

  //! Initializes #fEventsCorrelation and #fCoverageCorrelation 
  void InitGraph();

  //! Events (actual or simulated).
  const vector<TEvent> & fEvents;

  //! Number of events.
  unsigned int fNevents;
  
  //! Angular size of the step in the separation angle
  unsigned int fAlphaStep;

  //! Number of step in the maximum separation angle.
  unsigned int fNstep;

  //! Maximum separation angle for which we are looking for correlation.
  double fAlphaMax;

  //! Number of events pair with angular separation smaller than a given value.
  TGraph * fEventsCorrelation;

  //! Same that #fEventsCorrelation for the isotropic expectations.
  TGraphAsymmErrors * fCoverageCorrelation;

  //! \f$ \theta \f$ spherical angle (on the sky) of fEvents.
  double * fTheta;

  //! \f$ \phi \f$ spherical angle (on the sky) of fEvents.
  double * fPhi;

  //! \f$ \cos\theta\f$. Cosinus of the spherical angle.
  double * fCosTheta;

  //! \f$ \sin\theta\f$. Cosinus of the spherical angle.
  double * fSinTheta;

  //! Figure extension.
  string fExtension;
};

//! Computes \f$ C(\theta) \f$ from \f$ C_\ell \f$.
vector<double> Cl2Ctheta(int lmax, const vector<double> & cl, const vector<double> &  theta);

//! Computes \f$ C(\theta) \f$ from \f$ C_\ell \f$ at the given \f$ \theta \f$.
double GetCtheta(int lmax, const vector<double> & cl, double thetarad);

//! Using Cl2Ctheta, the histogram of the correlation angle is filled.
void ComputeCtheta(int lmax, TH1F ** histothetacov, const vector<double> & thecl);

#endif
