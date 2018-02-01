#ifndef _FITDIP_H_
#define _FITDIP_H_

#include <vector>

#include "healpixmap.h"

#include "TMinuit.h"
#include "TMatrixD.h"


/*!
  One directly fits the product of the coverage map with a general dipole (three degrees of freedom) superimposed 
  on a uniform background.
*/

//! Look for a dipolar modulation in the data set.
class TFitDipole
{
 public:
  //! Constructor.
  TFitDipole(const THealpixMap & evtMap, const THealpixMap & covMap, const double ampDipole = 0.1, const double RAdipole = 180., const vector<double> & lobe = vector<double>(), const vector<double> & thetaLobe = vector<double>());

  //! Destructor.
  ~TFitDipole();

  //! Fit procedure.
  void FitProcedure();

  //! Returns #fNParameters. 
  int GetNParameters() const {return fNParameters;}
  
  //! Returns #fParameters.
  vector<double> GetParameters() const {return fParameters;}
  
  //! Returns #fParametersErrors.
  vector<double> GetErrorParameters() const {return fParametersErrors;}
  
 private:
  //! Initializes the form factor and the coordinates of the pixels of my maps.
  void Init();

  //! Value of the beam.
  vector<double> fLobe;

  //! Position of the beam on the sphere.
  vector<double> fThetaLobe;

  //! Amplitude of the dipole. Input value.
  double fAmpDipole;

  //! Longitude of the dipole. Input value.
  double fRAdipole;
  
  //! Initializes the Minuit procedure.
  TMinuit* fMinuit;

  //! Number of parameters.
  int fNParameters;

  //! Parmeters of the fit.
  vector<double> fParameters;

  //! Errors on the parameters.
  vector<double> fParametersErrors;
};

//! \f$ \chi^2 \f$ function (minimization function)
void FCN(int& npar, double *gin, double& f, double *par, int iflag);

#endif
