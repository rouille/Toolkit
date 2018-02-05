#ifndef _USERFCN_H_
#define _USERFCN_H_

#include "common.h"
#include "TF1.h"
#include "TMath.h"

/*! 
  Definition of the parameters needed for the angular fit.
*/
enum eParameterStatus
  {
    eFree = 0,
    eFixed
  };


//! Parameters of the functions used for the fit of the angular distribution.
class TParameter
{
 public:
  //! Constructor.
  TParameter();

  //! Set #fStatus and #fValue.
  void SetParameter(eParameterStatus status, double value) {fStatus = status; fValue = value;}

  //! Status of the parameters of the fit. Either #eFree or #eFixed.
  eParameterStatus fStatus;

  //! Value taken by #fStatus.
  double fValue;
};



//! Fit of the angular distribution.
class TFitFunction
{
 public:
  //! Constructor.
  TFitFunction();

  //! Destructor.
  virtual ~TFitFunction();

  //! Function used for the angular fit.
  double (*fAngFitFunc)(double *t, double *par);

  //! Maximum degree the spline or the polynomial functions can reach.
  unsigned int fDegreeMax;

  //! \f$ \chi^2 \f$ of the angular fit.
  double fChi2;

  //! \f$ \chi^2 \f$ per NDF upper limit of the fit.
  double fChi2Limit;

  //! \f$ \chi^2 \f$ per NDF of the angular fit.
  double fChi2perNDF;

  //! Function used for the angular fit.
  TF1 * fFitFcn;

  //! Angular Distribution.
  TH1F * fDistribution;

  //! Set #fData.
  void SetData(const vector<double>& data) {fData = data;} 

  //! Set #fDataName.
  void SetDataName(string name) {fDataName = name;}

  //! Once the fit is done, this function compute #fDataFit in each #fDataBins.
  void ComputeFitFunction();

  //! Plot the distribution.
  void PlotDistribution() const;

  //! Set #fExtension.
  void SetExtension(string ext) {fExtension = ext;}

  //! Launch the angular fit.
  bool Run();

  //! Number of bins of the angular histogram which will be fitted.
  unsigned int fNBins;  

  //! Data. Could be either \f$ \theta \f$ or \f$ \phi \f$.
  vector<double> fData;

  //! First bin of the histogram.
  double fDataMin;

  //! Last bin of the histogram.
  double fDataMax;

  //! Name of the data.
  string fDataName;

  //! Figure extension.
  string fExtension;

  //! Value of angle regularly gridded.
  vector<double> fDataBins;

  //! Fit function evaluated at #fDataBins points.
  vector<double> fDataFit;

  //! Parameters of the fit.
  double * fPFitParameters;

  //! Parameters of the fit.
  vector<double> fFitParameters;

  //! Error on the parameters of the fit.
  vector<double> fFitParametersErrors;

  //! Set #fNPars, the number of parameters used for the angular fit.
  void SetNbParameters(unsigned int npar) {fNPars = npar; fParameters.resize(fNPars);}

  //! Number of parameters of the angular fit.
  unsigned int fNPars;

  //! Parameters of the fit initialized by the user.
  vector<TParameter> fParameters;
};


/*! 
  This function is called to establish the type of interpolating spline for a tabulated set of functional values 
  \f$ X_i, Y_i = f(X_i) \f$. This function is called only once to process an entire tabulated function in arrays 
  \f$ X_i \f$ and \f$ Y_i \f$. Once this has been done, values of the interpolated function for any value of x 
  are obtained by calls to a seperate routine spline_interp (Cubic Spline Interpolation). See Numerical Recipes in 
  C for details.
*/
vector<double> spline_init(const vector<double> & x, const vector<double> & y);

/*!
  Given the vectors xa and ya, which tabulate a function (with the xa's in order), and given the vector y2a, which 
  is the output from spline above and given the vector x, this routine returns a cubic spline interpolated value y. 
  See Numerical Recipes in C for details.
*/
vector<double> spline_interp(const vector<double> & xa, const vector<double> & ya, const vector<double> & y2a, const vector<double> & x);

/*!
  Given the vectors x and y, wich tabulate a function (with the x's in order), and given the vector u, this routine 
  returns a linear interpolation.
*/
vector<double> linear_interp(const vector<double> & x, const vector<double> & y, const vector<double> & u);

/*!
  Given the vectors x and y, wich tabulate a function (with the x's in order), this routine returns a linear 
  interpolation at point u.
*/
double linear_interp(const vector<double> & x, const vector<double> & y, double u);

//! Spline function.
double splineFunction(double * t, double * par);

//! Polynomial function.
double polyFunction(double * t, double * par);

//! Geometric function.
double geoFunction(double * t, double * par = 0x0);

//! Fermi-Dirac function.
double fdFunction(double * t, double * par);

//! Fermi Dirac - Geometric - Spline product function.
double fdsplFunction(double * t, double * par);

//! Fermi Dirac - Geometric - Polynomial product function.
double fdpolyFunction(double * t, double * par);

//! Geometric - Spline product function.
double geosplFunction(double * t, double * par);

//! Geometric - Polynomial product function.
double geopolyFunction(double * t, double * par);

//! Adapted to scintillator-like theta aperture (CODALEMA)
double scintillatorFunction(double * t, double * par);

//! Gaussian fit
double GaussFunction(double * t, double * par);

/*! 
  Usefull for the Surface Detector of the Pierre Auger Observatory. \f$ \phi \f$ modulation : large \f$ \theta \f$ 
  and low energies.
*/
double AugerPhiFunction(double * t, double * par);

// Gives the amplitude of the \f$ \phi \f$ modulation as a function of \f$ \theta \f$
double ModPhiThetaLaw(double * t, double * par);

//! The probability that a Poissonian distribution with a given MEAN will return a value \f$ \geq \f$ VALUE.
double PoissonFluctuation(double mean, double value);

//! Binomial
double Binomial(int NTot, int k, double p);

//! Cumulative binomial
double bigP(int Ttot, int k, double p);

//! 1\f$ \sigma \f$ Opening Angle based on detector errors dTheta and dPhi
double OneSigmaOpeningAngle(double theta, double dTheta, double dPhi);

/*! 
  This is a simple 5 points Newton-Cotes formula that allows numerical integration of the function given by x 
  and y. This (simple) algorithm requires the x values to be equally spaced and in increasing order. It is very 
  accurate as long as the function is well sampled (see Numerical Recipes in C for details).
*/
double integrate_nc5(const vector<double> & x, const vector<double> & y);

#endif
