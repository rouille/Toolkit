#ifndef _ANGDIST_H_
#define _ANGDIST_H_

#include <vector>

#include "TH1F.h"
#include "TF1.h"

#include "healpixmap.h"
#include "events.h"
#include "maptools.h"
#include "userfcn.h"

using namespace std;

/*! \file angdist.h
 This file is dedicated to the fit of the angular distribution of the events. Both azimuthal and zenithal 
 distributions can be fitted. Concerning the Surface Detector of the Pierre Auger Observatory, only the fit 
 of the zenith angle (\f$ \theta \f$) distribution is necessary since we know that this distribution is 
 (should) be uniform. The fit of the azimuthal angle (\f$ \phi \f$) is on the other hand essential for the 
 fluorescence detector and radio experiments.
 */

//! Angular distributions of the events.
class TAngularDistribution : public TFitFunction
{
public:
	//! Constructor.
	TAngularDistribution();
    
	//! Destructor.
	virtual ~TAngularDistribution()
	{
	}
	
	//! Set TFitFunction::fDataName.
	void SetAngleName(string name) {fDataName = name;}
	
	//! Set TFitFunction::fDataMin, TFitFunction::fDataMax, TFitFunction::fNBins and TFitFunction::fChi2Limit.
	void SetOptions(double xMin, double xMax, unsigned int nBins, double chi2);
	
	//! Returns the value of the fit  at a specific angle point.
	double GetAccAngle(double* angle) const;
	
	//! Returns TFitFunction::fDataMax.
	double GetMaxAngle() const {return fDataMax;}
	
	//! Returns TFitFunction::fDataMin.
	double GetMinAngle() const {return fDataMin;}
	
	//! Returns TFitFunction::fNBins.
	unsigned int GetNBins() const {return fNBins;}
	
	//! Returns TFitFunction::fChi2Limit.
	double GetChi2Limit() const {return fChi2Limit;}
	
	//! Returns TFitFunction::fDataBins.
	vector<double> GetAngleBins() const {return fDataBins;}
	
	//! Returns TFitFunction::fDataFit.
	vector<double> GetAngleFit() const {return fDataFit;}
	
	//! Returns TFitFunction::fData.
	vector<double> GetAngles() const {return fData;}
	
	//! Returns TFitFunction::fFitParameters.
	vector<double> GetFitParameters() const {return fFitParameters;}
	
	//! Returns TFitFunction::fFitParametersErrors.
	vector<double> GetFitParametersErrors() const {return fFitParametersErrors;}
	
protected:
	//! \f$ \phi \f$ angular modulation.
	bool fPhiModulation;
	
	//! Functions used for the fit.
	string fFitUsed;
	
	//! Value of angle regularly gridded.
	vector<double> fAngleBins;
	
	//! Fit function evaluated at #fAngleBins points.
	vector<double> fAngleFit;
	
	//! Set TFitFunction::fNBins.
	void SetNBins(unsigned int nBins) {fNBins = nBins;}
};



//! Fit of the azimuthal distribution of the events.
class TPhiModulation : public TAngularDistribution
{
public:
	//! Constructor
	TPhiModulation();
	
	//! Destructor
	virtual ~TPhiModulation();
	
	//! Set #fThetaMin
	void SetThetaMin(double x) {fThetaMin = x;}
	
	//! Set #fThetaMax
	void SetThetaMax(double x) {fThetaMax = x;}
	
	//! Set #fPhiLaw
	TF1 * SetPhiLaw(double (*f)(double *t, double *par), double thetamin, double thetamax, int npars);
	
	/*!
	 First step : for each \f$ \theta \f$ bin, we fit the \f$ \phi \f$ histogram with function 
	 TFitFunction::fAngFitFunc and collect the fitted parameters in each \f$ \theta \f$ bin. Second step : 
	 we fit the law followed by the parameters as a function of \f$ \theta \f$.
	 */
	void ComputePhiModulation(const vector<TEvent>&, unsigned int nthetabins);
	
	//! Returns the value of the fit angle point.
	double GetAccAngle(double *phi, double * pars);
	
	//! Function pointer describing the \f$ \phi \f$ modulation given by TFitFunction::fAngFitFunc.
	double (*fPhiLawFunc)(double *t, double *par);
	
	//! Number of parameters of the fit
	int fNParsPhiLaw;
	
private:
	//! Histogram of the \f$ \phi \f$ modulation.
	TH1F * fModPhi;
	
	//! Function used for the angular fit.
	TF1 * fPhiLaw;
	
	//! Minimum \f$ \theta \f$. 
	double fThetaMin;
	
	//! Maximum \f$ \theta \f$.
	double fThetaMax;
	
	//! Parameters of the fit of the \f$ \phi \f$ modulation.
	double * fModParams;
};

#endif
