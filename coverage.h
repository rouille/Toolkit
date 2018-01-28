#ifndef _COVERAGE_H_
#define _COVERAGE_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include "angdist.h"
#include "timemod.h"
#include "maptools.h"
#include "healpixmap.h"
#include "harmotools.h"
#include "events.h"
#include "STClibrary.h"
#include "common.h"

// ROOT
#include "TGraph.h"
#include "TRandom2.h"
#include "TMath.h"


/*!
 This class is dedicated to the coverage map computation. Both scrambling and Semi-Analytical methods are present. 
 It is possible to correct the map estimated with the Semi-Analytical technique from the acceptance effects using 
 an acceptance file or just fitting the UTC distribution of the events and interpolating the JD distribution.
 */

//! Coverage map.
class TCoverage
{
	public :
	//! Constructor.
	TCoverage(unsigned int nside = 0);
	
	//! Returns #fMap.
	THealpixMap  GetMap() const {return fMap;}
	
	//! Returns #fLatitude.
	double GetLatitude() const {return fLatitude;}
	
	//! Returns the nside of #fMap.
	unsigned int NSide() const {return fMap.NSide();}
	
	//! Set #fCoordSystem.
	void SetCoordSystem(char str) {fCoordSystem = str;}
	
	//! Returns #fCoordSystem.
	char GetCoordSystem() const {return fCoordSystem;}
	
	//! Set #fLatitude, #fcosLatitude and #fsinLatitude.
	void SetLatitude(double latitude) {fLatitude = latitude; fcosLatitude = cos(latitude*DTOR); fsinLatitude = sin(latitude*DTOR);}
	
	//! Returns #fLongitude.
	double GetLongitude() const {return fLongitude;}
	
	//! Set #fLongitude.
	void SetLongitude(double longitude) {fLongitude = longitude;}
	
	//! set filename extension for plots
	void SetExtension(string ext) {fExtension = ext;}
	
	//! Get filename extension for plots
	string GetExtension() const {return fExtension;}
	
	//! Returns TFitFunction::fDataMin using the TAngularDistribution::GetMinAngle() function of #fThetaDist.
	double GetThetaMin() const {return fThetaDist.GetMinAngle();}
	
	//! Returns TFitFunction::fDataMax using the TAngularDistribution::GetMaxAngle() function of #fThetaDist.
	double GetThetaMax() const {return fThetaDist.GetMaxAngle();}
	
	//! Returns #fDecMin.
	double GetDecMin() const {return fDecMin;}
	
	//! Returns #fDecMax.
	double GetDecMax() const {return fDecMax;}
	
	//! Look for pixels in the field of view. 
	void ComputeConstants() const;
	
	//! Get array of zenith Right Ascension as function of time.
	vector<double> GetZenithalRA(unsigned int nptsra = 1000) const;
	
	//! Returns true if #fDecMin \f$ \leq \f$ dec \f$ \leq \f$ #fDecMax.
	bool IsInFOV(double dec) const {return dec >= fDecMin-(fMap.GetPixSize()/2.) && dec <= fDecMax+(fMap.GetPixSize()/2.);}
	
	//! Compute the declination limits. To be called once #fLatitude and TFitFunction::fDataMax of #fThetaDist 
	//! are set.
	void ComputeDeclinationLimits();
	
	//! Geometric transformation to get the declination distribution of the events from the zenithal one.
	void ZenithToDeclination(unsigned int nb = 1000);
	
	/*! 
	 We use a linear interpolation to evaluate the value of the distribution at \b HEALPix declination. We then 
	 fill a THealpixMap vector with these values. The resulting coverage map is independent of the Right 
	 Ascension when the argument is FLAT (default).
	 */
	void ComputeCoverage(unsigned int npts = 1000);
	
	//! Declination.
	vector<double> fDeclination;
	
	//! Declination distribution.
	vector<double> fDecDistribution;
	
	//! Right ascension.
	vector<double> fRightAscension;
	
	//! A TAngularDistribution instance (\f$ \theta \f$).
	TAngularDistribution fThetaDist;
	
	//! A TPhiModulation instance (\f$ \phi \f$ modulated by \f$ \theta \f$).
	TPhiModulation fPhiMod;
	
	//! A TTimeModulation instance (UTC, UTC+JD or a file).
	TTimeModulation fTimeMod;
	
	//! Correct for angular (\f$ \theta \f$ or \f$ \phi \f$) modulation.
	void CorrectForAngularModulation(string angleName);
	
	//! Correct for time modulation
	void CorrectForTimeModulation(string timeModel);
	
	//! \f$ \theta \f$ angular modulation .
	bool fThetaModulation;
	
	//! \f$ \phi \f$ angular modulation.
	bool fPhiModulation;
	
	//! Time modulation.
	bool fTimeModulation;
	
	//! A THealpixMap instance.
	THealpixMap fMap;
	
	private :
	//! Coordinate system of the \b HEALPix map.
	char fCoordSystem;
	
	//! Latitude of the observatory.
	double fLatitude;
	
	//! Longitude of the observatory.
	double fLongitude;
	
	//! Cosinus of #fLatitude.
	double fcosLatitude;
	
	//! Sinus of #fLatitude.
	double fsinLatitude;
	
	//! Lower limit of the declination.
	double fDecMin;
	
	//! Upper limit of the declination.
	double fDecMax;
	
	//! filename extension for declination distribution plot
	string fExtension;
};


/*! 
 This function returns the coverage map for a detector with saturated acceptance up to thetamax. Meaning that 
 the zenith angle distribution of events is \f$ \sin \theta \cos \theta \f$ below thetamax and 0 above. One 
 needs to specify the nside of the \b HEALPix map and the latitude of the experiment. The calculation of the 
 coverage map in this case is described in : P. Sommers, Astropart.Phys., \b 14 (2001) 271-286.
 */
THealpixMap GetAnalyticalCoverage(int nside, double thetamax, double latsite);


//! The scrambling method to compute the coverage map.
THealpixMap ComputeCoverageScrambling(const vector<TEvent>& eventsInput, unsigned int nSide, double thetaMax, double latitude, double longitude, unsigned int nScramble, unsigned int nBins, string binning = "EVENTS", string variable = "UTC");

#endif
