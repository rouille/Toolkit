#ifndef _SIMUEVENTS_H
#define _SIMUEVENTS_H

#include <vector>

#include "events.h"
#include "healpixmap.h"

#include "TRandom2.h"

using namespace std;

/*! \file simuevents.h 
  This file contains the codes we developped for simulating events. The procedure is described in GAP2005-083. 
  We simulate the events according to a underlying anisotropic map provided by the user in input, this map can 
  contain large scale anisotropies, points sources or any sky pattern. It accounts for the zenith angle acceptance 
  and optionnaly a time variable acceptance of the  array. The events are returned as a vector<TEvent> that can be 
  directly used in our analysis tools.
*/


/*! 
  Event simulator: map is the input anisotropic map (the map one would get with an infinite number of events). 
  nEvents is the number of events to be simulated. thetaLaw and pThetaLaw are the shape of the zenith angle 
  acceptance. thetaLaw is in degrees and has to go up to \f$ 180^\circ \f$ (obviously it is zero between 
  \f$ 90^\circ \f$ and \f$ 180^\circ \f$). The optional arguments utcFile, jdFile and globalFile are files where 
  the relative acceptance is described as a function of time. The files have to be in ASCII format with two 
  columns : time (utcFile : hour, jdFile : julian day and global file: UTC seconds) and relative accetpance. 
  utcFile and jdFile may be used as the same time if your acceptance model depends separately on the time of the 
  day (utcFile) and date of the year (jdFile). The globalFile case is when your model is more general and depend 
  on the time in a complex way, you then specify the acceptance as a function of the time in general (in UTC 
  seconds - from January 1st, 1970). A detailed description of the algorithms used here can be found in 
  GAP2005-083.
*/
vector<TEvent> SimulateEvents(const THealpixMap & map, unsigned int nEvents, const vector<double> & thetaLaw, const vector<double> & pThetaLaw, double latSite, double lonSite, string utcFile = "", string jdFile = "", string globalFile = "");

//! Simple case of an isotropic distribution
vector<TEvent> SimulateIsoEvents(unsigned int nEvents, const vector<double> & thetaLaw, const vector<double> & pThetaLaw, double latSite, double lonSite);

//! It generates the UTC distribution. Uniform if utcFile is not supplied.
void UtcGenerator(string utcFile, vector<double> & utc, vector<double> & pUtc);

//! It generates the JD distribution. Uniform if jdFile is not supplied.
void JdGenerator(string jdFile, vector<double> & jd, vector<double> & pJd);

//! It generates the temporal distribution enclosed in globalFile.
void GlobalTimeGenerator(string globalFile, vector<double> & jd, vector<double> & pJd);


/*! 
  This function returns nb values distributed randomly following the distribution specified by x and y vectors. 
  The algorithm used here is the inverse transform method.
*/
vector<double> distri(const vector<double> & x, const vector<double> & y, unsigned int nb);


/*! 
  This function returns one value distributed randomly following the distribution specified by x and y vectors. 
  The algorithm used here is the inverse transform method.
*/
double distri(const vector<double> & x, const vector<double> & y);

/*!
  This function returns nb integer values distributed randomly following the distribution specified by x and y 
  vectors. The algorithm used here is the inverse transform method.
*/
vector<long> distriDiscrete(const vector<long> & x, const vector<double> & y, unsigned int nb);

/*!
  This function returns one integer value distributed randomly following the distribution specified by x and y 
  vectors. The algorithm used here is the inverse transform method.
*/
long distriDiscrete(const vector<long> & x, const vector<double> & y);

#endif
