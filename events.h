#ifndef _EVENT_H
#define _EVENT_H

#include <vector>
#include <string>
#include "common.h"

using namespace std;

/*! 
  This class only contains public attributes characterizing the events used in the Coverage & Anisotropy Toolkit. 
  Using the Herald reconstruction, the events file is given by : 
  \code
  awk -f Events/herald2events.awk heraldFile > eventsFile
  \endcode 
  where heraldFile is the file you can download directly on the 
  <a href="http://auger.colostate.edu/private/herald/">Herald</a> web page. Of course, The same procedure can be 
  used to produce an eventsFile based on the <a href="http://augerobserver.fzk.de">Auger Observer</a>.  
*/

//! Characterizes an event
class TEvent
{
 public :
  //! Constructor.
  TEvent() {}

  //! Event id.
  unsigned long long fId;
  
  //! \f$ \theta \f$.
  double fTheta;
  
  //! Error on \f$ \theta \f$.
  double fdTheta;

  //! \f$ \phi \f$ (\f$ \phi \f$ is measured northwards from East).
  double fPhi;
  
  //! Error on \f$ \phi \f$.
  double fdPhi;
  
  //! Galactic Longitude.
  double fL;
  
  //! Galactic Latitude.
  double fB;
  
  //! Right Ascension.
  double fRa;
  
  //! Declination.
  double fDec;
  
  //! UTC second time of the event : seconds since January 1, 1970 00h00 in UTC
  double fUTCs;
  
  //! Core time (in ns) taken from the UTC time of the event.
  string fTcore;
  
  //! Elevation (angle above the horizon).
  double fElevation;

  //! Error on elevation.
  double fdElevation;

  //! Azimuth (The azimuth is measured eastwards from North).
  double fAzimuth;

  //! Error on azimuth.
  double fdAzimuth;

  //! Year of the event.
  int fYear;

  //! Month of the event.
  int fMonth;

  //! Day of the event.
  int fDay;

  //! UTC of the event (decimal hour).
  double fUTCh;

  //! Energy.
  double fEnergy;
};

//! Reads an events file and returns a TEvent vector.
vector<TEvent> GetEvents(string fileNameInit, string cutFile = "");

//! This function allows one to choose any period of data. 
vector<TEvent> SelectPeriod(const vector<TEvent> &input, int yearMin, int monthMin, int dayMin, double UTChMin, int yearMax, int monthMax, int dayMax, double UTChMax = 24.);

//! This function keeps only the events falling in a period where the acceptance specified by the file 
//! "cutFileName" in non zero.
vector<TEvent> KillBadEvents(const vector<TEvent> &input, string cutFileName);

//! This function performs scrambling of the event dataset.
vector<TEvent> ScrambleData(const vector<TEvent>& eventsInput, unsigned int nBins, string binning, double latitude, double longitude, double thetaMax, string variable = "UTC");

//! This function returns a vector of vectors of TEvent binned in \f$ \theta \f$ for the scrambling.
vector<vector<TEvent> > BinEvents(const vector<TEvent>& eventsInput, unsigned int nBins, string binning, double thetaMax);

//! This is the function that actually does the scrambling in the bins of \f$ \theta \f$.
vector<TEvent> DoTheScrambling(const vector<vector<TEvent> >& eventsBinned, string variable, double latitude, double longitude);

//! Needed to get the events from the file.
char* SkipSpaces(char* input);

//! Display the local coordinates distribution and the user defined zenith angle law used to simulate events
void ShowLocalCoord(const vector<TEvent> & events, const vector<double> & theta, const vector<double> & pTheta);

//! Display the local coordinates distribution 
void ShowLocalCoord(const vector<TEvent> & events);

//! Display the local coordinates distribution
void ShowEqCoord(const vector<TEvent> & events);

//! Display the local coordinates distribution
void ShowArrivalTimesCoord(const vector<TEvent> & events);

#endif
