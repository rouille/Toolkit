#ifndef _SCANNER_H_
#define _SCANNER_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstring>
#include <iomanip>
#include <algorithm>

#include "events.h"
#include "healpixmap.h"
#include "maptools.h"
#include "agntools.h"

#include "TStopwatch.h"

//! Scan procedure
class TScanner
{
 public :
  //! Constructor
  TScanner(const vector<TEvent> & events, const vector<TAgn> & sources, const THealpixMap & cMap);
  
  //! Sets #fEmin and #fEmax
  void SetRangeEnergy(double eMin, double eMax);
  
  //! Sets #fNagn, #fZmin, #fZmax, #fZstep and #fNzStep
  void SetRangeRedshift(double zMin, double zMax, double zStep);
  
  //! Sets #fPsiMin, #fPsiMax, #fPsiStep and #fNpsi.
  void SetRangeAngle(double psiMin, double psiMax, double psiStep);
  
  //! Dump the ranges in energy, redshift and separation angle to the terminal.
  void PrintRanges();
  
  //! Calls the #bigN function. 
  int GetBigN(double eMinPar);
  
  //! Calls the #littleP function. 
  double GetLittleP(double zMaxPar, double psiMaxPar);
  
  //! Calls the #littleK function.
  int GetLittleK(double zMaxPar, double psiMaxPar, double eMinPar);
  
  //! Cheat, if you know your minimum-big-P parameters in advance
  void SetScanResults(double eMinAtMin, double zMaxAtMin, double psiMaxAtMin, double & minBigP, double & littlePAtMin, int & littleKAtMin, int & bigNAtMin); 
  
  //! Perform the scan
  void PerformScan(bool timerOn);
  
  //! Print the result of the scan
  void PrintScanResults(double & eMinAtMin, double & zMaxAtMin, double & psiMaxAtMin, double & minBigP, double & littlePAtMin, int & littleKAtMin, int & bigNAtMin); 
  
  //! Events that correlate with an AGN
  vector<TEvent> GetEventsThatCorrelate();
  
  //! Events that do not correlate with an AGN
  vector<TEvent> GetEventsThatDoNotCorrelate();
  
  //! AGNs within #fZmaxAtMin
  vector<TAgn> SourcesCloserZmax();
  
  //! Informations about events and AGNs at #fMinBigP
  void PrintEventsAndNearAGN();

  /*!
    Returns 1 if there is an AGN with z \f$ \leq \f$ zMaxPar with angular distance \f$ \psi \leq \f$ psiMaxPar 
    from test. Returns 0 if it's not the case.
  */
  int IsSourceNearTest(double lTest, double bTest, double zMaxPar, double psiMaxPar);
  
 private : 
  /*!
    Initializes #fNbPix, #fLpix, #fBpix, #fMinBigP, #fEminAtMin, #fZmaxAtMin, #fPsiMaxAtMin, #fLittlePAtMin, 
    #fLittleKAtMin and #fBigNAtMin.
  */
  void Init();

  //! Returns little-p for a given zMaxPar and psiMaxpar.
  double littleP(double zMaxPar, double psiMaxPar);
  
  //! Returns the number of events that correlate.
  int littleK(double zMaxPar, double psiMaxPar, double eMinPar);
  
  //! Returns the total number of events in set #fEmin \f$ \leq E \leq \f$ #fEmax with \f$ E \geq \f$ eMinPar.
  int bigN(double eMinPar);

  //! Events.
  vector<TEvent> fEvents;
  
  //! Number of events.
  unsigned int fNevents;
  
  //! Coverage Map.
  THealpixMap fCovMap;
  
  //! Number of pixels of the map.
  unsigned int fNbPix;

  //! Galactic longitude of th pixels.
  vector<double> fLpix;
  
  //! Galactic latitude of the pixels.
  vector<double> fBpix;
   
  //! Minimum energy of the events to be considered for the scan.
  double fEmin;
  
  //! Maximum energy of the events to be considered for the scan.
  double fEmax;
  
  //! Minimum redshift of the AGNs to be considered for the scan.
  double fZmin;
  
  //! maximum redshift of the AGNs of the sources to be considered.
  double fZmax;
  
  //! Step in redshift
  double fZstep; 
  
  //! Number of step in redshift
  unsigned int fNzStep;
   
  //! Minimum values of the redshift.
  vector<double> fZmins;
  
  //! AGNs.
  vector<TAgn> fAGNs;
  
  //! Number of AGNs.
  unsigned int fNagn;
  
  //! Minimum separation angle between AGNs and events.
  double fPsiMin;
  
  //! Maximum separation angle betxeen AGNs and events.
  double fPsiMax;
  
  //! Step in separation angles
  double fPsiStep;
  
  //! Minimum values of the separation angles
  vector<double> fPsis;
  
  //! Number of separation angles
  unsigned int fNpsi;
  
  //! Minimum Probability
  double fMinBigP;
  
  //! Lower value of the energy at #fMinBigP
  double fEminAtMin;
  
  //! Upper value of the redshift at #fMinBigP
  double fZmaxAtMin;
  
  //! Upper value of the separation angle at #fMinBigP
  double fPsiMaxAtMin;
  
  //! Little-P at #fMinBigP
  double fLittlePAtMin;
  
  //! Number of events that correlate at #fMinBigP
  int fLittleKAtMin;
  
  //! Total number of events at #fMinBigP
  int fBigNAtMin;
};

#endif
