#ifndef _AGNTOOLS_H_
#define _AGNTOOLS_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstring>
#include <iomanip>
#include <algorithm>

#include "STClibrary.h"
#include "events.h"


using namespace std;

/*!
  This class only contains public attributes characterizing the active galaxy nuclei in Catalogue/VCV.dat : 
  RA | Dec | Redshift | Name. This file is extracted from the \f$ 12^{th} \f$ edition of the catalog of quasars and 
  active galaxy nuclei by V&eacute;ron-Cetty & V&eacute;ron. The file contains 694 active galaxies with redshifts 
  \f$ z \leq 0.024 \f$, corresponding to distance \f$ D \leq 100 \f$ Mpc. You should refer to 
  M.-P. V&eacute;ron-Cetty, P. V&eacute;ron, Astron. Astrophys. 455, 773 (2006) for further explanations.
*/


//! AGNs in the V&eacute;ron-Cetty & V&eacute;ron catalog.
class TAgn
{
 public :
  TAgn() {}
  
  //! Redshift.
  double fZ;
  
  //! Right Ascension.
  double fRa;
  
  //! Declination.
  double fDec;
  
  //! Galactic longitude.
  double fL;
  
  //! Galactic Latitude.
  double fB;
  
  //! Name.
  string fName;
};

//! Get the AGN's list.
vector<TAgn> Get_AGNs(char* sourcedata, double zMax, double decLimit);

//! Dump AGN info to terminal.
void Print_AGN(vector<TAgn> sources, unsigned int IndexNum);

//! Sort the AGNs according to their redshift in the range \f$ [z_{min},z_{max}]\f$.
vector<TAgn> SortCutRedShift(vector<TAgn> zUnsorted, bool down, double zMin, double zMax);

//! Sort the events according to their energy in the range \f$ [E_{min},E_{max}]\f$.
vector<TEvent> SortCutEnergy(vector<TEvent> eUnsorted, bool down, double eMin, double eMax);

#endif
