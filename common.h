#ifndef _COMMON_H
#define _COMMON_H

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TColor.h>

#include "STCconstants.h"

using namespace std;

/*! \mainpage
  This is a set of C++ programs we are developing for the Pierre 
  Auger Collaboration in order to allow everyone to search for 
  small and large scale anisotropies on the sky. You will be able 
  to do the following things : 
  \arg Construct and manipulate easily \b HEALPix maps (Hierarchical Equal Area isoLatitude Pixelization of 
  the sphere).
  \arg Integrate maps through any beam window.
  \arg Draw sky maps using different coordinate systems and projections. Superimpose events or astrophysical 
  sources.
  \arg Simulate events.
  \arg Compute coverage maps using several methods
  \arg Compute the statistical significance of an excess in any given direction of the sky.
  \arg Estimate the \f$C_\ell\f$ angular power spectrum of the events.
  \arg Compute the angular correlation function of the events.
  \arg Look for a dipolar modulation in your data set.
  \arg Perform an harmonic analysis to the Right Ascension distribution of the events.
  \arg Do the scan procedure.
  \arg Perform a cross-correlation analysis of your data set with any catalog of galaxies.

  
  Our code relies on various libraries and software : 
  <a href="http://heasarc.nasa.gov/fitsio/">CFITSIO</a> (the FITS I/O library), 
  <a href="http://healpix.jpl.nasa.gov/">HEALPix</a>, \b STCoordinates (part of the 
  <a href="http://www.auger.org.ar/CDAS/">CDAS</a> software where coordinate system functions are coded) and the 
  well known <a href="http://root.cern.ch">ROOT</a>.

  <B>Authors :</B> J.-Ch. Hamilton, <a href="mailto:revenu@in2p3.fr">B. Revenu</a> and 
  <a href="mailto:rouille@kicp.uchicago.edu">B. Rouill&eacute; d'Orfeuil</a>.
  
  <B>Contributors :</B> J. D. Hague and M. C. Medina.
*/


/*! \file common.h
  You will find here tools that are used all along the code such as : coordinates transformations, a progress bar 
  showing the necessary time to complete a loop, the Pierre Auger Observatory coordinates and so on.
*/



/*! This is stricly identical to :
  \code
  vector<Type> NewObjectName(OriginalObject.size());
  for (unsigned int i=0; i<OriginalObject.size(); i++) NewObjectName(i) = OriginalObject(i).FieldName;
  \endcode
*/
#define DECLARE_VECTOR(Type,NewObjectName,OriginalObject,FieldName) std::vector<Type> NewObjectName;{size_t i=0;size_t size=OriginalObject.size();NewObjectName.reserve(size); while( i<size ){NewObjectName.push_back(OriginalObject[i].FieldName); i++;}}

/*! This is stricly identical to :
  \code
  Type* NewObjectName = new Type[OriginalObject.size()];
  for (unsigned int i=0; i<OriginalObject.size(); i++) NewObjectName(i) = OriginalObject(i).FieldName;
  \endcode
*/
#define DECLARE_POINTER(Type,NewObjectName,OriginalObject,FieldName) Type * NewObjectName=new Type[OriginalObject.size()]; {size_t i=0;size_t size=OriginalObject.size(); while( i<size ){NewObjectName[i]=OriginalObject[i].FieldName;i++;}}


/*! \namespace std
  All the files in the C++ standard library declare all of its entities within the std namespace.
*/
using namespace std;


/*! \namespace kSTC
  Uses the kSTC namespace defined in the STCoordinates package. It contains constants such as degrees to radians 
  conversion, usefull sky-transformations constants and so on.
*/
using namespace kSTC;


/*! \namespace kConstantsTK
  It defines the Piere Auger Observatory coordinates.
*/
namespace kConstantsTK
{
  //! Auger North latitude.
  static const double AugerSouthLatitude = -35.22;

  //! Auger North longitude.
  static const double AugerSouthLongitude = -69.28;

  //! Auger North latitude.
  static const double AugerNorthLatitude = 38.08;

  //! Auger North longitude.
  static const double AugerNorthLongitude = -102.62;
};


/*! 
  The TProgress bar class allows you to know the progression of your job. This can be usefull especially to have 
  an idea of the time it takes to complete a loop or to check whether the program is still running or not.
*/

//! Progress bar in the terminal. 
class TProgressBar
{
public:
  //! Constructor.
  TProgressBar();

  //! Initializes everything.
  void Zero();

  //! First value in the loop (0\%).
  unsigned int fBegin;

  //! Last value in the loop (100\%).
  unsigned int fEnd;

  //! Current value of the loop.
  unsigned int fReached;

  //! Current corresponding percentage value.
  unsigned int fCurrentPercentage;

  //! Initializes stderr.
  void InitPercent();

  //! Dump last stderr (full bar).
  void EndPercent() const;

  //! Print current percentage to stderr.
  void PrintPercent(unsigned int current);
};

/*!
  Given a input array of data containing good data and bad data (with 
  value ignore), detect contiguous points (not separated by the ignore 
  value) and return the list of the indexes of contiguous points. Example : 
  data = [10,11,12,13,-1,-1,14,15,16,-1,-1] with ignore = -1 returns [6,7,8] 
  and [0,1,2,3], index lists of good data.
*/
template<typename T> vector< vector<unsigned int> > GetContiguousPoints(unsigned int size, const T *input, T ignore_value);

//! Fills an histogram. The normalization is such as its integral equals to 1.
void GetHisto(TH1F*, const vector<double>& variable = vector<double>());

//! Plots an histogram.
void DrawHisto(TCanvas*, TH1F*, string Xaxis = "X", string Yaxis = "Y", string save = "histogram.png"); 

//! Plots a 2D graph.
void PlotXY(TCanvas*, TGraphErrors*, double xMin, double xMax, string name = "", string Xaxis = "X", string Yaxis = "Y");

//! Color table
void PaletteOrange(int nContours);

//! Color table
void PaletteBlue(int nContours);

//! Color table
void PaletteRGB(int nContours);

//! Color table
void PaletteBlueAndRed(int nContours);

//! Tells you whether the file is available or not.
bool CheckFile(string fileName);

//! Gives you the list of the fields that must have your input file of events.
void DumpFields();

#endif
