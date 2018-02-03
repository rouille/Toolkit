#ifndef _MAPTOOLS_H
#define _MAPTOOLS_H

#include <complex>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <string>

#include "healpixmap.h"
#include "events.h"
#include "agntools.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"

using namespace std;

/*! \file maptools.h
  This file contains general functions related to the map design and plot.
*/


//! Given a nside and galactic coordinates this function returns an events map.
THealpixMap map_events(unsigned int nside, const vector<double> & l, const vector<double> & b);

//! Get the events at less than the specified radius (degrees) from a given direction (l,b) in degrees.
vector<TEvent> GetEventsFromLB(const vector<TEvent> & events, double l, double b, double radius);

//! Get the pixels at less than the specified radius (degrees) from a given direction (l,b) in degrees.
vector<long> GetPixFromLB(const THealpixMap & map, double l, double b, double radius);

//! Get the pointing vector from the galactic coordinates (in degrees).
vector<double> ll2uv(double l, double b);

/*!
  Comuputes the angular disance between 2 points on the shpere : (lSky1,bSky1) and (lSky2,bSky2). 
  Takes angles in degrees and returns degrees. Works in BOTH Ra, Dec and L, B coord. systems : 
  \f$ -180^\circ < L \leq 180^\circ \f$ and \f$ -90^\circ \leq B \leq 90^\circ \f$ | 
  \f$ 0^\circ < Ra \leq 360^\circ \f$ and \f$ -90^\circ \leq Dec \leq 90^\circ \f$.
*/
double AngularDistance(double lSky1, double bSky1, double lSky2, double bSky2); 

/*!
  Comuputes the angular disance between 2 points on the shpere based on (\f$ \theta,\phi \f$) coordinates. 
  Takes angles in degrees and returns degrees. With \f$ 0 \leq \theta \leq \pi/2 \f$ measured from zenith and 
  \f$ -\pi < \phi \leq \pi \f$ measured from west.
*/
double TP_AngularDistance(double theta1, double phi1, double theta2, double phi2);

/*! 
  Map projections.
  KEEP THE SAME PROTOTYPE FOR THE FUNCTIONS IN ORDER TO USE A UNIQUE AND GENERIC static bool 
  (*gCoordFcnAng2XY)(double, double, double, double, double, double&, double&) 
  IN VARIOUS FUNCTIONS. \f$ \theta \f$ is in [\f$ -\pi/2,pi/2 \f$]. All angles in rad and (x,y) coordinates 
  scaled to [-1,1]. Converts \f$ \theta \f$, \f$ \phi \f$ to x, y in Lambert Azimuthal projection center of 
  the map is longituderefrad input : the \f$ \sin \f$ and \f$ \cos \f$ of the reference latitude.
*/
bool AngtoXYLambertAzimuthal(double thetarad, double phirad, double & xrad, double & yrad, double longituderefrad, double coslatref, double sinlatref);

/*! 
  Converts (x,y) to (\f$ \theta \f$,\f$ \phi \f$) in Lambert Azimuthal projectioni. Reciprocal of 
  AngtoXYLambertAzimuthal. Returns bool since the conversion is not always possible if \f$ \sqrt{x^2+y^2} \f$ 
  is too big !
*/
bool XYtoAngLambertAzimuthal(double x, double y, double & theta, double & phi, double longituderefrad = 0, double coslatref = 0, double sinlatref = 0);

//! Converts (\f$ \theta,\phi \f$) to (x,y) in Mollweide projection returns always true.
bool AngtoXYMollweide(double theta, double phi, double & x, double & y, double = 0, double = 0, double = 0);

/*! 
  Converts (x,y) to (\f$ \theta\,\phi \f$) in Mollweide projection. Returns true or false according to the value 
  of \f$ \sqrt{x^2+y^2} \f$.
*/
bool XYtoAngMollweide(double x, double y, double & theta, double & phi, double = 0, double = 0, double = 0);


//! Converts (x,y) to (\f$ \theta,\phi \f$) in radians
void xyz2thetaphi(double x, double y, double z, double & theta, double & phi);


/*!
  Returns the numbers of all pixels whose centers lie within radius of ipix in listpix.
  This method works only for the RING scheme 
  Adapted from HealpixBase.
*/
void query_disc(int nside, double theta, double phi, double radius, vector<int>& listpix);
int ring_above(int nside, double z);
void in_ring(int nside, int iz, double phi0, double dphi, vector<int>& listir);

#endif
