#ifndef _HARMOTOOLS_H
#define _HARMOTOOLS_H

#include <complex>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>

extern "C" 
{
#include "chealpix.h"
#include "fitsio.h"
}

#include "healpixmap.h"

using namespace std;


/*! \file harmotools.h
  Using the FITS (Flexible Image Transport System) I/O library, you will be able to read/write FITS standard which 
  is the natural binary file format of \b HEALPix. For example, you can write a 2D FITS image such as a map on the 
  disk. You will also find in this file functions related to the harmonic transform extensively used in our code 
  (smoothing, integration, upgrade/degrade the resolution of your map and so on). \b HEALPix binary (anafast and 
  synfast) are used in that case. 
*/ 


/*! 
  Writes a 2d FITS image on the disk. This image can be read and display with your favorite data analysis program 
  or more easily with the well known xv image viewer available on all Linux-like platforms.
*/
void write_fits_image(const vector<vector<double> > & image, char * filename);

//! Reads a 2D FITS image from a FITS file.
vector<vector<double> > read_fits_image(char * filename);

/*! 
  Writes a FITS file containing the \f$ a_{\ell m} \f$ that are input. The \f$ a_{\ell m} \f$ format is described 
  in the \b HEALPix documentation.
*/
void write_fits_alms(const vector<vector<complex<double> > > & alms, char * filename);

//! Writes a FITS file containing the \f$ C_\ell \f$ that is input.
void write_fits_cl(const vector<double> & ClT, char * filename);

//! Reads the \f$ a_{\ell m} \f$ contained in the file and returns them as a vector<vector<complex<double> > >
void read_fits_alms(char * filename, vector<vector<complex<double> > > & alms);

//! Reads the \f$ C_\ell \f$ contained in the file and returns them as an vector<double>.
vector<double>  read_fits_cl(char * filename);

//! Upgrade/Degrade an \b HEALPix map. The function returns an \b HEALPix map with an nside equal to nside_out.
THealpixMap ud_grade_healpix(const THealpixMap & map_in, int nside_out);

/*! 
  The anafast wrap: writes a temporary file with input map, launches anafast up to lmax on it, reads the output 
  power spectrum, removes the temporary files and returns the power spectrum.
*/
vector<double> anafast_healpix(const THealpixMap & map, int lmax);

/*! 
  The second anafast wrap: it does the same as the previous one but returns the \f$ a_{\ell m} \f$ of the input 
  map as vector<vector<complex<double> > >.
*/
vector<vector<complex<double> > > anafast_alm_healpix(const THealpixMap & map, int lmax);

/*! 
  The synfast wrap: writes a temporary file with the input \f$ C_\ell \f$, launches synfast which determines 
  \f$ a_{\ell m} \f$ with random phases and variance equal to the input \f$ C_\ell \f$. It also optionnaly 
  applies a beam transfer function given by beam (FWHM in arcminutes). The \f$ a_{\ell m} \f$ are then 
  transformed by synfast into a map that is read and returned as a THealpixMap. The temporary files are removed.
*/
THealpixMap synfast_healpix(int nside, int lmaxc, const vector<double> & ClT, double beam);

/*! 
  The second synfast wrap : this one just takes \f$ a_{\ell m} \f$ as an input, optionnaly applies the beam 
  (FWHM in arcminutes), writes the \f$ a_{\ell m} \f$ to a temporary file. Synfast is then invoked to transform 
  these \f$ a_{\ell m} \f$ into a FITS map (no random phases are drawn there of course) that is read and returned 
  after temporary files have been removed.
*/
THealpixMap synfast_alm_healpix(int nside, int lmaxc, const vector<vector<complex<double> > > & alms, double beam);

/*! 
  This function computes the legendre transform of the numerical input lobe whose abscissa and values are given by 
  theta and lobe up to lmax. This legendre transform is the one that has to be applied to the \f$ a_{\ell m} \f$ to 
  filter the map with the lobe. The legendre polynomials are determined within the function using the famous 
  recurrence formula. This allows to measure all in once all the legendre polynomials up to lmax (of course this is 
  not economic in terms of memory but it's fast !).
*/
vector<double> legendrelobe(const vector<double> & theta, const vector<double> & lobe, long lmax);

#endif
