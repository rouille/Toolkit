#ifndef CL_H
#define CL_H

#include <vector>
#include <string>

#include "events.h"
#include "healpixmap.h"

// ROOT
#include "TMatrixD.h"

using namespace std;

/*! \file Cl.h 
  This file gathers all the functions we need for computing the angular power spectrum from a set of events. 
  It has been shown that the power spectrum of a sky observed with a varying exposure is related to the underlying 
  one through a convolution operation. We use a method introduced in CMB data analysis to compute the convolution 
  kernel.
*/




//! Computation of the coupling matrix.
TMatrixD compute_Mll(vector<unsigned int>& wl, vector<double>& wCl, unsigned int lmax);

//! Computes the angular power spectrum with maps (partial exposure of the sky).
vector<double> compute_Cl(const vector<TEvent>& events, const THealpixMap& covmap, const THealpixMap& evtmap, unsigned int lmax, const vector<unsigned int>& lbins, vector<double>& ErrorCl);

//! Wigner 3j symbol initialization.
vector<double> wigner_init(unsigned int lmax);

//! Wigner 3j symbol computation.
vector<double> wigner_3j(unsigned int l1, unsigned int l2, const vector<unsigned int>& l3, const vector<double>& lngammaW);

//! Get the window function of the pixels.
vector<double> GetPixWindow(unsigned int nside, unsigned int lmax);

//! Binning of the power spectrum in l in order to reduce the correlations of the Cl-s induced by the cut sky, and 
//! also to reduce the errors on the resulting power spectrum estimator.
void makePQ(const vector<unsigned int>& lbins, unsigned int lmax, TMatrixD& matP, TMatrixD& matQ);

//! This function returns the center of the l bin and the error on l according to any binning used.
vector<vector<double> > lvalues(const vector<unsigned int>& lbins);

//! The variance of the angular power spectrum estimate.
TMatrixD covariance_Cltilde(const vector<unsigned int>& wl, const vector<double>& wCl, unsigned int lmax, double f1, unsigned int nevt);

//! The field we want to measure is obtained from the directions in the sky of a set of recorded events. We 
//! therefore only have access to a Poisson sampling of the field which induces an extra term in the power 
//! spectrum (GAP2004-007), the noise bias. This function computes the noise bias.
TMatrixD noise_bias(TMatrixD matP, double f1, unsigned int nevt, unsigned int lmax, unsigned int nbin);

#endif
