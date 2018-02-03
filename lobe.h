#ifndef _LOBE_H_
#define _LOBE_H_

#include "common.h"

#include <vector>
#include <string>


/*!
  This calss computes a beam window that will be used to integrate events and coverage maps.
*/

//! Beam window.
class TLobe
{
  public:
  //! Constructor
  TLobe(string type, double size, unsigned int nPoints);

  //! Computes the lobe
  void ComputeLobe();

  //! Returns #fTheta and #fLobe
  void GetLobe(vector<double> & theta, vector<double> & lobe) const {theta = fTheta; lobe = fLobe;}

  //! Returns #fType
  string GetType() const {return fType;}

  //! Returns #fSize
  double GetSize() const {return fSize;}
  
  //! Returns #fNpoints
  unsigned int GetNpoints() const {return fNpoints;}

  //! Write file. The file will contain 2 columns: #fTheta and the lobe value normalized to one at maximum. 
  void WriteFile(string fileName = "") const;

  private:
  //! Type of the lobe. Could be either Gaussian or circular
  string fType;

  //! Check the type of the lobe.
  void CheckType(string type);

  //! Size of the beam in degree. For a gaussian lobe the size corresponds to the standard deviation
  double fSize;

  //! Number of tabulation
  unsigned int fNpoints;

  //! Opening angle
  vector<double> fTheta;

  //! Lobe
  vector<double> fLobe;
};

#endif
