#ifndef _TIMEMOD_H_
#define _TIMEMOD_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include "STClibrary.h"
#include "angdist.h"

using namespace std;

/*! 
  This class is dedicated to the temporal distribution of the events. It is usefull when one wants to take this 
  modulation into account in the computation of the coverage map
*/

//! Temporal distribution of the events.
class TTimeModulation
{
 public:
  //! Constructor.
  TTimeModulation(); 

  //! Set #fAccTimeModel.
  void SetAccTimeModel(string model) {fAccTimeModel = model;}

  //! Returns #fAccTimeModel.
  string GetAccTimeModel() const {return fAccTimeModel;}

  //! Returns #fUTCs.
  vector<double> GetUTCs() const {return fUTCs;}

  //! Returns #fAccTime.
  vector<double> GetAccTime() const {return fAccTime;}

  /*! 
    This function computes the acceptance using distributions or files. If model = "UTC" then the acceptance is 
    computed using an interpolation of the UTC distribution. If model = "UTC+JD" then the JD variation are also 
    taken into account. If model = "filename" then the acceptance is read from a file with two columns : time 
    (UTC time) & acceptance (any relative value).
  */
  void ComputeAccTime(const vector<double> & UTCh, const vector<double> & UTCs, double timeStep, double chi2lim);

 private:
  //! Acceptance model. Either TIME_FLAT, UTC, UTC+JD or a file.
  string fAccTimeModel;

  //! UTCs regularly gridded.
  vector<double> fUTCs;

  //! Acceptance. Any relative value corresponding to #fUTCs.
  vector<double> fAccTime;
};

#endif
