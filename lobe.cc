#include <iostream>
#include <fstream>
#include <cmath>

#include "lobe.h"
using namespace std;


TLobe::TLobe(string type, double size, unsigned int nPoints)
{
  CheckType(type);
  fType = type;
  fSize = size;
  fNpoints = nPoints;
  fTheta.resize(nPoints);
  fLobe.resize(nPoints);
}



void TLobe::CheckType(string type)
{
  if( type != "gauss" && type != "circ")
  {
    cout << "type has to be set to either gauss (gaussian) or circ (top hat). Returning" << endl;
    exit(0);
  }
}



void TLobe::ComputeLobe()
{
  for(unsigned int i = 0; i < fNpoints; i++)
    {
      if(fType == "gauss")
        {
          fTheta[i] = i*90./(fNpoints-1);
          fLobe[i] = exp(-fTheta[i]*fTheta[i]*1./(2.*fSize*fSize));
        }
      else if(fType == "circ")
        {
          fTheta[i] = i*90./(fNpoints-1);
          if( fTheta[i] <= fSize ) fLobe[i] = 1.;
          else fLobe[i] = 0.;
        }
    }
}



void TLobe::WriteFile(string fileName) const
{
  char name[128];
  if( fileName == "") sprintf(name,"lobe_%s%3.1fdeg.txt", fType.c_str(), fSize);
  else strcpy(name,fileName.c_str());
  ofstream LobeFile(name);
  double lobeTmp = 0.;
  for(unsigned int i = 0; i < fNpoints; i++)
    {
      lobeTmp = fLobe[i];
      if(fLobe[i] < 1.e-20) lobeTmp = 0.;
      LobeFile << fTheta[i] << " " << lobeTmp << endl;
    }
  LobeFile.close();
}
