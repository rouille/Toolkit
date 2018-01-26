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
#include "agntools.h"

// ROOT
#include "TRint.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TMath.h"

using namespace std;


vector<TAgn> Get_AGNs(char* sourceData, double zMax, double decLimit)
{
  vector<TAgn> AGN;
  ifstream sourceFile(sourceData);
  if ( !sourceFile )
    {
      cout << "Cannot open file " << sourceData << ". Returning" <<endl;
      exit(0);
    }
  cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Using " << sourceData << " for AGN list" << endl;
  cout << "Keeping only z <= " << zMax << " and dec < " << decLimit << endl;
  unsigned int nKeep = 0;
  while( !sourceFile.eof() )
    {
      double raSource, decSource, zSource, lSource, bSource; 
      string nameSource;
      sourceFile >>  raSource >> decSource >> zSource >> nameSource;
      if( !sourceFile.eof() && zSource >= 0. && zSource <= zMax && decSource < decLimit )
        {
          TAgn AGNtmp;
          AGNtmp.fZ = zSource;
          radec2gal(raSource, decSource, &lSource, &bSource);
          AGNtmp.fRa = 15.0 * raSource; // in degrees
          AGNtmp.fDec = decSource;
          AGNtmp.fL = lSource;
          AGNtmp.fB = bSource;
          AGNtmp.fName = nameSource;
	
          AGN.push_back(AGNtmp);
          nKeep++;
        }
    }
  sourceFile.close(); //Close file
  cout << "There are " << AGN.size() << " AGNs in the list." << endl;
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  
  return(AGN);
}



void Print_AGN(vector<TAgn> sources, unsigned int indexNum)
{
  cout << "\n************************************************" << endl;
  cout << "  Info On Source # " << indexNum  << endl;
  cout << "-------------------------------------------------" << endl;
  cout << "              Z = " << sources[indexNum].fZ  << endl;
  cout << "       log_10(Z) = " << TMath::Log10(sources[indexNum].fZ)  << endl;
  cout << "              Ra = " << sources[indexNum].fRa  << endl;
  cout << "              Dec = " << sources[indexNum].fDec  << endl;
  cout << "              L = " << sources[indexNum].fL  << endl;
  cout << "              B = " << sources[indexNum].fB  << endl;
  cout << "            Name = " << sources[indexNum].fName  << endl;
  cout << "************************************************" << endl;
  
}



vector<TAgn> SortCutRedShift(vector<TAgn> zUnsorted, bool down, double zMin, double zMax)
{
  int nAll = zUnsorted.size();
  vector<TAgn> zSorted;
  
  double redShift[nAll];
  for(int i = 0; i < nAll; i++) { redShift[i] = zUnsorted[i].fZ; }

  int* index = new int[nAll];
  TMath::Sort(nAll, redShift, index, down);
  
  for(int i = 0; i < nAll; i++)
    {
      double redShiftTmp = zUnsorted[index[i]].fZ;
      if(redShiftTmp >= zMin && redShiftTmp <= zMax) zSorted.push_back(zUnsorted[index[i]]);
    }
  return(zSorted);
}



vector<TEvent> SortCutEnergy(vector<TEvent> eUnsorted, bool down, double eMin, double eMax)
{
  int nAll = eUnsorted.size();
  vector<TEvent> eSorted;
  
  double energy[nAll];
  for(int i = 0; i < nAll; i++) { energy[i] = eUnsorted[i].fEnergy; }

  int* index = new int[nAll]; 
  TMath::Sort(nAll, energy, index, down);

  for(int i = 0; i < nAll; i++)
    { 
      double energyTmp = eUnsorted[index[i]].fEnergy;
      if(energyTmp >= eMin && energyTmp <= eMax) eSorted.push_back(eUnsorted[index[i]]);
    }
  return(eSorted);
}

