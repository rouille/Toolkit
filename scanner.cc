#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <algorithm>

#include "agntools.h"
#include "scanner.h"
#include "userfcn.h"

// ROOT
#include "TRint.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TMath.h"

using namespace std;

static TProgressBar gProgS;
static TStopwatch gTimerS;


TScanner::TScanner(const vector<TEvent> & events, const vector<TAgn> & sources, const THealpixMap & cMap)
{
  fEvents = events;
  fNevents = fEvents.size();
  fAGNs = sources;
  fCovMap = cMap;
  Init();
}



void TScanner::Init()
{
  // Related to coverage map
  fNbPix = fCovMap.size();
  vector<long> ip;
  for(unsigned int i = 0; i < fNbPix-1; i++) ip.push_back(i);
  fCovMap.GiveLB(ip, fLpix, fBpix);
  ip.clear();
  
  // Related to the the scan
  fMinBigP = 1.; 
  fEminAtMin = 0.;
  fZmaxAtMin = 0.;
  fPsiMaxAtMin = 0.;
  fLittlePAtMin = 0.;
  fLittleKAtMin = 0;
  fBigNAtMin = 0;

  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "                  New TScanner Initialized                   " <<endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}



double TScanner::littleP(double zMaxPar, double psiMaxPar)
{
  cout << "Computing the probability that a CR isotropically distributed falls within " << psiMaxPar << " deg of an AGN with z < " << zMaxPar << endl;
  gProgS.Zero();
  gProgS.fBegin = 0;
  gProgS.fEnd = fNbPix;
  gProgS.InitPercent();
  gTimerS.Start();
  int stepNum = 0;
  double lp = 0.;
  for(unsigned int nPix = 0; nPix < fNbPix; nPix++)
    {
      if(fCovMap[nPix] > 0.) lp += fCovMap[nPix] * IsSourceNearTest(fLpix[nPix], fBpix[nPix], zMaxPar, psiMaxPar); 
      stepNum++;
      gProgS.PrintPercent(stepNum);
    }
  gProgS.EndPercent();
  gTimerS.Stop();
  gTimerS.Print();
  
  return lp;
}



int TScanner::littleK(double zMaxPar, double psiMaxPar, double eMinPar)
{
  int lk = 0;
  for(unsigned int nEvt = 0; nEvt < fNevents; nEvt++)
    if(fEvents[nEvt].fEnergy >= eMinPar) lk += IsSourceNearTest(fEvents[nEvt].fL, fEvents[nEvt].fB, zMaxPar, psiMaxPar);

  return lk;
}



int TScanner::bigN(double eMinPar)
{
  int lk = 0;
  for(unsigned int nEvt = 0 ; nEvt < fNevents; nEvt++)
    {
      double energy = fEvents[nEvt].fEnergy;
      if(energy >= eMinPar) lk++;
    }
  return lk;
}



// Public TScanner Members
void TScanner::SetRangeEnergy(double eMin, double eMax)
{
   bool grt2lst = true;
   vector<TEvent> eventsSorted = SortCutEnergy(fEvents, grt2lst, eMin, eMax);

   fEvents.resize(eventsSorted.size());
   for(unsigned int i = 0; i < eventsSorted.size(); i++) fEvents[i] = eventsSorted[i];
   fNevents = fEvents.size();
   fEmax = fEvents[0].fEnergy;
   fEmin = fEvents[fNevents-1].fEnergy;

   eventsSorted.clear();
}



void TScanner::SetRangeRedshift(double zMin, double zMax, double zStep)
{
  bool grt2lst = false;
  vector<TAgn> AGNsorted = SortCutRedShift(fAGNs, grt2lst, zMin, zMax);

  fAGNs.resize(AGNsorted.size());
  for(unsigned int i = 0;  i < AGNsorted.size(); i++) fAGNs[i] = AGNsorted[i];
  fNagn = fAGNs.size();
  fZmin = zMin;
  fZmax = zMax;
  fZstep = zStep;

  double zMinTmp = fZmin;
  while( zMinTmp <= fZmax )
    {
      zMinTmp += fZstep;
      fZmins.push_back(zMinTmp);
    }
  fNzStep = fZmins.size();
}



void TScanner::SetRangeAngle(double psiMin, double psiMax, double psiStep)
{
   fPsiMin = psiMin;
   fPsiMax = psiMax;
   fPsiStep = psiStep;

   double psiMinTmp = fPsiMin;
   while( psiMinTmp <= psiMax )
     {
       psiMinTmp += fPsiStep;
       fPsis.push_back(psiMinTmp);
     }
   fNpsi = fPsis.size();
}



void TScanner::PrintRanges()
{
  printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  printf("                  Scan Ranges                               \n");
  printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  printf("%3d events :: %2.3f <= E/EeV <= %2.3f :: %3d steps\n", fNevents, fEmin, fEmax, fNevents);
  printf(" %3d AGN   ::  %1.3f <=   Z   <=  %1.3f :: %3d steps of dZ = %1.3f\n", fNagn,  fZmin, fZmax, fNzStep, fZstep);
  printf("%3d steps  ::  %1.3f <=  Psi  <=  %1.3f :: %3d steps of dPsi = %1.3f\n", fNpsi, fPsiMin, fPsiMax, fNpsi, fPsiStep);
  printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  
}



void TScanner::SetScanResults(double eMinAtMin, double zMaxAtMin, double psiMaxAtMin, double &minBigP, double &littlePAtMin, int &littleKAtMin, int &bigNAtMin) 
{
  fEminAtMin = eMinAtMin;
  fZmaxAtMin = zMaxAtMin;
  fPsiMaxAtMin = psiMaxAtMin;
  fLittlePAtMin = littleP(fZmaxAtMin, fPsiMaxAtMin); //LittlePAtMin;
  fLittleKAtMin = littleK(fZmaxAtMin, fPsiMaxAtMin, fEminAtMin); // LittleKAtMin;
  fBigNAtMin = bigN(fEminAtMin); // BigNAtMin;
  fMinBigP = bigP(fBigNAtMin, fLittleKAtMin, fLittlePAtMin); // MinBigP;
  printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  printf("                 Set Scan Results                               \n");
  printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  printf(" big-P     |  Emin | Zmax  | Psimax | little-p | k  | N  \n");
  printf("------------------------------------------------------------\n");
  printf(" %1.3e | %2.2f | %1.3f | %2.3f  |   %2.2f   | %2d | %2d \n", fMinBigP, fEminAtMin, fZmaxAtMin, fPsiMaxAtMin, fLittlePAtMin, fLittleKAtMin, fBigNAtMin);
  printf("------------------------------------------------------------\n");
  littlePAtMin = fLittlePAtMin;
  littleKAtMin = fLittleKAtMin;
  bigNAtMin = fBigNAtMin;
  minBigP = fMinBigP;
}



void TScanner::PerformScan(bool timerOn)
{
  fMinBigP = 1.0;
  unsigned int TotNumSteps = fNevents*fNzStep*fNpsi;
  if(timerOn == true)
    {
      cout << "  Perfroming scan ... " << endl;
      gProgS.Zero();
      gProgS.fBegin = 0;
      gProgS.fEnd = TotNumSteps;
      gProgS.InitPercent();
      gTimerS.Start();
    }
  
  int stepNum = 0;
  for(unsigned int SrcN = 0; SrcN < fNzStep; SrcN++)
    {
      for(unsigned int PsiN = 0; PsiN < fNpsi; PsiN++)
        {
          double lP = littleP(fZmins[SrcN], fPsis[PsiN]);
          for(unsigned int EvN = 0; EvN < fNevents; EvN++)
            {
              int lK = littleK(fZmins[SrcN], fPsis[PsiN], fEvents[EvN].fEnergy);
              double bP = bigP(EvN+1, lK, lP);
              if(bP <= fMinBigP)
                {
                  fMinBigP = bP;
                  fEminAtMin = fEvents[EvN].fEnergy;
                  fZmaxAtMin = fZmins[SrcN];
                  fPsiMaxAtMin = fPsis[PsiN];
                  fLittlePAtMin = lP;
                  fLittleKAtMin = lK;
                  fBigNAtMin = EvN+1;
                }
              stepNum++;
              if(timerOn == true) gProgS.PrintPercent(stepNum);
            }
        }
    }
  if(timerOn == true)
    {
      gProgS.EndPercent();
      gTimerS.Stop();
      gTimerS.Print();
    }
}



void TScanner::PrintScanResults(double &eMinAtMin, double &zMaxAtMin, double &psiMaxAtMin, double &minBigP, double &littlePAtMin, int &littleKAtMin, int &bigNAtMin) 
{
  minBigP = fMinBigP;
  eMinAtMin = fEminAtMin;
  zMaxAtMin = fZmaxAtMin;
  psiMaxAtMin = fPsiMaxAtMin;
  littlePAtMin = fLittlePAtMin;
  littleKAtMin = fLittleKAtMin;
  bigNAtMin = fBigNAtMin;
  printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  printf("                 Scan Results                               \n");
  printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  printf(" big-P     | Emin | Zmax  | Psimax | little-p | k  | N  \n");
  printf("------------------------------------------------------------\n");
  printf(" %1.3e | %2.2f | %1.3f | %2.3f  |   %2.2f   | %2d | %2d \n", fMinBigP, fEminAtMin, fZmaxAtMin, fPsiMaxAtMin, fLittlePAtMin, fLittleKAtMin, fBigNAtMin);
  printf("------------------------------------------------------------\n");
}



vector<TEvent> TScanner::GetEventsThatCorrelate()
{
  vector<TEvent> eventsThatCorrelate;
  for(unsigned int i = 0; i < fNevents; i++)
    {
      if(fEvents[i].fEnergy >= fEminAtMin)
        {
          int isNear = IsSourceNearTest(fEvents[i].fL, fEvents[i].fB, fZmaxAtMin, fPsiMaxAtMin);
          if(isNear == 1) eventsThatCorrelate.push_back(fEvents[i]);
        }
    }
  return eventsThatCorrelate;
}



vector<TEvent> TScanner::GetEventsThatDoNotCorrelate()
{
  vector<TEvent> eventsThatDoNotCorrelate;
  for(unsigned int i=0; i<fNevents; i++)
    {
      if(fEvents[i].fEnergy >= fEminAtMin)
        {
          int isNear = IsSourceNearTest(fEvents[i].fL, fEvents[i].fB, fZmaxAtMin, fPsiMaxAtMin);
          if(isNear == 0) eventsThatDoNotCorrelate.push_back(fEvents[i]);
        }
    }
  return eventsThatDoNotCorrelate;
}



vector<TAgn> TScanner::SourcesCloserZmax()
{
  vector<double> lSources, bSources;
  for(unsigned int i = 0; i < fNagn; i++)
    {
      lSources.push_back(fAGNs[i].fL);
      bSources.push_back(fAGNs[i].fB);
    }  
  vector<double> cMapValAtAGNs = fCovMap.Values(lSources, bSources);   

  vector<TAgn> AGNs;
  for(unsigned int i = 0; i < fNagn; i++)
    {
      bool keep = (fAGNs[i].fZ <= fZmaxAtMin) && (cMapValAtAGNs[i] > 0.);
      if( keep ) AGNs.push_back(fAGNs[i]);
    }
  return AGNs;
}



void TScanner::PrintEventsAndNearAGN()
{
  cout << "-----------------------------------------------------------------------" << endl;
  cout << "    EventsID   | E/EeV |  galL  |  galB  ||  Psi  |   Z   | Type" << endl;
  cout << "-----------------------------------------------------------------------" << endl;

  for(unsigned int nEvt = 0; nEvt < fNevents; nEvt++)
    {
      if(fEvents[nEvt].fEnergy >= fEminAtMin)
        {
          double minAngDist = 180.0;
          unsigned int index = 0;
          for(unsigned int nSources = 0; nSources < fNagn; nSources++)
            {
              if(fAGNs[nSources].fZ <= fZmaxAtMin)
                {
                  double angDist = AngularDistance(fEvents[nEvt].fRa, fEvents[nEvt].fDec, fAGNs[nSources].fRa, fAGNs[nSources].fDec);
                  if(angDist <= minAngDist)
                    { 
                      minAngDist = angDist; 
                      index = nSources;
                    }
                }
            }
          if(minAngDist <= fPsiMaxAtMin)
            {
              cout.precision(4);
              cout << "DO CORRELATE" << endl;
              cout << "  " << fEvents[nEvt].fId << " | " << fEvents[nEvt].fEnergy << " | " << fEvents[nEvt].fL << " | " << fEvents[nEvt].fB << " || " << minAngDist << " | " << fAGNs[index].fZ << " | "<< fAGNs[index].fName << endl;
              cout << "-----------------------------------------------------------------------" << endl;
            }
          else
            {
              cout << "DO NOT CORRELATE" << endl;
              cout.precision(4);
              cout << "  " << fEvents[nEvt].fId << " | " << fEvents[nEvt].fEnergy << " | " << fEvents[nEvt].fL << " | " << fEvents[nEvt].fB << " || " << minAngDist << " | " << fAGNs[index].fZ << " | "<< fAGNs[index].fName << endl;
              cout << "-----------------------------------------------------------------------" << endl;
            }
        }
    }
}



int TScanner::IsSourceNearTest(double lTest, double bTest, double zMaxPar, double psiMaxPar)
{
  int isNear = 0;
  for(unsigned int i = 0; i < fNagn; i++)
    {
      if(fAGNs[i].fZ <= zMaxPar)
        {
          double angDist = AngularDistance(lTest, bTest, fAGNs[i].fL, fAGNs[i].fB);
          if(angDist <= psiMaxPar) {isNear = 1; break;}
        }
    }
  return isNear;
}



int TScanner::GetBigN(double eMinPar)
{
  return bigN(eMinPar);
}



double TScanner::GetLittleP(double zMaxPar, double psiMaxPar)
{
  return littleP(zMaxPar, psiMaxPar);
}



int TScanner::GetLittleK(double zMaxPar, double psiMaxPar, double eMinPar) 
{
  return littleK(zMaxPar, psiMaxPar, eMinPar);
}

