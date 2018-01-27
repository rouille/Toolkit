#include <vector>
#include <string>
#include <iostream>
#include <sys/stat.h>
#include "common.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"


bool CheckFile(string fileName)
{
  char fileNameCopy[100];
  struct stat fileStat;
  strcpy(fileNameCopy,fileName.c_str());
  int status = lstat(fileNameCopy,&fileStat);

  if( status == -1 )
    {
      cout << "-----------------------------------------------------" << endl;
      cout << "File " << fileName << " does not exist." << endl;
      cout << "-----------------------------------------------------" << endl;
      return false;
    } 
  cout << "-----------------------------------------------------" << endl;
  cout << "Checking " << fileName << " (" << fileStat.st_size << " bytes)." << endl;
  if( fileStat.st_size == 0 )
    {
      cout << "file size = 0 !!" << endl << endl;
      return false;
    }
  cout << "-----------------------------------------------------" << endl;
  return true;
}



const int kNFields = 12;
const string kEventFields[kNFields] =
  {
    "1 Event Id",
    "2 Theta",
    "3 Error on Theta",
    "4 Phi",
    "5 Error on Phi",
    "6 Galactic longitude",
    "7 Galactic latitude",
    "8 Right Ascension",
    "9 Declination",
    "10 UTC time",
    "11 TCore",
    "12 Energy"
  };



void DumpFields()
{
  cout << "------------------------------------------------------------------------------" << endl;
  cout << "Data fields in event file" << endl;
  cout << "------------------------------------------------------------------------------" << endl;
  int count = 0;
  for(int i = 0; i < kNFields; i++)
    {
      printf("%-30s,", kEventFields[i].c_str());
      count++;
      if( count == 2 ) {cout << endl; count = 0;}
    }
  cout << endl << "-------------------------------------------------------------------------------" << endl;
}



// Code from Xavier Bertou
TProgressBar::TProgressBar() {Zero();}



void TProgressBar::Zero()
{
  fCurrentPercentage = 0;
  fReached = 0;
  fBegin = 0;
  fEnd = 0;
}



void TProgressBar::InitPercent()
{
  fprintf(stderr,"[                                                  ]\r[");
  fCurrentPercentage = 0;
  fReached = 0;
}



void TProgressBar::EndPercent() const
{
  fprintf(stderr,"\r[##################################################]\n");
}



void TProgressBar::PrintPercent(unsigned int value)
{
  fReached=value;
  unsigned int newper=fReached*100/(fEnd-fBegin);
  // step of 2%
  while( newper>fCurrentPercentage+2 )
    {
      fprintf(stderr,"#");
      fflush(stderr);
      fCurrentPercentage+=2;
    }
}



template<typename T> vector< vector<unsigned int> > GetContiguousPoints(unsigned int size, const T *xx, T ignore)
{
  unsigned int xi(0), xstart;
  vector< vector<unsigned int> > good;
  // The idea is to find in xx contiguous portions of data.These portions are separated by the value ignore. 
  // For instance, for : 0 1 2 3 -1 4 5 6 -1 7 8 9
  // the function return 3 vectors of indices corresponding to the elements 0 1 2 3 then 4 5 6 and finally 7 8 9.
  // This function is used to correctly draw the isolongitude and isolatitude lines. Indeed, if you try to connect 
  // points lying inside your map with points lying outside then you will obtain large staright lines. This is what 
  // you will if you try to connect the -90 latitude with the 90 latitude
  // find first == ignore  
  while( xx[xi] != ignore && xi < size ) xi++;
  if( xi == size )
    {
      // only good data
      vector<unsigned int> tmp(size);
      for( unsigned int i = 0; i < size; i++ ) tmp[i] = i;
      // then return all indices as they are in the original vector
      good.push_back(tmp);
      return good;
    }
  xstart = xi;
  xi++; // skip detected ignore point
  while( xi != xstart ) // have to scan the complete vector until we are back to first ignore
    {
      // find first != ignore
      while( xx[xi] == ignore && xi != xstart ) xi=(xi+1+size)%size;
      // find next == ignore
      vector<unsigned int> contiguous;
      while( xx[xi] != ignore && xi != xstart )
        {
          contiguous.push_back(xi);
          xi = (xi+1+size)%size;
        }
      if( contiguous.size() ) good.push_back(contiguous);
    }
  return good;
}



void GetHisto(TH1F* Histo, const vector<double>& variable)
{
  // Fill with input values
  unsigned int sz = variable.size();
  if( !sz )
	{
		cout << "GetHisto: No data to make histogram. Returning." << endl;
		return;
	}

  for(unsigned int i = 0; i < sz; i++) Histo->Fill(variable[i]);
  Histo->Sumw2(); Histo->Scale(1./Histo->Integral());
}


void DrawHisto(TCanvas* cHisto, TH1F* Histo, string Xaxis, string Yaxis, string save)
{
  // X axis options
  Histo->SetXTitle(Xaxis.c_str());
  Histo->GetXaxis()->SetTitleFont(132);
  Histo->GetXaxis()->SetLabelSize(0.03);
  Histo->GetXaxis()->SetLabelFont(132);
  Histo->GetXaxis()->SetTitleSize(0.035);
  Histo->GetXaxis()->SetTitleOffset(1.1);
  // Y axis options
  Histo->SetYTitle(Yaxis.c_str());
  Histo->GetYaxis()->SetTitleFont(132);
  Histo->GetYaxis()->SetLabelSize(0.03);
  Histo->GetYaxis()->SetLabelFont(132);
  Histo->GetYaxis()->SetTitleSize(0.035);
  Histo->GetYaxis()->SetTitleOffset(1.4);
  // General
  Histo->SetStats(0);
  Histo->SetMarkerStyle(kFullCircle);
  Histo->SetMarkerSize(1.);

  if(save != "")
    { 
      Histo->Draw("e1p");
      cHisto->Update();
      cHisto->SaveAs(save.c_str());
    }
}


void PlotXY(TCanvas* cPlot, TGraphErrors* Plot, double xMin, double xMax, string name, string Xaxis, string Yaxis)
{
  Plot->SetTitle(name.c_str());
  // X axis options
  Plot->GetHistogram()->SetXTitle(Xaxis.c_str());
  Plot->GetHistogram()->SetAxisRange(xMin,xMax);
  Plot->GetXaxis()->SetTitleFont(132);
  Plot->GetXaxis()->SetLabelSize(0.03);
  Plot->GetXaxis()->SetLabelFont(132);
  Plot->GetXaxis()->SetTitleSize(0.035);
  Plot->GetXaxis()->SetTitleOffset(1.1);
  // Y axis options
  Plot->GetHistogram()->SetYTitle(Yaxis.c_str());
  Plot->GetYaxis()->SetTitleFont(132);
  Plot->GetYaxis()->SetLabelSize(0.03);
  Plot->GetYaxis()->SetLabelFont(132);
  Plot->GetYaxis()->SetTitleSize(0.035);
  Plot->GetYaxis()->SetTitleOffset(1.4);
  // General
  Plot->SetMarkerStyle(kFullCircle);
  Plot->SetMarkerColor(kBlack);
  Plot->SetMarkerSize(1.);
}


void PaletteOrange(int nContours)
{
  const int nStops = 3;

  double stops[nStops] = { 0.00, 0.5, 1.00 };
  double r[nStops]     = { 1.00, 1.00, 0.00 };
  double g[nStops]     = { 1.00, 0.50, 0.00 };
  double b[nStops]     = { 1.00, 0.00, 0.00 };
  TColor::CreateGradientColorTable(nStops, stops, r, g, b, nContours);
  gStyle->SetNumberContours(nContours);
}


void PaletteBlue(int nContours)
{
  const int nStops = 3;

  double stops[nStops] = { 0.00, 0.5, 1.00 };
  double r[nStops]     = { 1.00, 0.25, 0.00 };
  double g[nStops]     = { 1.00, 0.25, 0.00 };
  double b[nStops]     = { 1.00, 1.00, 0.00 };
  TColor::CreateGradientColorTable(nStops, stops, r, g, b, nContours);
  gStyle->SetNumberContours(nContours);
}


void PaletteRGB(int nContours)
{
  const int nStops = 5;

  double stops[nStops] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  double r[nStops]     = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  double g[nStops]     = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  double b[nStops]     = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(nStops, stops, r, g, b, nContours);
  gStyle->SetNumberContours(nContours);
}


void PaletteBlueAndRed(int nContours)
{
  const int nStops = 3;

  double stops[nStops] = { 0.00, 0.5, 1.00 };
  double r[nStops]     = { 0.00, 1.00, 1.00 };
  double g[nStops]     = { 0.00, 1.00, 0.00 };
  double b[nStops]     = { 1.00, 1.00, 0.00 };
  TColor::CreateGradientColorTable(nStops, stops, r, g, b, nContours);
  gStyle->SetNumberContours(nContours);
}


//! Instanciation of templates
template vector< vector<unsigned int> > GetContiguousPoints(unsigned int size, const double * xx, double ignore_value);
template vector< vector<unsigned int> > GetContiguousPoints(unsigned int size, const float * xx, float ignore_value);
template vector< vector<unsigned int> > GetContiguousPoints(unsigned int size, const int * xx, int ignore_value);


