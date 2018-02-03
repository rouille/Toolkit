#ifndef _PROJMAP_H_
#define _PROJMAP_H_

#include <vector>
#include <string>

#include "healpixmap.h"
#include "events.h"

#include "TList.h"
#include "TCanvas.h"


/*! \file projMap.h
 You will find in this file everything related to the projection of the sky map.
 */

//! Defines the style of the lines that will be used to plot objects on a map. 
struct TDrawLine
{
	//! Number of points in the line.
	unsigned int fNPts;
	
	//! Style of the line.
	int fStyle;
	
	//! Color of the line.
	int fColor;
	
	//! Width of the line.
	int fWidth;
};

//! Projection of the sky map.
class TProjMap
{
public:
	//! Constructor for Mollweide projection.
	TProjMap(const THealpixMap& map, int sizeX, int sizeY, double decLimit = 90.);
	
	//! Constructor for Lambert azimutal projection.
	TProjMap(const THealpixMap& map, int sizeX, int sizeY, double longCenter, double latCenter, double radius, double decLimit = 90.);
	
	//! Destructor.
	~TProjMap();
	
	//! Color table
	void SetPalette(void (*PaletteFunc)(int), int nContours);
	
	//! Plot the map.
	void SkyMap(string mapTitle="", double vmin = 0, double vmax = 0);
	
	//! Logarithmic bar
	void SetLogz() {fLogZ=true;}

	//! Enable the grid representing the isolatitude and isolongitude lines.
	void ShowGrid(double longStep, double latStep, int style = 2, int color = kBlack, int width = 1);
	
	//! Enable the field of view of the experiments.
	void ShowFOV(double decLimit, int style = 9, int color = kBlack, int width = 2);
	
	//! Enable the Super Galactic Plane.
	void ShowSGP(int style = 1, int color = kBlue, int width = 2);
	
	//! Plot Astrophysical objects using markers.
	void PlotSources(const vector<double> & lSources, const vector<double> & bSources, const vector<string> & nameSources, bool writeNameSources = false, int style = kOpenStar, int color = kBlack, double size = 1.);
	
	//! Plot Astrophysical objects using markers.
	void PlotSource(double lSource, double bSource, string nameSource, bool writeNameSource = false, int style = kOpenStar, int color = kBlack, double size = 1.);
	
	//! Plot events using markers.
	void PlotEvents(const vector<double> & lEvents, const vector<double> & bEvents, int style = kDot, int color = kBlack, double size = 1.);
	
	//! Plot events using markers of size proportional to their energy.
	void PlotEvents(const vector<TEvent>& events, int style = kFullCircle, int color = kBlack);
	
	//! Plot events using markers of size proportional to their energy.
	void PlotEvents(const vector<double> & lEvents, const vector<double> & bEvents, const vector<double> & eEvents, int style = kFullCircle, int color = kBlack);
	
	//! Plot a circle of radius \f$\alpha\f$ centered on each direction (l,b).
	void PlotCircle(const vector<double> & l, const vector<double> & b, double alpha, int style = 1, int color = kBlack, int width = 1);
	
	//! Plot a circle of radius \f$\alpha\f$ centered on direction (l,b).
	void PlotCircle(const double  l, const double b, double alpha, int style = 1, int color = kBlack, int width = 1);

	//! Plot Maxima. Relevant when doing a blind search.
	void PlotMaxima(const vector<long> & ipAboveThreshold, int style = kOpenCircle, int color = kBlack, double size = 2.);
	
	//! Save the map.
	void Save(string fileName);
	
	
private:
	//! Build the isolongitude and isolatitude lines.
	void MakeLongitudesLatitudes(double longStep, double latStep, const TDrawLine & drawLine);
	
	//! Draw the lines.
	void DrawRegion(const vector<double> & longitudes, const vector<double> & latitudes, const TDrawLine & drawLine);
	
	//! Is enclosed in the projected map.
	bool IsInWindow(double xVal, double yVal);
	
	//! Map.
	const THealpixMap & fMap;
	
	//! Canvas.
	TCanvas * fCanvas;
	
	//! X size of the canvas.
	int fSizeX;
	
	//! Y size of the canvas.
	int fSizeY;
	
	//! Galactic longitude of the center of the map. Relevant for the Lambert projection.
	double fLongCenter;
	
	//! Galactic latitude of the center of the map. Relevant for the Lambert projection.
	double fLatCenter;
	
	//! Projection of the map.
	vector<vector<double> > fImage;
	
	//! X minimum value of the projected map.
	double fXMin;
	
	//! X maximum value of the projected map.
	double fXMax;
	
	//! Y minimum value of the projected map.
	double fYMin;
	
	//! Y maximum value of the projected map.
	double fYMax;
	
	//! \f$ \cos \f$ of the reference latitude for the projection
	double fCosLatRef;
	
	//! \f$ \sin \f$ of the reference latitude for the projection.
	double fSinLatRef;

	bool fLogZ;
};

#endif
