#include <iostream>
#include <fstream>

#include "harmotools.h"
#include "maptools.h"
#include "angdist.h"
#include "common.h"
#include "STClibrary.h"

#ifdef gcc323
char* operator+( std::streampos&, char* );
#endif

using namespace std;



THealpixMap map_events(unsigned int nside, const vector<double> & l, const vector<double> & b)
{
  THealpixMap map(nside);
  vector<long> ipix = map.Ip(l,b);
  for(unsigned int i = 0; i < l.size(); i++) map[ipix[i]]++;
  return map;
}



vector<TEvent> GetEventsFromLB(const vector<TEvent> & events, double l, double b, double radius)
{
  vector<TEvent> eventsKeep;
  vector<double> uv0 = ll2uv(l,b), uvEvents;
  double ang;
  for(unsigned int i = 0; i < events.size(); i++)
    {
      uvEvents = ll2uv(events[i].fL,events[i].fB);
      ang = acos(uvEvents[0]*uv0[0]+uvEvents[1]*uv0[1]+uvEvents[2]*uv0[2]);
      if( ang <= radius*DTOR ) eventsKeep.push_back(events[i]);
    }
  return eventsKeep;
}


vector<long> GetPixFromLB(const THealpixMap & map, double l, double b, double radius)
{
  vector<long> pixKeep;
  vector<long> allPix(map.NPix());
  for(unsigned int i = 0; i < map.NPix(); i++) allPix[i] = i;
  vector<double> lAllPix, bAllPix;
  map.GiveLB(allPix,lAllPix,bAllPix);

  vector<double> uv0 = ll2uv(l,b), uvPix;
  double ang;
  for(unsigned int i = 0; i < map.NPix(); i++)
    {
      uvPix = ll2uv(lAllPix[i],bAllPix[i]);
      ang = acos(uvPix[0]*uv0[0]+uvPix[1]*uv0[1]+uvPix[2]*uv0[2]);
      if( ang <= radius*DTOR ) pixKeep.push_back(allPix[i]);
    }
  return pixKeep;
}



vector<double> ll2uv(double l, double b)
{
  double lon = l*DTOR;
  double lat = b*DTOR;
  double cos_lat = cos(lat);
  double sin_lat = sin(lat);
  double cos_lon = cos(lon);
  double sin_lon = sin(lon);
  vector<double> uv(3);
  uv[0] = cos_lat * cos_lon;
  uv[1] = cos_lat * sin_lon;
  uv[2] = sin_lat;
  return uv;
}



double AngularDistance(double lSky1, double bSky1, double lSky2, double bSky2)
{
  vector<double> uv1 = ll2uv(lSky1, bSky1);
  vector<double> uv2 = ll2uv(lSky2, bSky2);
  double cosAngle, angle;
  cosAngle = uv1[0]*uv2[0]+uv1[1]*uv2[1]+uv1[2]*uv2[2];
  if( cosAngle >= 1 ) cosAngle = 1;
  if( cosAngle <= -1 ) cosAngle = -1;
  angle = acos(cosAngle)*1./DTOR;

  return angle;
}



double TP_AngularDistance(double theta1, double phi1, double theta2, double phi2)
{
  return ( (1./DTOR)*acos( sin(theta1*DTOR)*cos(phi1*DTOR)*sin(theta2*DTOR)*cos(phi2*DTOR)
			   + sin(theta1*DTOR)*sin(phi1*DTOR)*sin(theta2*DTOR)*sin(phi2*DTOR)
			   + cos(theta1*DTOR)*cos(theta2*DTOR) ) );
}



bool XYtoAngLambertAzimuthal(double x, double y, double & theta, double & phi, double phi0, double cl, double sl)
{
  double rho, c, sc, cc, val;
  bool status = true;
  double thex2=x*x;
  rho=sqrt(thex2+y*y);
  c=2*asin(rho/2);
  cc=cos(c);
  sc=sin(c);
  phi=mod(phi0+atan2(-x*sc,rho*cl*cc-y*sl*sc),TwoPi);
  val=cc*sl+y*sc*cl/rho;
  if( fabs(val)<=1 ) theta=asin(val);
  else status=false;
  return status;
}



bool AngtoXYLambertAzimuthal(double theta, double phi, double & xrad, double & yrad, double lambda0, double cp1, double sp1)
{  
  double cp=cos(theta), sp=sin(theta);
  double dphi=phi-lambda0, sdphi=sin(dphi), cdphi=cos(dphi);
  double k=sqrt(2./(1+sp1*sp+cp1*cp*cdphi));  
  xrad=-k*cp*sdphi;
  yrad=k*(cp1*sp-sp1*cp*cdphi);
  return true; // the transformation is always possible
}



bool XYtoAngMollweide(double x, double y, double & theta, double & phi, double, double, double)
{
  static double prevy = y; // initialization
  static double factor = TwoPi/sin(acos(2.*y));
  if( y != prevy )// recompute the factor
    {
      factor = TwoPi/sin(acos(2.*y));
      prevy=y;
    }
  theta = (0.5-y)*M_PI;
  phi = -x*factor;
  bool status;
  if( phi <= M_PI && phi >=( -M_PI) ) status = true;
  else status = false;
  return status;
}



bool AngtoXYMollweide(double theta, double phi, double & x, double & y, double, double, double)
{
  double factor = sin(acos(theta/PiOver2));
  if( phi>M_PI ) phi -= TwoPi;
  x = -phi*factor;
  y = theta;
  return true;
}



void xyz2thetaphi(double x, double y, double z, double & theta, double & phi)
{
  theta = 0., phi = 0.;
  double norm = sqrt(x*x+y*y+z*z);
  if( norm != 0 )
    {
      theta = acos(z/norm);
      if( x != 0. ) phi = atan2(y,x);
      if( x == 0. ) phi = M_PI/2.;
    }
  if( phi < 0.) phi = phi+TwoPi;
}



void query_disc(int nside, double theta, double phi, double radius, vector<int>& listpix)
{
  listpix.clear();
 
  unsigned int nPix = nside2npix(nside); 
  double dth1 = 4./nPix;
  double dth2 = 2*nside*dth1;
  double cosang = cos(radius);

  double z0 = cos(theta);
  double xa = 1./sqrt((1-z0)*(1+z0));

  double rlat1  = theta - radius;
  double zmax = cos(rlat1);
  int irmin = ring_above(nside,zmax)+1;

  if( rlat1 <= 0 ) // north pole in the disc
    for(int m = 1; m < irmin; ++m) // rings completely in the disc
      in_ring(nside,m,0,M_PI,listpix);
  
  double rlat2 = theta + radius;
  double zmin = cos(rlat2);
  int irmax = ring_above(nside,zmin);

  // ------------- loop on ring number ---------------------
  for(int iz = irmin; iz <= irmax; ++iz) // rings partially in the disc
    {
      double z;
      if( iz < nside ) // north polar cap
        z = 1.0 - iz*iz*dth1;
      else if(iz <= (3*nside)) // tropical band + equat.
        z = (2*nside-iz) * dth2;
      else
        z = -1.0 + (4*nside-iz)*(4*nside-iz)*dth1;
      
      // --------- phi range in the disc for each z ---------
      double x = (cosang-z*z0)*xa;
      double ysq = 1-z*z-x*x;
      double dphi = atan2(sqrt(ysq),x);
      in_ring(nside,iz,phi,dphi,listpix);
    }

  if (rlat2 >= M_PI ) // south pole in the disc
    for(int m = irmax+1; m < (4*nside); ++m)  // rings completely in the disc
      in_ring(nside,m,0,M_PI,listpix);
}


int ring_above(int nside, double z)
{
  double az = abs(z);
  double twothird = 0.6666666666666666;
  if( az > twothird ) // polar caps
    {
      int iring = int(nside*sqrt(3*(1-az)));
      return (z>0) ? iring : 4*nside-iring-1;
    }
  else // ----- equatorial region ---------
    return int(nside*(2-1.5*z));
}



void in_ring(int nside, int iz, double phi0, double dphi, vector<int> &listir) 
{
  int nr, ir, ipix1;
  double shift=0.5;
  unsigned int nPix = nside2npix(nside); 
  
  if( iz < nside ) // north pole
    {
      ir = iz;
      nr = ir*4;
      ipix1 = 2*ir*(ir-1); // lowest pixel number in the ring
    }
  else if( iz > (3*nside) ) // south pole
    {
      ir = 4*nside - iz;
      nr = ir*4;
      ipix1 = nPix - 2*ir*(ir+1); // lowest pixel number in the ring
    }
  else // equatorial region
    {
      ir = iz - nside + 1; // within {1, 2*nside + 1}
      nr = nside*4;
      if ( (ir&1) == 0 ) shift = 0;
      ipix1 = 2*(nside*nside-nside) + (ir-1)*nr; // lowest pixel number in the ring
    }
  
  int ipix2 = ipix1 + nr - 1; // highest pixel number in the ring
  
  // ----------- constructs the pixel list --------------
  if( dphi > (M_PI-1e-7) )
    for(int i = ipix1; i <= ipix2; ++i) listir.push_back(i);
  else
    {
      int ip_lo = (int) floor(nr*(1./TwoPi)*(phi0-dphi) - shift)+1;
      int ip_hi = (int) floor(nr*(1./TwoPi)*(phi0+dphi) - shift);
      int pixnum = ip_lo+ipix1;
      if( pixnum < ipix1) pixnum += nr;
      for(int i = ip_lo; i <= ip_hi; ++i, ++pixnum)
        {
          if(pixnum > ipix2) pixnum -= nr;
          listir.push_back(pixnum);
        }
    }
}

