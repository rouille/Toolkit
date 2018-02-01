#include "harmotools.h"
#include "fitsio.h"
#include <fstream>

vector<vector<double> >  read_fits_image(char * filename)
{
  vector<vector<double> > img;
  
  fitsfile* fptr; /* pointer to the FITS file, defined in fitsio.h */
  int status,  nfound, anynull;
  long naxes[2], fpixel, npixels;
  
  float datamin, datamax, nullval;

  status = 0;

  if ( fits_open_file(&fptr, filename, READONLY, &status) ) printerror( status );

  /* read the NAXIS1 and NAXIS2 keyword to get image size */
  char naxis[] = "NAXIS";
  if ( fits_read_keys_lng(fptr, naxis, 1, 2, naxes, &nfound, &status) ) printerror( status );

  npixels  = naxes[0] * naxes[1]; /* number of pixels in the image */
  fpixel   = 1;
  nullval  = 0; /* don't check for null values in the image */
  datamin  = 1.0E30;
  datamax  = -1.0E30;

  double* data;
  data=new double[npixels];

  if ( fits_read_img(fptr, TDOUBLE, fpixel, npixels, &nullval, data, &anynull, &status) ) printerror( status );
  
  cout << "Data read" << endl;

  img.resize(naxes[0]);
  for(long i = 0; i < naxes[0]; i++)
    {
      img[i].resize(naxes[1]);
      for(long j = 0; j < naxes[1]; j++) img[i][j]=data[i+j*naxes[0]];
    }
  delete [] data;
  return(img);
}


void write_fits_image(const vector<vector<double> > & image, char * filename)
{
  char newfilename[1000];
  sprintf(newfilename,"!%s",filename);

  fitsfile* fptr; /* pointer to the FITS file, defined in fitsio.h */
  int status;
  long  fpixel, nelements;
  
  /* initialize FITS image parameters */
  int bitpix = DOUBLE_IMG;
  long naxis = 2;  /* 2-dimensional image */  
  long naxes[2];
  naxes[0] = image.size();
  naxes[1] = image[0].size();
  
  nelements = naxes[0]*naxes[1];
  double* array[naxes[1]];
  array[0] = (double *)malloc( naxes[0] * naxes[1] * sizeof(double ) );
  for(long ii = 1; ii < naxes[1]; ii++) array[ii] = array[ii-1] + naxes[0];
  

  for (long jj = 0; jj < naxes[1]; jj++) for (long ii = 0; ii < naxes[0]; ii++) array[jj][ii] = image[ii][jj];

  status=0;
  /* create new FITS file */
  if ( fits_create_file(&fptr, newfilename, &status) ) printerror( status );
  if ( fits_create_img(fptr, bitpix, naxis, naxes, &status) ) printerror( status );          
  fpixel = 1; /* first pixel to write */
  
  /* write the array of unsigned integers to the FITS file */
  if ( fits_write_img(fptr, TDOUBLE, fpixel, nelements, array[0], &status) ) printerror( status );
  free( array[0] );
  if ( fits_close_file(fptr, &status) ) printerror( status );           
  
  return;
  
}

void write_fits_alms(const vector<vector<complex<double> > > & alms, char * filename)
 {
   fitsfile* fptr;
   int status = 0;
   char newfilename[1000];
   sprintf(newfilename,"!%s",filename);
   
   char* coltype[3]={(char*)"INDEX",(char*)"REAL",(char*)"IMAG"};
   char* colform[3]={(char*)"1J",(char*)"1E",(char*)"1E"};
   
   long sizeofl=alms.size();
   long sizeofall= sizeofl*(sizeofl+1)/2;
   
   long* index;
   double *repart,*impart;
   index = new long[sizeofall];
   repart = new double[sizeofall];
   impart = new double[sizeofall];
   
   long u = 0;
   for(long i = 0; i < sizeofl; i++)
     {
       for(long j = 0; j < i+1; j++)
         {
           long ind = i*i+i+j+1;
           index[u] = ind;
           repart[u] = alms[i][j].real();
           impart[u] = alms[i][j].imag();
           u++;
         }
     }
   
   ffinit(&fptr,newfilename,&status);
   ffcrtb(fptr,BINARY_TBL,sizeofall,0,NULL,NULL,NULL,NULL,&status);
   ffmahd(fptr,2,NULL,&status);
   fficls(fptr,1,3,coltype,colform,&status);
   
   ffpcl(fptr,TLONG,1,1,1,sizeofall,index,&status);
   ffpcl(fptr,TDOUBLE,2,1,1,sizeofall,repart,&status);
   ffpcl(fptr,TDOUBLE,3,1,1,sizeofall,impart,&status);
   
   ffclos(fptr,&status);
   
   delete [] index;
   delete [] repart;
   delete [] impart;   
}

void write_fits_cl(const vector<double> & ClT, char * filename)
{
  fitsfile *fptr;
  int status = 0;
  double * transfert;
  transfert = new double [ClT.size()];
  char* coltype[4]={(char*)"ClT",(char*)"ClE",(char*)"ClB",(char*)"ClExT"};
  char* colform[4]={(char*)"E15.7",(char*)"E15.7",(char*)"E15.7",(char*)"E15.7"};
  
  char newfilename[1000];
  sprintf(newfilename,"!%s",filename);

  ffinit(&fptr,newfilename,&status);
  ffcrtb(fptr,ASCII_TBL,ClT.size(),0,NULL,NULL,NULL,NULL,&status);
  ffmahd(fptr,2,NULL,&status);
  fficls(fptr,1,4,coltype,colform,&status);
  for (unsigned int j = 0; j < ClT.size(); j++) transfert[j] = ClT[j];
  ffpcl(fptr,TDOUBLE,1,1,1,ClT.size(),transfert,&status);
  for (unsigned int j = 0; j < ClT.size(); j++) transfert[j] = ClT[j]*0.;
  ffpcl(fptr,TDOUBLE,2,1,1,ClT.size(),transfert,&status);
  for (unsigned int j = 0; j < ClT.size(); j++) transfert[j]=ClT[j]*0.;
  ffpcl(fptr,TDOUBLE,3,1,1,ClT.size(),transfert,&status);
  for (unsigned int j = 0; j < ClT.size(); j++) transfert[j]=ClT[j]*0.;
  ffpcl(fptr,TDOUBLE,4,1,1,ClT.size(),transfert,&status);

  delete [] transfert;
  ffclos(fptr,&status);
}

void read_fits_alms(char * filename, vector<vector<complex<double> > > & alms)
{
  char comment[81];
  fitsfile *fptr;
  int status = 0, hdutyp, anynul;
  long _Nele;
  double * _donneesr;
  double * _donneesi;
 
  ffopen(&fptr,filename, 0, &status);
  ffmahd(fptr, 2, &hdutyp, &status);
  ffgky(fptr, TLONG, (char*)"NAXIS2", &_Nele, comment, &status);
  _donneesr = new double [_Nele];
  _donneesi = new double [_Nele];
  ffgcvd(fptr, 2, 1 , 1, _Nele, -9999, _donneesr, &anynul, &status);
  ffgcvd(fptr, 3, 1 , 1, _Nele, -9999, _donneesi, &anynul, &status);
 
  ffclos(fptr, &status);
  
  long sizeofl;

  sizeofl = long((sqrt(8.*_Nele+1)-1)/2);
  alms.resize(sizeofl);
  for (short i = 0; i < sizeofl; i++) alms[i].resize(i+1);
  for (short i = 0; i < sizeofl; i++) 
    for (short j = 0;j < i+1; j++) alms[i][j] = complex<double>(_donneesr[i*(i+1)/2+j],_donneesi[i*(i+1)/2+j]);
  
  delete [] _donneesr;
  delete [] _donneesi;
}

vector<double> read_fits_cl(char * filename)
{
  char comment[81];
  fitsfile *fptr;
  int status = 0, hdutyp, anynul;
  long _Nele;
  double * _donnees;
  vector<double> ClT;
 
  ffopen(&fptr,filename, 0, &status);
  ffmahd(fptr, 2, &hdutyp, &status);
  ffgky(fptr, TLONG, (char*)"NAXIS2", &_Nele, comment, &status);
  _donnees = new double [_Nele];
  ffgcvd(fptr, 2, 1 , 1, _Nele,-9999, _donnees, &anynul, &status);
  
  ffclos(fptr, &status);

  ClT.resize(_Nele);
  for (long j = 0;j < _Nele; j++) ClT[j] = _donnees[j];
  
  delete [] _donnees;
  return(ClT);
}

THealpixMap ud_grade_healpix(const THealpixMap & map_in, int nside_out)
{
  unsigned int seed;
  string healpixdir(getenv("HEALPIX_DIR"));
  struct timeval mytimeval;
  struct timezone mytimezone;
  gettimeofday(&mytimeval,&mytimezone);
  seed = (unsigned int)(mytimeval.tv_usec+(mytimeval.tv_sec % 1000)*1000000);
  srandom(seed);
  char random_str[256];
  sprintf(random_str, "%d", (int)random());
  
  // Defining filenames
  char fitsmap_in[256], fitsmap_out[256], txtfile[256];
  char* tmpbatch;
  tmpbatch = getenv("TMPBATCH");
  if (tmpbatch == NULL)
    {
      tmpbatch = new char [256];
      sprintf(tmpbatch,"/tmp");
    }
  sprintf(fitsmap_in, "%s/tmpmap_in%s.fits", tmpbatch, random_str);
  sprintf(fitsmap_out, "%s/tmpmap_out%s.fits", tmpbatch, random_str);
  sprintf(txtfile, "%s/tmptxt%s.fits", tmpbatch, random_str);

  // Writing temporary map
  map_in.WriteFits(fitsmap_in);

  // writing temporary file
  ofstream oftxt(txtfile);
  oftxt << "infile = " << fitsmap_in << endl;
  oftxt << "outfile = " << fitsmap_out << endl;
  oftxt << "nside = " << nside_out << endl;
  oftxt.close();

  // Preparing command to be launched
  char commande[1000];
  sprintf(commande, "%s/bin/udgrade_cxx %s", healpixdir.c_str(), txtfile);
  system(commande);

  // Reading map
  THealpixMap map_out(fitsmap_out);

  // removing temporary files
  sprintf(commande,"rm -f %s %s %s", fitsmap_in, fitsmap_out, txtfile);
  system(commande);
  
  // returning map
  return(map_out);
}

vector<double> anafast_healpix(const THealpixMap & map, int lmax)
{
  vector<double> thecl(lmax+1);
  string healpixdir(getenv("HEALPIX_DIR"));
  
  unsigned int seed;
  struct timeval mytimeval;
  struct timezone mytimezone;
  gettimeofday(&mytimeval,&mytimezone);
  seed=(unsigned int) (mytimeval.tv_usec + (mytimeval.tv_sec % 1000)*1000000);
  // Ini_Ranf_Quick(seed,0);
  srandom(seed);
  char random_str[256];
  sprintf(random_str,"%d",(int)random());
  // defining filenames
  char fitsmap[256],fitscl[256],txtfile[256];
  char * tmpbatch;
  tmpbatch = getenv("TMPBATCH");
  if (tmpbatch == NULL)
    {
      tmpbatch = new char [256];
      sprintf(tmpbatch,"/tmp");
    }
  sprintf(fitsmap,"%s/tmpmap%s.fits",tmpbatch,random_str);
  sprintf(fitscl,"%s/tmpcl%s.fits",tmpbatch,random_str);
  sprintf(txtfile,"%s/tmptxt%s.fits",tmpbatch,random_str);
  // writing temporary map
  map.WriteFits(fitsmap);
  int order = 3;
  // writing temporary file
  ofstream oftxt(txtfile);
  oftxt << "nlmax = " << lmax << endl;
  oftxt << "infile = " << fitsmap << endl;
  oftxt << "outfile = " << fitscl << endl;
  oftxt << "polarisation = false" << endl;
  oftxt << "iter_order = " << order << endl;
  oftxt.close();
  //preparing and launching command
  char commande[1000];
  sprintf(commande,"%s/bin/anafast_cxx %s >& /dev/null",healpixdir.c_str(),txtfile);
  system(commande);
  thecl=read_fits_cl(fitscl);
  // removing temporary filenames
  sprintf(commande,"rm -f %s %s %s",fitsmap,fitscl,txtfile);
  system(commande);  
  return(thecl);
}

vector<vector<complex<double> > > anafast_alm_healpix(const THealpixMap & map, int lmax)
{
  vector<vector<complex<double> > > thealm;
  string healpixdir(getenv("HEALPIX_DIR"));
  
  unsigned int seed;
  struct timeval mytimeval;
  struct timezone mytimezone;
  gettimeofday(&mytimeval,&mytimezone);
  seed=(unsigned int) (mytimeval.tv_usec + (mytimeval.tv_sec % 1000)*1000000);
  //Ini_Ranf_Quick(seed,0);
  srandom(seed);
  char random_str[256];
  sprintf(random_str,"%d",(int)random());
  // defining filenames
  char fitsmap[256],fitscl[256],txtfile[256],fitsalm[256];
  char * tmpbatch;
  tmpbatch = getenv("TMPBATCH");
  if (tmpbatch == NULL)
    {
      tmpbatch = new char [256];
      sprintf(tmpbatch,"/tmp");
    }
  sprintf(fitsmap,"%s/tmpmap%s.fits",tmpbatch,random_str);
  sprintf(fitscl,"%s/tmpcl%s.fits",tmpbatch,random_str);
  sprintf(txtfile,"%s/tmptxt%s.fits",tmpbatch,random_str);
  sprintf(fitsalm,"%s/tmpalm%s.fits",tmpbatch,random_str);
  // writing temporary map
	map.WriteFits(fitsmap);
	int order = 3;
  // writing temporary file
  ofstream oftxt(txtfile);
  oftxt << "nlmax = " << lmax << endl;
  oftxt << "infile = " << fitsmap << endl;
  oftxt << "outfile_alms = " << fitsalm << endl;
  oftxt << "outfile = " << fitscl << endl;
  oftxt << "polarisation = false" << endl;
  oftxt << "iter_order = " << order << endl;
  oftxt.close();
  // preparing and launching the command
  char commande[1000];
  sprintf(commande,"%s/bin/anafast_cxx %s > /dev/null",healpixdir.c_str(),txtfile);
	system(commande);
	read_fits_alms(fitsalm,thealm);
	// removing temporary files
  sprintf(commande,"rm -f %s %s %s %s",fitsmap,fitscl,txtfile,fitsalm);
	system(commande);
	
  return(thealm);
}

THealpixMap synfast_healpix(int nside, int lmaxc, const vector<double> & ClT, double beam)
{
  string healpixdir(getenv("HEALPIX_DIR"));
  unsigned int seed;
  struct timeval mytimeval;
  struct timezone mytimezone;
  gettimeofday(&mytimeval,&mytimezone);
  seed = (unsigned int)(mytimeval.tv_usec+(mytimeval.tv_sec % 1000)*1000000);
  srandom(seed);
  char random_str[256];
  sprintf(random_str,"%d",(int)random());

  // defining filenames
  char fitsmap[256],fitscl[256],txtfile[256],fitsalm[256];
  char* tmpbatch;
  tmpbatch = getenv("TMPBATCH");
  if (tmpbatch == NULL)
    {
      tmpbatch = new char [256];
      sprintf(tmpbatch,"/tmp");
    }
  sprintf(fitsmap, "%s/tmpmap%s.fits", tmpbatch, random_str);
  sprintf(fitsalm, "%s/tmpalm%s.fits", tmpbatch, random_str);
  sprintf(fitscl, "%s/tmpcl%s.fits", tmpbatch, random_str);
  sprintf(txtfile, "%s/tmptxt%s.txt", tmpbatch, random_str);
  // writing temporary Cl file
  write_fits_cl(ClT, fitscl);
  // writing temporary text file
  ofstream oftxt(txtfile);
  oftxt << "fwhm_arcmin = " << beam << endl;
  oftxt << "infile = " << fitscl << endl;
  oftxt << "rand_seed = " << 1234 << endl;
  oftxt << "nlmax = " << lmaxc << endl;
  oftxt << "outfile = " << fitsalm << endl;//!test.alm
  oftxt << "polarisation = false" << endl;
  oftxt.close();
  char commande[1000];
  sprintf(commande, "%s/bin/syn_alm_cxx %s", healpixdir.c_str(), txtfile);
  system(commande);
  // creating alms
  sprintf(commande, "rm -f %s", txtfile);
  system(commande);
  oftxt.open(txtfile);
  oftxt << "nlmax = " << lmaxc << endl;
  oftxt << "infile = " << fitsalm << endl;
  oftxt << "outfile = " << fitsmap << endl;
  oftxt << "nside = " << nside << endl;
  oftxt << "polarisation = false" << endl;
  oftxt.close();
  sprintf(commande, "%s/bin/alm2map_cxx %s > /dev/null", healpixdir.c_str(), txtfile);
  system(commande);
  // reading map
  THealpixMap map(fitsmap);
  // removing temporary files
  sprintf(commande, "rm -f %s %s %s %s", fitsmap, fitscl, txtfile, fitsalm);
  system(commande);
  // returning the map
  return(map);
}

THealpixMap synfast_alm_healpix(int nside, int lmaxc, const vector<vector<complex<double> > > & alms, double beam)
{
	string healpixdir(getenv("HEALPIX_DIR"));
  unsigned int seed;
  struct timeval mytimeval;
  struct timezone mytimezone;
  gettimeofday(&mytimeval,&mytimezone);
  seed = (unsigned int)(mytimeval.tv_usec+(mytimeval.tv_sec % 1000)*1000000);
  srandom(seed);
  char random_str[256];
  sprintf(random_str,"%d",(int)random());

  // Defining filenames
  char fitsmap[256],fitscl[256],txtfile[256],fitsalm[256];
  char * tmpbatch;
  tmpbatch = getenv("TMPBATCH");
  if (tmpbatch == NULL)
    {
      tmpbatch = new char [256];
      sprintf(tmpbatch,"/tmp");
    }
  sprintf(fitsmap, "%s/tmpmap%s.fits", tmpbatch, random_str);
  sprintf(fitsalm, "%s/tmpalm%s.fits", tmpbatch, random_str);
  sprintf(fitscl, "%s/tmpcl%s.fits", tmpbatch, random_str);
  sprintf(txtfile, "%s/tmptxt%s.txt", tmpbatch, random_str);
  // Writing temporary Cl
	write_fits_alms(alms, fitsalm);
	char commande[1000];
  // Writing temporary text file
  ofstream oftxt(txtfile);
  oftxt << "nlmax = " << lmaxc << endl;
  oftxt << "infile = " << fitsalm << endl;
  oftxt << "outfile = " << fitsmap << endl;
  oftxt << "nside = " << nside << endl;
  oftxt << "polarisation = false" << endl;
  oftxt.close();
	
  sprintf(commande, "%s/bin/alm2map_cxx %s > /dev/null", healpixdir.c_str(), txtfile);
	system(commande);

  THealpixMap map(fitsmap);
	// Removing temporary files
  sprintf(commande,"rm -f %s %s %s %s", fitsmap, fitscl, txtfile, fitsalm);
  system(commande);
  // returning the map
	return(map);
}

vector<double> legendrelobe(const vector<double> & theta, const vector<double> & lobe, long lmax)
{
  long sz = theta.size();
  vector<double> thetarad;
  thetarad.resize(sz);
  for(long i=0; i<sz; i++) thetarad[i] = theta[i]*M_PI/180;

  vector<vector<double> > poln;
  poln.resize(lmax+1);
  for(long i=0; i<lmax+1; i++) poln[i].resize(sz);

  for(long i=0; i<sz; i++) poln[0][i] = 1.;
  
  if(lmax >=1 ) for(long i=0; i<sz; i++) poln[1][i] = cos(thetarad[i]);
  if(lmax >= 2)
    {
      for(long k=2; k<lmax+1; k++)
        {
          for(long i=0; i<sz; i++)
            {
              poln[k][i] = (1./k)*((2.*k-1.)*poln[k-1][i]*cos(thetarad[i])-(k-1)*poln[k-2][i]);
            }
        }
    }

  vector<double> bl;
  bl.resize(lmax+1);
  vector<double> y,y0;
  y.resize(sz);
  y0.resize(sz);
  for(long j=0;j<sz;j++) y0[j] = lobe[j]*sin(thetarad[j]);

  for(long i=0; i<lmax+1; i++)
    {
      for(long j=0; j<sz; j++) y[j] = y0[j]*poln[i][j];
      bl[i] = integrate_nc5(thetarad,y);
    }
  
  double renorm=bl[0];
  for(long i=0; i<lmax+1; i++) bl[i] = bl[i]/renorm;
  
  return(bl);
}

