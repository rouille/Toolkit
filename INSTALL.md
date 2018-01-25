# COVERAGE & ANISOTROPY TOOLKIT INSTALLATION NOTES
February 25th, 2010

You need various softwares to be able to compile our code :
- CFITSIO: the FITS I/O library (C language). This is a world standard for astronomy binary files. It basically consists in an ASCII header describing the rest of the file that is binary. You therefore just have to use the FITSIO routines to access files, without a prior knowledge of the file content, which is cool. HEALPix uses this standard.
- HEALPix: this is the sky pixelization that we use in our *COVERAGE & ANISOTROPY TOOLKIT*. It is the standard in the CMB community and has many advantages : iso-surface pixels regularly spaced at a given latitude allowing Harmonic transform to be almost replaced by Fast Fourier Transforms.
- ROOT: well you know it ...
- STCoordinates: part of the CDAS software where coordinate system functions are coded.
- Our *COVERAGE & ANISOTROPY TOOLKIT*


## Install HEALPix and CFITSIO
You can get the files from our web page but we recommend that you get it directly from HEALPix people by registering at the following URL: https://sourceforge.net/projects/healpix. Once you have uncompressed HEALPix (`tar -xzvf Healpix_2.13a.tar.gz`), define an environment variable in your .cshrc, .bashrc or .whatevershrc named HEALPIX_DIR giving the full path of the HEALPix directory. For instance for bash:
```
export HEALPIX_DIR=/Users/yourself/.../Healpix_2.13a
```
You need to copy the file cfitsio3240.tar.gz in `$HEALPIX_DIR/src/cxx/libcfitsio/`. This file can be downloaded on: http://heasarc.gsfc.nasa.gov/fitsio/ and must be placed in `$HEALPIX_DIR/src/cxx/libcfitsio/` as a tar.gz. Then, create the CFITSIO_DIR environment variable in your .whatevershrc. For instance for bash:
```
export CFITSIO_DIR=$HEALPIX_DIR
```
In `$HEALPIX_DIR`, type :
```
./configure
```
This scripts takes care of the configuration and compilation of the C, C++, F90 and IDL packages of the HEALPix distribution. You want to configure the HEALPix C++ package, and edit a Makefile so enter:
```
4
```
You are asked to set the configuration you want to use for the C++ compilation. If a recent version of gcc and g++ is installed on the system, `generic_gcc` should always work. Enter:
```
2
```
For the compilation, HEALPix needs to modify your .profile file, located in your home directory. Don't worry, just one line will be added at the end of your .profile file without any other modifications. Note that a copy of of your actual .profile file will be saved in your home directory as .profile.save. Enter:
```
y
```
The C++ configuration is done. You can exit. Enter:
```
0
```
You can go now in `$HEALPIX_DIR/src/cxx` and do:
```
make
```
In order to compile our *COVERAGE & ANISOTROPY TOOLKIT* you will need the C part as a library. In the
`$HEALPIX_DIR/src/C` directory use:
```
./doinstall
```
The following message should appear:
```
Warning: The following directories could not be found:
$HEALPIX_DIR/include
$HEALPIX_DIR/lib
Should I attempt to create these directories (Y|n)?
```
Answer no. The directory in question already exist but not in `$HEALPIX_DIR`. You need to create in `$HEALPIX_DIR` the following symbolic links:
```
ln -s src/cxx/generic_gcc/bin bin
ln -s src/cxx/generic_gcc/lib lib
ln -s src/cxx/generic_gcc/include include
```
Return in `$HEALPIX_DIR/src/C` and type `./doinstall`.

You are then asked to enter the C compiler you want to use. `gcc` is the default one and is fine. Press `return`.  
To choose the options for the compilation. `-O2 -Wall` is fine. Press `return`.  
To enter the archive and indexing creation command. `ar -rsv is fine`. Press `return`.  
To enter full name of CFITSIO library (`libcfitsio.a`). Press `return`.  
To enter location of CFITSIO libray (`$HEALPIX_DIR/lib`). Press `return`.  
To enter location of CFITSIO header cfitsio.h (`usr/local/include`). Enter `$HEALPIX_DIR/include`.  
A static library is produced by default. You are asked if you want a shared library too. Press `n`.


## Install ROOT (available at https://root.cern.ch/)
Don't ask me how to do it ...
Don't forget to define your `ROOTSYS` environment variable. And also don't forget to add the `$ROOTSYS/lib` directory in your `LD_LIBRARY_PATH` environment variable, in your .bashrc file for instance, by adding in it:
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib)
```


## Install and compile STCoordinates
Download from our web page STCoordinates.tar.gz and uncompress it. Define the environment variable `STC_DIR` in your .whatevershrc file pointing towards the directory where it has been uncompressed.
```
cd $STC_DIR
make
```
it should create a file called libSTCoordinates.a.


## Install the "COVERAGE & ANISOTROPY TOOLKIT"
Go to our web page and download the file Toolkit-v3.0.tar.gz. If you have correctly set up the environment variables `$CFITSIO_DIR`, `$HEALPIX_DIR`, `$ROOTSYS` and `$STC_DIR` then you should just go into the directory where you uncompressed Toolkit-v3.0.tar.gz and type:
```
make
```


## Get started
We joined some exe all called `example_*.exe` to show what can be done with the *COVERAGE & ANISOTROPY TOOLKIT*. This programs are thoroughly documented. In addition, the information on the functions used in our code are available on the web site through a DOXYGEN documentation. We also provide the DOXYGEN configuration file (doxyconf.txt) so you can easily produce this documentation in your *COVERAGE & ANISOTROPY TOOLKIT* directory. If you have DOXYGEN installed on your machine just do:
```
doxygen doxyconf.txt
```
and a `doc/` directory will be created with the html documentation. To see it, you simply need to open the index.html file in `doc/html/` with your web browser.


NOTE : during the execution of the code, we have to create some temporary files in `/tmp`. So please, be sure you have the appropriate permissions in this directory and enough space.
