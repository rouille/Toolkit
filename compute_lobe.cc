#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include <cstdio>
#include <cstdlib>

using namespace std;

void Usage(string myname)
{
  cout << endl;
  cout << " Synopsis : " << endl;
  cout << myname << " -h, --help to obtain this message" << endl;
  cout << myname << " <lobe type> <size> <number>" << endl;
  cout << "<lobe type> : gauss (gaussian) or circ (top hat)" << endl;
  cout << "<size> : size of the beam in deg. For <gauss>, the size corresponds to the standard deviation" << endl;
  cout << "<number> : the function is tabulated n times" << endl << endl;
  
  cout << " Description :" << endl;
  cout << myname << " computes a beam window that will be used to integrate events and coverage maps. The file "
       << "will contain two columns : theta (in deg.) and lobe value normalized to one at maximum." <<  endl;
  
  exit(0);
}



int main(int argc, char* argv[])
{
  // Command line
  if(argc != 4) Usage(argv[0]);
  if(strncmp(argv[1],"-",1) == 0) Usage(argv[0]);

  string type = argv[1];
  if(type != "gauss" && type != "circ")
    {
      cout << "Only gaussian and top hat beams can be made with " << argv[0] << endl;
      exit(0);
    }
  float sig = atof(argv[2]);
  unsigned int nb = atoi(argv[3]);
  vector<double> theta(nb);
  vector<double> value(nb);

  for(unsigned int i = 0; i < nb; i++)
    {
      value[i] = 0;
      if(type == "gauss")
        {
          theta[i] = i*90./(nb-1);
          value[i] = exp(-theta[i]*theta[i]*1./(2.*sig*sig));
        } 
      else if(type == "circ")
        {
          theta[i] = i*90./(nb-1);
          if( theta[i] <= sig ) value[i] = 1.;
        }
    }

  char tmpName[256];
  sprintf(tmpName, "lobe_%s%4.2fdeg.txt", type.c_str(), sig);
  string name = (string) tmpName;
  ofstream lobeFile(name.c_str());
  double theval;
  for(unsigned int i = 0; i < nb; i++)
    {
      theval=value[i];
      if(value[i] < 1.e-20) theval=0.;
      lobeFile << theta[i] << " " << theval << endl;
      cout << theta[i] << " " << theval << endl;
    }
  lobeFile.close();
}
