#include <fstream>
#include <iostream>
#include <cmath>
#include "util.hpp"
#include "linearmodel.hpp"

int main(int argc, char *argv[]){
  // Argument check:
  if(argc!=2){
    std::cout << "Wrong number of inputs. Usage";
    std::cout << "$ ./main plasma_data_file" << std::endl;
    return 0;
  }
  // Variable definition:
  std::string plasma_data_file;
  double mi, t, Te, B0, n0, E, E0, R0, LB, Ln, De0, vix, kz, kx, 
         my_max, my_min, ne_over_n0;
  std::pair<double,double> voltage_scale;
  double voltage_increment, my_increment;
  unsigned int N_voltages, N_modes;
  
  plasma_data_file = argv[1];
  std::ifstream f(plasma_data_file);

  ReadDatafromFile(f, mi, Te, B0, n0, E0, R0, LB, Ln, De0, vix, kz, kx, 
                    my_max, my_min, t, voltage_scale, voltage_increment, 
                    my_increment, N_voltages, N_modes, ne_over_n0);
  
  std::ofstream s("solution.txt");  
  PrintInputHeader(s);
  PrintInputData(s,mi,Te,B0,n0,R0,t,LB,Ln,kz,kx,t,ne_over_n0,N_voltages,N_modes);
  PrintOutputHeader(s);
  s.setf(std::ios::scientific, std::ios::floatfield);
  s.precision(4);  
  LinearModel discharge(mi, Te, B0, E0, R0, LB, Ln, De0, vix, kz, kx, n0, t);
  
  for(unsigned int i = 0; i < N_voltages; i++){
    E = E0*(voltage_scale.first+voltage_increment*i);
    discharge.set_E(E);
    discharge.Save2File(s,ne_over_n0);
  }
   
  std::ofstream g("w_solution.txt");

  // Read data from file
  if (g.is_open()){

    g.setf(std::ios::scientific, std::ios::floatfield);
    g.precision(4);
    discharge.set_E(E0);

    for( double my = 2.0; my < my_max; my+=my_increment){
    double ky = my/R0;
    
    g  << my << "   ";	
    g  << discharge.wR(ky) << "   ";	
    g  << discharge.A() << "   ";	
    g  << discharge.B(ky)	   << "   ";	
    g  << discharge.C(ky)	   << "   ";	     
    g  << discharge.wI(ky) << "   ";
    g  << discharge.alpha(ky) << "   ";
    g  << discharge.beta(ky) << "   ";
    g  << discharge.gamma(ky) << "\n";	      

    }
  }
  return 0;
}
  
