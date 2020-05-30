#include <fstream>
#include <iostream>
#include <cmath>
#include "util.hpp"
#include "linearmodel.hpp"
#include "nonlinearmodel.hpp"


int main(int argc, char *argv[]){
  // Argument check:
  if(argc!=2){
    std::cout << "Wrong number of inputs. Usage";
    std::cout << "$ ./main plasma_data_file" << std::endl;
    return 0;
  }
  // Variable definition:
  std::string plasma_data_file;
  double mi, t, Te, B0, n0, E0, R0, LB, Ln, De0, vix, kz, kx, 
         my_max, my_min, ne_over_n0;
  std::pair<double,double> voltage_scale;
  double voltage_increment, my_increment;
  unsigned int N_voltages, N_modes;
  
  // Additional arguments for nonlinear analysis
  int Ni, Nj;
  double tf, n10_n0;
  
  plasma_data_file = argv[1];
  std::ifstream f(plasma_data_file);


  ReadNonlinearDatafromFile(f, mi, Te, B0, n0, E0, R0, LB, Ln, De0, vix, kz, kx, 
                    my_max, my_min, t, voltage_scale, voltage_increment, 
                    my_increment, N_voltages, N_modes, ne_over_n0, 
                    Ni, Nj, tf, n10_n0);

  
  int Nj_val = 2000000;
  NonlinearModel NLdischarge_val(mi, Te, B0, E0, R0, LB, Ln, De0, vix, kz, kx, n0, t,
                             Ni, Nj_val, tf, n10_n0);
  for (int N_val = 1000; N_val < Nj_val; N_val*=2){
    NLdischarge_val.resize(1,N_val);
    NLdischarge_val.solver_validation();
  }

  return 0;
}
  
