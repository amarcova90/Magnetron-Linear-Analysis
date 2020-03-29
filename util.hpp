#include <fstream>
#include <iostream>
#include <cmath>
#include <utility>
#include <string>

void PrintInputHeader(std::ofstream& s){
  s << "###########################################\n";
	s << "#             Input parameters            #\n";
	s << "###########################################\n";
}
void PrintOutputHeader(std::ofstream& s){
	s << "###########################################\n";
	s << "#             Output parameters           #\n";
	s << "###########################################\n";	
  // mode independent values
	s <<  "#          Voltage [V]";
  s <<  "             vi0 [m/s]";
  s <<  "            De [m^2/s]";
  s <<  "               Te [eV]";
  s <<  "          Jxi0 [A/m^2]";
  // decreasing modes values
  s <<  "             m decr []";
  s <<  "           f decr [Hz]";
  s <<  "         dphi decr [V]";
  s <<  "       vix1 decr [m/s]";
  s <<  "       vex1 decr [m/s]";
  s <<  "         dJi1 decr [A]";
  s <<  "         dJe1 decr [A]";
  s <<  "    dIix ampl decr [A]";
  s <<  "    dIex ampl decr [A]";
  s <<  "       Ix tot decr [A]";  
  s <<  "        Lbmin decr [m]";
  s <<  "        A decr [rad/s]";  
  s <<  "        B decr [rad/s]";  
  s <<  "        C decr [rad/s]";
  s <<  "  alpha decr [rad/s]^2";
  s <<  "   beta decr [rad/s]^2";
  s <<  "  gamma decr [rad/s]^2";  
  // increasing modes values
  s <<  "             m incr []";
	s <<  "           f incr [Hz]";
  s <<  "          dphi incr []";
  s <<  "       vix1 incr [m/s]";
  s <<  "       vex1 incr [m/s]";
  s <<  "         dJi1 incr [A]";
  s <<  "         dJe1 incr [A]";
  s <<  "    dIix ampl decr [A]";
  s <<  "    dIex ampl decr [A]";
  s <<  "       Ix tot decr [A]"; 
  s <<  "        Lbmin incr [m]";
  s <<  "        A incr [rad/s]";  
  s <<  "        B incr [rad/s]";  
  s <<  "        C incr [rad/s]";
  s <<  "  alpha incr [rad/s]^2";
  s <<  "   beta incr [rad/s]^2";
  s <<  "  gamma incr [rad/s]^2";  
  s <<  "\n";

}

void ReadDatafromFile(std::ifstream& f, double& mi,     double& Te,
                                        double& B0,     double& n0,     
                                        double& E0,     double& R0,  
                                        double& LB,     double& Ln,
                                        double& De0,    double& vix,    
                                        double& kz,     double& kx, 
                                        double& my_max, double& my_min,
                                        double& t,      std::pair<double,double>& voltage_scale, 
                                        double& voltage_increment, double& my_increment, 
                                        unsigned int& N_voltages, unsigned int& N_modes){
  double ua_to_kg = 1.66053892e-27;                                      
  if (f.is_open()){
    std::string name;
    while (f >> name){
      if (name == "ion_mass_u"){ f >> mi; mi=mi*ua_to_kg;}
      else if (name == "Te") f >> Te;
      else if (name == "B0") f >> B0;
      else if (name == "E0") f >> E0;
      else if (name == "R0") f >> R0;
      else if (name == "LB") f >> LB;
      else if (name == "Ln") f >> Ln;
      else if (name == "n0") f >> n0;
      else if (name == "plasma_thickness") f >> t;
      else if (name == "electron_diffusion") f >> De0;
      else if (name == "vix") f >> vix;
      else if (name == "lambda_z"){ f >> kz; kz=2*M_PI/kz;}
      else if (name == "lambda_x"){ f >> kx; kx=2*M_PI/kx;}
      else if (name == "N_voltages") f >> N_voltages;
      else if (name == "N_modes") f >> N_modes;
      else if (name == "my_max") f >> my_max;
      else if (name == "my_min") f >> my_min;
      else if (name == "voltage_scale") {f >> voltage_scale.first >> voltage_scale.second;
      }      
    }
    voltage_increment=(voltage_scale.second-voltage_scale.first)/(N_voltages-1);
    my_increment = (my_max - my_min) / (N_modes - 1);
    f.close();
  }
  
}

void PrintInputData(std::ofstream& s, double& mi, double& Te, double& B0, 
               double& R0, double& LB, double& Ln, double& kz, double& kx,
               unsigned int& N_voltages, unsigned int& N_modes){
  s << "  ion_mass_u                     " << mi / 1.66053892e-27 << "\n";
	s << "  Te                             " << Te << "\n";
  s << "  B0                             " << B0 << "\n";
	s << "  R0                             " << R0 << "\n";
	s << "  LB                             " << LB << "\n";
	s << "  Ln                             " << Ln << "\n";
	s << "  kz                             " << kz << "\n";
	s << "  kx                             " << kx << "\n";
	s << "  N_voltages                     " << N_voltages << "\n";
	s << "  N_modes                        " << N_modes << "\n\n";
}

