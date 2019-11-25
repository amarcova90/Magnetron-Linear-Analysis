#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
//#include <Python.h>

int main(int argc, char *argv[]){
  if(argc!=2){
    std::cout << "Wrong number of inputs. Usage";
    std::cout << "$ ./main plasma_data_file" << std::endl;
    return 0;
  }
  std::string plasma_data_file;
  double mi, Te, B0, ne, E0, R0, LB, Ln, De0, vi0, kz, kx, my_max, my_min;
  double voltage_scale[2];
  double voltage_increment, my_increment;
  unsigned int N_voltages, N_modes;
  //double me = 9.10938291e-31;
  double e = 1.60217657e-19;
  //double kb = 1.3806488e-23;
  double ua_to_kg = 1.66053892e-27;

  
  plasma_data_file = argv[1];
  std::ifstream f(plasma_data_file);

  // Read data from file
  if (f.is_open()){
      std::string name;

      while (f >> name){
	if (name == "ion_mass_u"){
	  f >> mi;
          mi=mi*ua_to_kg;
	  }
	else if (name == "Te")
	  f >> Te;
	else if (name == "B0")
	  f >> B0;
	else if (name == "E0")
	  f >> E0;
	else if (name == "R0")
	  f >> R0;
	else if (name == "LB")
	  f >> LB;
	else if (name == "Ln")
	  f >> Ln;
    else if (name == "ne")
	  f >> ne;
	else if (name == "electron_diffusion")
	  f >> De0;
	else if (name == "vi0")
	  f >> vi0;
	else if (name == "lambda_z"){
	  f >> kz;
	  kz=2*M_PI/kz;
	  }
	else if (name == "lambda_x"){
	  f >> kx;
	  kx=2*M_PI/kx;
	  }
	else if (name == "N_voltages")
	  f >> N_voltages;
	else if (name == "N_modes")
	  f >> N_modes;
	else if (name == "my_max")
	  f >> my_max;
	else if (name == "my_min")
	  f >> my_min;
	else if (name == "voltage_scale")
	  f >> voltage_scale[0] >> voltage_scale[1]; 
      }
      voltage_increment=(voltage_scale[1]-voltage_scale[0])/(N_voltages-1);
      my_increment = (my_max - my_min) / (N_modes - 1);
      f.close();
  }

  std::ofstream s("solution.txt");
  if(s.is_open()){
	s << "###########################################\n";
	s << "#             Input parameters            #\n";
	s << "###########################################\n";
	s << "  ion_mass_u                     " << mi / ua_to_kg << "\n";
	s << "  Te                             " << Te << "\n";
    s << "  B0                             " << B0 << "\n";
	s << "  R0                             " << R0 << "\n";
	s << "  LB                             " << LB << "\n";
	s << "  Ln                             " << Ln << "\n";
	s << "  kz                             " << kz << "\n";
	s << "  kx                             " << kx << "\n";
	s << "  N_voltages                     " << N_voltages << "\n";
	s << "  N_modes                        " << N_modes << "\n\n";
	s << "###########################################\n";
	s << "#             Output parameters           #\n";
	s << "###########################################\n";	
	s << "#  Voltage [V]   m decr []      f decr [Hz]   Ia decr [A]   m incr []     "	
	      " f incr [Hz]   Ia decr [A]     Ion velocity [m/s]   De [m^2/s]    Te [eV]\n";
    double E, T, cs, Ix, nupara, scale_factor, De, vd, vD, v0, vix; 
	double wi0, my, ky, ky_incr, ky_decr, kp, wd, wD, w0, wR, wR_incr, wR_decr;
	double  wI_new, wI_old, wI_previous;

	s.setf(std::ios::scientific, std::ios::floatfield);
    s.precision(4);
	for(unsigned int i = 0; i < N_voltages; i++){
      
	  E = E0*(voltage_scale[0]+voltage_increment*i);
      scale_factor = E/E0;
      T = Te*scale_factor*e;
      De = De0*scale_factor;
      cs = sqrt(T/mi);
      nupara = kz*kz*De;
      vd = -T/e/B0/Ln;
      vD = -2*T/e/B0/LB;
      v0 = -E/B0;
      vix = vi0*sqrt(scale_factor);
      wi0 = kx*vix;

	  wI_new = -1.0; wI_old = -1.0; wI_previous = -1.0;
	  bool flag = false;
      for( my = my_min; my < my_max; my+=my_increment){
		  wI_old=wI_previous;
	      wI_previous=wI_new;
	      ky = my/R0;
	      kp = sqrt(kx*kx+ky*ky);
	      wd = vd*ky;
	      wD = vD*ky;
	      w0 = v0*ky;
		  
	      wI_new = 0.5/sqrt(2)*sqrt(sqrt(pow(cs,8)*pow(kp,8)/pow((wd-wD),4) +
                   16*pow(cs,4)*pow(kp,4)/pow((wd-wD),2)*(pow((wi0-w0-wD),2)+nupara*nupara)+ 
	               8*pow(cs,6)*pow(kp,6)/pow((wd-wD),3)*(wi0-w0-wD))-
                   pow(cs,4)*pow(kp,4)/pow((wd-wD),2)-4*cs*cs*kp*kp/(wd-wD)* 
	               (wi0-w0-wD));

	      wR = 0.5*(2*wi0+kp*kp*cs*cs/(wd-wD))+
               0.5/sqrt(2)*sqrt(sqrt(pow(cs,8)*pow(kp,8)/pow((wd-wD),4) +
               16*pow(cs,4)*pow(kp,4)/pow((wd-wD),2)*(pow((wi0-w0-wD),2)+nupara*nupara)+
               8*pow(cs,6)*pow(kp,6)/pow((wd-wD),3)*(wi0-w0-wD))+
               pow(cs,4)*pow(kp,4)/pow((wd-wD),2)+4*cs*cs*kp*kp/(wd-wD)*
               (wi0-w0-wD));
		
	      if (wI_old > 0.0 and wI_previous > 0.0 and wI_new > 0.0 and 
	          wI_previous > wI_old and wI_previous > wI_new){

			  s << "  " << E*Ln << "       " ;
			  
			  ky_decr = floor(my)/R0;
			  kp = sqrt(kx*kx+ky_decr*ky_decr);
			  wd = vd*ky_decr;
	          wD = vD*ky_decr;
	          w0 = v0*ky_decr;
			  wR_decr = 0.5*(2*wi0+kp*kp*cs*cs/(wd-wD))+
                   0.5/sqrt(2)*sqrt(sqrt(pow(cs,8)*pow(kp,8)/pow((wd-wD),4) +
                   16*pow(cs,4)*pow(kp,4)/pow((wd-wD),2)*(pow((wi0-w0-wD),2)+nupara*nupara)+
                   8*pow(cs,6)*pow(kp,6)/pow((wd-wD),3)*(wi0-w0-wD))+
                   pow(cs,4)*pow(kp,4)/pow((wd-wD),2)+4*cs*cs*kp*kp/(wd-wD)*
                   (wi0-w0-wD));	

			  s << (int)floor(my) << "           " << wR_decr/M_PI/2;
			  
			  
			  double dphi = E*Ln/10;//5;
			  Ix = sqrt(pow((ky*ne*(wd-wD)*(wI_previous+kz*kz*De)*e*e*dphi*dphi/T/B0)/
			              (pow(wR-w0-wD,2)+pow(wI_previous+kz*kz*De,2)),2) +
				          pow((ky*ne*(wd-wD)*(wR-w0-wD)*e*e*dphi*dphi/T/B0)/
			              (pow(wR-w0-wD,2)+pow(wI_previous+kz*kz*De,2)),2))*M_PI*R0*R0/8;
			  s << "    " << Ix ;
			 
			  
			  ky_incr = ceil(my)/R0;
			  kp = sqrt(kx*kx+ky_incr*ky_incr);
			  wd = vd*ky_incr;
	          wD = vD*ky_incr;
	          w0 = v0*ky_incr;
			  wR_incr = 0.5*(2*wi0+kp*kp*cs*cs/(wd-wD))+
                   0.5/sqrt(2)*sqrt(sqrt(pow(cs,8)*pow(kp,8)/pow((wd-wD),4) +
                   16*pow(cs,4)*pow(kp,4)/pow((wd-wD),2)*(pow((wi0-w0-wD),2)+nupara*nupara)+
                   8*pow(cs,6)*pow(kp,6)/pow((wd-wD),3)*(wi0-w0-wD))+
                   pow(cs,4)*pow(kp,4)/pow((wd-wD),2)+4*cs*cs*kp*kp/(wd-wD)*
                   (wi0-w0-wD));
			  
			  s << "       " << (int)ceil(my)<< "           " << wR_incr/M_PI/2;
			  

			  Ix = sqrt(pow((ky*ne*(wd-wD)*(wI_previous+kz*kz*De)*e*e*dphi*dphi/T/B0)/
			              (pow(wR-w0-wD,2)+pow(wI_previous+kz*kz*De,2)),2) +
				          pow((ky*ne*(wd-wD)*(wR-w0-wD)*e*e*dphi*dphi/T/B0)/
			              (pow(wR-w0-wD,2)+pow(wI_previous+kz*kz*De,2)),2))*M_PI*R0*R0/8;
		      s << "    " << Ix ;

			  
			  
			  

			  flag = true;

			  
			  s << "      " << vix << "          " << De << "    " << T/e << "\n"; 
			  

	  
	        }
	
	
	
      }
      if (flag == false){
		  std::cout << "Maximum growth rate not found for V = " << E*Ln << " V. ";
		  std::cout << "Consider increasing lowest voltage or highest mode number" << std::endl;
		  return 0;
        }
	}
    
  }
  return 0;
}
  
