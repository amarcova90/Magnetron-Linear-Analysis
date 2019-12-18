#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>


int main(int argc, char *argv[]){
  //using namespace std::complex_literals; // use this only with C++17 (then you can use 1i)
  if(argc!=2){
    std::cout << "Wrong number of inputs. Usage";
    std::cout << "$ ./main plasma_data_file" << std::endl;
    return 0;
  }
  std::string plasma_data_file;
  double mi, Te, B0, n0, E0, R0, LB, Ln, De0, vix, kz, kx, my_max, my_min, Aspoke = 0;
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
    else if (name == "n0")
	  f >> n0;
    else if (name == "Aspoke")
	  f >> Aspoke;
	else if (name == "electron_diffusion")
	  f >> De0;
	else if (name == "vix")
	  f >> vix;
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
	s << "#  Voltage [V]      vi0 [m/s]        De [m^2/s]    Te [eV]       m decr []      f decr [Hz]   ne/n0 decr []   vi1 decr [m/s]   ve1 decr [m/s]   Ji0 decr [A]    Ji1 decr [A]    Je1 decr [A]    Ix decr [A]    m incr []     "	
	      " f incr [Hz]   ne/n0 incr []   vi1 incr [m/s]   ve1 incr [m/s]   Ji0 incr [A]    Ji1 incr [A]    Je1 incr [A]    Ix incr [A]\n";
    double E, T, cs, ne0, Ix, nupara, scale_factor, De, vd, vD, v0, vi0; 
	double wi0, my, my_decr, my_incr, ky, ky_incr, ky_decr, kp, wd, wD, w0, wR_incr, wR_decr;
	double  wI_new, wI_old, wI_previous, wI_decr, wI_incr;
	std::complex<double> vi1_complex, ve1_complex ,ne_over_n0_complex, Jxi0, Jxi1, Jxe1;

	s.setf(std::ios::scientific, std::ios::floatfield);
    s.precision(4);
	for(unsigned int i = 0; i < N_voltages; i++){
      
	  E = E0*(voltage_scale[0]+voltage_increment*i);
      scale_factor = E/E0;
      T = Te*scale_factor*e;
      De = De0*scale_factor;
	  //ne0=n0*scale_factor;
	  ne0=n0;
      cs = sqrt(T/mi);
      nupara = kz*kz*De;
      vd = -T/e/B0/Ln;
      vD = -2*T/e/B0/LB;
      v0 = -E/B0;
      vi0 = vix*sqrt(scale_factor);
      wi0 = kx*vi0;

	  wI_new = -1.0; wI_old = -1.0; wI_previous = -1.0;
	  bool flag = false;
      for( my = my_min; my < my_max; my+=my_increment){
		  
		  //std::cout << wI_old << "   " << wI_previous << "   " << wI_new << std::endl;
		  
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

/* 	      wR = 0.5*(2*wi0+kp*kp*cs*cs/(wd-wD))+
               0.5/sqrt(2)*sqrt(sqrt(pow(cs,8)*pow(kp,8)/pow((wd-wD),4) +
               16*pow(cs,4)*pow(kp,4)/pow((wd-wD),2)*(pow((wi0-w0-wD),2)+nupara*nupara)+
               8*pow(cs,6)*pow(kp,6)/pow((wd-wD),3)*(wi0-w0-wD))+
               pow(cs,4)*pow(kp,4)/pow((wd-wD),2)+4*cs*cs*kp*kp/(wd-wD)*
               (wi0-w0-wD)); */
		
	      if (wI_old > 0.0 and wI_previous > 0.0 and wI_new > 0.0 and 
	          wI_previous > wI_old and wI_previous > wI_new){

			  s << "  " << E*Ln << "       " ;
			  
			  s << vi0 << "      " << De << "    " << T/e ;
			  
			  // DECREASING VOLTAGES
			  
			  my_decr = floor(my);
			  ky_decr = my_decr/R0;
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
				   
			  wI_decr = 0.5/sqrt(2)*sqrt(sqrt(pow(cs,8)*pow(kp,8)/pow((wd-wD),4) +
                   16*pow(cs,4)*pow(kp,4)/pow((wd-wD),2)*(pow((wi0-w0-wD),2)+nupara*nupara)+ 
	               8*pow(cs,6)*pow(kp,6)/pow((wd-wD),3)*(wi0-w0-wD))-
                   pow(cs,4)*pow(kp,4)/pow((wd-wD),2)-4*cs*cs*kp*kp/(wd-wD)* 
	               (wi0-w0-wD));

			  s << "       " << (int)my_decr << "           " << wR_decr/M_PI/2;
			  
			  
			  double dphi = E*Ln/10;//5;
			  
			  // density flucuation
			  
			  ne_over_n0_complex.real( (wd-wD)*e*dphi/T*(wR_decr-w0-wD)/(pow((wR_decr-w0-wD),2)+pow(wI_decr+kz*kz*De,2)) );
			  ne_over_n0_complex.imag( -(wd-wD)*e*dphi/T*(wI_decr+kz*kz*De)/(pow((wR_decr-w0-wD),2)+pow(wI_decr+kz*kz*De,2)) );
  
			  s << "    " << abs(ne_over_n0_complex)  ;			  

			  //fact=(n0*(wd-wD)*e*e*dphi/T)/(pow(wR-w0-wD,2)+pow(wI_previous+kz*kz*De,2))
			  //vire= vi0*(wR_decr)
			  
			  // Ion velocity fluctuation
			  
			  //vi1 = dphi*(e/mi*kx*(wR_decr-kx*vi0))/(pow(wR_decr-kx*vi0,2)+wI_decr*wI_decr)*sqrt(1+wI_decr*wI_decr/pow(wR_decr-kx*vi0,2));
			  
			  vi1_complex.real ( dphi*(e/mi*kx*(wR_decr-kx*vi0))/(pow(wR_decr-kx*vi0,2)+wI_decr*wI_decr) );
			  vi1_complex.imag ( -dphi*(e/mi*kx*wI_decr)/(pow(wR_decr-kx*vi0,2)+wI_decr*wI_decr) );
			  
			  s << "      " << -abs(vi1_complex) ;
			  
			  ve1_complex.real (0.0);
			  ve1_complex.imag (ky_decr*dphi/B0);
			  
			  s << "      " << -abs(ve1_complex) ;
			  
			  // Current density fluctuation
			  
			  Jxi0 =  e*ne0*(1.0+ne_over_n0_complex)*vi0;
			  Jxi1 =  e*ne0*(1.0+ne_over_n0_complex)*vi1_complex;
			  Jxe1 = -e*ne0*(1.0+ne_over_n0_complex)*ve1_complex;

			  
			  s << "      " << abs(Jxi0) << "      " << abs(Jxi1) << "      " << abs(Jxe1) ;
			  
			  Ix =  -std::real(Jxi0 + Jxi1 + Jxe1) * Aspoke * my_decr;
			  //Ix =  std::real(Jxi0) * Aspoke * my_decr;
			  
			  //std::cout << Ix << std::endl;
			  
			  s << "      " << Ix ;
			 

			  
			  // INCREASING VOLTAGES
			  my_incr = ceil(my);
			  ky_incr = my_incr/R0;
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
			  wI_incr = 0.5/sqrt(2)*sqrt(sqrt(pow(cs,8)*pow(kp,8)/pow((wd-wD),4) +
                   16*pow(cs,4)*pow(kp,4)/pow((wd-wD),2)*(pow((wi0-w0-wD),2)+nupara*nupara)+ 
	               8*pow(cs,6)*pow(kp,6)/pow((wd-wD),3)*(wi0-w0-wD))-
                   pow(cs,4)*pow(kp,4)/pow((wd-wD),2)-4*cs*cs*kp*kp/(wd-wD)* 
	               (wi0-w0-wD));
			  s << "        " << (int)my_incr << "           " << wR_incr/M_PI/2;

			  ne_over_n0_complex.real( (wd-wD)*e*dphi/T*(wR_incr-w0-wD)/(pow((wR_incr-w0-wD),2)+pow(wI_incr+kz*kz*De,2)) );
			  ne_over_n0_complex.imag( -(wd-wD)*e*dphi/T*(wI_incr+kz*kz*De)/(pow((wR_incr-w0-wD),2)+pow(wI_incr+kz*kz*De,2)) );

			  s << "    " << abs(ne_over_n0_complex) ;
		  
			  // Ion velocity fluctuation increasing voltages 
			  
			  vi1_complex.real ( dphi*(e/mi*kx*(wR_incr-kx*vi0))/(pow(wR_incr-kx*vi0,2)+wI_incr*wI_incr) );
			  vi1_complex.imag ( -dphi*(e/mi*kx*wI_incr)/(pow(wR_incr-kx*vi0,2)+wI_incr*wI_incr) );
			    
			  
			  //s << "    " << vi1 ;
			  s << "      " << -abs(vi1_complex) ;
			  
			  ve1_complex.real (0.0);
			  ve1_complex.imag (ky_incr*dphi/B0);
			  
			  s << "      " << -abs(ve1_complex) ;
			  
			  // Current density fluctuation
			  
			  Jxi0 =  e*ne0*(1.0+ne_over_n0_complex)*vi0;
			  Jxi1 =  e*ne0*(1.0+ne_over_n0_complex)*vi1_complex;
			  Jxe1 = -e*ne0*(1.0+ne_over_n0_complex)*ve1_complex;


			  s << "      " << abs(Jxi0) << "      " << abs(Jxi1) << "      " << abs(Jxe1) ;
			  
			  

			  // Current fluctuation increasing voltages 
			  
			  Ix =  -std::real(Jxi0 + Jxi1 + Jxe1) * Aspoke * my_incr;
			  //Ix =  std::real(Jxi0) * Aspoke * my_incr;
			  
			  //std::cout << Ix << std::endl;
			  
			  s << "      " << Ix ;

			  // density flucuation

			  

			  flag = true;

			  
			  s << "\n"; 
			  

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
  
