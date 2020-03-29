#include <fstream>
#include <complex>
#include <iomanip>

std::complex<double>i1(0.0,1.0); // imaginary number defined as 
                                 // global variable ensures compatibility
                                 // with older versions of C++

class LinearModel{
  public:
  LinearModel(double& mi_,  double& Te_, double& B0_, double& E0_ ,
              double& R0_,  double& LB_, double& Ln_, double& De0_,
              double& vix_, double& kz_, double& kx_, double& n0_,
              double& t_):
              mi{mi_}, Te{Te_}, B0{B0_}, E0{E0_}, R0{R0_}, LB{LB_}, 
              Ln{Ln_}, De0{De0_}, vix{vix_}, kz{kz_}, kx{kx_},
              n0{n0_}, t{t_}{ E = E0; scale_factor = 1.0; update_all(); }
              
  void set_E(const double& newE){E = newE; scale_factor = E/E0; update_all();}

  double wI(const double ky){
    double kp = sqrt(kx*kx+ky*ky);
    double wd = vd*ky;
    double wD = vD*ky;
    double w0 = v0*ky;      
    return 0.5/sqrt(2)*sqrt(sqrt(pow(cs,8)*pow(kp,8)/pow((wd-wD),4) +
           16*pow(cs,4)*pow(kp,4)/pow((wd-wD),2)*(pow((wi0-w0-wD),2)+nupara*nupara)+ 
           8*pow(cs,6)*pow(kp,6)/pow((wd-wD),3)*(wi0-w0-wD))-
           pow(cs,4)*pow(kp,4)/pow((wd-wD),2)-4*cs*cs*kp*kp/(wd-wD)* 
           (wi0-w0-wD));
  }
  
  double alpha(const double ky){
    return wI(ky)-beta(ky)-gamma(ky);
  }
  
  double beta(const double ky){
    double kp = sqrt(kx*kx+ky*ky);
    double wd = vd*ky;
    double wD = vD*ky;
    return 1.0/8.0*pow(cs,4)*pow(kp,4)/pow((wd-wD),2);
  } 
 
  double gamma(const double ky){
    double kp = sqrt(kx*kx+ky*ky);
    double wd = vd*ky;
    double wD = vD*ky;
    double w0 = v0*ky;
    return 0.5*cs*cs*kp*kp/(wd-wD)*(wi0-w0-wD);
  } 

  double wR(const double ky){
    double kp = sqrt(kx*kx+ky*ky);
    double wd = vd*ky;
    double wD = vD*ky;
    double w0 = v0*ky;      
    return  0.5*(2*wi0+kp*kp*cs*cs/(wd-wD))+
            0.5/sqrt(2)*sqrt(sqrt(pow(cs,8)*pow(kp,8)/pow((wd-wD),4) +
            16*pow(cs,4)*pow(kp,4)/pow((wd-wD),2)*(pow((wi0-w0-wD),2)+nupara*nupara)+
            8*pow(cs,6)*pow(kp,6)/pow((wd-wD),3)*(wi0-w0-wD))+
            pow(cs,4)*pow(kp,4)/pow((wd-wD),2)+4*cs*cs*kp*kp/(wd-wD)*
            (wi0-w0-wD));
  }
  
  double A(){
    return wi0;
  }
  
  double B(const double ky){
    double kp = sqrt(kx*kx+ky*ky);
    double wd = vd*ky;
    double wD = vD*ky;   
    return 0.5*(kp*kp*cs*cs/(wd-wD));    
  }
  
  double C(double ky){      
    return wR(ky)-A()-B(ky);    
  }
    
  std::complex<double> dphi_decr(const std::complex<double>& dne_over_n0){
    double wd = vd*my_decr/R0;
    double wD = vD*my_decr/R0;
    double w0 = v0*my_decr/R0;     
    return dne_over_n0*T*((w_decr-w0-wD)+i1*(kz*kz*De))/(wd-wD) ;
  }

  std::complex<double> dphi_incr(const std::complex<double>& dne_over_n0){
    double wd = vd*my_incr/R0;
    double wD = vD*my_incr/R0;
    double w0 = v0*my_incr/R0;     
    return dne_over_n0*T*((w_incr-w0-wD)+i1*(kz*kz*De))/(wd-wD) ;
  }  
  
  std::complex<double> dvix1_decr(const std::complex<double>& dne_over_n0){     
    return  e/mi*kx/(w_decr-kx*vi0)*dphi_decr(dne_over_n0);
  }   

  std::complex<double> dvix1_incr(const std::complex<double>& dne_over_n0){     
    return  e/mi*kx/(w_incr-kx*vi0)*dphi_incr(dne_over_n0);
  }   

  std::complex<double> dviy1_decr(const std::complex<double>& dne_over_n0){     
    return  e/mi*my_decr/R0/(w_decr-kx*vi0)*dphi_decr(dne_over_n0);
  } 

  std::complex<double> dviy1_incr(const std::complex<double>& dne_over_n0){     
    return  e/mi*my_incr/R0/(w_incr-kx*vi0)*dphi_incr(dne_over_n0);
  }

  std::complex<double> dvex1_decr(const std::complex<double>& dne_over_n0){     
    return  -my_decr/R0*dphi_decr(dne_over_n0)/B0;
  }  

  std::complex<double> dvex1_incr(const std::complex<double>& dne_over_n0){     
    return  -my_incr/R0*dphi_incr(dne_over_n0)/B0;
  }  

  double Jxi0() const{ return e*ne0*vi0;}
  double Ixi0() const{ return Jxi0()*A_plasma;}
  
  double Jxi1_decr(const std::complex<double>& dne_over_n0){ 
    return 0.5*e*ne0*std::real(dne_over_n0*std::conj(dvix1_decr(dne_over_n0)));
  }
  double Jxi1_incr(const std::complex<double>& dne_over_n0){ 
    return 0.5*e*ne0*std::real(dne_over_n0*std::conj(dvix1_incr(dne_over_n0)));
  }

  double Ixi1_decr(const std::complex<double>& dne_over_n0){ 
    return Jxi1_decr(dne_over_n0)*A_plasma;
  }
  double Ixi1_incr(const std::complex<double>& dne_over_n0){ 
    return Jxi1_incr(dne_over_n0)*A_plasma;
  }
  
  double Jxe1_decr(const std::complex<double>& dne_over_n0){
//  The following commented code compares the return value to the equation
//  in the paper and leads to the same result. The equation in the paper is 
//  for real dne_over_n0, so this method is more general.

 /*    double wd = vd*my_decr/R0;
    double wD = vD*my_decr/R0;
    double w0 = v0*my_decr/R0;    
    double eta = (w_decr.real()-w0-wD)/(wd-wD);
    double expl = 0.5*(ne0*my_decr*T*e*eta)/(B0*R0)*dne_over_n0.real()*dne_over_n0.real();
    double nonexpl = -0.5*e*ne0*std::real(dne_over_n0*std::conj(dvex1_decr(dne_over_n0)));
    double paper_estim = 0.5*(ne0*my_decr*T*e*eta)/(B0*R0)*dne_over_n0.real()*dne_over_n0.real()*A_plasma;
    std::cout << paper_estim << " " << nonexpl*A_plasma << " " << expl*A_plasma << std::endl;
 */ 
 
    return -0.5*e*ne0*std::real(dne_over_n0*std::conj(dvex1_decr(dne_over_n0)));
  }
  double Jxe1_incr(const std::complex<double>& dne_over_n0){
    return -0.5*e*ne0*std::real(dne_over_n0*std::conj(dvex1_incr(dne_over_n0)));
  }

  double Ixe1_decr(const std::complex<double>& dne_over_n0){
    return Jxe1_decr(dne_over_n0)*A_plasma;
  }
  double Ixe1_incr(const std::complex<double>& dne_over_n0){
    return Jxe1_incr(dne_over_n0)*A_plasma;
  }
  
  double Ix_decr(const std::complex<double>& dne_over_n0){
    return Ixi0()+Ixi1_decr(dne_over_n0)+Ixe1_decr(dne_over_n0);
  }
  double Ix_incr(const std::complex<double>& dne_over_n0){
    return Ixi0()+Ixi1_incr(dne_over_n0)+Ixe1_incr(dne_over_n0);
  }  
 
  double LB_limit_decr(){
    return 2*T/( std::abs(vi0)*kx/(my_decr/R0)*B0 + std::abs(E));
  }
  
  double LB_limit_incr(){
    return 2*T/( std::abs(vi0)*kx/(my_incr/R0)*B0 + std::abs(E));
  }
  
  void Save2File(std::ofstream& s, const std::complex<double>& dne_over_n0){
    unsigned int width{22};
    // mode independent values
    s << std::setw(width);
    s << -get_plasma_potential() << std::setw(width);
    s << get_vi0() << std::setw(width);
    s << get_De() << std::setw(width);
    s << get_T() << std::setw(width);
    s << Jxi0() << std::setw(width);
  // decreasing modes values    
    s << get_my_decr() << std::setw(width);
    s << get_wR_decr()/2/M_PI << std::setw(width);
    s << std::abs(dphi_decr(dne_over_n0)) << std::setw(width);
    s << std::abs(dvix1_decr(dne_over_n0)) << std::setw(width); 
    s << std::abs(dvex1_decr(dne_over_n0)) << std::setw(width);
    s << Jxi1_decr(dne_over_n0) << std::setw(width);
    s << Jxe1_decr(dne_over_n0) << std::setw(width);
    s << Ixi1_decr(dne_over_n0) << std::setw(width);
    s << Ixe1_decr(dne_over_n0) << std::setw(width);
    s << Ix_decr(dne_over_n0) << std::setw(width);
    s << LB_limit_decr() << std::setw(width);
    s << A() << std::setw(width);
    s << B(my_decr/R0) << std::setw(width);
    s << C(my_decr/R0) << std::setw(width);
    s << alpha(my_decr/R0) << std::setw(width);
    s << beta(my_decr/R0) << std::setw(width);
    s << gamma(my_decr/R0) << std::setw(width);
  // increasing modes values    
    s << get_my_incr() << std::setw(width);
    s << get_wR_incr()/2/M_PI << std::setw(width);
    s << std::abs(dphi_incr(dne_over_n0)) << std::setw(width);
    s << std::abs(dvix1_incr(dne_over_n0)) << std::setw(width); 
    s << std::abs(dvex1_incr(dne_over_n0)) << std::setw(width);    
    s << Jxi1_incr(dne_over_n0) << std::setw(width);
    s << Jxe1_incr(dne_over_n0) << std::setw(width);
    s << Ixi1_incr(dne_over_n0) << std::setw(width);
    s << Ixe1_incr(dne_over_n0) << std::setw(width);
    s << Ix_incr(dne_over_n0) << std::setw(width);
    s << LB_limit_incr() << std::setw(width);
    s << A() << std::setw(width);
    s << B(my_incr/R0) << std::setw(width);
    s << C(my_incr/R0) << std::setw(width);
    s << alpha(my_incr/R0) << std::setw(width);
    s << beta(my_incr/R0) << std::setw(width);
    s << gamma(my_incr/R0);
    s << "\n";
  }
  
  inline double get_wImaxgrowth() const{ return wI_maxgrowth;}
  inline double get_T() const{ return T;}
  inline double get_De() const{ return De;}
  inline double get_ne0() const{ return ne0;}
  inline double get_cs() const{ return cs;}
  inline double get_plasma_potential() const{ return E*Ln;}
  inline double get_vi0() const{ return vi0;}
  inline double get_E() const{ return E;}
  inline double get_my_incr() const{ return my_incr;}
  inline double get_my_decr() const{ return my_decr;}
  inline double get_wI_incr() const{ return w_incr.imag();}
  inline double get_wI_decr() const{ return w_decr.imag();}  
  inline double get_wR_incr() const{ return w_incr.real();}
  inline double get_wR_decr() const{ return w_decr.real();} 
  
  private:
  
  void update_coherent_modes(){
    my_decr = floor(my_maxgrowth);
    w_decr = std::complex<double>(wR(my_decr/R0),wI(my_decr/R0));
    my_incr = ceil(my_maxgrowth);
    w_incr = std::complex<double>(wR(my_incr/R0),wI(my_incr/R0));
  }
  
  void update_wI_maxgrowth(){
    double my_min{0.01};
    double my_max{14};
    double dmy = (my_max-my_min)/100;
    double wI_new;
    wI_maxgrowth = -10;
    for( double my = my_min; my <= my_max; my+=dmy){
      wI_new = wI(my/R0);
      if (wI_new > wI_maxgrowth){
        wI_maxgrowth = wI_new; 
        my_maxgrowth = my;
      }        
    }
    if (wI_maxgrowth < 0){
      throw "negative_wI";
    }
    if (my_maxgrowth == my_min or my_maxgrowth == my_max){
      std::cout << "out" << std::endl;
      throw "maximum_not_found";
    }
  }
  
  void update_all(){
    T = Te*scale_factor;
    De = De0*scale_factor;
    ne0=n0*scale_factor;
    A_plasma = 2*M_PI*R0*t;
    //ne0=n0;
    cs = sqrt(T*e/mi);
    nupara = kz*kz*De;
    vd = -T/B0/Ln;
    vD = -2*T/B0/LB;
    v0 = -E/B0;
    vi0 = vix*sqrt(scale_factor);
    wi0 = kx*vi0;
    update_wI_maxgrowth();
    update_coherent_modes();
  }
  double& mi, Te, B0, E0, R0, LB, Ln, De0, vix, kz, kx, n0, t;
  double e{1.60217657e-19};
  double scale_factor, E, T, ne0, De, cs, nupara, vd, vD, v0, vi0, wi0,
         my_maxgrowth, wI_maxgrowth, A_plasma;
  double my_decr, my_incr;
  std::complex<double> w_incr, w_decr;
};