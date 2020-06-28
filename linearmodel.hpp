#pragma once

#include <fstream>
#include <complex>
#include <iomanip>


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

  // assumes 1/Ln << kx 
  double wI(const double& ky){
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
  
  double alpha(const double& ky){
    return wI(ky)-beta(ky)-gamma(ky);
  }
  
  double beta(const double& ky){
    double kp = sqrt(kx*kx+ky*ky);
    double wd = vd*ky;
    double wD = vD*ky;
    return 1.0/8.0*pow(cs,4)*pow(kp,4)/pow((wd-wD),2);
  } 
 
  double gamma(const double& ky){
    double kp = sqrt(kx*kx+ky*ky);
    double wd = vd*ky;
    double wD = vD*ky;
    double w0 = v0*ky;
    return 0.5*cs*cs*kp*kp/(wd-wD)*(wi0-w0-wD);
  } 

  // assumes 1/Ln << kx 
  double wR(const double& ky){
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
  
  double B(const double& ky){
    double kp = sqrt(kx*kx+ky*ky);
    double wd = vd*ky;
    double wD = vD*ky;   
    return 0.5*(kp*kp*cs*cs/(wd-wD));    
  }
  
  double C(double& ky){
    return wR(ky)-A()-B(ky);    
  }
  
  double wR_alter(const double& ky){
    std::complex<double> comp = w_alter(ky);
    return comp.real();
  }

  double wI_alter(const double& ky){
    std::complex<double> comp = w_alter(ky);
    return comp.imag();
  }
    
  std::complex<double> dphi_decr(const std::complex<double>& dne_over_n0){
    double kp = sqrt(kx*kx+ky_decr*ky_decr);
    //double wd = vd*ky_decr;
    //double wD = vD*ky_decr;
    //double w0 = v0*ky_decr;     
    //return dne_over_n0*T*(w_decr-w0-wD+i1*nupara)/(wd-wD) ;
    // The following is equivalent to the form on top
    return dne_over_n0*(w_decr-kx*vi0)*(w_decr-kx*vi0)/(e/mi*(kp*kp - i1*kx/Ln));
  }
  
  std::complex<double> dphi_incr(const std::complex<double>& dne_over_n0){
    double wd = vd*ky_incr;
    double wD = vD*ky_incr;
    double w0 = v0*ky_incr;     
    return dne_over_n0*T*(w_incr-w0-wD+i1*nupara)/(wd-wD) ;
  }  

  
  double dphi_phase_decr(const std::complex<double>& dne_over_n0){
    return std::arg( dphi_decr(dne_over_n0) );
  }
  double dphi_phase_incr(const std::complex<double>& dne_over_n0){
    return std::arg( dphi_incr(dne_over_n0) );
  }
  
  std::complex<double> dvxi1_decr(const std::complex<double>& dne_over_n0){
    return  e/mi*kx/(w_decr-kx*vi0)*dphi_decr(dne_over_n0);
  }   
  std::complex<double> dvxi1_incr(const std::complex<double>& dne_over_n0){
    return  e/mi*kx/(w_incr-kx*vi0)*dphi_incr(dne_over_n0);
  }   

  double dvxi1_phase_decr(const std::complex<double>& dne_over_n0){
    return std::arg( dvxi1_decr(dne_over_n0) );
  }
  double dvxi1_phase_incr(const std::complex<double>& dne_over_n0){
    return std::arg( dvxi1_incr(dne_over_n0) );
  }

  std::complex<double> dvyi1_decr(const std::complex<double>& dne_over_n0){
    return  e/mi*ky_decr/(w_decr-kx*vi0)*dphi_decr(dne_over_n0);
  } 
  std::complex<double> dvyi1_incr(const std::complex<double>& dne_over_n0){
    return  e/mi*ky_incr/(w_incr-kx*vi0)*dphi_incr(dne_over_n0);
  }

  double dvyi1_phase_decr(const std::complex<double>& dne_over_n0){
    return std::arg( dvyi1_decr(dne_over_n0) );
  }
  double dvyi1_phase_incr(const std::complex<double>& dne_over_n0){
    return std::arg( dvyi1_incr(dne_over_n0) );
  }

  std::complex<double> dvxe1_EXB_decr(const std::complex<double>& dne_over_n0){
    return  -i1*ky_decr*dphi_decr(dne_over_n0)/B0;
  }  
  std::complex<double> dvxe1_EXB_incr(const std::complex<double>& dne_over_n0){
    return  -i1*ky_incr*dphi_incr(dne_over_n0)/B0;
  }  
  
  double dvxe1_EXB_phase_decr(const std::complex<double>& dne_over_n0){
    return std::arg( dvxe1_EXB_decr(dne_over_n0) );
  }  
  double dvxe1_EXB_phase_incr(const std::complex<double>& dne_over_n0){
    return std::arg( dvxe1_EXB_incr(dne_over_n0) );
  }  
  
  std::complex<double> dvxe1_D_decr(const std::complex<double>& dne_over_n0){
    return  i1*T/B0*ky_decr*dne_over_n0;
  }
  std::complex<double> dvxe1_D_incr(const std::complex<double>& dne_over_n0){
    return  i1*T/B0*ky_incr*dne_over_n0;
  }

  double dvxe1_D_phase_decr(const std::complex<double>& dne_over_n0){
    return std::arg( dvxe1_D_decr(dne_over_n0) );
  }  
  double dvxe1_D_phase_incr(const std::complex<double>& dne_over_n0){
    return std::arg( dvxe1_D_incr(dne_over_n0) );
  }  
  
  double Jxi0() const{ return e*ne0*vi0;}
  double Ixi0() const{ return Jxi0()*A_plasma;}
  
  double Jxi1_ave_decr(const std::complex<double>& dne_over_n0){
   /*  // comparison with equations in the paper (leads tto same result)
    // equation in paper is speciaol case when deltan is real
    double Ji = 0.5*my_decr*e*ne0*kx/(kx*kx+ky_decr*ky_decr)*(w_decr.real()-kx*vi0);
    std::cout << Ji << std::endl;
    std::cout << my_decr*0.5*e*ne0*std::real(dne_over_n0*std::conj(dvxi1_decr(dne_over_n0))) << "\n" << std::endl;
 */    
    return 0.5*e*ne0*std::real(dne_over_n0*std::conj(dvxi1_decr(dne_over_n0)));
  }
  double Jxi1_ave_incr(const std::complex<double>& dne_over_n0){
    return 0.5*e*ne0*std::real(dne_over_n0*std::conj(dvxi1_incr(dne_over_n0)));
  }
  
  double Ixi1_ave_decr(const std::complex<double>& dne_over_n0){
    return Jxi1_ave_decr(dne_over_n0)*A_plasma;
  }
  double Ixi1_ave_incr(const std::complex<double>& dne_over_n0){
    return Jxi1_ave_incr(dne_over_n0)*A_plasma;
  }
  
  double Jxe1_EXB_ave_decr(const std::complex<double>& dne_over_n0){
    return -0.5*e*ne0*std::real(dne_over_n0*std::conj(dvxe1_EXB_decr(dne_over_n0)));
  }
  double Jxe1_EXB_ave_incr(const std::complex<double>& dne_over_n0){
    return -0.5*e*ne0*std::real(dne_over_n0*std::conj(dvxe1_EXB_incr(dne_over_n0)));
  }
  
  double Ixe1_EXB_ave_decr(const std::complex<double>& dne_over_n0){
    return Jxe1_EXB_ave_decr(dne_over_n0)*A_plasma;
  }
  double Ixe1_EXB_ave_incr(const std::complex<double>& dne_over_n0){
    return Jxe1_EXB_ave_incr(dne_over_n0)*A_plasma;
  }
  
  double Jxe1_D_ave_decr(const std::complex<double>& dne_over_n0){
    return -0.5*e*ne0*std::real(dne_over_n0*std::conj(dvxe1_D_decr(dne_over_n0)));
  }
  double Jxe1_D_ave_incr(const std::complex<double>& dne_over_n0){
    return -0.5*e*ne0*std::real(dne_over_n0*std::conj(dvxe1_D_incr(dne_over_n0)));
  }  
  
  double Ixe1_D_ave_decr(const std::complex<double>& dne_over_n0){
    return Jxe1_D_ave_decr(dne_over_n0)*A_plasma;
  }
  double Ixe1_D_ave_incr(const std::complex<double>& dne_over_n0){
    return Jxe1_D_ave_incr(dne_over_n0)*A_plasma;
  }
  
  double Ix_ave_decr(const std::complex<double>& dne_over_n0){
    return Ixi0()+Ixe1_D_ave_decr(dne_over_n0)+Ixi1_ave_decr(dne_over_n0)+Ixe1_EXB_ave_decr(dne_over_n0);
  }
  double Ix_ave_incr(const std::complex<double>& dne_over_n0){
    return Ixi0()+Ixe1_D_ave_incr(dne_over_n0)+Ixi1_ave_incr(dne_over_n0)+Ixe1_EXB_ave_incr(dne_over_n0);
  } 
  
  double LB_limit_decr(){
    return 2*T/( std::abs(vi0)*kx/(ky_decr)*B0 + std::abs(E));
  }
  
  double LB_limit_incr(){
    return 2*T/( std::abs(vi0)*kx/(ky_incr)*B0 + std::abs(E));
  }
  
  void Save2File(std::ofstream& s, const std::complex<double>& dne_over_n0){
    unsigned int width{22};
    // mode independent values
    s << std::setw(width);
    s << -get_plasma_potential() << std::setw(width);
    s << ne0 << std::setw(width);
    s << get_vi0() << std::setw(width);
    s << get_De() << std::setw(width);
    s << get_T() << std::setw(width);
    s << Jxi0() << std::setw(width);
    s << Ixi0() << std::setw(width);
  // decreasing modes values    
    s << get_my_decr() << std::setw(width);
    s << get_wR_decr()/2/M_PI << std::setw(width);
    s << std::abs(dphi_decr(dne_over_n0)) << std::setw(width);
    s << std::abs(dvxi1_decr(dne_over_n0)) << std::setw(width); 
    s << std::abs(dvxe1_EXB_decr(dne_over_n0)) << std::setw(width);
    s << std::abs(dvxe1_D_decr(dne_over_n0)) << std::setw(width);
    
    s << dphi_phase_decr(dne_over_n0) << std::setw(width);
    s << dvxi1_phase_decr(dne_over_n0) << std::setw(width);
    s << dvxe1_EXB_phase_decr(dne_over_n0) << std::setw(width);
    s << dvxe1_D_phase_decr(dne_over_n0) << std::setw(width);
    
    s << Jxi1_ave_decr(dne_over_n0) << std::setw(width);
    s << Jxe1_EXB_ave_decr(dne_over_n0) << std::setw(width);
    s << Jxe1_D_ave_decr(dne_over_n0) << std::setw(width);     
    s << Ixi1_ave_decr(dne_over_n0) << std::setw(width);
    s << Ixe1_EXB_ave_decr(dne_over_n0) << std::setw(width);
    s << Ixe1_D_ave_decr(dne_over_n0) << std::setw(width);    
    s << Ix_ave_decr(dne_over_n0) << std::setw(width);
    s << LB_limit_decr() << std::setw(width);
    s << A() << std::setw(width);
    s << B(ky_decr) << std::setw(width);
    s << C(ky_decr) << std::setw(width);
    s << alpha(ky_decr) << std::setw(width);
    s << beta(ky_decr) << std::setw(width);
    s << gamma(ky_decr) << std::setw(width);
  // increasing modes values    
    s << get_my_incr() << std::setw(width);
    s << get_wR_incr()/2/M_PI << std::setw(width);
    s << std::abs(dphi_incr(dne_over_n0)) << std::setw(width);
    s << std::abs(dvxi1_incr(dne_over_n0)) << std::setw(width); 
    s << std::abs(dvxe1_EXB_incr(dne_over_n0)) << std::setw(width);  
    s << std::abs(dvxe1_D_incr(dne_over_n0)) << std::setw(width);  
    
    s << dphi_phase_incr(dne_over_n0) << std::setw(width);
    s << dvxi1_phase_incr(dne_over_n0) << std::setw(width);
    s << dvxe1_EXB_phase_incr(dne_over_n0) << std::setw(width);
    s << dvxe1_D_phase_incr(dne_over_n0) << std::setw(width);
    
    s << Jxi1_ave_incr(dne_over_n0) << std::setw(width);
    s << Jxe1_EXB_ave_incr(dne_over_n0) << std::setw(width);
    s << Jxe1_D_ave_incr(dne_over_n0) << std::setw(width);     
    s << Ixi1_ave_incr(dne_over_n0) << std::setw(width);
    s << Ixe1_EXB_ave_incr(dne_over_n0) << std::setw(width);
    s << Ixe1_D_ave_incr(dne_over_n0) << std::setw(width);   
    s << Ix_ave_incr(dne_over_n0) << std::setw(width);
    s << LB_limit_incr() << std::setw(width);
    s << A() << std::setw(width);
    s << B(ky_incr) << std::setw(width);
    s << C(ky_incr) << std::setw(width);
    s << alpha(ky_incr) << std::setw(width);
    s << beta(ky_incr) << std::setw(width);
    s << gamma(ky_incr);
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
  
  // without assumption 1/Ln << kx 
  std::complex<double> w_alter(const double& ky){
    double kp = sqrt(kx*kx+ky*ky);
    double wd = vd*ky;
    double wD = vD*ky;
    double w0 = v0*ky;  
    std::complex<double> b,c;
    b = -2.0*wi0 - (cs*cs*(kp*kp-i1*kx/Ln))/(wd-wD);
    c = wi0*wi0 + (cs*cs*(kp*kp-i1*kx/Ln))/(wd-wD)*(w0+wD-i1*nupara);
    return -b/2.0+sqrt(b*b/4.0-c);
  }
  
  void update_coherent_modes(){
    my_decr = floor(my_maxgrowth);
    ky_decr = my_decr/R0;
    w_decr = std::complex<double>(wR(ky_decr),wI(ky_decr));
    my_incr = ceil(my_maxgrowth);
    ky_incr = my_incr/R0;
    w_incr = std::complex<double>(wR(ky_incr),wI(ky_incr));
  }

  void update_coherent_modes_alter(){
    my_decr = floor(my_maxgrowth);
    ky_decr = my_decr/R0;
    w_decr = w_alter(ky_decr);
    my_incr = ceil(my_maxgrowth);
    ky_incr = my_incr/R0;
    w_incr = w_alter(ky_incr);
  }
  
  void update_wI_maxgrowth(){
    double my_min{0.01};
    double my_max{14};
    double dmy = (my_max-my_min)/100.0;
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
      std::cout << "maximum_not_found" << std::endl;
      throw "maximum_not_found";
    }
  }

  void update_wI_maxgrowth_alter(){
    double my_min{0.01};
    double my_max{14};
    double dmy = (my_max-my_min)/100.0;
    double wI_new;
    wI_maxgrowth = -10;
    for( double my = my_min; my <= my_max; my+=dmy){
      wI_new = wI_alter(my/R0);
      if (wI_new > wI_maxgrowth){
        wI_maxgrowth = wI_new; 
        my_maxgrowth = my;
      }        
    }
    if (wI_maxgrowth < 0){
      throw "negative_wI";
    }
    if (my_maxgrowth == my_min or my_maxgrowth == my_max){
      std::cout << "maximum_not_found" << std::endl;
      throw "maximum_not_found";
    }
  }
  
  void update_all(){
    T = Te*scale_factor;
    De = De0*scale_factor;
    //ne0=n0*scale_factor;
    A_plasma = 2*M_PI*R0*t;
    ne0=n0;
    cs = sqrt(T*e/mi);
    nupara = kz*kz*De;
    vd = -T/B0/Ln;
    vD = -2*T/B0/LB;
    v0 = -E/B0;
    vi0 = vix*sqrt(scale_factor);
    wi0 = kx*vi0;
    //update_wI_maxgrowth();
    //update_coherent_modes();
    update_wI_maxgrowth_alter();
    update_coherent_modes_alter();
  }
  
  protected:
  
  double mi, Te, B0, E0, R0, LB, Ln, De0, vix, kz, kx, n0, t;
  double e{1.60217657e-19};
  double scale_factor, E, T, ne0, De, cs, nupara, vd, vD, v0, vi0, wi0,
         my_maxgrowth, wI_maxgrowth, A_plasma;
  double my_decr, ky_decr, my_incr, ky_incr;
  std::complex<double> w_incr, w_decr;
  
  std::complex<double> i1{0.0,1.0};
  
};