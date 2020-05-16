#include <fstream>
#include <complex>
#include <iomanip>
#include <vector>

class NonlinearModel: public LinearModel {
  public:
  NonlinearModel(double& mi_,  double& Te_, double& B0_, double& E0_ ,
              double& R0_,  double& LB_, double& Ln_, double& De0_,
              double& vix_, double& kz_, double& kx_, double& n0_,
              double& t_, unsigned& Ni_, unsigned& Nj_, double& tf_,
              double& n10_n0_, unsigned& m_n10_) 
              : LinearModel(mi_, Te_, B0_, E0_, R0_, LB_, 
              Ln_, De0_, vix_, kz_, kx_, n0_, t_), 
              tf{tf_}, n10{n10_n0_*n0}, Ni{Ni_}, Nj{Nj_}, m_n10{m_n10_}{
                
                // initialize the memory for solutions
                n1.resize(Ni, std::vector<std::complex<double>>(Nj));
                vx1.resize(Ni, std::vector<std::complex<double>>(Nj));
                vy1.resize(Ni, std::vector<std::complex<double>>(Nj));
                phi1.resize(Ni, std::vector<std::complex<double>>(Nj));
                sigma_phi.resize(Nj);
                r_phi.resize(Nj);
                
                Q0 = n0/B0/Ln; Z0 = E0/B0 + 2*Te/B0/LB;
                t = 0; // simulation time initialized to 0
                dy = 2*M_PI*R0/(Nj - 1);
                dt = tf/(Ni - 1);
                
                BuildInitialConditions();
                
                for (auto it = n1[0].begin(); it != n1[0].end(); ++it)
                  std::cout << (*it).real() << std::endl;                
                
                //std::cout << "dy = " << dy << "; dt =" << dt <<std::endl;
/*              print to check
                for (auto it = n1.begin(); it != n1.end(); ++it){
                  for (auto it2 = (*it).begin(); it2 != (*it).end(); ++it2)
                    std::cout << (*it2).real() << " ";
                  std::cout << std::endl;
                }*/
              } 
  void BuildInitialConditions(){
    for( unsigned j = 0; j < Nj; j++){
      double y = dy*j;
      n1[0][j].real( n10 * sin(m_n10/R0*y) );
    }
  }

  private:
  
  double tf, t, dt, dy;
  double n10;
  double Q0, Z0;
  unsigned Ni, Nj, m_n10;
  std::vector<std::vector<std::complex<double>>> n1, vx1, vy1, phi1;
  std::vector<std::complex<double>> sigma_phi, r_phi;
};