#include <fstream>
#include <complex>
#include <iomanip>
#include <vector>
#include <string> 
#include "nonlinearmodel.hpp"


namespace mtl{
  template <>
  struct Collection<A_phi>
  {
      typedef double value_type;
      typedef int    size_type;
  };

  namespace ashape {
      template <> struct ashape_aux<A_phi> 
      {       typedef nonscal type;    };
  }
}

namespace mtl{
  template <>
  struct Collection<A_chi>
  {
      typedef double value_type;
      typedef int    size_type;
  };

  namespace ashape {
      template <> struct ashape_aux<A_chi> 
      {       typedef nonscal type;    };
  }
}

NonlinearModel::NonlinearModel(double& mi_,  double& Te_, double& B0_, double& E0_ ,
              double& R0_,  double& LB_, double& Ln_, double& De0_,
              double& vix_, double& kz_, double& kx_, double& n0_,
              double& t_, int& Ni_, int& Nj_, double& tf_,
              double& n10_n0_): LinearModel(mi_, Te_, B0_, E0_, R0_, LB_, 
              Ln_, De0_, vix_, kz_, kx_, n0_, t_), 
              tf{tf_}, n0{n0_}, n10{n10_n0_*n0} {
              
              // initialize the memory for solutions
              resize(Ni_,Nj_);              
              Z0 = n0/B0/Ln - 2.0*n0/B0/LB; Q0 = E0/B0 + 2.0*Te/B0/LB;
              wce = e*B0/me;
              nu = 0.0*1.0e5;  
              rL2 = Te/B0/wce;
              std::cout << "wce = " << wce << std::endl; 
              std::cout << "rL2 = " << rL2 << std::endl; 
              std::cout << "nupara = " << nupara << std::endl; 
}

void NonlinearModel::BuildInitialConditions(){
  double y;
  //double dy_ = (2*M_PI*R0 - dy)/ Nj;
  for( int j = 0; j < Nj+1; j++){
    y = dy*j;
    n1[0][j]   = n10                 * ( cos(my_decr/R0*y) +i1*sin(my_decr/R0*y) );
    vx1[0][j]  = dvxi1_decr(n10/n0)  * ( cos(my_decr/R0*y) +i1*sin(my_decr/R0*y) );
    vy1[0][j]  = dvyi1_decr(n10/n0)  * ( cos(my_decr/R0*y) +i1*sin(my_decr/R0*y) );
    vz10[j]    = dvzi1_decr(n10/n0)  * ( cos(my_decr/R0*y) +i1*sin(my_decr/R0*y) );
  }
  
  for( int j = 0; j < Nj+1; j++){
    y = dy*j;
    lap1[0][j]    = - i1*kx*vx1[0][j] - d_y(vy1[0],j) - i1*kz*vz10[j];
  }
  
/*   dense_vector<complex_number> test3(Nj+1);
   */
  
/*   double y1{-dy};
  
  for( int j = 0; j < Nj+1; j++){
    y1 += dy;
    y = dy*j;
    std::cout << "y = " << y << "; y1 = " << y1 << std::endl;
    std::cout << "y^4 = " << y*y*y*y << "; y1 = " << y1*y1*y1*y1 << std::endl;
    test3[j]    = cos(y) ;
  }  
  
  for( int j = 10; j < 20; j++){
    std::cout << "d4 = "<< d4_y(test3,j) << std::endl;
  }  */

}
  
void NonlinearModel::solve(){
  solve_chi1(0);
  solve_phi1(0);
  for(int i = 1; i < Ni+1; i++){
    march_RK4(i);
    std::cout << "simulation time [sec]: " << dt*i  << std::endl;      
  }
  SaveSolution();
}

  
complex_number NonlinearModel::lambda(int i, int j){
  return fact_lin*i1*kx/B0*d_y(n1[i], j) / (Z0 + (i1*kx/B0 - 2.0/B0/LB)*n1[i][j] );

}

/* complex_number NonlinearModel::g(int& i, int j){
  return -( (Q0-fact_lin*d_y(chi1[i],j))*d_y(n1[i],j) + i1*kx*vi0*n1[i][j] +
            (kx*kx + kz*kz)*chi1[i][j]*n1[i][j] - i1*n0/Ln*kx*chi1[i][j] - 
            nupara*n1[i][j] - (n0 + n1[i][j])*lap1[i][j]  ) /
         ( Z0 + fact_lin*( i1*kx/B0 - 2.0/B0/LB )*n1[i][j] );
} */

complex_number NonlinearModel::g(int& i, int j){
  return -( (Q0-fact_lin*d_y(chi1[i],j))*d_y(n1[i],j) + i1*kx*vi0*n1[i][j] +
            (kx*kx )*chi1[i][j]*n1[i][j] - i1*n0/Ln*kx*chi1[i][j] - 
            nupara*n1[i][j] - (n0 + n1[i][j])*(d2_y(chi1[i],j)-kx*kx*chi1[i][j]) ) /
         ( Z0 + fact_lin*( i1*kx/B0 - 2.0/B0/LB )*n1[i][j] );
}
  
/* complex_number NonlinearModel::g(int& i, int j){
  return -( (Q0+fact_lin*vy1[i][j])*d_y(n1[i],j) + (n0+fact_lin*n1[i][j])*(i1*kx*vx1[i][j]
              + d_y(vy1[i],j))  + nu*rL2*d2_y(n1[i],j)
            + (vi0+fact_lin*vx1[i][j])*i1*kx*n1[i][j] - (nupara + 0.0*nu*rL2*kx*kx)*n1[i][j] 
            + vx1[i][j]*n0/Ln ) /
         ( Z0 + fact_lin*( i1*kx/B0 - 2.0/B0/LB )*n1[i][j] ); 
} */

void NonlinearModel::build_b(int& i){
  for(int j = 0; j < spacesteps_number(); j++){
                b[j] = 12.0*dy*g(i,j);
      }
}

void NonlinearModel::build_b_chi(int& i){
  for(int j = 0; j < spacesteps_number(); j++){
                b_chi[j] = 12.0*dy*dy*lap1[i][j];
      }
}

void NonlinearModel::solve_phi1(int i){
  A_phi A(this,i);
  build_b(i);                       
  pc::identity<A_phi>        P(A);  
  
  if (i == 0){
    r0 = dense_vector<complex_number>(spacesteps_number(),0.0);
  }else{
    r0 = b - A*phi1[i-1];
  }
  
  basic_iteration<double> iter{r0, Nj, err_rel, err_abs};    
  //noisy_iteration<double> iter{r0, 1000, err_rel, err_abs};
  bicgstab(A, phi1[i], b, P, iter);

}

void NonlinearModel::solve_chi1(int i){
  A_chi A_c(this);
  build_b_chi(i);                       
  pc::identity<A_chi>        P(A_c);  
  
  if (i == 0){
    r0 = dense_vector<complex_number>(spacesteps_number(),0.0);
  }else{
    r0 = b_chi - A_c*chi1[i-1];
  }
    
  basic_iteration<double> iter{r0, Nj, err_rel, err_abs};
  bicgstab(A_c, chi1[i], b_chi, P, iter);

}
  
complex_number NonlinearModel::f_n1(int i, int j){
/*   return  (Q0 - fact_lin* i1*kx/B0*phi1[i][j])*d_y(n1[i],j) + 
          (Z0 + fact_lin* i1*kx/B0*n1[i][j] - fact_lin* 2*n1[i][j]/B0/LB)*d_y(phi1[i],j) -
           nupara*n1[i][j] ; */
  return  (Q0 - fact_lin* i1*kx/B0*phi1[i][j])*d_y(n1[i],j) + 
          (Z0 + fact_lin* i1*kx/B0*n1[i][j] - fact_lin* 2.0*n1[i][j]/B0/LB)*d_y(phi1[i],j)
          - nu*n0*rL2/Te*d2_y(phi1[i],j) + nu*n0*rL2/Te*kx*kx*phi1[i][j]
          + nu*rL2*d2_y(n1[i],j) - (nupara+0.0*nu*rL2*kx*kx)*n1[i][j] -
          D_*(kx*kx*kx*kx + kz*kz*kz*kz)*n1[i][j] - D_*d4_y(n1[i],j);
}

complex_number NonlinearModel::f_vx1(int i, int j){
    return  - fact_lin* vy1[i][j]*d_y(vx1[i],j) - vi0*i1*kx*vx1[i][j] - 
             e*i1*kx/mi*phi1[i][j];
}  

complex_number NonlinearModel::f_vy1(int i, int j){
  return  -vi0*i1*kx*vy1[i][j] - fact_lin* vx1[i][j]*i1*kx*vy1[i][j] - 
           e/mi*d_y(phi1[i],j);               
}  

complex_number NonlinearModel::f_lap1(int i, int j){
  return  -vi0*i1*kx*lap1[i][j] + e/mi*( d2_y(phi1[i],j) - (kx*kx + kz*kz)*phi1[i][j] ) + 
           d2_y(chi1[i],j)*d2_y(chi1[i],j) + d_y(chi1[i],j)*d3_y(chi1[i],j) - 
           3.0*(kx*kx + kz*kz)*d_y(chi1[i],j)*d_y(chi1[i],j) - 
           (kx*kx + kz*kz)*chi1[i][j]*d2_y(chi1[i],j) + 
           2.0*(kx*kx + kz*kz)*(kx*kx + kz*kz)*chi1[i][j]*chi1[i][j] -
           D_*(kx*kx*kx*kx + kz*kz*kz*kz)*lap1[i][j] - D_*d4_y(lap1[i],j);               
}  

void NonlinearModel::march_RK4(int i){
  
  // Update k1
  for(int j = 0; j < spacesteps_number(); j++){
    k1_n1[j] =  f_n1(i-1,j)  *dt;
    k1_vx1[j] = f_vx1(i-1,j) *dt;
    k1_vy1[j] = f_vy1(i-1,j)*dt;
    k1_lap1[j] = f_lap1(i-1,j)*dt;
  }    
  
  for(int j = 0; j < spacesteps_number(); j++){
    n1[i][j] = n1[i-1][j] + 0.5*k1_n1[j]; // temp for k2
    vx1[i][j] = vx1[i-1][j] + 0.5*k1_vx1[j]; // temp for k2
    vy1[i][j] = vy1[i-1][j] + 0.5*k1_vy1[j]; // temp for k2 
    lap1[i][j] = lap1[i-1][j] + 0.5*k1_lap1[j];
  }    
  
  // Update phi1 and place temporarely in i
  solve_chi1(i);
  solve_phi1(i);     
  
  // Update k2
  for(int j = 0; j < spacesteps_number(); j++){
    k2_n1[j] =  f_n1(i,j)*dt;
    k2_vx1[j] = f_vx1(i,j)*dt;    
    k2_vy1[j] = f_vy1(i,j)*dt;    
    k2_lap1[j] = f_lap1(i,j)*dt;    
  }

  for(int j = 0; j < spacesteps_number(); j++){           
    n1[i][j] = n1[i-1][j] + 0.5*k2_n1[j]; // temp for k3 
    vx1[i][j] = vx1[i-1][j] + 0.5*k2_vx1[j]; // temp for k3 
    vy1[i][j] = vy1[i-1][j] + 0.5*k2_vy1[j]; // temp for k3 
    lap1[i][j] = lap1[i-1][j] + 0.5*k2_lap1[j];
  }
  
  // Update phi1 and place temporarely in i
  solve_chi1(i);
  solve_phi1(i);
  // Update k3
  for(int j = 0; j < spacesteps_number(); j++){
    k3_n1[j] = f_n1(i,j)*dt;
    k3_vx1[j] = f_vx1(i,j)*dt;
    k3_vy1[j] = f_vy1(i,j)*dt;
    k3_lap1[j] = f_lap1(i,j)*dt;
  }

  for(int j = 0; j < spacesteps_number(); j++){  
    n1[i][j] = n1[i-1][j] + k3_n1[j]; // temp for k3 
    vx1[i][j] = vx1[i-1][j] + k3_vx1[j]; // temp  
    vy1[i][j] = vy1[i-1][j] + k3_vy1[j]; // temp for k4
    lap1[i][j] = lap1[i-1][j] + k3_lap1[j];
  }
  
  // Update phi1 and place temporarely in i
  solve_chi1(i);
  solve_phi1(i);
  
  for(int j = 0; j < spacesteps_number(); j++){
    k4_n1[j] = f_n1(i,j)*dt;
    k4_vx1[j] = f_vx1(i,j)*dt;
    k4_vy1[j] = f_vy1(i,j)*dt;
    k4_lap1[j] = f_lap1(i,j)*dt;
  }
  
  // Compute next step
  for(int j = 0; j < spacesteps_number(); j++){
    n1[i][j] =  n1[i-1][j]  + (k1_n1[j]  + 2*k2_n1[j]  + 2*k3_n1[j]  + k4_n1[j])/6;
    vx1[i][j] = vx1[i-1][j] + (k1_vx1[j] + 2*k2_vx1[j] + 2*k3_vx1[j] + k4_vx1[j])/6;
    vy1[i][j] = vy1[i-1][j] + (k1_vy1[j] + 2*k2_vy1[j] + 2*k3_vy1[j] + k4_vy1[j])/6;
    lap1[i][j] = lap1[i-1][j] + (k1_lap1[j] + 2*k2_lap1[j] + 2*k3_lap1[j] + k4_lap1[j])/6;
   }    
  solve_chi1(i);
  solve_phi1(i);
  
}
  
void NonlinearModel::resize(int Ni_new,int Nj_new){

  Nj = Nj_new; Ni = Ni_new;
  n1.resize(Ni+1, dense_vector<complex_number>(Nj+1));
  vx1.resize(Ni+1, dense_vector<complex_number>(Nj+1));
  vy1.resize(Ni+1, dense_vector<complex_number>(Nj+1));
  chi1.resize(Ni+1, dense_vector<complex_number>(Nj+1));
  lap1.resize(Ni+1, dense_vector<complex_number>(Nj+1));
  vz10.resize(Nj+1);
  phi1.resize(Ni+1, dense_vector<complex_number>(Nj+1));
  
  // initialize temporary containers
  sigma_phi.resize(Nj+1);
  r_phi.resize(Nj+1);
  
  k1_n1.resize(Nj+1);
  k1_vx1.resize(Nj+1);
  k1_vy1.resize(Nj+1);
  k1_lap1.resize(Nj+1);
  
  k2_n1.resize(Nj+1);
  k2_vx1.resize(Nj+1);
  k2_vy1.resize(Nj+1);
  k2_lap1.resize(Nj+1);

  k3_n1.resize(Nj+1);
  k3_vx1.resize(Nj+1);
  k3_vy1.resize(Nj+1);
  k3_lap1.resize(Nj+1);

  k4_n1.resize(Nj+1);
  k4_vx1.resize(Nj+1);
  k4_vy1.resize(Nj+1);    
  k4_lap1.resize(Nj+1);

  dy = 2*M_PI*R0/(Nj+1.0);
  dt = tf/Ni;
  b.change_dim(spacesteps_number());
  b_chi.change_dim(spacesteps_number());
  r0.change_dim(spacesteps_number());
  BuildInitialConditions();
}
  
void NonlinearModel::solver_validation(){
  // linearize potential equation
  
  solution_container space(1,dense_vector<complex_number>(Nj+1));
  solution_container exact_linear(1,dense_vector<complex_number>(Nj+1));
  solution_container solution_linear(1,dense_vector<complex_number>(Nj+1));
  solution_container solution_nonlinear(1,dense_vector<complex_number>(Nj+1));
  
  std::ofstream space_file("validation_files/val_space_Nj" + std::to_string (Nj) + ".txt");
  std::ofstream exact_linear_file("validation_files/val_exact_linear_Nj" + std::to_string (Nj) + ".txt");
  std::ofstream solution_linear_file("validation_files/val_solution_linear_Nj" + std::to_string (Nj) + ".txt");
  std::ofstream solution_nonlinear_file("validation_files/val_solution_nonlinear_Nj" + std::to_string (Nj) + ".txt");
  
  // compute exact solution from linear analysis
  double y;
  //double dy_;
  complex_number phi1_amp;
/*     phi1_amp = - R0*n0/Z0/my_decr * ( ( Q0*my_decr/R0 + vi0*kx + i1*De0*kz*kz) * n10/n0 + \
                                     kx*dvxi1_decr(n10/n0) + my_decr*dvyi1_decr(n10/n0)/R0 - \
                                     i1*dvxi1_decr(n10/n0)/Ln ); */
  
  phi1_amp = dphi_decr(n10/n0);
  
  for( int j = 0; j < Nj+1; j++){
    y = dy*j;
    space[0][j] = dy*j;
    exact_linear[0][j] = phi1_amp  * ( cos(my_decr/R0*y) +i1*sin(my_decr/R0*y) );
  }
  
  // compute linear solution using solver
  fact_lin = 0.0; // setting nonlinear terms to 0
  int i = 0; // first timestep (t = 0)
  A_phi A(this,i);
  build_b(i);                       
  pc::identity<A_phi>        P(A);   
  
  //r0 = dense_vector<complex_number>(spacesteps_number(),1.0);
  
  dense_vector<complex_number> x0 = exact_linear[0];
  
  r0 = b - A*x0;
  
  noisy_iteration<double> iter_lin{r0, 10000, err_rel, err_abs};
  
  bicgstab(A, solution_linear[0], b, P, iter_lin);
  //bicg(A, solution_linear[0], b, P, iter_lin);
  
  // compute nonlinear solution using solver
  fact_lin = 1.0; // turning on nonlinear terms    
  
  //r0 = dense_vector<complex_number>(spacesteps_number(),1.0);
  r0 = b - A*x0;

  
  noisy_iteration<double> iter_nonlin{r0, 10000, err_rel, err_abs};
  build_b(i);                          
  bicgstab(A, solution_nonlinear[0], b, P, iter_nonlin);    
  //bicg(A, solution_nonlinear[0], b, P, iter_nonlin);   
  
  // save solutions
  Save2File(space,space_file);
  Save2File(exact_linear,exact_linear_file);
  Save2File(solution_linear,solution_linear_file);
  Save2File(solution_nonlinear,solution_nonlinear_file);
  
}
  
void NonlinearModel::RK4_validation(){
  solution_container time(1,dense_vector<complex_number>(Ni+1));
  solution_container exact_linear(1,dense_vector<complex_number>(Ni+1));
  solution_container solution_linear(1,dense_vector<complex_number>(Ni+1));
  solution_container solution_nonlinear(1,dense_vector<complex_number>(Ni+1));

  std::ofstream time_file("validation_files/valRK4_time_Ni" + std::to_string (Ni) + ".txt");
  std::ofstream exact_linear_file("validation_files/valRK4_exact_linear_Ni" + std::to_string (Ni) + ".txt");
  std::ofstream solution_linear_file("validation_files/valRK4_solution_linear_Ni" + std::to_string (Ni) + ".txt");
  std::ofstream solution_nonlinear_file("validation_files/valRK4_solution_nonlinear_Ni" + std::to_string (Ni) + ".txt");
 
  fact_lin = 1.0;
  solve();
  
  for( int i = 0; i < Ni+1; i++){
    t = dt*i;
    time[0][i] = dt*i;
    exact_linear[0][i] = n10/n0 * exp(get_wI_decr()*t) * ( cos(get_wR_decr()*t) + i1*sin(get_wI_decr()*t) );
    solution_nonlinear[0][i] = n1[i][0]/n0;
  }    
  


  fact_lin = 0.0;
  solve();
  for( int i = 0; i < Ni+1; i++){
    solution_linear[0][i] = n1[i][0]/n0;
  }      
  //std::cout << "wI = " << get_wI_decr() << std::endl;

  // save solutions
  Save2File(time,time_file);
  Save2File(exact_linear,exact_linear_file);
  Save2File(solution_linear,solution_linear_file);        
  Save2File(solution_nonlinear,solution_nonlinear_file);    
    
}  

void NonlinearModel::SaveSolution(){
 std::ofstream s_phi1("phi1.txt");
 Save2File(phi1,s_phi1);
 std::ofstream s_n1("n1.txt");
 Save2File(n1,s_n1);
 std::ofstream s_vx1("vx1.txt");
 Save2File(vx1,s_vx1);  
 std::ofstream s_vy1("vy1.txt");
 Save2File(vy1,s_vy1);    
 std::ofstream s_chi1("chi1.txt");
 Save2File(chi1,s_chi1); 
 std::ofstream s_lap1("lap1.txt");
 Save2File(lap1,s_lap1); 
}  
  
void NonlinearModel::Save2File(solution_container& v, std::ofstream& file){
  unsigned int width{30};
  file.setf(std::ios::scientific, std::ios::floatfield);
  file.precision(17); 
  for (unsigned int i = 0; i< v.size(); i++){
    for (unsigned int j = 0; j < size(v[0]); j++){
      file << v[i][j].real() << std::setw(width);
    }
    file << "\n";
  }
  file.close(); 
}

int NonlinearModel::spacesteps_number() const{ return Nj + 1;}
int NonlinearModel::timesteps_number() const{ return Ni + 1;}
  
int NonlinearModel::index(int j){
  if (j < 0)
    return Nj + 1 + j;
  else if (j > Nj)
    return j - Nj - 1;
  else
    return j;
}
  
// compute and return space derivatives of variable f
complex_number NonlinearModel::d_y(const dense_vector<complex_number>& f, int& j){
  return (f[index(j-2)]-8.0*f[index(j-1)]+8.0*f[index(j+1)]-f[index(j+2)])/12.0/dy;
}

complex_number NonlinearModel::d2_y(const dense_vector<complex_number>& f, int& j){
  return (-f[index(j+2)]+16.0*f[index(j+1)]-30.0*f[index(j)]+16.0*f[index(j-1)]-f[index(j-2)])/12.0/dy/dy;  
}

complex_number NonlinearModel::d3_y(const dense_vector<complex_number>& f, int& j){
  return (-f[index(j+3)]+8.0*f[index(j+2)]-13.0*f[index(j+1)]+13.0*f[index(j-1)]-8.0*f[index(j-2)]+f[index(j-3)])/8.0/dy/dy/dy;  
}

complex_number NonlinearModel::d4_y(const dense_vector<complex_number>& f, int& j){
  double dy4 = dy*dy*dy*dy; 
  return (-f[index(j+3)]+12.0*f[index(j+2)]-39.0*f[index(j+1)]+56.0*f[index(j)]-39.0*f[index(j-1)]+12.0*f[index(j-2)]-f[index(j-3)])/6.0/dy4;  
}  

/* complex_number NonlinearModel::d4_y(const dense_vector<complex_number>& f, int& j){
  std::complex<double> dy4 = 6.0*dy*dy*dy*dy;
  std::complex<double> f_p3 =  -f[index(j+3)];
  std::complex<double> f_p2 =  12.0*f[index(j+2)];
  std::complex<double> f_p1 =  -39.0*f[index(j+1)];
  std::complex<double> f_ =  56.0*f[index(j)];
  std::complex<double> f_m1 =  -39.0*f[index(j-1)];
  std::complex<double> f_m2 =  12.0*f[index(j-2)];
  std::complex<double> f_m3 =  -f[index(j-3)];
  //std::cout << (-f_p3+f_p2-f_p1+f_-f_m1+f_m2-f_m3) << std::endl;
  //std::cout << dy4 << std::endl;
  std::complex<double> out_long = (f_p3+f_p2+f_p1+f_+f_m1+f_m2+f_m3)/dy4;  
  
  return (complex_number)out_long;
}   */

/* complex_number NonlinearModel::d4_y2(const dense_vector<complex_number>& f, int& j){
  long double dy4 = dy*dy*dy*dy;
  return (1.0L*(std::complex<long double>)f[index(j+2)]-4.0L*(std::complex<long double>)f[index(j+1)]+6.0L*(std::complex<long double>)f[index(j)]-4.0L*(std::complex<long double>)f[index(j-1)]+1.0L*(std::complex<long double>)f[index(j-2)])/dy4;  
}   */

