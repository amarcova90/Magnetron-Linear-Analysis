#include <fstream>
#include <complex>
#include <iomanip>
#include <vector>
#include <string> 
#include "nonlinearmodel.hpp"
#include <iomanip>      // std::setprecision

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
              
              
              wR_decr = get_wR_decr();
              vp = wR_decr/ky_decr;
              chix =kx/ky_decr; chiz = kz/ky_decr;
              chix2 = chix*chix;
              chiz2 = chiz*chiz;
              chix4 = chix2*chix2;
              chiz4 = chiz2*chiz2;
              
              
              // initialize the memory for solutions
              resize(Ni_,Nj_);              
              Z0 = n0/B0/Ln - 2.0*n0/B0/LB; Q0 = E0/B0 + 2.0*Te/B0/LB;
              wce = e*B0/me;
              nu = 0.0*1.0e5;  
              rL2 = Te/B0/wce;
              
              

              
              cn0 = Q0/vp;
              cn1 = i1*kx/B0/vp*Te;
              cn2 = Z0/n0*Te/vp;
              cn3 = (i1*kx/B0-2.0/B0/Ln)*Te/vp;
              cn4 = -(nupara/wR_decr + D_*(chix4 + chiz4));
              
              cl0 = -(i1*kx*vi0/wR_decr + D_*(chix4 + chiz4));
              cl1 = e/mi*Te/vp/vp;
              cl2 = -e/mi*Te*(kx*kx/wR_decr/wR_decr + kz*kz/wR_decr/wR_decr);
              cl3 = chix2 + chiz2;
              cl4 = 2.0*cl3*cl3;
              
              cp0 = nupara/wR_decr - i1*kx*vi0/wR_decr;
              cp1 = i1*chix/Ln/ky_decr;

              std::cout << "vp = " << vp << std::endl; 
              std::cout << "cn0 = " << cn0 << std::endl; 
              std::cout << "cn1 = " << cn1 << std::endl; 
              std::cout << "cn2 = " << cn2 << std::endl;
              std::cout << "cn3 = " << cn3 << std::endl;
              std::cout << "cn4 = " << cn4 << std::endl;
              std::cout << "cl0 = " << cl0 << std::endl;
              std::cout << "cl1 = " << cl1 << std::endl;
              std::cout << "cl2 = " << cl2 << std::endl;
              std::cout << "cl3 = " << cl3 << std::endl;
              std::cout << "cl4 = " << cl4 << std::endl;
              std::cout << "cp0 = " << cp0 << std::endl;
              std::cout << "cp1 = " << cp1 << std::endl;
              //getchar();
}

void NonlinearModel::BuildInitialConditions(){
  double y;
  //double dy_ = (2*M_PI*R0 - dy)/ Nj;
  for( int j = 0; j < Nj+1; j++){
    y = dy*j;
    n1[0][j]   = n10/n0              * ( cos(y) + i1*sin(y) );
    vx1[0][j]  = dvxi1_decr(n10/n0)  * ( cos(y) + i1*sin(y) );
    vy1[0][j]  = dvyi1_decr(n10/n0)  * ( cos(y) + i1*sin(y) );
    vz10[j]    = dvzi1_decr(n10/n0)  * ( cos(y) + i1*sin(y) );
  }
  
  for( int j = 0; j < Nj+1; j++){
    lap1[0][j]    = - i1*kx*vx1[0][j]/wR_decr - d_y(vy1[0],j)/vp - i1*kz*vz10[j]/wR_decr;
    //std::cout << "lap1[0][j] = "<< lap1[0][j] << std::endl;
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
  clock_t t = clock();
  solve_chi1(0);
  solve_phi1(0);
  for(int i = 1; i < Ni+1; i++){
    march_RK4(i);
/*     std::cout << std::setprecision(5);
    std::cout <<  "simulation time [sec]: " << i*dt/wR_decr; */
    std::cout << std::setprecision(4);
    std::cout << "Progress [%] = "<< (double)i/Ni*100 << std::endl;      
  }
  t = clock() - t;
  print_simulation_time(((double)t)/CLOCKS_PER_SEC);
  //std::cout << "simulation time: " 
  //            << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()*1e-9
  //            << " s" << std::endl;
  SaveSolution();
}

  
complex_number NonlinearModel::lambda(int i, int j){
  return fact_lin*cn1*d_y(n1[i], j) / (cn2+cn3*n1[i][j] );

}

/* complex_number NonlinearModel::g(int& i, int j){
  return -( (Q0-fact_lin*d_y(chi1[i],j))*d_y(n1[i],j) + i1*kx*vi0*n1[i][j] +
            (kx*kx + kz*kz)*chi1[i][j]*n1[i][j] - i1*n0/Ln*kx*chi1[i][j] - 
            nupara*n1[i][j] - (n0 + n1[i][j])*lap1[i][j]  ) /
         ( Z0 + fact_lin*( i1*kx/B0 - 2.0/B0/LB )*n1[i][j] );
} */

complex_number NonlinearModel::g(int& i, int j){
  return (-cn0*d_y(n1[i],j) - cl3*chi1[i][j]*n1[i][j] + cp0*n1[i][j] + cp1*chi1[i][j] + 
          d_y(chi1[i],j)*d_y(n1[i],j) + (1.0+n1[i][j])*lap1[i][j])/(cn2+cn3*n1[i][j]);
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
  return  (cn0-cn1*phi1[i][j])*d_y(n1[i],j) + (cn2 + cn3*n1[i][j])*d_y(phi1[i],j) +
           cn4*n1[i][j] - D_*d4_y(n1[i],j);
} 

complex_number NonlinearModel::f_lap1(int i, int j){
  return  cl0*lap1[i][j] + cl1*d2_y(phi1[i],j) + cl2*phi1[i][j] + d2_y(chi1[i],j)*d2_y(chi1[i],j) +
          d_y(chi1[i],j)*d3_y(chi1[i],j) - cl3*(3.0*d_y(chi1[i],j)*d_y(chi1[i],j) +
          chi1[i][j]*d2_y(chi1[i],j)) + cl4*chi1[i][j]*chi1[i][j] - D_*d4_y(lap1[i],j);               
}  

void NonlinearModel::march_RK4(int i){
  
  // Update k1
  //#pragma omp parallel for
  for(int j = 0; j < spacesteps_number(); j++){
    k1_n1[j] =  f_n1(i-1,j)  *dt;
    k1_lap1[j] = f_lap1(i-1,j)*dt;
  }    
  
  //#pragma omp parallel for
  for(int j = 0; j < spacesteps_number(); j++){
    n1[i][j] = n1[i-1][j] + 0.5*k1_n1[j]; // temp for k2
    lap1[i][j] = lap1[i-1][j] + 0.5*k1_lap1[j];
  }    
  
  // Update phi1 and place temporarely in i
  solve_chi1(i);
  solve_phi1(i);     
  
  // Update k2
  //#pragma omp parallel for
  for(int j = 0; j < spacesteps_number(); j++){
    k2_n1[j] =  f_n1(i,j)*dt;
    k2_lap1[j] = f_lap1(i,j)*dt;    
  }

  //#pragma omp parallel for
  for(int j = 0; j < spacesteps_number(); j++){           
    n1[i][j] = n1[i-1][j] + 0.5*k2_n1[j]; // temp for k3 
    lap1[i][j] = lap1[i-1][j] + 0.5*k2_lap1[j];
  }
  
  // Update phi1 and place temporarely in i
  solve_chi1(i);
  solve_phi1(i);
  // Update k3
  //#pragma omp parallel for
  for(int j = 0; j < spacesteps_number(); j++){
    k3_n1[j] = f_n1(i,j)*dt;
    k3_lap1[j] = f_lap1(i,j)*dt;
  }

  //#pragma omp parallel for
  for(int j = 0; j < spacesteps_number(); j++){  
    n1[i][j] = n1[i-1][j] + k3_n1[j]; // temp for k3 
    lap1[i][j] = lap1[i-1][j] + k3_lap1[j];
  }
  
  // Update phi1 and place temporarely in i
  solve_chi1(i);
  solve_phi1(i);
  
  //#pragma omp parallel for
  for(int j = 0; j < spacesteps_number(); j++){
    k4_n1[j] = f_n1(i,j)*dt;
    k4_lap1[j] = f_lap1(i,j)*dt;
  }
  
  // Compute next step
  //#pragma omp parallel for
  for(int j = 0; j < spacesteps_number(); j++){
    n1[i][j] =  n1[i-1][j]  + (k1_n1[j]  + 2*k2_n1[j]  + 2*k3_n1[j]  + k4_n1[j])/6.0;
    lap1[i][j] = lap1[i-1][j] + (k1_lap1[j] + 2*k2_lap1[j] + 2*k3_lap1[j] + k4_lap1[j])/6.0;
   }    
  solve_chi1(i);
  solve_phi1(i);

// updating velocities
  //#pragma omp parallel for
  for(int j = 0; j < spacesteps_number(); j++){
    vx1[i][j] =  -i1*chix*vp*chi1[i][j];
    vy1[i][j] =  -vp*d_y(chi1[i],j);
   }    
  
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

  dy = 2*M_PI*my_decr/(Nj+1.0);
  dt = tf/Ni*wR_decr;
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
  
  phi1_amp = dphi_decr(n10/n0)/Te;
  
  for( int j = 0; j < Nj+1; j++){
    y = dy*j;
    space[0][j] = dy*j;
    exact_linear[0][j] = phi1_amp*( cos(y) + i1*sin(y) );
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
  
  noisy_iteration<double> iter_lin{r0, Nj, err_rel, err_abs};
  
  bicgstab(A, solution_linear[0], b, P, iter_lin);
  //bicg(A, solution_linear[0], b, P, iter_lin);
  
  // compute nonlinear solution using solver
  fact_lin = 1.0; // turning on nonlinear terms    
  
  //r0 = dense_vector<complex_number>(spacesteps_number(),1.0);
  r0 = b - A*x0;

  
  noisy_iteration<double> iter_nonlin{r0, Nj, err_rel, err_abs};
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
    exact_linear[0][i] = n10/n0 * exp(get_wI_decr()/wR_decr*t) * ( cos(wR_decr*t) + i1*sin(get_wI_decr()*t) );
    solution_nonlinear[0][i] = n1[i][0];
  }    
  


  fact_lin = 0.0;
  solve();
  for( int i = 0; i < Ni+1; i++){
    solution_linear[0][i] = n1[i][0];
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


void print_simulation_time(double ts){
  //ts = ts*1e-9; //ns to s
  int h,m,s;
  h = (int)(ts/3600.0);
  m = (int)((ts-h*3600.0)/60.0);
  s = (int)(ts-h*3600.0-m*60.0);
  std::cout << "simulation time = " << h << " hours " << m << " minutes " << s << " seconds." << std::endl;
}