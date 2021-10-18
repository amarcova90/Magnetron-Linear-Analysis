#pragma once
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <boost/thread.hpp>
#include "boost/filesystem/operations.hpp"
#include <complex>
#include <vector>
#include "linearmodel.hpp"
#include <time.h>

using namespace mtl;
using namespace itl;

using complex_number = std::complex<double>;
using solution_container = std::vector<dense_vector<complex_number>>;
//typedef std::chrono::high_resolution_clock Clock;

class NonlinearModel: public LinearModel {
  
  public:
  
  NonlinearModel(double& mi_,  double& Te_, double& B0_, double& E0_ ,
              double& R0_,  double& LB_, double& Ln_, double& De0_,
              double& vix_, double& kz_, double& kx_, double& n0_,
              double& t_, int& Ni_, int& Nj_, double& tf_,
              double& n10_n0_);
               
  void solve();
  
  void solver_validation();

  void RK4_validation();
  
  void resize(int Ni_new,int Nj_new);

  private:

  void BuildInitialConditions();
  
  complex_number f_n1(int i, int j);
  
  complex_number f_lap1(int i, int j);
  
  complex_number lambda(int i, int j);
  
  complex_number g(int& i, int j);
  
  void build_b(int& i);
  
  void build_b_chi(int& i);

  void solve_phi1(int i);
  
  void solve_chi1(int i);
  
  void march_RK4(int i);
  
  void SaveSolution();
  
  void Save2File(solution_container& f, std::ofstream& file);

  int spacesteps_number() const;
  int timesteps_number() const;
  
  int index(int j);
  
  // compute and return space derivative of variable f
  complex_number d_y(const dense_vector<complex_number>& f, int& j);
  complex_number d2_y(const dense_vector<complex_number>& f, int& j);
  complex_number d3_y(const dense_vector<complex_number>& f, int& j);
  complex_number d4_y(const dense_vector<complex_number>& f, int& j);
  //complex_number d4_y2(const dense_vector<complex_number>& f, int& j);
  
  private:
  friend class A_phi;
  friend class A_chi;
  double tf, dt, dy;
  double n0, n10;
  double Q0, Z0;
  double wce, rL2, nu;
  double D_ = 5.0e-3;
  double err_rel{1e-3};
  double err_abs{1e-6};
  double vp, wR_decr; //phase velocity
  double chix, chix2, chix4, chiz, chiz2, chiz4;
  complex_number cn0, cn1, cn2, cn3, cn4;
  complex_number cl0, cl1, cl2, cl3, cl4;
  complex_number cp0, cp1;
  int i, Ni, Nj;
  solution_container n1, vx1, vy1, phi1, chi1, lap1;
  
  std::vector<complex_number> vz10;
  
  std::vector<complex_number> sigma_phi, r_phi;
  
  std::vector<complex_number> k1_n1, k1_vx1, k1_vy1, k1_lap1;
  std::vector<complex_number> k2_n1, k2_vx1, k2_vy1, k2_lap1;
  std::vector<complex_number> k3_n1, k3_vx1, k3_vy1, k3_lap1;
  std::vector<complex_number> k4_n1, k4_vx1, k4_vy1, k4_lap1;
  
  dense_vector<complex_number> b, b_chi;
  dense_vector<complex_number> r0;
  //A_phi* A;
 
  double fact_lin{1.0};
  
  
};



class A_phi{
  
  public:
  
  // default constructor
  A_phi(){
    ptr = nullptr;
    i = -1; 
  }
  
  A_phi(class NonlinearModel* ptr_, int i_):ptr{ptr_}, i{i_}{}
  

  
  /** Helper function to perform multiplication . Allows for delayed
  * evaluation of results .
  * Assign::apply(a, b) resolves to an assignment operation such as
  * a += b, a -= b, or a = b.
  * @pre @a size (v) == size (w) */
  
  template <typename VectorIn , typename VectorOut , typename Assign >
  void mult ( const VectorIn& v, VectorOut& w, Assign ) const{
    //double temp = 0.0;
    double dy = ptr->dy;
/*     
    for(int j = 0; j < number_rows(); j++){
        Assign::apply(w[j],                        v[ptr->index(j-2)]
                                             - 8.0*v[ptr->index(j-1)]
               - (ptr->dy)*12.0*(ptr->lambda(i,j))*v[j] 
                                             + 8.0*v[ptr->index(j+1)]  
                                                -  v[ptr->index(j+2)]); 
    } */
    //#pragma omp parallel for
    for(int j = 0; j < number_rows(); j++){
        Assign::apply(w[j],                        v[ptr->index(j-2)]
                                             - 8.0*v[ptr->index(j-1)]
               - 12.0*(ptr->lambda(i,j))*v[j]*dy
                                             + 8.0*v[ptr->index(j+1)]
                                                -  v[ptr->index(j+2)]); 
    }    
    
/*     double dy = ptr->dy;
    complex_number f;
    for(int j = 0; j < number_rows(); j++){
        f = dy*12.0*(ptr->lambda(i,j));
        Assign::apply(w[j],                        v[ptr->index(j-2)]/f
                                             - 8.0*v[ptr->index(j-1)]/f
               - v[j] 
                                             + 8.0*v[ptr->index(j+1)]/f  
                                                -  v[ptr->index(j+2)]/f ); 
    }   */  
    
    }
  
  /** Matvec forwards to MTL 's lazy mat_cvec_multiplier operator */
  template <typename Vector>
  mtl::vec::mat_cvec_multiplier<A_phi, Vector>
  operator *( const Vector & v) const{
    return mtl::vec::mat_cvec_multiplier<A_phi, Vector>(*this,v); 
  }

  dense_vector<double>  operator [](int j) const{
    dense_vector<double> out(number_columns(),0.0);
    out[ptr->index(j-2)] = 1.0;
    out[ptr->index(j-1)] = -8.0;
    out[ptr->index(j)] = - std::abs( (ptr->dy)*12.0*(ptr->lambda(i,j)) );
    out[ptr->index(j+1)] = 8.0;
    out[ptr->index(j+2)] = -1.0;
    
    return out; 
  }
  
  void print_address(){
    std::cout << (*this).number_rows() << std::endl;
    
  }
  int number_rows() const{return ptr->spacesteps_number();}
  int number_columns() const{return ptr->spacesteps_number();}
  int size() const{return number_rows()*number_columns();}
  
  private:
  NonlinearModel* ptr;
  int  i;
};


class A_chi{
  
  public:
  
  // default constructor
  A_chi(){
    ptr = nullptr;
  }
  
  A_chi(class NonlinearModel* ptr_):ptr{ptr_}{}
  

  
  /** Helper function to perform multiplication . Allows for delayed
  * evaluation of results .
  * Assign::apply(a, b) resolves to an assignment operation such as
  * a += b, a -= b, or a = b.
  * @pre @a size (v) == size (w) */
  
  template <typename VectorIn , typename VectorOut , typename Assign >
  void mult ( const VectorIn& v, VectorOut& w, Assign ) const{
    //double temp = 0.0;
    double dy = ptr->dy;
    complex_number cl3 = ptr->cl3;
    
      //#pragma omp parallel for
      for(int j = 0; j < number_rows(); j++){
          Assign::apply(w[j],                        -v[ptr->index(j-2)]
                                               + 16.0*v[ptr->index(j-1)]
                              - (30.0+12.0*dy*dy*cl3)*v[j]
                                               + 16.0*v[ptr->index(j+1)]
                                                     -v[ptr->index(j+2)]); 
      }    
      
    }
  
  /** Matvec forwards to MTL 's lazy mat_cvec_multiplier operator */
  template <typename Vector>
  mtl::vec::mat_cvec_multiplier<A_chi, Vector>
  operator *( const Vector & v) const{
    return mtl::vec::mat_cvec_multiplier<A_chi, Vector>(*this,v); 
  }

/*   dense_vector<double>  operator [](int j) const{
    double dy = ptr->dy;
    dense_vector<double> out(number_columns(),0.0);
    out[ptr->index(j-2)] = -1.0;
    out[ptr->index(j-1)] = 16.0;
    out[ptr->index(j)] = - (30.0+12.0*dy*dy*(ptr->kx*ptr->kx+ptr->kz*ptr->kz));
    out[ptr->index(j+1)] = 16.0;
    out[ptr->index(j+2)] = -1.0;
    
    return out; 
  } */
  
  void print_address(){
    std::cout << (*this).number_rows() << std::endl;
    
  }
  int number_rows() const{return ptr->spacesteps_number();}
  int number_columns() const{return ptr->spacesteps_number();}
  int size() const{return number_rows()*number_columns();}
  
  private:
  NonlinearModel* ptr;
};


/** The number of elements in the matrix . */
template <typename mat>
inline std::size_t size (const mat & A){ return A.size();}

/** The number of rows in the matrix . */
template <typename mat>
inline std::size_t num_rows (const mat & A){ return A.number_rows();}


/** The number of columns in the matrix . */
template <typename mat>
inline std::size_t num_cols (const mat & A){ return A.number_columns();}

void print_simulation_time(double ts);