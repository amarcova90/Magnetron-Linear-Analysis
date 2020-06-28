#pragma once
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <complex>
#include <vector>
#include "linearmodel.hpp"

using namespace mtl;
using namespace itl;

using complex_number = std::complex<double>;
using solution_container = std::vector<dense_vector<complex_number>>;

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
  
  complex_number f_vx1(int i, int j);
  
  complex_number f_vy1(int i, int j);
  
  complex_number lambda(int i, int j);
  
  complex_number g(int& i, int j);
  
  void build_b(int& i);

  void solve_phi1(int i);
  
  void march_RK4(int i);
  
  void SaveSolution();
  
  void Save2File(solution_container& f, std::ofstream& file);

  int spacesteps_number() const;
  int timesteps_number() const;
  
  int index(int j);
  
  // compute and return space derivative of variable f
  complex_number d_y(const dense_vector<complex_number>& f, int& j);
  
  private:
  friend class A_phi;
  double tf, dt, dy;
  double n0, n10;
  double Q0, Z0;
  double err_rel{1e-3};
  double err_abs{1e-6};
  int i, Ni, Nj;
  solution_container n1, vx1, vy1, phi1;
  std::vector<complex_number> sigma_phi, r_phi;
  
  std::vector<complex_number> k1_n1, k1_vx1, k1_vy1;
  std::vector<complex_number> k2_n1, k2_vx1, k2_vy1;
  std::vector<complex_number> k3_n1, k3_vx1, k3_vy1;
  std::vector<complex_number> k4_n1, k4_vx1, k4_vy1;
  
  dense_vector<complex_number> b;
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

    for(int j = 0; j < number_rows(); j++){
        Assign::apply(w[j],                        v[ptr->index(j-2)]/dy/12
                                             - 8.0*v[ptr->index(j-1)]/dy/12
               - (ptr->lambda(i,j))*v[j] 
                                             + 8.0*v[ptr->index(j+1)]/dy/12
                                                -  v[ptr->index(j+2)]/dy/12); 
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
  int i;
};


/** The number of elements in the matrix . */
inline std::size_t size (const A_phi & A){ return A.size();}

/** The number of rows in the matrix . */
inline std::size_t num_rows (const A_phi & A){ return A.number_rows();}


/** The number of columns in the matrix . */
inline std::size_t num_cols (const A_phi & A){ return A.number_columns();}