#pragma once
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <complex>
#include <vector>
#include "linearmodel.hpp"
//#include "matrix_phi_OLDD.hpp"

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
  
  void resize(int Ni_new,int Nj_new);
  
  private:

  void BuildInitialConditions();
  
  complex_number lambda(int i, int j);
  
  complex_number g(int& i, int j);
  
  void build_b(int& i);

  void solve_phi1(int& i);
  
  void SaveSolution();
  
  void Save2File(solution_container& f, std::ofstream& file);

  int spacesteps_number() const;
  int timesteps_number() const;
  
  int index(int j);
  
  // compute and return space derivative of variable f
  complex_number d_y(const dense_vector<complex_number>& f, int& j);
  
  private:
  friend class A_phi;
  double tf, t, dt, dy;
  double n0, n10;
  double Q0, Z0;
  int i, Ni, Nj;
  solution_container n1, vx1, vy1, phi1;
  std::vector<complex_number> sigma_phi, r_phi;
  
  dense_vector<complex_number> b;
  dense_vector<complex_number> r0;
  //A_phi* A;
 
  double fact_phi{1.0};
  
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
    
    for(int j = 0; j < number_rows(); j++){
        Assign::apply(w[j],                        v[ptr->index(j-2)]
                                             - 8.0*v[ptr->index(j-1)]
               - (ptr->dy)*12.0*(ptr->lambda(i,j))*v[j] 
                                             + 8.0*v[ptr->index(j+1)]  
                                                -  v[ptr->index(j+2)]); 
    }
    
    
    }
  
  /** Matvec forwards to MTL 's lazy mat_cvec_multiplier operator */
  template <typename Vector>
  mtl::vec::mat_cvec_multiplier<A_phi, Vector>
  operator *( const Vector & v) const{
    return mtl::vec::mat_cvec_multiplier<A_phi, Vector>(*this,v); 
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