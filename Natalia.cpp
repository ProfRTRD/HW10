#include <iostream>
#include <valarray>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <vector>
#include <algorithm>
#include <fstream>

typedef std::valarray<double> state_t;

void initial_conditions(state_t & y);

struct Functor_deriv_t {
  double W;
  double operator()(double time,double time_1, state_t & dydt,const state_t & y) {
      return dydt[0] = y[1];
  dydt[1] = -W*W*y[0];
    }
  };

template <class deriv_t, class system_t>
void solve_heun(Functor_deriv_t deriv, system_t & y, double 
tinit, double tend, double dt);

template <class deriv_t, class system_t>
void  integrate_euler(deriv_t deriv, system_t & y, double tinit, double tend, double dt);

int main(int argc, char **argv) {
 const double dt = std::atof(argv[1]);
 const double tinit = std::atof(argv[2]);
 const double tend = std::atof(argv[3]);
 const double W = std::atof(argv[4]);
 const int N = 2;
 state_t y(N);
 initial_conditions(y);
 Functor_deriv_t deriv;
 //solve_heun( deriv, y,tinit,tend, dt  );

 integrate_euler( deriv, y, tinit , tend, dt);
  
  return 0;
  }

void initial_conditions(state_t & y)
{
  y[0] = 0.9876; // x
  y[1] = 0.0; // v
}

template <class deriv_t, class system_t>
void solve_heun(Functor_deriv_t deriv, system_t & y, double tinit, double tend, double dt) 
{
 int N = y.size();
 system_t dydt(N, 0.0);
 std::ofstream fout {"data.txt"};
 fout.precision(15); fout.setf (std::ios::scientific);

 std::vector<double> K1;
 int ii = 0;   
 auto faux = [&ii, dydt](double & x){ x = dydt[ii] ; ii++; }; 
 std::for_each(K1.begin(), K1.end(), faux);

  
 std::vector<double> K2;
 int ff = 0;
  auto fau = [&ff, y, dt, dydt](double & x){ x = y[ff]+ dt*dydt[ff] ; ff++; }; //dydt[ff]=k1
 std::for_each(K2.begin(), K2.end(), fau);
  
 double nsteps = (tend - tinit)/dt;
 double time = 0;
 double time_1 = 0;
 for(int ii = 0; ii < nsteps; ++ii) {
  time = tinit + ii*dt;
  time_1 = tinit + dt + ii*dt;
  for (int ii = 0; ii < N; ++ii) {
      deriv(y, dydt, time,0);
      for (int ff = 0; ff < N; ++ff){
        deriv(y,dydt,time_1,0);
       y[ii] = y[ii] + (dt/2) * (dydt[ii]+ y[ff]+ (dt*dydt[ff])); 
        }
     fout << time_1 << "\t" << fderiv(y,dydt,time_1) <<  "\n";
    }
	
  }
 fout.close();

}


template <class deriv_t, class system_t>
void integrate_euler(Functor_deriv_t deriv, system_t & y, double tinit, double tend, double dt)
{
  int N = y.size();
  system_t dydt(N, 0.0);
  std::ofstream fout {"data.txt"};
  fout.precision(15); fout.setf (std::ios::scientific);
  double time = 0;
  int nsteps = (tend - tinit)/dt;
  for(int ii = 0; ii < nsteps; ++ii) {
	  time = tinit + ii*dt;
    deriv(y, dydt, time,0);
    for (int ii = 0; ii < N; ++ii) {
      y[ii] += dt*dydt[ii]; //y[ii] = y[ii] + dt*dydt[ii]; // EULER
   }
	fout << time << "\t" << fderiv(y,dydt,time) <<  "\n";
  }
  fout.close();
  
}



