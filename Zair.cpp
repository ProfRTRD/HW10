#include <iostream>
#include <cstdlib>
#include <vector>
#include <numeric>
#include <cmath>
#include <valarray>


typedef std::valarray<double> state_t; 

void initial_conditions(state_t & y);
void print(const state_t & y, double time);

struct functor_E{
  double W=0;
  void operator()(const state_t & y, state_t & dydt, double t){
  dydt[0] = y[1];
  dydt[1] = -W*W*y[0];
  }
};
template <class deriv_t, class system_t, class printer_t>
void integrate_euler(deriv_t deriv, system_t & y, double tinit, double tend, double dt, printer_t writer);
template <class deriv_t, class system_t>
void integrate_huell(deriv_t deriv, system_t & y, double por, double sas, double d_time);

int main()
{
  functor_E fderiv;
  fderiv.W = 3.1237;
  double DT = 0.001;
  double TF = 4*2*M_PI/fderiv.W;
  int N = 2;
  double W = 3.1237;
  state_t y(N);
  initial_conditions(y);
  integrate_euler(fderiv, y, 0.0, TF, DT, print);
  integrate_huell(fderiv, y, 0.0, TF, DT);
  return 0;
}

void initial_conditions(state_t & y)
{
  y[0] = 0.9876; // x
  y[1] = 0.0; // v
}

void print(const state_t & y, double time)
{
  std::cout << time << "\t" << y[0] << "\t" << y[1] << std::endl;
}
template <class deriv_t, class system_t, class printer_t>
void integrate_euler(deriv_t deriv, system_t & y, double tinit, double tend, double dt, printer_t writer)
{
  int N = y.size();
  system_t dydt(N, 0.0);
  double time = 0;
  int nsteps = (tend - tinit)/dt;
  for(int ii = 0; ii < nsteps; ++ii) {
	time = tinit + ii*dt;
    deriv(y, dydt, time);
    for (int ii = 0; ii < N; ++ii) {
      y[ii] += dt*dydt[ii]; 
    }

  }
}

template <class deriv_t, class system_t>
void integrate_huell(deriv_t deriv, system_t & y, double por, double sas, double d_time)
{
  functor_E fderiv;
  int Z = y.size();
  system_t dydt(Z,0.0);
  double tempo= 0.0;
  int num_pas= (sas-por)/d_time;
  for(int ii=0; ii <num_pas;ii++){
    tempo= por + ii*d_time;
    fderiv(y,dydt,tempo);
    for(int i= 0; i<Z; i++){
      y[i] += d_time*dydt[i];
    }
  }
  
}

