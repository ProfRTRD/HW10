#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <valarray>

typedef std::valarray<double> state_t; // alias for state type

void initial_conditions(state_t & y);

template <class deriv_t, class system_t>
void integrate_euler(deriv_t deriv, system_t & y, double tinit, double tend, double dt);

template <class deriv_t, class system_t>
void solve_heun(deriv_t deriv, system_t & y, double tinit, double tend, double dt);


int main(int argc, char **argv)
{
  
  const double DT = std::atof(argv[1]);
  const double T0 = std::atof(argv[2]);
  const double TF = std::atof(argv[3]);
  const double W = std::atof(argv[4]);
  const int N = 2;
  state_t y(N);
  initial_conditions(y);
  
  auto fderiv = [W](const state_t & y, state_t & dydt, double t){
    dydt[0] = y[1];
    dydt[1] = -W*W*y[0];
  };
    
  integrate_euler(fderiv, y, T0, TF, DT);
  solve_heun(fderiv, y, T0, TF, DT);
  
  return 0;
}

void initial_conditions(state_t & y)
{
  y[0] = 0.9876; // x
  y[1] = 0.0; // v
}

template <class deriv_t, class system_t>
void integrate_euler(deriv_t deriv, system_t & y, double tinit, double tend, double dt)
{
  std::ofstream fout("output-euler.txt");
  int N = y.size();
  system_t dydt(N, 0.0);
  double time = 0;
  int nsteps = (tend - tinit)/dt;
  for(int ii = 0; ii < nsteps; ++ii) {
      time = tinit + ii*dt;
      deriv(y, dydt, time);
      y += dt*dydt;
      fout << time << "\t" << y[1] << "\n";
  }
  fout.close();
}

template <class deriv_t, class system_t>
void solve_heun(deriv_t deriv, system_t & y, double tinit, double tend, double dt) {
  std::ofstream fout("output-heun.txt");
  int N = y.size();
  system_t k1(N, 0.0);
  system_t k2(N, 0.0);
  system_t in(N, 0.0);
  double time = 0;
  double time1 = 0;
  int nsteps = (tend - tinit)/dt;
  for(int ii = 0; ii < nsteps; ++ii) {
      time = tinit + ii*dt;
      time1 = tinit + (ii+1)*dt;
      deriv(y, k1, time);
      in = y + k1*dt;
      deriv(in, k2, time1);
      for(int ii = 0; ii < N; ++ii) {
	y[ii] += (dt/2)*(k1[ii]+k2[ii]);
      }
      fout << time << "\t" << y[1] << std::endl;
  }
  fout.close();
}
