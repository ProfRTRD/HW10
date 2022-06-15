#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <valarray>
#include <fstream>

typedef std::valarray<double> state_t;

void initial_conditions(state_t & y);
template <class deriv_t, class system_t>
void integrate_euler(deriv_t dlambda, system_t & y, double tinit, double tend, double dt, double w);
template <class deriv_t, class system_t>
void solve_heun(deriv_t dlambda, system_t & y, double tinit, double tend, double dt, double w);

int main(int argc, char *argv[])
{
  const double DT = std::atof(argv[1]);
  const double TO = std::atof(argv[2]);
  const double TF = std::atof(argv[3]);
  const double W = std::atof(argv[4]);
  const int N = 2;
  state_t y(N);
  initial_conditions(y);
  auto dlambda = [W](const state_t &y,state_t &dydt, double time){
    dydt[0] = y[1];
    dydt[1] = -W*W*y[0];
  };
  
  integrate_euler(dlambda,y, TO, TF, DT, W);
  solve_heun(dlambda,y, TO, TF, DT, W);

  return 0;
}

void initial_conditions(state_t & y)
{
  y[0] = 0.9876; // x
  y[1] = 0.0; // v
}

template <class deriv_t, class system_t>
void integrate_euler (deriv_t dlambda, system_t & y, double tinit, double tend, double dt, double w)
{ std::ofstream fout("output-euler.txt");
  
  int N = y.size();
  system_t dydt(N, 0.0);
  double time = 0;
  int nsteps = (tend - tinit)/dt;
  for(int ii = 0; ii < nsteps; ++ii) {
    time = tinit + ii*dt;
    dlambda(y,dydt,time);
    for (int ii = 0; ii < N; ++ii) {
      y[ii] += dt*dydt[ii];
    } ;
    fout << time <<"\t"<< y[1] <<"\n";
  };
  fout.close();
}

template <class deriv_t, class system_t>
void solve_heun (deriv_t dlambda, system_t & y, double tinit, double tend, double dt, double w)
{ std::ofstream fout("output-heun.txt");
  
  int N = y.size();
  system_t k1(N, 0.0);
  system_t k2(N, 0.0);
  double time = 0;
  double timeplus1 = 0;
  int nsteps = (tend - tinit)/dt;
  for(int ii = 0; ii < nsteps; ++ii) {
    time = tinit + ii*dt;
    timeplus1 = tinit + (ii+1)*dt;
    dlambda(y,k1,time); //k1
    system_t ydtk1 = y + dt*k1;
    dlambda(ydtk1,k2,timeplus1); //k2
    for (int ii = 0; ii < N; ++ii) {
      y[ii] += (dt/2)*(k1[ii]+k2[ii]);
    } ;
    fout << time <<"\t"<< y[1] <<"\n";
  };
  fout.close();
}
