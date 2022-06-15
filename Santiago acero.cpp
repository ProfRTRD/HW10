
#include <iostream>
#include <valarray>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>

typedef std::valarray<double> state_t; // alias for state type

// constants


void initial_conditions(state_t & y);
void print(const state_t & y, double time);
void fderiv(const state_t & y, state_t & dydt, double t);
template <class deriv_t, class system_t>
void integrate_euler(deriv_t deriv, system_t  y, double tinit, double tend, double dt);
template <class deriv_t, class system_t>
void solve_heun(deriv_t deriv, system_t  y, double tinit, double tend, double dt);

double W=0;

int main(int argc, char *argv[])
{
    const double DT =std::atof( argv[0]);
    const double TI = std::atof( argv[1]);
    const double TF = std::atof( argv[2]);
    W = std::atof( argv[3]);

    const int N = 2;
    std::valarray<double> y(N);
    initial_conditions(y);
    integrate_euler(fderiv, y, TI, TF, DT);
    solve_heun(fderiv, y, TI, TF, DT);

    return 0;
}

void initial_conditions(state_t & y)
{
  y[0] = 0.9876; // x
  y[1] = 0.0; // v
}

void print(const state_t & y, double time, char name)
{

  std::cout << time << "\t" << y[0] << "\t" << y[1] << std::endl;
}

void fderiv(const state_t & y, state_t & dydt, double t)
{
  dydt[0] = y[1];
  dydt[1] = -W*W*y[0];
}



template <class deriv_t, class system_t>
void integrate_euler(deriv_t deriv, system_t  y, double tinit, double tend, double dt)
{

  int N = y.size();
  system_t dydt(N, 0.0);
  double time = 0;
  int nsteps = (tend - tinit)/dt;
  for(int ii = 0; ii < nsteps; ++ii) {
	time = tinit + ii*dt;
    deriv(y, dydt, time);
    for (int ii = 0; ii < N; ++ii) {
        y[ii] += dt*dydt[ii]; //y[ii] = y[ii] + dt*dydt[ii]; // EULER
    }
    std::ofstream fout("euler.txt");

    for(auto n:y)
    {
        fout<<time<<" "<<n<<std::endl;
    }
  }
}

template <class deriv_t, class system_t>
void solve_heun(deriv_t deriv, system_t  y, double tinit, double tend, double dt)
{
  int N = y.size();
  system_t k1(N, 0.0);
  system_t k2(N, 0.0);
  double time = 0;
  double time1 = 0;
  int nsteps = (tend - tinit)/dt;


  for(int ii = 0; ii < nsteps; ++ii) {
	time1 = tinit + ii*dt;
    deriv(y, k1, time);
    deriv(y, k2, time1);

    for (int ii = 0; ii < N; ++ii) {
      y[ii] += (dt/2)*(k1[ii]+k2[ii]);
    }
    time=time1;
    std::ofstream fout("heun.txt");
    for(auto n:y)
    {
        fout<<n<<" "<<time<<std::endl;
    }

  }

}
