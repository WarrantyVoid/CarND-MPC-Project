#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;


class FGEval
{
public:

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

  /**
   * @brief FGEval contructor.
   * @param coeffs Fitted polynomial coefficients
   */
  FGEval(const MPCParams &params, Eigen::VectorXd coeffs)
    : mP(params)
    , mCoeffs(coeffs)
  {

  }

  /**
   * @brief operator ()
   * @param fg Vector containing the cost and constraints
   * @param vars Vector containing the variable values (state & actuators)
   */
  void operator()(ADvector& fg, const ADvector& vars)
  {
    // The cost
    fg[0] = 0;

    // The part of the cost based on the reference state.
    for (int t = 0; t < mP.N; ++t)
    {
      fg[0] += CppAD::pow(vars[mP.cte + t], 2);
      fg[0] += CppAD::pow(vars[mP.epsi + t], 2);
      fg[0] += CppAD::pow(vars[mP.v + t] - mP.refV, 2);
    }

    // Minimize the use of actuators.
    for (int t = 0; t < mP.N - 1; ++t)
    {
      fg[0] += CppAD::pow(vars[mP.steer + t], 2);
      fg[0] += CppAD::pow(vars[mP.throttle + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int t = 1; t < mP.N - 1; ++t)
    {
      fg[0] += CppAD::pow(vars[mP.steer + t] - vars[mP.steer + t - 1], 2);
      fg[0] += CppAD::pow(vars[mP.throttle + t] - vars[mP.throttle + t - 1], 2);
    }

    // Initial constraints
    fg[mP.x + 1] = vars[mP.x];
    fg[mP.y + 1] = vars[mP.y];
    fg[mP.psi + 1] = vars[mP.psi];
    fg[mP.v + 1] = vars[mP.v];
    fg[mP.cte + 1] = vars[mP.cte];
    fg[mP.epsi + 1] = vars[mP.epsi];
    //std::cout << "[Init] x=" << fg[mP.x + 1] << " y=" << fg[mP.y + 1] << " psi=" << fg[mP.psi + 1]
    //          << " v=" << fg[mP.v + 1] << " cte=" << fg[mP.cte + 1] << " epsi=" << fg[mP.epsi + 1] << std::endl;

    // The rest of the constraints
    for (int t = 1; t < mP.N; ++t)
    {
      AD<double> x0 = vars[mP.x + t - 1];
      AD<double> y0 = vars[mP.y + t - 1];
      AD<double> psi0 = vars[mP.psi + t - 1];
      AD<double> v0 = vars[mP.v + t - 1];
      AD<double> cte0 = vars[mP.cte + t - 1];
      AD<double> epsi0 = vars[mP.epsi + t - 1];
      AD<double> steer = vars[mP.steer + t - 1];
      AD<double> throttle = vars[mP.throttle + t - 1];

      AD<double> x1 = vars[mP.x + t];
      AD<double> y1 = vars[mP.y + t];
      AD<double> psi1 = vars[mP.psi + t];
      AD<double> v1 = vars[mP.v + t];
      AD<double> cte1 = vars[mP.cte + t];
      AD<double> epsi1 = vars[mP.epsi + t];

      AD<double> f0(0.0);
      for (int d = 0; d < mCoeffs.size(); ++d)
      {
        f0 += mCoeffs(d) * CppAD::pow(x0, d);
      }
      AD<double> psides0(0.0);
      for (int d = 1; d < mCoeffs.size(); ++d)
      {
        psides0 += d * mCoeffs(d) * CppAD::pow(x0, d - 1);
      }
      psides0 = CppAD::atan(psides0);

      fg[mP.x + t + 1] = x1 - (x0 + v0 * CppAD::cos(psi0) * mP.deltaT);
      fg[mP.y + t + 1] = y1 - (y0 + v0 * CppAD::sin(psi0) * mP.deltaT);
      fg[mP.psi + t + 1] = psi1 - (psi0 + v0 * steer / mP.Lf * mP.deltaT);
      fg[mP.v + t + 1] = v1 - (v0 + throttle * mP.deltaT);
      fg[mP.cte + t + 1] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * mP.deltaT));
      fg[mP.epsi + t + 1] = epsi1 - ((psi0 - psides0) + v0 * steer / mP.Lf * mP.deltaT);

      //std::cout << "x=" << fg[mP.x + t + 1] << " y=" << fg[mP.y + t + 1] << " psi=" << fg[mP.psi + t + 1]
      //          << " v=" << fg[mP.v + t + 1] << " cte=" << fg[mP.cte + t + 1] << " epsi=" << fg[mP.epsi + t + 1] << std::endl;
    }
  }

private:
  /** Helper indices */
  const MPCParams mP;

  /** Fitted polynomial coefficients */
  Eigen::VectorXd mCoeffs;
};


//=================================================================================


MPCParams::MPCParams(double refV, double Lf, int N, double deltaT)
{
  this->refV = refV;
  this->Lf = Lf;
  this->N = N;
  this->deltaT = deltaT;

  // Calculate helper ndices
  x = 0;
  y = x + N;
  psi = y + N;
  v = psi + N;
  cte = v + N;
  epsi = cte + N;
  steer = epsi + N;
  throttle = steer + N - 1;
}


MPC::MPC(const MPCParams &params)
: mP(params)
{

}


MPC::~MPC()
{

}


bool MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs, std::vector<Eigen::VectorXd> &actuations)
{
  typedef CPPAD_TESTVECTOR(double) Dvector;
  const int nVars = state.size() * mP.N  + 2 * (mP.N - 1);
  const int nConstraints = state.size() * mP.N;

  // Values for variables
  Dvector vars(nVars);
  for (int i = 0; i < nVars; ++i)
  {
    vars[i] = 0;
  }
  vars[mP.x] = state(0);
  vars[mP.y] = state(1);
  vars[mP.psi] = state(2);
  vars[mP.v] = state(3);
  vars[mP.cte] = state(4);
  vars[mP.epsi] = state(5);

  // Lower and upper limits for variables
  Dvector varsLowBound(nVars);
  Dvector varsUppBound(nVars);
  for (int i = 0; i < mP.steer; ++i)
  {
    varsLowBound[i] = std::numeric_limits<double>::lowest();
    varsUppBound[i] = std::numeric_limits<double>::max();
  }
  for (int i = mP.steer; i < mP.throttle; ++i)
  {
    varsLowBound[i] = deg2rad(-25.0);
    varsUppBound[i] = deg2rad(25.0);
  }
  for (int i = mP.throttle; i < nVars; ++i)
  {
    varsLowBound[i] = -1.0;
    varsUppBound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  Dvector constraintsLowBound(nConstraints);
  Dvector constraintsUppBound(nConstraints);
  for (int i = 0; i < nConstraints; ++i)
  {
    constraintsLowBound[i] = 0;
    constraintsUppBound[i] = 0;
  }
  constraintsLowBound[mP.x] = state[0];
  constraintsUppBound[mP.x] = state[0];
  constraintsLowBound[mP.y] = state[1];
  constraintsUppBound[mP.y] = state[1];
  constraintsLowBound[mP.psi] = state[2];
  constraintsUppBound[mP.psi] = state[2];
  constraintsLowBound[mP.v] = state[3];
  constraintsUppBound[mP.v] = state[3];
  constraintsLowBound[mP.cte] = state[4];
  constraintsUppBound[mP.cte] = state[4];
  constraintsLowBound[mP.epsi] = state[5];
  constraintsUppBound[mP.epsi] = state[5];

  // object that computes objective and constraints
  FGEval fgEval(mP, coeffs);

  // options for IPOPT solver
  std::string options;

  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FGEval>(
        options, vars, varsLowBound, varsUppBound,
        constraintsLowBound, constraintsUppBound,
        fgEval, solution);

  actuations.clear();
  if (solution.status == CppAD::ipopt::solve_result<Dvector>::success)
  {
    // Cost
    auto cost = solution.obj_value;
    std::cout << "Cost " << cost << std::endl;

    // Set actuation result
    for(int t = 0; t < mP.N; ++t)
    {
      Eigen::VectorXd actuation(8);
      actuation << solution.x[mP.x + t]
                 , solution.x[mP.y + t]
                 , solution.x[mP.psi + t]
                 , solution.x[mP.v + t]
                 , solution.x[mP.cte + t]
                 , solution.x[mP.epsi + t]
                 , (t < mP.N - 1) ? rad2deg(solution.x[mP.steer + t]) / 25.0 : 0.0
                 , (t < mP.N - 1) ? solution.x[mP.throttle + t] : 0.0;
      actuations.push_back(actuation);
    }
  }
  else
  {
    std::cout << "Error." << std::endl;
  }
  return solution.status == CppAD::ipopt::solve_result<Dvector>::success;
}

