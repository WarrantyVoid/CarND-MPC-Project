#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

namespace
{
  constexpr double pi()
  {
    return M_PI;
  }

  double deg2rad(double x)
  {
    return x * pi() / 180;
  }

  double rad2deg(double x)
  {
    return x * 180 / pi();
  }
}


//=================================================================================


class FGEval
{
public:

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

  /**
   * @brief FGEval contructor.
   * @param coeffs Fitted polynomial coefficients
   */
  FGEval(const CPPADIndices &idx, Eigen::VectorXd coeffs)
    : mIdx(idx)
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
    for (int t = 0; t < mIdx.N; ++t)
    {
      fg[0] += CppAD::pow(vars[mIdx.cte + t], 2);
      fg[0] += CppAD::pow(vars[mIdx.epsi + t], 2);
      fg[0] += CppAD::pow(vars[mIdx.v + t] - CPPADIndices::REF_V, 2);
    }

    // Minimize the use of actuators.
    for (int t = 0; t < mIdx.N - 1; ++t)
    {
      fg[0] += CppAD::pow(vars[mIdx.steer + t], 2);
      fg[0] += CppAD::pow(vars[mIdx.throttle + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int t = 1; t < mIdx.N - 1; ++t)
    {
      fg[0] += CppAD::pow(vars[mIdx.steer + t] - vars[mIdx.steer + t - 1], 2);
      fg[0] += CppAD::pow(vars[mIdx.throttle + t] - vars[mIdx.throttle + t - 1], 2);
    }

    // Initial constraints
    fg[mIdx.x + 1] = vars[mIdx.x];
    fg[mIdx.y + 1] = vars[mIdx.y];
    fg[mIdx.psi + 1] = vars[mIdx.psi];
    fg[mIdx.v + 1] = vars[mIdx.v];
    fg[mIdx.cte + 1] = vars[mIdx.cte];
    fg[mIdx.epsi + 1] = vars[mIdx.epsi];

    // The rest of the constraints
    for (int t = 1; t < mIdx.N; ++t)
    {
      AD<double> x0 = vars[mIdx.x + t - 1];
      AD<double> y0 = vars[mIdx.y + t - 1];
      AD<double> psi0 = vars[mIdx.psi + t - 1];
      AD<double> v0 = vars[mIdx.v + t - 1];
      AD<double> cte0 = vars[mIdx.cte + t - 1];
      AD<double> epsi0 = vars[mIdx.epsi + t - 1];
      AD<double> steer = vars[mIdx.steer + t - 1];
      AD<double> throttle = vars[mIdx.throttle + t - 1];

      AD<double> x1 = vars[mIdx.x + t];
      AD<double> y1 = vars[mIdx.y + t];
      AD<double> psi1 = vars[mIdx.psi + t];
      AD<double> v1 = vars[mIdx.v + t];
      AD<double> cte1 = vars[mIdx.cte + t];
      AD<double> epsi1 = vars[mIdx.epsi + t];

      AD<double> f0 = mCoeffs[0] + mCoeffs[1] * x0;// + mCoeffs[2] * CppAD::pow(x0, 2);
      AD<double> psides0 = CppAD::atan(mCoeffs[1]/* + 2.0 * mCoeffs[2] * x0*/);

      fg[mIdx.x + t + 1] = x1 - (x0 + v0 * CppAD::cos(psi0) * mIdx.deltaT);
      fg[mIdx.y + t + 1] = y1 - (y0 + v0 * CppAD::sin(psi0) * mIdx.deltaT);
      fg[mIdx.psi + t + 1] = psi1 - (psi0 + v0 * mIdx.steer / CPPADIndices::LF * mIdx.deltaT);
      fg[mIdx.v + t + 1] = v1 - (v0 + throttle * mIdx.deltaT);
      fg[mIdx.cte + t + 1] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * mIdx.deltaT));
      fg[mIdx.epsi + t + 1] = epsi1 - ((psi0 - psides0) + v0 * steer / CPPADIndices::LF * mIdx.deltaT);
    }
  }

private:
  /** Helper indices */
  CPPADIndices mIdx;

  /** Fitted polynomial coefficients */
  Eigen::VectorXd mCoeffs;
};


//=================================================================================


const double CPPADIndices::LF = 2.67;
const double CPPADIndices::REF_V = 40;

CPPADIndices::CPPADIndices(int N, double deltaT)
{
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


MPC::MPC(int N, double deltaT)
: mIdx(N, deltaT)
{

}


MPC::~MPC()
{

}


bool MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs, std::vector<Eigen::VectorXd> &actuations)
{
  typedef CPPAD_TESTVECTOR(double) Dvector;
  const int nVars = state.size() * mIdx.N  + 2 * (mIdx.N - 1);
  const int nConstraints = state.size() * mIdx.N;

  // Values for variables
  Dvector vars(nVars);
  for (int i = 0; i < nVars; ++i)
  {
    vars[i] = 0;
  }
  vars[mIdx.x] = state[0];
  vars[mIdx.y] = state[1];
  vars[mIdx.psi] = state[2];
  vars[mIdx.v] = state[3];
  vars[mIdx.cte] = state[4];
  vars[mIdx.epsi] = state[5];

  // Lower and upper limits for variables
  Dvector varsLowBound(nVars);
  Dvector varsUppBound(nVars);
  for (int i = 0; i < mIdx.steer; ++i)
  {
    varsLowBound[i] = std::numeric_limits<double>::min();
    varsUppBound[i] = std::numeric_limits<double>::max();
  }
  for (int i = mIdx.steer; i < mIdx.throttle; ++i)
  {
    varsLowBound[i] = deg2rad(-25.0);
    varsUppBound[i] = deg2rad(25.0);
  }
  for (int i = mIdx.throttle; i < nVars; ++i)
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
  constraintsLowBound[mIdx.x] = state[0];
  constraintsUppBound[mIdx.x] = state[0];
  constraintsLowBound[mIdx.y] = state[1];
  constraintsUppBound[mIdx.y] = state[1];
  constraintsLowBound[mIdx.psi] = state[2];
  constraintsUppBound[mIdx.psi] = state[2];
  constraintsLowBound[mIdx.v] = state[3];
  constraintsUppBound[mIdx.v] = state[3];
  constraintsLowBound[mIdx.cte] = state[4];
  constraintsUppBound[mIdx.cte] = state[4];
  constraintsLowBound[mIdx.epsi] = state[5];
  constraintsUppBound[mIdx.epsi] = state[5];

  // object that computes objective and constraints
  FGEval fgEval(mIdx, coeffs);

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

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  actuations.clear();
  if (solution.status == CppAD::ipopt::solve_result<Dvector>::success)
  {
    // Return actuation result
    for(int t = 0; t < mIdx.N; ++t)
    {
      Eigen::VectorXd actuation(2);
      actuation << solution.x[mIdx.x + t]
                 , solution.x[mIdx.y + t]
                 , solution.x[mIdx.psi + t]
                 , solution.x[mIdx.v + t]
                 , solution.x[mIdx.cte + t]
                 , solution.x[mIdx.epsi + t]
                 , (t < mIdx.N - 1) ? rad2deg(solution.x[mIdx.steer + t]) / 25.0 : 0.0
                 , (t < mIdx.N - 1) ? solution.x[mIdx.throttle + t] : 0.0;
    }
  }
  return solution.status == CppAD::ipopt::solve_result<Dvector>::success;
}

