#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"


static constexpr double pi()
{
  return M_PI;
}

static double deg2rad(double x)
{
  return x * pi() / 180;
}

static double rad2deg(double x)
{
  return x * 180 / pi();
}


/**
 * @brief Structure containing helper indices.
 */
struct MPCParams
{
  MPCParams(double refV, double Lf, int N, double deltaT);

  /** The number of states to evaluate */
  int N;

  /** Delta time between two successive actuations */
  double deltaT;

  /** Length from front to center of gravity */
  double Lf;

  /** Reference speed not to cause costs */
  double refV;

  // Helper indices
  int x;
  int y;
  int psi;
  int v;
  int cte;
  int epsi;
  int steer;
  int throttle;
};


class MPC
{

public:

  MPC(const MPCParams &params);

  virtual ~MPC();

public:

  /**
   * @brief Solve the model given an initial state and polynomial coefficients.
   * @param state Current car state (x, y, psi, v, cte, epsi)
   * @param coeffs Coefficients of planned trajectory path
   * @param actuation The calculated actuation (car states + steer, throttle)
   * @return True, if a solution was found, false otherwise
   */
  bool Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs, std::vector<Eigen::VectorXd> &actuations);

private:

  /** Parametrization */
  const MPCParams mP;
};

#endif /* MPC_H */
