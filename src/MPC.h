#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"


/**
 * @brief Structure containing helper indices.
 */
struct CPPADIndices
{
  CPPADIndices(int N, double deltaT);

  /** The number of states to evaluate */
  int N;

  /** Delta time between two successive actuations */
  int deltaT;

  int x;
  int y;
  int psi;
  int v;
  int cte;
  int epsi;
  int steer;
  int throttle;

  /** Length from front to center of gravity */
  static const double LF;

  /** Reference speed not to cause costs */
  static const double REF_V;
};


class MPC
{

public:

  MPC(int N, double deltaT);

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

  /** Helper indices */
  CPPADIndices mIdx;
};

#endif /* MPC_H */
