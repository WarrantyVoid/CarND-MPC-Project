#ifndef TOOLS_H
#define TOOLS_H

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

#endif /* TOOLS_H */
