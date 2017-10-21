#ifndef TOOLS_H
#define TOOLS_H

#include <math.h>

class Tools
{
public:
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

  template<typename TAngle>
  static TAngle calculateAngleDelta(TAngle  a1, TAngle a2)
  {
    return Tools::normalizeAngle(Tools::normalizeAngle(a2) - Tools::normalizeAngle(a1));
  }

  template<typename TAngle>
  static TAngle normalizeAngle(TAngle a)
  {
    while (a > M_PI)
    {
      a -= 2.0f * M_PI;
    }
    while (a < -M_PI)
    {
      a += 2.0f * M_PI;
    }
    return a;
  }

  template<typename TValue>
  static bool isZero(TValue f)
  {
    return fabs(f) < std::numeric_limits<double>::epsilon();
  }
};

#endif /* TOOLS_H */
