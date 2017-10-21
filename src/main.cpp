#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <list>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"
#include "tools.h"

// for convenience
using json = nlohmann::json;

namespace
{
  // Checks if the SocketIO event has JSON data.
  // If there is data the JSON object in string format will be returned,
  // else the empty string "" will be returned.
  std::string hasData(std::string s)
  {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.find_last_of("]");
    if (found_null != std::string::npos)
    {
      return "";
    }
    else if (b1 != std::string::npos && b2 != std::string::npos)
    {
      return s.substr(b1, b2 - b1 + 1);
    }
    return "";
  }

  // Evaluate a polynomial.
  double polyEval(const Eigen::VectorXd &coeffs, double x)
  {
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); i++)
    {
      result += coeffs[i] * pow(x, i);
    }
    return result;
  }

  // Fit a polynomial.
  // Adapted from
  // https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
  Eigen::VectorXd polyFit(const Eigen::VectorXd &xvals, const Eigen::VectorXd &yvals, int order)
  {
    assert(xvals.size() == yvals.size());
    assert(order >= 1 && order <= xvals.size() - 1);
    Eigen::MatrixXd A(xvals.size(), order + 1);

    for (int i = 0; i < xvals.size(); i++)
    {
      A(i, 0) = 1.0;
    }

    for (int j = 0; j < xvals.size(); j++)
    {
      for (int i = 0; i < order; i++)
      {
        A(j, i + 1) = A(j, i) * xvals(j);
      }
    }

    Eigen::HouseholderQR<Eigen::MatrixXd> Q = A.householderQr();
    Eigen::VectorXd result = Q.solve(yvals);
    return result;
  }
}


//=================================================================================


class TimeTracker
{
 public:

  TimeTracker(unsigned N)
    : mN(N)
    , mBeginPoint(std::chrono::high_resolution_clock::now())
    , mMeasurements()
  {

  }

  void begin()
  {
    mBeginPoint = std::chrono::high_resolution_clock::now();
  }

  void end()
  {
    std::chrono::time_point<std::chrono::system_clock> endPoint= std::chrono::high_resolution_clock::now();
    unsigned elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(endPoint - mBeginPoint).count();
    mMeasurements.push_back(elapsed);
    if (mMeasurements.size() > mN)
    {
      mMeasurements.pop_front();
    }
  }

  double getAvgElapsedMs() const
  {
    int n = mMeasurements.size();
    if (n == 0)
    {
      return 0;
    }
    else
    {
      return std::accumulate(mMeasurements.begin(), mMeasurements.end(), 0.0) / 1000.0 / n;
    }
  }

private:
  unsigned mN;
  std::chrono::time_point<std::chrono::system_clock> mBeginPoint;
  std::list<unsigned> mMeasurements;
};

//=================================================================================


int main()
{
  static const double targetSpeed = 80.0;
  uWS::Hub h;
  TimeTracker tracker(3);
  MPCParams params(targetSpeed, 2.67, 10.0, 8.0 / (targetSpeed + 4.0));
  MPC mpc(params);

  h.onMessage([&tracker, &params, &mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode)
  {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    std::string sdata = std::string(data).substr(0, length);
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2')
    {
      std::string s = hasData(sdata);
      if (s != "")
      {
        auto j = json::parse(s);
        std::string event = j[0].get<std::string>();
        if (event == "telemetry")
        {
          json msgJson;
          tracker.begin();

          // Fetch car state from msg
          double x = j[1]["x"];
          double y = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double steer = j[1]["steering_angle"];
          double throttle = j[1]["throttle"];

          // Fetch planned path from message
          std::vector<double> px = j[1]["ptsx"];
          std::vector<double> py = j[1]["ptsy"];

          // Translate path to vehicle coord system
          assert(px.size() == py.size());
          Eigen::VectorXd refXPlanned(px.size());
          Eigen::VectorXd refYPlanned(py.size());
          for (int i = 0; i < static_cast<int>(px.size()); ++i)
          {
            refXPlanned[i] = cos(-psi) * (px[i] - x) - sin(-psi) * (py[i] - y);
            refYPlanned[i] = sin(-psi) * (px[i] - x) + cos(-psi) * (py[i] - y);
          }

          // Match 3rd order polynom
          Eigen::VectorXd refCoeffs = polyFit(refXPlanned, refYPlanned, 3);

          // Sample waypoints from reference polynom for display
          std::vector<double> refWaypointsX;
          std::vector<double> refWaypointsY;
          for (int i = 0; i < 125; i += 5)
          {
            refWaypointsX.push_back(i);
            refWaypointsY.push_back(polyEval(refCoeffs, i));
          }
          msgJson["next_x"] = refWaypointsX;
          msgJson["next_y"] = refWaypointsY;

          // Create car state with latency and relative errors
          double latency = tracker.getAvgElapsedMs();
          double yawRate = v * -steer / params.Lf;
          double dX = Tools::isZero(yawRate) ? v * latency : v / yawRate * sin(yawRate * latency);
          double dY = Tools::isZero(yawRate) ? 0.0 : v / yawRate * (-cos(yawRate * latency) + 1.0);
          double dPsi = yawRate * latency;
          double dV = throttle * latency;
          double cte = polyEval(refCoeffs, dX) - dY;
          double epsi = Tools::calculateAngleDelta(atan(refCoeffs[1] + 2 * refCoeffs[2] * dX), dPsi);
          Eigen::VectorXd state(6);
          state << dX
                 , dY
                 , dPsi
                 , v + dV
                 , cte
                 , epsi;

          // Calculate actuations using mpc
          std::vector<Eigen::VectorXd> actuations;
          if (mpc.Solve(state, refCoeffs, actuations))
          {
            Eigen::VectorXd actuation = actuations.front();
            msgJson["steering_angle"] = -Tools::rad2deg(actuation(6)) / 25.0 / params.Lf;
            msgJson["throttle"] = actuation(7);

            //Display the MPC predicted trajectory
            std::vector<double> mpcWaypointsX;
            std::vector<double> mpcWaypointsY;
            for (int t = 0; t < actuation.size(); ++t)
            {
              mpcWaypointsX.push_back(actuations[t](0));
              mpcWaypointsY.push_back(actuations[t](1));
              msgJson["mpc_x"] = mpcWaypointsX;
              msgJson["mpc_y"] = mpcWaypointsY;
            }
          }
          else
          {
            msgJson["steering_angle"] = 0.0;
            msgJson["throttle"] = 0.1;
            msgJson["mpc_x"] = refWaypointsX;
            msgJson["mpc_y"] = refWaypointsY;
          }


          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          std::this_thread::sleep_for(std::chrono::milliseconds(100));
          tracker.end();

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      }
      else
      {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t)
  {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1)
    {
      res->end(s.data(), s.length());
    }
    else
    {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req)
  {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length)
  {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port))
  {
    std::cout << "Listening to port " << port << std::endl;
  }
  else
  {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
