
# Reflection

## Model Predictive Control Project

The goals / steps of this project are the following:

* Your code must compile without errors with cmake and make.
* Student describes their model in detail. This includes the state, actuators and update equations.
* Student discusses the reasoning behind the chosen N (timestep length) and dt (elapsed duration between timesteps) values. Additionally the student details the previous values tried.
* A polynomial is fitted to waypoints.
* The student implements Model Predictive Control that handles a 100 millisecond latency. Student provides details on how they deal with latency.

[//]: # (Image References)
[image1]: ./img/error.png
[image2]: ./img/error2.png
[image3]: ./img/mcp.png

## Rubric Points
Here I will consider the reflection points of the [rubric](https://review.udacity.com/#!/rubrics/896/view) individually and describe how I addressed each point in my implementation.

---

### 1. Describe the model in detail. This includes the state, actuators and update equations.

The model used for prediction is a kinematic model, which nicely corresponds to the non-dynamic model implemented within the used udacity simulator.

#### Model state

The state of the model is described by Vector of dimension 6, just as taught in the lessons:

| Dim. | Value | Description  |
|------|-------|--------------|
| 0    |   x   | x pos in vehicle coordinate system based on initial state |
| 1    |   y   | y pos in vehicle coordinate system based on initial state  |
| 2    |   ψ   | Angle of vehicle in vehicle coordinate system based on initial state |
| 3    |   v   | Velocity of vehicle in simulator speed unit |
| 4    |  cte  | Cross track error (distance to reference path) |
| 5    |  eψ   | Vehicle angle error (angle difference to reference path |

The initial MPC state is directly derived from current vehicle state. Thus, in contrast to the PID controller, the prediction of the MPC is not depending on any previously made prediction.

#### Model actuators

The model is able to control two actuators in order to modify the predicted trajectory of the vehicle:

| Dim. |  Value   | Description  |
|------|----------|--------------|
| 0    |  steer   | Steering from -1.0 to 1.0, influencing the vehicles raw rate |
| 1    | throttle | Acceleration from -1.0 to 1.0 as factor of current speed |

Actuator values are at no point passed to the MCP, these are purely derived using
* The initial state
* The update equations
* The cost function which quantifies the difference to the reference trajectory

#### Model update equations

I have been reusing the update equations from the *Scentless Kalman Filter* project here.
The advantage was, that it integrates the yaw rate over time which more closely matches the simulators implementation. Here's the codes in C++:

```c++
  if (Tools::isZero(yawRate))
  {
    // Straight case
    x1 = x0 + v0 * CppAD::cos(psi0) * mP.deltaT;
    y1 = y0 + v0 * CppAD::sin(psi0) * mP.deltaT;
  }
  else
  {
    // Cyclic case
    x1 = x0 + v0 / yawRate * ( CppAD::sin(psi0 + yawRate * mP.deltaT) - CppAD::sin(psi0));
    y1 = y0 + v0 / yawRate * (-CppAD::cos(psi0 + yawRate * mP.deltaT) + CppAD::cos(psi0));
  }
```

The formula for the `yawRate` in this case was given by the lesson as:
```c++
AD<double> yawRate = v0 * steer / mP.Lf;
```
The constant factor `Lf` hereby has been derived from testing with a constant steer value against against the actual (simulated) vehicle and comparing turn radius. It has to depend on speed, because the final angle turned depends on the distance traveled during delta time.

The deltas for angle and speed update are straight-forward calculated by just multiplying their associated actuator values by the time passed:
```c++
psi1 = psi0 + yawRate * mP.deltaT;
v1 = v0 + throttle * mP.deltaT;
```

As we calculate the new cross track error, we are performing two approximations:
 1) We evaluate the polynomial function which describes the reference trajectory at `x0` in order to re-calculate `cte0`
 2) We evaluate the distance that was moved orthogonally to the angle of the reference polynomial at `x0` in order to approximate the delta to `cte0` in next state

The new error between vehicle psi and polynomial psi, is calculated similar:
 1) We re-calculate `psi0` by evaluating the angle of the reference polynomial at `x0`
 2) We simply add the already known delta between `psi1` and `psi0` to it.

The formulas in C++ are:
```c++
cte1 = (yPoly(x0) - y0) + (v0 * CppAD::sin(epsi0) * mP.deltaT);
epsi1 = (psi0 - psiPoly(x0) + yawRate * mP.deltaT;
```

To be suitable for `CppAd` we further transform all above formulas, so that the left side of the equation has to evaluate to 0.


### 2. Discuss the reasoning behind the chosen N and dt values

I did settle with a constant `N` of 10 and a variable `dt` based on target speed of controller:
```c++
deltaT = 8.0 / (targetSpeed + 4.0)
```

While I have chosen N my manually evaluating the controllers driving quality on simulator, I evaluated the effect of `dt` using an extended version of the `mpc_to_line` quizzle.

The plot visualizations of the quizzle clear showed, how the amplitude of state delta is dependant on the product of vehicle velocity and delta time between two successive states. The effect of large state delta is, that the solver fails to converge the trajectory with the reference polynomial and instead rushes the car into a loop trajectory. So as a counter-measure, I'd decrease the `dt` with growing target speed.

Another important effect is caused by the product of `N` and `dt`. If this "future horizon" of the solver is chosen too small, the controller reacts to curves and corrections too late causing the vehicle to slide on an unsteady wave along the reference trajectory.

Another effect I monitored is that if 'N' is chosen larger and larger, it contributes more and more to the latency, thus amplifying all errors between motion model and actual (simulated) vehicle.

Also I found that the tendency for cyclic trajectories can further be reduced by tweaking the used cost model.
Here is an example erratic trajectory (enforced by high `dt` compared to `v`):

![Error case plot][image1]

The first important step is to decrease the impact of the velocity reference cost, as we need the solver to be able to vary driving speed as needed. I did so by simply applying a constant factor of 0.01 to the reference speed cost and 0.1 to the absolute throttle  cost. I did not touch the delta throttle cost, as I still do want the car to de/accelerate consistently.

Furthermore I wanted to encourage the optimizer to slow the vehicle down when the distance to the trajectory grows too high. I tried to achieve this by multiplying each cte value in the trajectory with the corresponding velocity value.

Last but not least, I want the optimizer to approach the reference path more smoothly and avoid ending in a loop. As a solution, I am  multiplying the psi error cost with the according time spent inside the trajectory.

Changes to the default cost formulas:
```c++
    for (int t = 0; t < mP.N; ++t)
    {
      fg[0] += 1.0  * CppAD::abs(vars[mP.v]) * CppAD::pow(vars[mP.cte + t], 2);
      fg[0] += 1.0  * CppAD::pow(t * vars[mP.epsi + t], 2);
      fg[0] += 0.01 * CppAD::pow(vars[mP.v + t] - mP.refV, 2);
    }

    // Minimize the use of actuators.
    for (int t = 0; t < mP.N - 1; ++t)
    {
      fg[0] += 1.0 * CppAD::pow(vars[mP.steer + t], 2);
      fg[0] += 0.1 * CppAD::pow(vars[mP.throttle + t], 2);
    }
```

Resulting trajectories for same case as above:

![Error case plot with adjusted cost model][image2]

We can see that steering and throttle are much more smooth now. If we would have chosen a smaller `dt` for the given velocity now, the vehicle would have smoothly approached the reference path and consistently followed it.

### 3. Provide details on how you've dealed with latency.

I haven taken the hint from the lessons and applied the update equasions once on the initial state before feeding it into the MCP. In this step, the equasion can be simplified quite alot as many factory e.g. cos/sin are either zero or one. As `Dt` I have chosen the latency as projected from the mean of the last 5 algorithm invocations.

With this setup, the car drives safely and smoothly along the course with a speed between 70 on straight sections and 40 inside curves. The algorithm has potential for higher speeds, but I played it safe.
The following video shows one lap:


[Vehicle driving one save lap (mcp.ogv)](./img/mcp.ogv)

![Screenshot of video][image3]