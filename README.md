# CarND-Controls-MPC
Self-Driving Car Engineer Nanodegree Program

---
Please click on this image to see video of PID controller driving the car for over a lap.
[![MPC](https://img.youtube.com/vi/e74l1EP15mI/0.jpg)](https://youtu.be/e74l1EP15mI)

## Reflections
### Build And Run Issues
* The simulator was excruciatingly slow on Ubuntu running on a VM hosted on Windows.
* Build on Windows with Visual Studio was becoming very complex, esp. the Fortran dependencies for BLAS. Tried [OpenBLAS](http://www.openblas.net/) and [NetLib](http://www.netlib.org/lapack/index.html)'s [LAPACK for windows](http://icl.cs.utk.edu/lapack-for-windows/lapack/)
* Eventually used cygwin toolchain and gcc, gfortran to compile on Windows.
* Same uWS version issues encountered as was the case with PID project. Windows build required changes and uWS version available is v0.14.2 which is source incompatible with v0.13. **NB Reviewers: please use v0.14+ of uWS**
* I acknowledge the information provided in help forums and slack channels.

### The Model
The model of the vehicle comprises of following:
  * State: 
    * Vehicle's location coordinates (`x`,`y`), orientation (`psi`), and velocity (`v`)
    * Two further augmented parameters to model error: cross trace error (`cte`) and orientation (`epsi`)
  * Actuators:  There are two actuators, both of are constrained in range [-1,1]
    * Steering angle (`delta`)
    * Throttle/acceleration (`a`), negative values imply breaking

The inputs from simulator have the vehicle’s location and orientation and the path to follow in global coordinates.  To make the processing simpler, these are converted to car's coordinates.

The state update equations follow those from the lesson, assuming elapsed duration between timesteps is `dt`
  ```
    x[t+1]   = x[t]   + v * cos(psi) * dt
    y[t+1]   = y[t]   + v * sin(psi) * dt
    psi[t+1] = psi[t] + v/Lf * delta * dt
    a[t+1]   = v      + a[t]         * dt
  ```
Where `Lf` is the distance between centre of front wheels and car's centre of gravity.  This parameter's value was supplied for the model.

I have further constrained the model's `delta` actuator to be dynamically adjusted based on speed of the car.  For higher speeds, the `delta` range is reduced, as small steering angle changes still result in large turns for the car at high speeds.

### Timestep Length (`N`) & Elapsed Duration (`dt`)
To help experiment with these (and other parameter's) values, I allowed them to be (optionally) specified on the command line.  After a several experiments, I have selected default values for these parameters such that it runs out of box (at least on my PC).

Since the processing in the model was quite performant, I had initial thinking that higher values `N` and smaller values of `dt` will give me very good and stable results.  However, as was pointed in lesson, I experienced that the further out I projected (`N*dt`), the higher the error was.  On the other hand, very low values of `N*dt` also didn't provide sufficiently accurate tracking of the path.  This observation also implied that the optimal values of `N` and `dt` will differ for different target velocities.  At higher velocity, same `N` and `dt` will project out much further and be unstable.

After a few experiments, for default target velocity of 50mph, I chose `N=8` and `dt=0.1`.

### Polynomial Fitting and MPC Preprocessing
The waypoint path from simulator is in global coordinates.  It is converted to car's coordinates before passing into the model.  This implies, the vehicles position becomes (0,0), orientation becomes 0.  Rest of the processing within model is then done in car's coordinates.  It could have been done in global coordinates too, but many other transforms would be required, whilst thinking from car's point of view makes the processing slightly simpler.

Another point to note is that for the simulator steering values are [-25,25] degrees *normalised* to range of [-1,1], however the model works in radians.  Hence both while feeding the model and using its outputs, the angles are (un)normalised and converted to radians.

### Model Tuning
I added several hyper-parameters to tune the weights of various parameters of the model.  During my experiments, I found stabilising `delta` had a very significant effect on car's stability, so its weight is orders of magnitude higher than even those of `cte`.  Velocity has a very low weight, such that is always penalised in favour of stability and tracking the path. Just like `N` and `dt`, these too need experimentation and optimal values differ for different target speeds.

### Model Predictive Control with Latency
Actuator latency is handled by determining where the vehicle will be after latency before feeding to the model.  Since the processing is from car's point of view and `psi` is already converted to 0, we can predict (x,y) coordinates of where car will end-up after latency (100ms) using steering angle and velocity.  I had to un-normalise the steering angle into radians and also convert the speed from mph to meters-per-second.  I pass-in the current values of steering angle and throttle (which I use as proxy for acceleration) to the model and set the initial few steps (until the latency time is passed) of model constraints to those current values.  I also look further out beyond latency time steps in output array and use those values for actuator controls.

## Command Line Usage
  ```
    Usage: mpc [N] [dt] [ref_v] [sleep_ms] [delay] [MPC_pts_offsets] [w_v w_cte w_epsi w_delta w_a w_delta2 w_a2 w_exp]
    Eg:
      # no arguments, uses built-in defaults
      mpc

      # N=10, dt=0.1, ref_v=50, latency_ms=0, delay=0.3
      mpc 10 0.1 50 0 0.3
  ```


---
## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets 
    cd uWebSockets
    git checkout e94b6e1
    ```
    Some function signatures have changed in v0.14.x. See [this PR](https://github.com/udacity/CarND-MPC-Project/pull/3) for more details.
* Fortran Compiler
  * Mac: `brew install gcc` (might not be required)
  * Linux: `sudo apt-get install gfortran`. Additionall you have also have to install gcc and g++, `sudo apt-get install gcc g++`. Look in [this Dockerfile](https://github.com/udacity/CarND-MPC-Quizzes/blob/master/Dockerfile) for more info.
* [Ipopt](https://projects.coin-or.org/Ipopt)
  * Mac: `brew install ipopt`
  * Linux
    * You will need a version of Ipopt 3.12.1 or higher. The version available through `apt-get` is 3.11.x. If you can get that version to work great but if not there's a script `install_ipopt.sh` that will install Ipopt. You just need to download the source from the Ipopt [releases page](https://www.coin-or.org/download/source/Ipopt/) or the [Github releases](https://github.com/coin-or/Ipopt/releases) page.
    * Then call `install_ipopt.sh` with the source directory as the first argument, ex: `bash install_ipopt.sh Ipopt-3.12.1`. 
  * Windows: TODO. If you can use the Linux subsystem and follow the Linux instructions.
* [CppAD](https://www.coin-or.org/CppAD/)
  * Mac: `brew install cppad`
  * Linux `sudo apt-get install cppad` or equivalent.
  * Windows: TODO. If you can use the Linux subsystem and follow the Linux instructions.
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). This is already part of the repo so you shouldn't have to worry about it.
* Simulator. You can download these from the [releases tab](https://github.com/udacity/self-driving-car-sim/releases).
* Not a dependency but read the [DATA.md](./DATA.md) for a description of the data sent back from the simulator.


## Basic Build Instructions


1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./mpc`.

## Tips

1. It's recommended to test the MPC on basic examples to see if your implementation behaves as desired. One possible example
is the vehicle starting offset of a straight line (reference). If the MPC implementation is correct, after some number of timesteps
(not too many) it should find and track the reference line.
2. The `lake_track_waypoints.csv` file has the waypoints of the lake track. You could use this to fit polynomials and points and see of how well your model tracks curve. NOTE: This file might be not completely in sync with the simulator so your solution should NOT depend on it.
3. For visualization this C++ [matplotlib wrapper](https://github.com/lava/matplotlib-cpp) could be helpful.

## Editor Settings

We've purposefully kept editor configuration files out of this repo in order to
keep it as simple and environment agnostic as possible. However, we recommend
using the following settings:

* indent using spaces
* set tab width to 2 spaces (keeps the matrices in source code aligned)

## Code Style

Please (do your best to) stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html).

## Project Instructions and Rubric

Note: regardless of the changes you make, your project must be buildable using
cmake and make!

More information is only accessible by people who are already enrolled in Term 2
of CarND. If you are enrolled, see [the project page](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/f1820894-8322-4bb3-81aa-b26b3c6dcbaf/lessons/b1ff3be0-c904-438e-aad3-2b5379f0e0c3/concepts/1a2255a0-e23c-44cf-8d41-39b8a3c8264a)
for instructions and the project rubric.

## Hints!

* You don't have to follow this directory structure, but if you do, your work
  will span all of the .cpp files here. Keep an eye out for TODOs.

## Call for IDE Profiles Pull Requests

Help your fellow students!

We decided to create Makefiles with cmake to keep this project as platform
agnostic as possible. Similarly, we omitted IDE profiles in order to we ensure
that students don't feel pressured to use one IDE or another.

However! I'd love to help people get up and running with their IDEs of choice.
If you've created a profile for an IDE that you think other students would
appreciate, we'd love to have you add the requisite profile files and
instructions to ide_profiles/. For example if you wanted to add a VS Code
profile, you'd add:

* /ide_profiles/vscode/.vscode
* /ide_profiles/vscode/README.md

The README should explain what the profile does, how to take advantage of it,
and how to install it.

Frankly, I've never been involved in a project with multiple IDE profiles
before. I believe the best way to handle this would be to keep them out of the
repo root to avoid clutter. My expectation is that most profiles will include
instructions to copy files to a new location to get picked up by the IDE, but
that's just a guess.

One last note here: regardless of the IDE used, every submitted project must
still be compilable with cmake and make./
