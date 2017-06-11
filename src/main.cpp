#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <cmath>
#include <iomanip>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Fit a polynomial. Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

template <typename T>
void dump_on_cerr(int iter, string prefix, const vector<T>& v) {
  cerr << iter << "," << prefix;
  for (auto& e : v) cerr << "," << e;
  cerr << endl;
}

void dump_on_cerr(int iter, string prefix, const Eigen::VectorXd& v) {
  cerr << iter << "," << prefix;
  for (size_t i = 0; i < v.size(); ++i) cerr << "," << v(i);
  cerr << endl;
}

Params read_args(int argc, char* argv[], int& N, int& latency_ms) {
  static const string usage = R"USAGE(
    Usage: mpc [N] [dt] [ref_v] [sleep_ms] [delay] [MPC_pts_offsets] [w_v w_cte w_epsi w_delta w_a w_delta2 w_a2 w_exp]
    Eg:
      # no arguments, uses built-in defaults
      mpc

      # N=10, dt=0.1, ref_v=50, latency_ms=0, delay=0.3
      mpc 10 0.1 50 0 0.3
)USAGE";

  // Set the (default) timestep length and duration, ref velocity
  // All of these can be overridden from command line
  N = 8; latency_ms = 100;
  double dt = 0.1, ref_v = 50, delay = 0.2; 
  int MPC_pts_offsets = 1;
  
  try {
    auto i = 1;
    if (argc > i) N = stoi(argv[i++]);
    if (argc > i) dt = stod(argv[i++]);
    if (argc > i) ref_v = stod(argv[i++]);
    if (argc > i) latency_ms = stoi(argv[i++]);
    if (argc > i) delay = stod(argv[i++]);
    if (argc > i) MPC_pts_offsets = stoi(argv[i++]);
    cerr << "-1,params,N=," << N << ",dt=," << dt << ",ref_v=," << ref_v << ",latency_ms=" << latency_ms << ",delay=," << delay 
         << ",MPC_pts_offsets=," << MPC_pts_offsets << endl;

    auto p = Params { MPC_pts_offsets, dt, ref_v, delay };
    if (argc > i) p.w_v = stod(argv[i++]);
    if (argc > i) p.w_cte = stod(argv[i++]);
    if (argc > i) p.w_epsi = stod(argv[i++]);
    if (argc > i) p.w_delta = stod(argv[i++]);
    if (argc > i) p.w_a = stod(argv[i++]);
    if (argc > i) p.w_delta2 = stod(argv[i++]);
    if (argc > i) p.w_a2 = stod(argv[i++]);
    if (argc > i) p.w_exp = stod(argv[i++]);
    cerr << "-1,weights,w_v=," << p.w_v << ",w_cte=," << p.w_cte << ",w_epsi=," << p.w_epsi << ",w_delta=" << p.w_delta 
         << ",w_a=," << p.w_a << ",w_delta2=," << p.w_delta2 << ",w_a2=," << p.w_a2 << ",w_exp=" << p.w_exp 
         << ",projected=," << p.projected << endl;

    return p;
  } catch (...) {
    cerr << usage << endl;
    exit(1);
  }
}

inline void to_car_coordinates(const double psi, const double px, const double py, double& x, double& y) {
  double tx = (x - px),     ty = (y - py);
  //double c = cos(psi),      s = sin(psi);
  //x = tx*c + ty*s;          y = ty*c - tx*s;
  double c = cos(0 - psi),  s = sin(0 - psi);  // rotate such that psi=0
  x = tx*c - ty*s;          y = tx*s + ty*c;
}

inline void to_car_coordinates(const double psi, const double px, const double py, vector<double>& X, vector<double>& Y) {
  for (size_t i = 0; i < X.size(); ++i) {
    to_car_coordinates(psi, px, py, X[i], Y[i]);
  }
}

int main(int argc, char* argv[]) {
  uWS::Hub h;
  int N, latency_ms;
  Params p = read_args(argc, argv, N, latency_ms);

  // MPC is initialized here!
  MPC mpc(N, p);
  cerr.setf(ios::fixed, ios::floatfield); cerr.precision(6);
  auto iter = 0ul;

  h.onMessage([&mpc, &iter, &latency_ms](uWS::WebSocket<uWS::SERVER>* ws, char *data, size_t length, uWS::OpCode opCode) {
    static const auto STEERING_NORMALISER = deg2rad(25); //*Lf;
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double steering_angle = j[1]["steering_angle"];
          double throttle = j[1]["throttle"];
          cerr << iter << ",px=," << px << ",py=," << py << ",psi=," << psi << ",speed=," << v  << ",angle=," << steering_angle << ",throttle=," << throttle << endl;
          dump_on_cerr(iter, "ptsx", ptsx);
          dump_on_cerr(iter, "ptsy", ptsy);
          
          /*
          * Calculate steeering angle and throttle using MPC.
          * Both are in between [-1, 1].
          */

          //Display the waypoints/reference line
          to_car_coordinates(psi, px, py, ptsx, ptsy);
          dump_on_cerr(iter, "NEW ptsx", ptsx);
          dump_on_cerr(iter, "NEW ptsy", ptsy);
          //auto coeffs = polyfit(xvals, yvals, 3);
          auto coeffs = polyfit(Eigen::Map<Eigen::VectorXd>(&ptsx[0], 6), Eigen::Map<Eigen::VectorXd>(&ptsy[0], 6), 3);
          if (std::isnan(coeffs(0))) return;
          dump_on_cerr(iter, "coeffs", coeffs);

          // kinematic model predict after latency_ms
          auto v_dt = (latency_ms / 1000.0) * (v * 0.447038889); // convert from mph to m/s (mph * 1.60934 * 1000) / (60 * 60);
          steering_angle *= STEERING_NORMALISER;
          px = v_dt * cos(steering_angle);
          py = v_dt * sin(steering_angle);
          psi = 0; //-v_ * steering_angle / Lf;
          v += throttle * latency_ms / 1000.0;
          cerr << iter << ",NEW px=," << px << ",py=," << py << ",psi=," << psi << ",speed=," << v  << ",angle=," << steering_angle << ",throttle=," << throttle << endl;

          auto cte = polyeval(coeffs, px) - py;
          auto deriv = polyevalderiv(coeffs, px);
          auto epsi = -atan(deriv); //psi - atan(deriv);

          Eigen::VectorXd state(6);
          //state << 0.0 /*px*/, 0.0 /*py*/, 0.0 /*psi*/, v, cte, epsi;
          state << px, py, psi, v, cte, epsi;
          //state << px, 0.0, 0.0, v, cte, epsi;
          //state << v_ /*px*/, 0.0 /*py*/, -v_ * steering_angle / 2.67 /*psi*/, v + throttle*latency_ms, cte, epsi; // predicted after latency_ms
          
          dump_on_cerr(iter, "state", state);
          auto r = mpc.Solve(state, coeffs, steering_angle, throttle);
          dump_on_cerr(iter, "r", r);
          double steer_value = r[0];
          double throttle_value = r[1]; //0.3;

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value / STEERING_NORMALISER;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory 
          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line
          msgJson["mpc_x"] = mpc.mpc_x_vals;
          msgJson["mpc_y"] = mpc.mpc_y_vals;
          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line
          msgJson["next_x"] = ptsx;
          msgJson["next_y"] = ptsy;

          dump_on_cerr(iter, "ptsx", ptsx);
          dump_on_cerr(iter, "ptsy", ptsy);
          dump_on_cerr(iter, "mpc_x", mpc.mpc_x_vals);
          dump_on_cerr(iter, "mpc_y", mpc.mpc_y_vals);
          iter++;

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(latency_ms)); // project Rubric requires this to be 100, which is the default value
          ws->send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws->send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER>* ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER>* ws, int code,
                         char *message, size_t length) {
    ws->close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen("0.0.0.0", port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
