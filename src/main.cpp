#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"

#include "Helper.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  Vehicle ego_car = Vehicle();

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&ego_car](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

			/*Intepolate waypoint to the close area of ego_car*/
			int num_waypoints = map_waypoints_s.size();
			int next_waypoint = NextWaypoint(car_x, car_y, car_yaw, map_waypoints_x, map_waypoints_y);

			vector<double> close_waypoint_x, close_waypoint_y, close_waypoint_s, close_waypoint_dx, close_waypoint_dy;
			for(int i = -NUM_WAYPOINT_BEFORE; i<NUM_WAYPOINT_AFTER; i++)
			{
				int idx = next_waypoint + i;
				if(idx < 0)
				{
					idx += num_waypoints;

				}
				/*to detect the case vehicle is clsoe to the end of the track, and the close_waypoints has cover the end of the track, need to make it continus*/
				s_point_idx = map_waypoints_s[idx];
				s_next_waypoint = map_waypoints_s[map_waypoints_s];
				if(i < 0 && s_point_idx > s_next_waypoint)
				{
					s_point_idx -= TRACK_LENGTH; 
				}
				if(i>0 && s_point_idx<s_next_waypoint)
				{
					s_point_idx += TRACK_LENGTH;
				}
				
				close_waypoint_x.push_back(map_waypoints_x[idx]);
				close_waypoint_y.push_back(map_waypoints_y[idx]);
				close_waypoint_s.push_back(map_waypoints_s[idx]);
				close_waypoint_dx.push_back(map_waypoints_dx[idx]);
				close_waypoint_dy.push_back(map_waypoints_dy[idx]);
				
			}

			double dist_s_inc = 0.5;
			int num_intepolate_points = (close_waypoint_s[close_waypoint_s.size() - 1] - close_waypoint_s[0])/dist_s_inc;
			vector<double> intepolate_waypoint_x, intepolate_waypoint_y, intepolate_waypoint_s, intepolate_waypoint_dx, intepolate_waypoint_dy;
			for(int i =0; i<num_intepolate_points;i++ )
			{
				intepolate_waypoint_s.push_back(close_waypoint_s[0] + i * dist_s_inc);
			}
			intepolate_waypoint_x = interpolate_points(close_waypoint_s, close_waypoint_x, dist_s_inc, num_intepolate_points);
			intepolate_waypoint_y = interpolate_points(close_waypoint_s, close_waypoint_y, dist_s_inc, num_intepolate_points);
			intepolate_waypoint_dx = interpolate_points(close_waypoint_s, close_waypoint_dx, dist_s_inc, num_intepolate_points);
			intepolate_waypoint_dy = interpolate_points(close_waypoint_s, close_waypoint_dy, dist_s_inc, num_intepolate_points);


			/* calculate ego car parameter*/
			//ego car has following parameters(s,s_d,s_dd,d,d_d,d_dd)
			double pos_s, pos_s_d, pos_s_dd, pos_d, pos_d_d, pos_d_dd,pos_s2,pos_d2,pos_s3,pos_d3,pos_s_d2, pos_d_d2;
			double pos_x, pos_y, angle, pos_x2,pos_y2, pos_x3,pos_y3, pos_x4, pos_y4;

			int pre_closepath_size = min(NUM_PREVIOUSPATH_KEEP, int previous_path_x.size());

			/* calculate ego parameter based on the current vehicle signals if the previous path points is not enough*/
			if(pre_closepath_size < 5)
			{
				pos_x = car_x;
				pos_y = car_y;
				angle = deg2rad(car_yaw);
				pos_s = car_s;
				pos_s_d = car_speed;
				pos_s_dd = 0;
				pos_d = car_d;
				pos_d_d = 0;
				pos_d_dd = 0;
			}
			else
			{
				pos_x = previous_path_x[pre_closepath_size -1];
				pos_y = previous_path_y[pre_closepath_size -1];
				pos_x2 = previous_path_x[pre_closepath_size -2];
				pos_y2 = previous_path_y[pre_closepath_size -2];
				pos_x3 = previous_path_x[pre_closepath_size -3];
				pos_y3 = previous_path_y[pre_closepath_size -3];
				pos_x4 = previous_path_x[pre_closepath_size -4];
				pos_y4 = previous_path_y[pre_closepath_size -4];
				angle = atan2(pos_y-pos_y2, pos_x-pos_x2);
				angle2 = atan2(pos_y2-pos_y3, pos_x2-pos_x3);
				angle3 = atan2(pos_y3-pos_y4, pos_x3-pos_x4);

				vector<double> Frenet_s_d = getFrenet(pos_x, pos_y, angle, intepolate_waypoint_x, intepolate_waypoint_y);
				vector<double> Frenet_s_d2 = getFrenet(pos_x2, pos_y2, angle2, intepolate_waypoint_x, intepolate_waypoint_y);
				vector<double> Frenet_s_d3 = getFrenet(pos_x3, pos_y3, angle3, intepolate_waypoint_x, intepolate_waypoint_y);
				
				pos_s = Frenet_s_d[0];
				pos_d = Frenet_s_d[1];
				s2 = Frenet_s_d[0];
				d2 = Frenet_s_d[1];
				s3 = Frenet_s_d[0];
				d3 = Frenet_s_d[1];

				pos_s_d = (pos_s- pos_s2)/PATH_DT;
				pos_d_d = (pos_d- pos_d2)/PATH_DT;
				pos_s_d2 = (pos_s2- pos_s3)/PATH_DT;
				pos_d_d2 = (pos_d2- pos_d3)/PATH_DT;

				pos_s_dd = (pos_s_d - pos_s_d2)/PATH_DT;
				pos_d_dd = (pos_d_d - pos_d_d2)/PATH_DT;
			}

				ego_car.s = pos_s;
				ego_car.s_d = pos_s_d;
				ego_car.s_dd = pos_s_dd;
				ego_car.d = pos_d;
				ego_car.d_d = pos_s_d;
				ego_car.d_dd = pos_d_dd;

			/* Predictions */
			double prediction_time = PREDICTION_SAMPLES * PREDICTION_DT;

			vector<Vehicle> other_cars;
			map<int, vector<vector<double>>> predictions;
			//data format for sensor_fusion: [ID, x, y, dx, dy, s, d]
			for(auto sf:sensor_fusion)
			{
				double othercar_vel = sqrt(pow((double)sf[3],2)+pow((double)sf[4],2));
				Vehicle othercar = Vehicle(sf[5], othercar_vel, 0 , sf[6],0,0)
				other_cars.push_back(othercar);
				vector<vector<double>> prediction;
				for(i=0;i<PREDICTION_SAMPLES;i++)
				{
					double predict_s = this->s + this->s_d *PREDICTION_DT;
					double predict_d = this ->d;
					vector<double> s_and_d = {predict_s, predict_d};
					prediction.push_back(s_and_d);
				}
				predictions[sf[0]] = prediction;
				
			}

			/* calculate best trajectory */
			//calculate whether there is car in ego lane, left lane and right lane
			bool car_at_left = false, car_at_right = false, car_in_ahead = false;
			for(Vehicle othercar: other_cars)
			{
				if(fabs(car_s - other_car.s) < SAFETY_DISTANCE)
				{
					double d_diff = car_d - othercar.d; 
					if(2 < d_diff < 6)
					{
						car_at_left = true;
					}
					else if(-6 < d_diff < -2)
					{
						car_at_right = true;
					}
					else if( -2<d_diff<2)
					{
						car_in_ahead = true;
					}
				}
			}

			ego_car.possible_state = {"KL"};
			if(ego_car.d > 4 && !car_at_left)
			{
				ego_car.possible_state.push_back("LCL");
			}
			if(ego_car.d < 8 && !car_at_right)
			{
				ego_car.possible_state.push_back("LCR");
			}
			
			// calculate the vehicle target state and the end of the prediction time according to each state 
			double target_s, target_s_d, target_s_dd, target_d, target_d_d, target_d_dd;
			target_s_d = SPEED_LIMIT;
			target_s_dd = 0;
			target_d_d = 0;
			target_d_dd = 0;
			target_s = ego_car.s + (ego_car.s_d + target_s_d)/2*prediction_time;
			int target_lane, current_lane = ego_car.d/4;

			vector<vector<double>> best_frenet_traj, best_target;
			double best_cost = 999999;
			string best_traj_state = "";
			
			for(string state: ego_car.possible_state)
			{
				if(state.compare("KL") == 0)
				{
					target_d = current_lane*4 + 2;
				}
				else if(state.compare("LCL") == 0)
				{
					target_d = (current_lane -1)*4 + 2;
				}
				else if(state.compare("LCR") == 0)
				{
					target_d = (current_lane +1)*4 + 2;
				}
				
			
				target_lane = target_d/4;

				// calculate predicted leading vehicle state
				for(auto prediction: predictions)
				{
					vector<vector<double>> pre_traj = prediction.second;
					int prediction_lane = pre_traj[0][1]/4;
					if(prediction_lane == target_lane)
					{
						prediction_start_s = pre_traj[0][0];
						prediction_end_s = pre_traj[pre_traj.size() -1][0];
						prediction_next2end_s = pre_traj[pre_traj.size() -2][0];
						prediction_end_s_d = (prediction_end_s - prediction_next2end_s)/PREDICTION_DT;
						if(target_s > prediction_end_s && prediction_start_s > ego_car.s)
						{
							target_s = prediction_end_s;
							target_s_d = prediction_end_s_d;
						}
						
					}
				}
				vector<double> target_state_s = {target_s, target_s_d, target_s_dd};
				vector<double> target_state_d = {target_d, target_d_d, target_d_dd};
				vector<double> current_state_s = {ego_car.s, ego_car.s_d, ego_car.s_dd};
				vector<double> current_state_d = {ego_car.d, ego_car.d_d, ego_car.d_dd};
				
				

				// caculate the possible target trajectory
				ego_car.s_trajectory_coeffs = get_traj_coeffs(current_state_s, target_state_s, prediction_time);
				ego_car.d_trajectory_coeffs = get_traj_coeffs(current_state_d, target_state_d, prediction_time);

				  // populate s and t trajectories at each time step
				vector<double> s_traj;
  				vector<double> d_traj;
				for (int i = 0; i < PREDICTION_SAMPLES; i++) {
				    double t = i * PREDICTION_DT;
				    double s_val = 0, d_val = 0;
				    for (int j = 0; j < ego_car.s_traj_coeffs.size(); j++) {
				      s_val += ego_car.s_traj_coeffs[j] * pow(t, j);
				      d_val += ego_car.d_traj_coeffs[j] * pow(t, j);
				    }
				    s_traj.push_back(s_val);
				    d_traj.push_back(d_val);
				  }

				double current_state_cost = calculate_total_cost(s_traj, d_traj, predictions);



				if (current_state_cost < best_cost) {
					best_cost = current_state_cost;
					best_frenet_traj = {s_traj, d_traj};
					best_traj_state = state;
					best_target = {target_state_s, target_state_d};
				}
	
				
			}

			// ********************* PRODUCE NEW PATH ***********************
			// begin by pushing the last and next-to-last point from the previous path for setting the 
			// spline the last point should be the first point in the returned trajectory, but because of 
			// imprecision, also add that point manually

			vector<double> coarse_s_traj, coarse_x_traj, coarse_y_traj, interpolated_s_traj, 
												 interpolated_x_traj, interpolated_y_traj;

			double prev_s = pos_s - pos_s_d * PATH_DT;
					
			// first two points of coarse trajectory, to ensure spline begins smoothly
			if (pre_closepath_size >= 2) {
				coarse_s_traj.push_back(prev_s);
				coarse_x_traj.push_back(previous_path_x[subpath_size-2]);
				coarse_y_traj.push_back(previous_path_y[subpath_size-2]);
				coarse_s_traj.push_back(pos_s);
				coarse_x_traj.push_back(previous_path_x[subpath_size-1]);
				coarse_y_traj.push_back(previous_path_y[subpath_size-1]);
				} else {
						double prev_s = pos_s - 1;
						double prev_x = pos_x - cos(angle);
						double prev_y = pos_y - sin(angle);
						coarse_s_traj.push_back(prev_s);
						coarse_x_traj.push_back(prev_x);
						coarse_y_traj.push_back(prev_y);
						coarse_s_traj.push_back(pos_s);
						coarse_x_traj.push_back(pos_x);
						coarse_y_traj.push_back(pos_y);
				}

			// last two points of coarse trajectory, use target_d and current s + 30,60
			double target_s1 = pos_s + 30;
			double target_d1 = best_target[1][0];
			vector<double> target_xy1 = getXY(target_s1, target_d1, intepolate_waypoint_s, intepolate_waypoint_x, intepolate_waypoint_y);
			double target_x1 = target_xy1[0];
			double target_y1 = target_xy1[1];
			coarse_s_traj.push_back(target_s1);
			coarse_x_traj.push_back(target_x1);
			coarse_y_traj.push_back(target_y1);
			double target_s2 = target_s1 + 30;
			double target_d2 = target_d1;
			vector<double> target_xy2 = getXY(target_s2, target_d2, intepolate_waypoint_s, intepolate_waypoint_x, intepolate_waypoint_y);
			double target_x2 = target_xy2[0];
			double target_y2 = target_xy2[1];
			coarse_s_traj.push_back(target_s2);
			coarse_x_traj.push_back(target_x2);
			coarse_y_traj.push_back(target_y2);


			// next s values
			double target_s_dot = best_target[0][1];
			double current_s = pos_s;
			double current_v = pos_s_d;
			double current_a = pos_s_dd;
			for (int i = 0; i < (NUM_PATH_POINTS - pre_closepath_size); i++) {
				double v_incr, a_incr;
				if (fabs(target_s_dot - current_v) < 2 * VELOCITY_INCREMENT_LIMIT) {
					v_incr = 0;
				} else {

					v_incr = (target_s_dot - current_v)/(fabs(target_s_dot - current_v)) * VELOCITY_INCREMENT_LIMIT;
				}
				current_v += v_incr;
				current_s += current_v * PATH_DT;
				interpolated_s_traj.push_back(current_s);


			interpolated_x_traj = interpolate_points(coarse_s_traj, coarse_x_traj, interpolated_s_traj);
			interpolated_y_traj = interpolate_points(coarse_s_traj, coarse_y_traj, interpolated_s_traj);


			for(int i = 0; i < pre_closepath_size; i++) {
				next_x_vals.push_back(previous_path_x[i]);
				next_y_vals.push_back(previous_path_y[i]);
			} 
			// add xy points from newly generated path
			for (int i = 0; i < interpolated_x_traj.size(); i++) {
				//if (subpath_size == 0 && i == 0) continue; // maybe skip start position as a path point?
				next_x_vals.push_back(interpolated_x_traj[i]);
				next_y_vals.push_back(interpolated_y_traj[i]);
			} 

				
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
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

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
