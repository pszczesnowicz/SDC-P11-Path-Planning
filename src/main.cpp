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
#include "spline.h"

using namespace std;
using json = nlohmann::json;

constexpr double pi() { return M_PI; }

/**
 Degrees to radians conversion

 @param x angle in degrees
 @return angle in radians
 */
double deg2rad(double x) { return x * pi() / 180; }

/**
 Radians to degrees conversion

 @param x angle in radians
 @return angle in degrees
 */
double rad2deg(double x) { return x * 180 / pi(); }

// Checking if the SocketIO event has JSON data
// If there is data, the JSON object in string format will be returned,
// else the empty string "" will be returned
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  }
  else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

/**
 Distance between two coordinate points

 @param x1 first x coordinate
 @param y1 first y coordinate
 @param x2 second x coordinate
 @param y2 second y coordinate
 @return distance
 */
double distance(double x1, double y1, double x2, double y2) {
  
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

/**
 <#Description#>

 @param x <#x description#>
 @param y <#y description#>
 @param maps_x <#maps_x description#>
 @param maps_y <#maps_y description#>
 @return <#return value description#>
 */
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y) {

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++) {
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
    
		if(dist < closestLen) {
			closestLen = dist;
			closestWaypoint = i;
		}
	}
	return closestWaypoint;
}


/**
 <#Description#>

 @param x <#x description#>
 @param y <#y description#>
 @param theta <#theta description#>
 @param maps_x <#maps_x description#>
 @param maps_y <#maps_y description#>
 @return <#return value description#>
 */
int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y) {

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];
	double heading = atan2( (map_y-y),(map_x-x) );
	double angle = abs(theta-heading);

	if(angle > pi()/4) {
		closestWaypoint++;
	}
	return closestWaypoint;
}

/**
 Transform from Cartesian x,y coordinates to Frenet s,d coordinates

 @param x <#x description#>
 @param y <#y description#>
 @param theta <#theta description#>
 @param maps_x <#maps_x description#>
 @param maps_y <#maps_y description#>
 @return <#return value description#>
 */
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y) {
	
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);
	int prev_wp;
	prev_wp = next_wp-1;
  
	if(next_wp == 0) {
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

	if(centerToPos <= centerToRef) {
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++) {
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};
}

/**
 Transform from Frenet s,d coordinates to Cartesian x,y

 @param s <#s description#>
 @param d <#d description#>
 @param maps_s <#maps_s description#>
 @param maps_x <#maps_x description#>
 @param maps_y <#maps_y description#>
 @return <#return value description#>
 */
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y) {
	
  int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1))){
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

  // Waypoint map to read from
//  string map_file_ = "../data/highway_map.csv";
  string map_file_ = "/Users/skunkworks/programming/SDCND/Projects/SDCND-P11-PathPlanning/data/highway_map.csv";
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
  
  int target_lane = 1;
  double target_speed = 0;
  double max_speed = 49.0 * 0.447;
  
  h.onMessage([&map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy,
               &max_s, &target_lane, &target_speed, &max_speed]
              (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Ego's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
//          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"]; // Degrees
          double car_speed = j[1]["speed"]; // Miles per hour
          car_speed *= 0.447; // Meters per second
          
          // Previous path data given to the simulator
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          
          // Previous path's end s and d values
//          double end_path_s = j[1]["end_path_s"];
//          double end_path_d = j[1]["end_path_d"];
          
          // ###################################
          // # Sensor fusion data sorting step #
          // ###################################
          
          // Sensor fusion data: a list of all other cars on the same side of the road
          // 1 = x, 2 = y, 3 = vx, 4 = vy, 5 = s, 6 = d
          auto sensor_fusion = j[1]["sensor_fusion"];
          
          vector<vector<double>> sensor_fusion_left;
          vector<vector<double>> sensor_fusion_center;
          vector<vector<double>> sensor_fusion_right;
          
          // Looping through the list of all other cars on the same side of the road as ego
          for (int i = 0; i < sensor_fusion.size(); ++i) {
            vector<double> sensor_car = sensor_fusion[i];
            // Checking if the car is within a zone of 100 meters ahead and 20 meters behind ego
            if ((sensor_car[5] < (car_s + 100.0)) && (sensor_car[5] > (car_s - 20.0))) {
              // Sorting the car into lane 0
              if ((sensor_car[6] > 0) && (sensor_car[6] < 4.0)) {
                sensor_fusion_left.push_back(sensor_car);
              }
              // Sorting the car into lane 1
              else if ((sensor_car[6] > 4.0) && (sensor_car[6] < 8.0)) {
                sensor_fusion_center.push_back(sensor_car);
              }
              // Sorting the car into lane 2
              else if ((sensor_car[6] > 8.0) && (sensor_car[6] < 12.0)) {
                sensor_fusion_right.push_back(sensor_car);
              }
            }
          }
          
          vector<vector<vector<double>>> sensor_fusion_sorted;
          
          sensor_fusion_sorted.push_back(sensor_fusion_left);
          sensor_fusion_sorted.push_back(sensor_fusion_center);
          sensor_fusion_sorted.push_back(sensor_fusion_right);
          
          // ##################
          // # Lane cost step #
          // ##################
          
          vector<double> car_dists;
          car_dists.assign(3, max_s);
          
          vector<double> car_speeds;
          car_speeds.assign(3, max_speed);
          
          vector<double> lane_costs;
          lane_costs.assign(3, 0.0);
          
          // Time between each path point in seconds
          double dt = 0.02;
          
          // Buffer forward of ego
          double forward_buffer = car_s + car_speed * dt + 15.0;
          
          double follow_dist = 30.0;
          
          // Looping through the lanes
          for (int i = 0; i < 3; ++i) {
            
            // Adding a lane change penalty to the lanes other than ego's
            lane_costs[i] += abs(i - target_lane);
            
            // Looping through the sorted list of all other cars on the same side of the road as ego
            for (int j = 0; j < sensor_fusion_sorted[i].size(); ++j) {
              
              double sensor_s = sensor_fusion_sorted[i][j][5];
              double sensor_dist = sensor_s - car_s;
              double sensor_speed = sqrt(pow(sensor_fusion_sorted[i][j][3], 2.0) +
                                         pow(sensor_fusion_sorted[i][j][4], 2.0));
              
              // Buffer forward of other car; rearward of ego
              double rearward_buffer = sensor_s + sensor_speed * dt + 5.0;
              
              // Checking if the other car is the nearest one in its lane ahead of ego
              if ((sensor_dist > 0) && (sensor_dist < car_dists[i])) {
                car_dists[i] = sensor_dist;
                car_speeds[i] = sensor_speed;
                
                // Adding a distance and speed penalty to the lane in which the other car is in
                lane_costs[i] += (follow_dist / sensor_dist) * 2.0 + (max_speed / sensor_speed);
              }
              
              // Checking if the car is in between ego's forward and rearward buffers
              if ((sensor_s < forward_buffer) && (car_s < rearward_buffer)) {
                // Adding a large value to the lane in which the other car is in to prevent a lane change
                lane_costs[i] += 100.0;
                
                // Checking if ego is in one of the outer lanes
                if (target_lane != 1) {
                  // Adding a large value to the lane adjacent to the one the other car is in to prevent a lane change
                  lane_costs[i - target_lane + 1] += 100.0;
                }
              }
            }
          }
          
          cout << lane_costs[0] << " " << lane_costs[1] << " " << lane_costs[2] << endl;
          
          // ##############################
          // # Lane and speed change step #
          // ##############################
          
          int current_lane = target_lane;
          double lowest_cost = 100.0;
          
          // Looping through the lanes to find the one with the lowest cost
          for (int i = 0; i < 3; ++i) {
            if (lane_costs[i] < lowest_cost) {
              target_lane = i;
              lowest_cost = lane_costs[i];
            }
          }
          
          // Checking if the new target lane is 2 lane changes away
          if (abs(target_lane - current_lane) == 2) {
            // Setting the new target lane to the center lane to prevent 2 simultaneous lane changes
            target_lane = 1;
          }
          
          // Checking if ego is stuck behind a slower car
          if ((target_lane == current_lane) && (car_dists[target_lane] <= follow_dist)) {
            // Decelerating to match the slower car
            target_speed -= (car_speed - car_speeds[target_lane]) / 100.0;
          }
          else {
            // Accelerating to match the speed limit
            target_speed += (max_speed - car_speed) / 200.0;
          }
          
          // ######################
          // # Path creation step #
          // ######################
          
          // Number of points in the previous path
          int prev_path_size = previous_path_x.size();
          
          // Reference state
          double ref_x;
          double ref_y;
          double ref_yaw;
          
          double prev_x;
          double prev_y;
          
          // Checking if the previous path has enough points to create a tangent spline
          if (prev_path_size < 2) {
            // Using the car's coordinates as the second spline point
            ref_x = car_x;
            ref_y = car_y;
            
            ref_yaw = deg2rad(car_yaw);
            // Creating a point behind the car to be the first spline point
            prev_x = ref_x - cos(ref_yaw);
            prev_y = ref_y - sin(ref_yaw);
          }
          else {
            // Using the previous path's last point as the second spline point
            ref_x = previous_path_x[prev_path_size - 1];
            ref_y = previous_path_y[prev_path_size - 1];
            // Using the previous path's second last point as the first spline point
            prev_x = previous_path_x[prev_path_size - 2];
            prev_y = previous_path_y[prev_path_size - 2];
            // Calculating the car's yaw angle
            ref_yaw = atan2((ref_y - prev_y), (ref_x - prev_x));
          }
          
          // Spline points vectors
          vector<double> spline_ptsx;
          vector<double> spline_ptsy;
          
          // First two spline points
          spline_ptsx.push_back(prev_x);
          spline_ptsx.push_back(ref_x);
          spline_ptsy.push_back(prev_y);
          spline_ptsy.push_back(ref_y);
          
          // Next three spline points
          for (int i = 1; i <= 3; ++i) {
            // Creating evenly spaced points in Frenet coordinates and then transforming them into x and y coordinates
            vector<double> next_pt = getXY((car_s + 50.0 * i), (2.0 + 4.0 * target_lane),
                                           map_waypoints_s, map_waypoints_x, map_waypoints_y);
            
            spline_ptsx.push_back(next_pt[0]);
            spline_ptsy.push_back(next_pt[1]);
          }
          
          // Transforming the spline points from map to car coordinates
          // x, y, and yaw = 0 in car coordinates
          for (int i = 0; i < spline_ptsx.size(); ++i) {
            
            // Shifting the spline points to the car's origin
            double shift_x = spline_ptsx[i] - ref_x;
            double shift_y = spline_ptsy[i] - ref_y;
            
            // Rotating the spline points to the car's x-axis (heading)
            spline_ptsx[i] = shift_x * cos(-ref_yaw) - shift_y * sin(-ref_yaw);
            spline_ptsy[i] = shift_y * cos(-ref_yaw) + shift_x * sin(-ref_yaw);
            
          }
          
          // Initializing the spline with the points created above
          tk::spline s;
          s.set_points(spline_ptsx, spline_ptsy);
          
          vector<double> next_ptsx;
          vector<double> next_ptsy;
          
          // Adding the remaining points in the previous path to the next path
          for (int i = 0; i < prev_path_size; ++i) {
            
            next_ptsx.push_back(previous_path_x[i]);
            next_ptsy.push_back(previous_path_y[i]);
            
          }
          
          // Ego's horizon
          double target_x = 30;
          double target_y = s(target_x);
          
          // Linear distance to the horizon
          double target_dist = sqrt(pow(target_x, 2.0) + pow(target_y, 2.0));
          
          // Number of points required to travel the linear distance at the target velocity
          int num_pts = target_dist / (target_speed * dt);
          
          // Distance increment along the x axis that results in the target velocity
          double x_inc = target_x / num_pts;
          
          // Creating more Pac-Man nuggets :P
          for (int i = 1; i <= (50 - prev_path_size); ++i) {
            
            double next_x = x_inc * i;
            double next_y = s(next_x);
            
            double car_coor_x = next_x;
            double car_coor_y = next_y;
            
            // Rotating the spline points back to the map's coordinates
            next_x = car_coor_x * cos(ref_yaw) - car_coor_y * sin(ref_yaw);
            next_y = car_coor_y * cos(ref_yaw) + car_coor_x * sin(ref_yaw);
            
            // Shifting the spline points back to the map's coordinates
            next_x += ref_x;
            next_y += ref_y;
            
            next_ptsx.push_back(next_x);
            next_ptsy.push_back(next_y);
            
          }
          
          // Sending the next path points to the simulator
          json msgJson;
          msgJson["next_x"] = next_ptsx;
          msgJson["next_y"] = next_ptsy;
          auto msg = "42[\"control\","+ msgJson.dump()+"]";
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      }
      
      // Manual driving mode
      else {
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });
  
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    }
    else {
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  }
  else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
