This is my submission for the Udacity Self-Driving Car Nanodegree Path Planning Project. You can find my C++ source code [here](https://github.com/pszczesnowicz/SDC-P11-Path-Planning/tree/master/src). The goal of this project was to implement a path planner that can navigate a car on a three lane highway while respecting the speed, acceleration, and jerk limits.

# Summary

The essential bits of my code are simple logical statements and costs, yet they can navigate a car (ego car from here on) through traffic by taking the most efficient path.

# Sensor Fusion Data

The sensor fusion data of all the other cars on the same side of the road as the ego car are first filtered based on the distance of the other cars to the ego car (within 100 m ahead or 20 m behind) and then sorted into the appropriate lane. This was done to reduce computational time by ignoring cars that are far away and simplify the lane cost logic.

# Lane Costs

The filtered and sorted sensor fusion data are then looped through lane by lane and car by car to assign costs to the three lanes. The lane costs are stored in a vector with an element for each lane.

A lane change penalty is the first to be added to the lane costs. Its purpose is to prevent unnecessary lane changes when there is no traffic ahead.

```c
lane_costs[i] += abs(i - target_lane);
```

While looping through the lanes and the cars in them, the car to ego car distances, like the costs, are stored in a vector. These distances are used to determine the closest cars ahead of the ego car in each lane. The distances and speeds of the closest cars are used to calculate a penalty for their lanes.

```c
lane_costs[i] += (follow_dist / sensor_dist) * 2.0 + (max_speed / sensor_speed);
```

A buffer forward of the ego car is calculated using the current s-coordinate, speed, time between path points, and safety margin. The rearward buffer is calculated the same way except using the s-coordinate and speed of the other car.

```c
forward_buffer = car_s + car_speed * dt + 15.0;

rearward_buffer = sensor_s + sensor_speed * dt + 5.0;
```

If a car is within the two buffer brackets&mdash;like the black and red cars below&mdash;then a large value is added to its lane to prevent a lane change.

<img src="https://raw.githubusercontent.com/pszczesnowicz/SDC-P11-Path-Planning/master/readme_images/example1.jpg" width="640">

In the case that the ego car is in one of the outside lanes, there is a car adjacent to it&mdash;like the white car below&mdash;and within its buffer then a large value is also added to the empty lane adjacent to the car to prevent a lane change.

```c
if (target_lane != 1) {
  lane_costs[i - target_lane + 1] += 100.0;
}
```

<img src="https://raw.githubusercontent.com/pszczesnowicz/SDC-P11-Path-Planning/master/readme_images/example4.jpg" width="640">

Finally, the lane with the lowest cost is chosen as the next target lane. There is also a check to make sure the ego car does not make two simultaneous lane changes.

```c
if (abs(target_lane - current_lane) == 2) {
  target_lane = 1;
}
```

# Path Creation

I used the cubic spline interpolation library to create the path for the ego car. A spline was used to create a smooth trajectory in Frenet (s and d) coordinates which makes specifying the next target lane easier than in Cartesian (x and y) coordinates.

<img src="https://d17h27t6h515a5.cloudfront.net/topher/2017/July/595e74e6_frenet-5/frenet-5.png" width="640">

If the previous path has enough points (at least two) to create a spline tangent to it then those are used as the first two spline points. Otherwise, the ego car's coordinates and a point behind and in line with its heading are used as the first two spline points. The next three equally spaced spline points are created in Frenet coordinates and then transformed into Cartesian coordinates. This is also where the next target lane plays a role in defining the path. These five spline points (in vectors of Cartesian coordinates) are first shifted and rotated to the ego car's origin and heading and then used to fit a spline.

If the previous path has any remaining points they are added to the next path; this ensures a smooth transition between the two paths. The number of points required to travel a linear distance at the target speed is calculated and used to come up with the distance increment along the ego car's heading (x-axis). The difference between the number of remaining path points and desired number of path points is made up of points along the previously defined spline. These points are created by adding the distance increment to the x-coordinate and its corresponding y-coordinate (using the previously defined spline). Finally these new path points are rotated back to the map's heading and origin and sent to the simulator.

# Results

The ego car was able to drive 5 miles in just over 6 mins on the simulated track without any incident.

[<img src="https://raw.githubusercontent.com/pszczesnowicz/SDC-P11-Path-Planning/master/readme_images/path_planning.jpg" width="640">](https://youtu.be/zbxIljDAO5o "Click to watch")

# Conclusion

The combination of if-statements and tuned costs yields a surprisingly good path planner that makes safe driving decisions and stays within set limits.

# References

[Udacity Self-Driving Car ND](http://www.udacity.com/drive)

[Udacity Self-Driving Car ND - Path Planning Project](https://github.com/udacity/CarND-Path-Planning-Project)

[Cubic Spline Interpolation in C++](http://kluge.in-chemnitz.de/opensource/spline/)
