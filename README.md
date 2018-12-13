To run Gazebo simulation:
 - command line:
   #+BEGIN_SRC sh
     roslaunch infinite_interaction gazebo_interface.launch
   #+END_SRC
 - now both joint positions and force torque sensor reading should be
   available! You can check by listing the topics with

      rostopic list

   and check that the topics are shown.


Notes on loading the system parameters:
 - load parameters with:
   #+BEGIN_SRC sh
     roscd infinite_interaction/config && rosparam load admittance_params.yaml
   #+END_SRC


Notes on running the actual system:
- Using the "right Denso" now, the left one is moved to Eureka
- commands
  #+BEGIN_SRC sh
    # robot driver and joint position control interface
    roslaunch denso_control rc8_ros_driver.launch rate:=125 ip:=192.168.0.21
    roslaunch denso_control joint_position_controllers.launch

    # ftsensor driver and ros interface
    roslaunch netft_control netft_ros_driver.launch ip:=192.168.0.22
    roslaunch netft_control ros_controllers.launch

    # or, instead of the above step: launch netft via a node
	roslaunch infinite_interaction load_netft.launch

  #+END_SRC

Cartesian control experiment

roslaunch infinite_interaction admittance_controller.yaml

# Requirements:

## resource at crigroup
  - robotic_description

## Eigen 3.3

Best is to built Eigen from source
  
  ```bash
  wget http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2
  tar xvf eigen-eigen-5a0156e40feb.tar.bz2
  cd eigen-eigen-5a0156e40feb
  mkdir build && cd build
  cmake .. && make install
  ```
## qpOASES

  ``` bash
  git clone https://github.com/hungpham2511/qpoases
  cd qpoases && mkdir build -p && cd build
  cmake .. 
  make && make install
  ```
	
