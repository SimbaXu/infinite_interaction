# To run force control experiment

``` shell
# sliding experiment (Jan 2019)
roslaunch infinite_interaction 2019_Jan_sliding_experiment.launch

# teaching experiment (Jan 2019)
2019_jan_teaching_experiment.launch
```

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
	
