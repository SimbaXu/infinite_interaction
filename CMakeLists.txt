cmake_minimum_required(VERSION 2.8.3)
project(infinite_interaction)

set(CMAKE_CXX_STANDARD 11)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
)

# openrave dep
find_package(OpenRAVE REQUIRED)
find_package(Boost ${OpenRAVE_Boost_VERSION} EXACT COMPONENTS iostreams python thread system)
if( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
    add_definitions("-fno-strict-aliasing -Wall")
endif( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
# end openrave dep

if( Boost_INCLUDE_DIRS )
    include_directories(${Boost_INCLUDE_DIRS})
endif()
include_directories(
  ${OpenRAVE_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} 
  ${catkin_INCLUDE_DIRS}
)

catkin_package()
add_executable(admittance_controller src/admittance_controller.cpp)
target_include_directories(admittance_controller PUBLIC include)
target_link_libraries(admittance_controller
  ${catkin_LIBRARIES}
  # openrave depencidenc
  ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES}
  ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY}
  )
set_target_properties(admittance_controller PROPERTIES COMPILE_FLAGS "${OpenRAVE_CXX_FLAGS}")
set_target_properties(admittance_controller PROPERTIES LINK_FLAGS "${OpenRAVE_LINK_FLAGS}")


