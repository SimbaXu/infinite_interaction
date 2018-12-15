cmake_minimum_required(VERSION 2.8.3)
project(infinite_interaction)

set(CMAKE_CXX_STANDARD 11)
# to prevent optimizing out variable during debugging, uncomment the
# below line
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

find_package(catkin REQUIRED COMPONENTS
		roscpp
		rospy
		std_msgs
		denso_control
)

find_package(Eigen3 3.3 REQUIRED)

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
		include
		${OpenRAVE_INCLUDE_DIRS}
		${CMAKE_CURRENT_SOURCE_DIR}
		${catkin_INCLUDE_DIRS}
)
catkin_package()

add_library(infinite_interaction
		src/libinfinite_interaction/infinite_interaction.cpp
		src/libinfinite_interaction/filters.cpp
		src/libinfinite_interaction/rt_utils.cpp
		)
target_link_libraries(infinite_interaction
		${catkin_LIBRARIES}
		Eigen3::Eigen
		)
target_include_directories(infinite_interaction PUBLIC include)


# admittance controller: (main program)
add_executable(joint_admittance_controller src/admittance_controller.cpp )
target_include_directories(joint_admittance_controller PUBLIC include)
target_link_libraries(joint_admittance_controller
		${catkin_LIBRARIES}
		infinite_interaction qpOASES
		# openrave dependencies
		${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES}
		${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY})

set_target_properties(joint_admittance_controller PROPERTIES COMPILE_FLAGS "${OpenRAVE_CXX_FLAGS}")
set_target_properties(joint_admittance_controller PROPERTIES LINK_FLAGS "${OpenRAVE_LINK_FLAGS}")

#add_executable(sinuisoidal_motion_generator src/sinuisoidal_motion_generator.cpp)
#target_include_directories(sinuisoidal_motion_generator PUBLIC include)
#target_link_libraries(sinuisoidal_motion_generator
#        ${catkin_LIBRARIES}
#		infinite_interaction qpOASES
#		# openrave depencidencies
#        ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES}
#        ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY}
#  )
#set_target_properties(sinuisoidal_motion_generator PROPERTIES COMPILE_FLAGS "${OpenRAVE_CXX_FLAGS}")
#set_target_properties(sinuisoidal_motion_generator PROPERTIES LINK_FLAGS "${OpenRAVE_LINK_FLAGS}")

add_executable(test_discrete_time_filter src/tests/test_discrete_time_filter.cpp
		src/libinfinite_interaction/infinite_interaction.cpp
		src/libinfinite_interaction/filters.cpp
		)
target_link_libraries(test_discrete_time_filter gtest_main
        ${catkin_LIBRARIES}
		qpOASES Eigen3::Eigen
        )
target_include_directories(test_discrete_time_filter PUBLIC include)


add_executable(test_fir_controller src/tests/test_fir_controller.cpp
		src/libinfinite_interaction/infinite_interaction.cpp
		src/libinfinite_interaction/filters.cpp
		)
target_link_libraries(test_fir_controller gtest_main
		${catkin_LIBRARIES}
		qpOASES Eigen3::Eigen
        )
target_include_directories(test_fir_controller PUBLIC include)

add_executable(test_cart_position_tracker src/tests/test_cart_position_tracker.cpp
		src/libinfinite_interaction/infinite_interaction.cpp
		src/libinfinite_interaction/filters.cpp
		)
target_link_libraries(test_cart_position_tracker gtest_main
        ${catkin_LIBRARIES}
        qpOASES Eigen3::Eigen
        # openrave depencidencies
        ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES}
        ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY}
        )
target_include_directories(test_cart_position_tracker PUBLIC include)

add_executable(test_topic_debugger src/tests/test_topic_debugger.cpp
		src/libinfinite_interaction/infinite_interaction.cpp
		)
target_include_directories(test_topic_debugger PUBLIC include)
target_link_libraries(test_topic_debugger gtest
        ${catkin_LIBRARIES} )
