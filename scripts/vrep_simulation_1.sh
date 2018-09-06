#!/bin/bash          

# Parameters
VREP_DIR="/home/hung/apps/V-REP_PRO_EDU_V3_5_0_Linux"
VREP_SCENARIO="data/vrep_scene/test_kuka.ttt"

# ref: http://www.coppeliarobotics.com/helpFiles/en/commandLine.htm
for dir in ${ROS_PACKAGE_PATH//:/ }
do
    if [ -f "$dir/infinite_interaction/$VREP_SCENARIO" ]; then
	echo "-- Scenario $VREP_SCENARIO found! Starting Vrep."
	/bin/bash $VREP_DIR/vrep.sh $dir/infinite_interaction/$VREP_SCENARIO -gREMOTEAPISERVERSERVICE_19999_FALSE_TRUE
	echo "-- Terminating vrep binary."
	exit 1
    fi
done

exit 0
