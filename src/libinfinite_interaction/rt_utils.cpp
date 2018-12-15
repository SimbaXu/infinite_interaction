//
// Created by hung on 14/12/18.
//

#include "infinite_interaction/infinite_interaction_lib.h"

namespace RTUtils {

    void increment_timespec(timespec &tspec, const int &inc_nsec) {
        tspec.tv_nsec += inc_nsec;
        if (tspec.tv_nsec > NANOSECS_IN_SEC) {
            tspec.tv_sec += 1;
            tspec.tv_nsec -= NANOSECS_IN_SEC;
        }
    }


    inline void diff_timespec(int &diff_nsec, timespec &tspec1, timespec &tspec2) {
        diff_nsec = NANOSECS_IN_SEC * (tspec2.tv_sec - tspec1.tv_sec) + (tspec2.tv_nsec - tspec1.tv_nsec);
    }


    int set_policy_fifo() {
        // set scheduling policy and priority
        struct sched_param param;
        int ret;
        param.sched_priority = 80;
        ret = sched_setscheduler(0, SCHED_FIFO, &param);
        if (ret == -1) {
            ROS_ERROR("Unable to switch policty to FIFO");
            return 1;
        }
        return 0;
    }
}

