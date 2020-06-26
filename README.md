# ERKF
Source code and sample data for estimating IMU orientation using an Error State Kalman Filter (ESKF).

orientation_ErSKF.m is the main source code file.

validation1.mat contains IMU data in the data structure format the code is expecting.

Note that sensor characteristics in the 'define noise matrices' section on lines 63-82 correspond to APDM Opal IMUs (https://www.apdm.com/wearable-sensors/).

Developed with Matlab Version 2019b. Dependencies on Statistics and Machine Learning Toolbox, Navigation Toolbox, Robotics System Toolbox, and ROS Toolbox.
-motionDynamics.m    -> chi2inv function depends on the Statistics and Machine Learning Toolbox
-magneticFiled.m     -> chi2inv function depends on the Statistics and Machine Learning Toolbox
                        quat2eul function depends on Navigation Toolbox, Robotics System Toolbox, and ROS Toolbox
-orientation_ErSKF.m -> rotm2quat function depends on Navigation Toolbox, Robotics System Toolbox, and ROS Toolbox
                        quat2eul function depends on Navigation Toolbox, Robotics System Toolbox, and ROS Toolbox
