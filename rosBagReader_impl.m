% (C) Rahul Bhadani
% Useage of Bagreader class

% Windows Example
%B = ROSBagReader('C:\Users\Ivory - rahulbhadani\Documents\VersionControl\ROSBagReader\2019-09-11-14-59-34.bag');

% Linux/Mac Example
B = ROSBagReader('/home/ivory/VersionControl/catvehicle_ws/src/sparkle/src/Circle_Test_n_6_updateRate_20.0_max_update_rate_50.0_time_step_0.01_logtime_35.0_laser_true_2020-01-10-15-23-23.bag');

% Extract Laser Scan Data
B.extractLaserData();

% Extract Velocity Data
B.extractVelData();
% plot velocity data
B.ts_plotVelData();

B.extractStandard();
% plot standard data
B.ts_plotStandard();
 
% Extract Camera Data
B.extractCompressedImages();

% Extract Odometry Data
B.extractOdometryData();
% plot odometry data
B.ts_plotOdometryData( );

% Extract Wrench Data
B.extractWrenchData();