% (C) Rahul Bhadani
% Useage of Bagreader class

B = ROSBagReader('/home/ivory/2019-09-11-14-59-34.bag');

% Extract Laser Scan Data
B.extractLaserData();

% Extract Velocity Data
B.extractVelData();
 
% Extract Camera Data
B.extractCompressedImages();

% Extract Odometry Data
Data = B.extractOdometryData();

