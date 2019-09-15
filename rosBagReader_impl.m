% (C) Rahul Bhadani
% Useage of Bagreader class

% Windows Example
%B = ROSBagReader('C:\Users\Ivory - rahulbhadani\Documents\VersionControl\ROSBagReader\2019-09-11-14-59-34.bag');

% Linux/Mac Example
B = ROSBagReader('/home/ivory/2019-09-11-14-59-34.bag');

% Extract Laser Scan Data
B.extractLaserData();

% Extract Velocity Data
B.extractVelData();
 
% Extract Camera Data
B.extractCompressedImages();

% Extract Odometry Data
B.extractOdometryData();

% Extract Wrench Data
B.extractWrenchData();