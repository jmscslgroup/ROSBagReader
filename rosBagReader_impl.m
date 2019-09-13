% (C) Rahul Bhadani
% Useage of Bagreader class

B = ROSBagReader('/home/ivory/2019-09-12-16-20-25.bag');

% Extract Laser Scan Data
B.extractLaserData();

% Extract Velocity Data
B.extractVelData();

% Extract Camera Data
%B.extractCameraFrames();