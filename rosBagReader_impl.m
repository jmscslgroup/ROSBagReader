% (C) Rahul Bhadani
% Useage of Bagreader class

% Windows Example
%B = ROSBagReader('C:\Users\Ivory - rahulbhadani\Documents\VersionControl\ROSBagReader\2019-09-11-14-59-34.bag');

% Linux/Mac Example
B = ROSBagReader('/home/ivory/CyverseData/JmscslgroupData/ARED/2016-07-28/data_by_test/Bag Files/test3_09-23-59.bag');

% Extract Laser Scan Data
% B.extractLaserData();

% Extract Velocity Data
% V = B.extractVelData();
% plot velocity data
% B.ts_plotVelData();

% B.extractStandard();
% % plot standard data
% B.ts_plotStandard();
%  
% % Extract Camera Data
% B.extractCompressedImages();
% 
% % Extract Odometry Data
% B.extractOdometryData();
% % plot odometry data
% B.ts_plotOdometryData( );
% 
% % Extract Wrench Data
B.extractWrenchData();