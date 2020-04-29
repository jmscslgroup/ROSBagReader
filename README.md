# ROSBagReader
__A MATLAB Class to facilitate the reading of a rosbag file based on semantic datatypes.__

__`ROSBagReader`__ is a wrapper class written in MATLAB that provides an easy to use interface for reading [bag files](http://wiki.ros.org/Bags) recorded by `rosbag record` command. This wrapper class uses MATLAB's API `rosbag()`
internally to perform all operations. One of the interesting features about using __`ROSBagReader`__ is that a user doesn't need to
supply [rostopic name](http://wiki.ros.org/rostopic) to extract relevant data. One can extract data based on the type of data the
user is seeking.

## Usage principle
The philosophy behind developing this project is to make everything as simple and less confusing as possible. As a result, there are
not too many options or freedom of usage being provided by __`ROSBagReader`__ class. If you need a wide variety of options, users
can directly use Robotics System Toolbox APIs provided by MATLAB. However, if you are looking for simplicity, __`ROSBagReader`__ is an
elegant choice.

An object/instance of __`ROSBagReader`__ class takes the absolute path of the bag file as an argument and provides several member functions
to extract data of choice. Data are extracted in a subfolder with the same name as of the provided bag file in the directory where bag
file is located. Further, data are extracted in two formats: csv and mat file. To further process the data in MATLAB, mat file is a 
recommended choice but CSV file format is a good start for amateurs or working with other languages such as Python or C++. However, 
image data are an exception. Image frames are extracted as png files in a subsubfolder of the subfolder with the same name as the topic name.
Image's timestamp is used as the file name for image frames.

## System Requirements and Dependencies
- Basic knowledge of Object-Oriented Program (OOP) is recommended but not required. Usages are simple but knowledge of
OOP can boost the understanding of the core principle behind this package.
- Operating System: Windows 8 or higher, Ubuntu 16.04 LTS or higher, Mac (not tested on Mac)
- MATLAB 2018a or higher (MATLAB 2019a or higher recommended)
- MATLAB's Robotic System Toolbox
- MATLAB's Parallel Computing Toolbox (recommended but not required)
- ROS Indigo or higher

## ROS Data type Supported
Currently, __`ROSBagReader`__ does not support every ROS datatype. However, __`ROSBagReader`__ supports only those ROS data types
that are important when working with Autonomous Vehicle Testbed [CAT Vehicle Testbed](https://github.com/jmscslgroup/catvehicle)
developed at the University of Arizona. The testbed publishes and subscribes to topics on laser data, camera data, velocity data,
steering data, brake and acceleration data. As a result, developers provide support for only the following data types:

- sensor_msgs/LaserScan (For laser data)
- geometry_msgs/Twist (For input and output velocity data)
- sensor_msgs/CompressedImage (For compressed images captured from cameras)
- nav_msgs/Odometry (For GPS odometry data)
- geometry_msgs/Wrench (For output steering, output brake, and output accelerator data)

## Installation and Usage
1. Clone this git repository 

`git clone https://github.com/jmscslgroup/ROSBagReader ROSBagReader`

2. Launch MATLAB. Add the directory where you clone this repository to MATLAB's search path or change the MATLAB's current path to that path.

3. Create an object of class __`ROSBagReader`__ in MATLAB with the absolute path of the bag file as argument as follows:

Linux Example:

`>> B = ROSBagReader('/home/wildcat/2019-09-11-14-59-34.bag');`

Windows Example:

`>> B = ROSBagReader('C:\Users\WildCat\Documents\VersionControl\2019-09-11-14-59-34.bag');`

4. Now we are ready to extract data from bag file:

Extract Laser Scan Data:

`>>B.extractLaserData();`

In this case, a mat file and a csv file is created in the folder `C:\Users\WildCat\Documents\VersionControl\2019-09-11-14-59-34\` 
with the same name as the topic name. Mat file is a struct datatype, while csv file has 182 columns: 1st column is a timestamp in POSIX format while column 2-182 are radial distance/range data.


Extract Velocity Data

`>>B.extractVelData();`
In this case, a mat file and a csv file is created in the folder `C:\Users\WildCat\Documents\VersionControl\2019-09-11-14-59-34\` 
with the same name as the topic name. Mat file is a struct datatype, while csv file has following columns:

*[timestamp in POSIX format, Linear.X, Linear.Y, Linear.Z, Angular.X, Angular.Y, Angular.Z]*

Extract Camera Data

`>>B.extractCompressedImages();`

In this case, a subfolder is created in the folder `C:\Users\WildCat\Documents\VersionControl\2019-09-11-14-59-34\` 
with the same name as the topic name. This subfolder contains all the extracted image frame. Each image frame has file name that denotes
its timestamp in POSIX format.

Extract Odometry Data

`>>B.extractOdometryData();`

In this case, a mat file and a csv file is created in the folder `C:\Users\WildCat\Documents\VersionControl\2019-09-11-14-59-34\` 
with the same name as the topic name. Mat file is a struct datatype, while csv file has following columns:

*[timestamp in POSIX format, sequence id, frameId, childFrameId, PoseX, PoseY, PoseZ, QueternionX, QueternionY, QueternionZ, QueternionW, Linear.X, Linear.Y, Linear.Z, Angular.X, Angular.Y, Angular.Z, Position Covariance (36 columns), Twist Covariance (36 Columns)]*

Extract Wrench Data

`>>B.extractWrenchData();`

In this case, a mat file and a csv file is created in the folder `C:\Users\WildCat\Documents\VersionControl\2019-09-11-14-59-34\` 
with the same name as the topic name. Mat file is a struct datatype, while csv file has following columns:

*[timestamp in POSIX format, ForceX, ForceY, ForceZ, TorqueX, TorqueY, TorqueZ]*

Note: For brake and accelerator data captured in bag files from CAT Vehicle, only timestamp and Force.X is populated and remaining columns stay zero, whereas for steering data, only Torque.Z is populated.

## Getting help
Once you have added ROSBagReader project folder to MATLAB's path using `addpath` command, you can tye `doc ROSBagReader` to read the documentation.

## Issues
If you run into any issues, please use the issue feature of the GitHub to log your issues. I will try my best to address any issue as soon as
possible.

## Contributing to this project
If you like to contribute to this project, please fork this repository to your GitHub account, create a new branch for yourself and
send a pull request for the merge. After reviewing the changes, we will decide if this is a good place to add your changes.
If you like to see new data types being supported, please raise an enhancement issue and provide a relevant bag file that contains the 
ros message of desired data types.

## Authors and Contributors
- Rahul Bhadani ( rahulbhadani@email.arizona.edu)

## Licensing

    License: MIT License 
    Copyright 2019-2020 Rahul Bhadani
    Initial Date: Sept 12, 2019
    Permission is hereby granted, free of charge, to any person obtaining 
    a copy of this software and associated documentation files 
    (the "Software"), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to 
    permit persons to whom the Software is furnished to do so, subject 
    to the following conditions:

    The above copyright notice and this permission notice shall be 
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF 
    ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED 
    TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT 
    SHALL THE AUTHORS, COPYRIGHT HOLDERS OR ARIZONA BOARD OF REGENTS
    BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
    AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE 
    OR OTHER DEALINGS IN THE SOFTWARE.
