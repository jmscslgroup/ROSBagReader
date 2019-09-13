classdef ROSBagReader < matlab.mixin.Copyable
    % Author: Rahul Bhadani
    % Maintainer Email: rahulbhadani@email.arizona.edu
    % License: MIT License 
    % Copyright 2019-2020 Rahul Bhadani
    % Permission is hereby granted, free of charge, to any person obtaining 
    % a copy of this software and associated documentation files 
    % (the "Software"), to deal in the Software without restriction, including
    % without limitation the rights to use, copy, modify, merge, publish,
    % distribute, sublicense, and/or sell copies of the Software, and to 
    % permit persons to whom the Software is furnished to do so, subject 
    % to the following conditions:

    % The above copyright notice and this permission notice shall be 
    % included in all copies or substantial portions of the Software.

    % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF 
    % ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED 
    % TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
    % PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT 
    % SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR 
    % ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
    % ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
    % OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE 
    % OR OTHER DEALINGS IN THE SOFTWARE.
    
    
    % Description 
    % This is a matlab class called as Bagreader that reads a bag file
    % and 
    
    %% Properties
    properties
        bagfile; % Full path of the bag  file, e.g /home/ece446/2019-08-21-22-00-00.bag
        filename; % Name of the bag file 2019-08-21-22-00-00.bag
        dir; % Directory where bag file is located
        bagReader;
        topicTable;
        availableTopics;
        numMessages;
        messageType;
        scan_select;
        datafolder;
    end
    
    %% Methods
    methods
        %% Constructors
        function obj = ROSBagReader(file)
            if nargin == 1
                obj.bagfile = file;
                
                file = file;
                slashIndices = strfind(file, '/');
                obj.filename = file(slashIndices(end)+1: end);
                

                if ~isempty(slashIndices)
                    obj.dir = file(slashIndices(1): slashIndices(end));
                else
                    obj.dir = '';
                end


                fprintf('\nReading the bag file %s\n\n', file);
                obj.bagReader = rosbag(file);
                obj.topicTable = obj.bagReader.AvailableTopics; % Get the topic table
                obj.availableTopics = obj.bagReader.AvailableTopics.Row; % Store the available Topic
                obj.numMessages = obj.bagReader.AvailableTopics.NumMessages; % store number of messages for each topic
                obj.messageType = obj.bagReader.AvailableTopics.MessageType; % Store message type for each topic
                
                obj.datafolder = strcat(obj.dir ,obj.filename(1:end-4),'/');
                
                [~, ~, msgID] =  mkdir(obj.datafolder );
                if( strcmp(msgID, 'MATLAB:MKDIR:DirectoryExists') )
                    fprintf('\nData folder for %s already exists.\n\n', file);
                end
                
            else
                aac = matlab.lang.correction.AppendArgumentsCorrection('"bagfilefullpath"');
                error(aac, 'MATLAB:notEnoughInputs', 'Too many input arguments.')  
            end
        end
        
        % Function to extract  Laser Scan Data and save in the file where
        % rosbag file is located
        function extractLaserData(obj)
            
            % Note down the index of laserscan messages from table topic
            % Note that there can be multiple topics of
            % sensor_msgs/LaserScan type. extractLaserData function
            % retrieves all such topics of type sensor_msgs/LaserScan
            index_scan = obj.messageType == 'sensor_msgs/LaserScan';
            
            % find all the topics of type sensor_msgs/LaserScan
            topic_of_scan = obj.availableTopics(index_scan);
            num_msgs = obj.numMessages(index_scan);
            % Now we will iterate over all the topics and retrieve data
            % from every topic of message type sensor_msgs/LaserScan
            % and save them in mat format as well as csv format
            for i = 1:length(topic_of_scan)
                topic_to_read = topic_of_scan{i};
                fprintf('\nReading the Laser Scan messages on the topic %s\n\n', topic_to_read);
                
                % Select object that selects particular topic from rosbag object
                 obj.scan_select = select(obj.bagReader, 'Topic', topic_to_read);

                % Save all messages from this topic in an struct format
                msgStructs = readMessages(obj.scan_select,'DataFormat','struct');
                
                % since ROS topic names contain slashes, we will replace
                % every slash with a dash. This topic name will be used to
                % create the mat file, meaning mat file that saves message
                % has same name as the corresponding topic
                
                topic_to_save = strrep(topic_to_read,'/','-');
                
                % Now save the retrieved data in the datafolder
                matfile = strcat(obj.datafolder, topic_to_save,'.mat');
                
                save(matfile,'msgStructs');

                fprintf('Writing Laser Scan Data  to file %s from topic %s completed!!', matfile, topic_to_read);

                                
                num_messages = num_msgs(i);
                
                Data = zeros(num_messages, length(msgStructs{1}.Ranges) + 1);
                
                for i = 1:num_messages
                    
                    tempLen = length(msgStructs{i}.Ranges);
                    % TODO: How to node hardcode 181.
                    % We are doing this because for a particular timestep,
                    % there are less than 181 datapoints
                    if(tempLen < 181)
                         NaNPOSTFIX = zeros( 181 - tempLen, 1)*NaN;
                         msgStructs{i}.Ranges = [msgStructs{i}.Ranges; NaNPOSTFIX];
                    end
                    
                    secondtime = double(msgStructs{i}.Header.Stamp.Sec)+double(msgStructs{i}.Header.Stamp.Nsec)*10^-9;
                    Data(i, 1) = secondtime;
                    Data(i, 2:end) = msgStructs{i}.Ranges;
                end
                
                csvfile = strcat(obj.datafolder, topic_to_save,'.csv');
                writematrix(Data,csvfile,'Delimiter',',');
                
                fprintf('Writing Laser Scan Data  to file %s from topic %s completed!!\n\n', csvfile, topic_to_read);
               
                
            end
        end
    end % end of methods
    
            
end % of classdef

    
