function save_Thermal_Images(dirin,dirout,start,stop,undistort,metadata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function extracts the data from the DJI .JPG files. For the visible
% images, it will save the image into a matlab file. For the thermal
% images, it will save the decoded thermal image as well as the thermal
% calibration. Both visible and thermal images will also save the metadata
% as a table, conveniently called table.
%
% NOTE: You should have this code in a directory one step lower than your
% directory for the JPG images. For example, I have this in a directory
% called 'DJI_Mavic2EA_IR_test/Codes/' and my images (specified by dirin)
% are in a directory 'DJI_Mavic2EA_IR_test/dirin'
%
% Helper codes you'll need: 
% - Exiftool (installed in your bin/ directory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters:
% dirin - the path to the directory containing the image files
%
% dirout - the name and location of the directory where images will be
%           saved (make sure to end both with '/')
%
% Note: If having trouble with dirin, dirout, use the full path for each
%
% start - (optional) the starting index of images (defaults to 1)
%
% stop - (optional) the ending index of images (defaults to last one)
%
% undistort - (optional) if set to 1, will also undistort the images from
%                        the camera parameters given from Matlab's camera
%                        calibration saved in DJI_Mavic2EA_IR_test in
%                        Camera_Calibration and Camera_Calibration_Vis.
%                        Defaults to 1.
%
% metadata (optional) - if 1, will also calculate the metadata and save to
%                       folder Metadata/ within dirout. If you already have
%                       the metadata you can leave this blank or 0
%
% cameraParams_vis, cameraParams_therm (optional)
%  - if left blank, it will load from camera_calibration folders
%  - can load in yourself
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Aidan Blaser (ablaser@ucsd.edu)
% Last Edited: 09/07/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
    metadata = 1;
end

if nargin < 5
    undistort = 1;
end

if nargin < 3
    start = 1;
    stop = [];
end

% First, create the outgoing directory
mkdir(dirout);
% Create directory for metadata
mkdir(strcat(dirin,'/Metadata'))
% Load in camera parameters (for undistorting)
if undistort
    load('camera_calibration.mat','cameraParams_vis','cameraParams_therm')
end

% Turn image data into metadata
if metadata
    current_folder = pwd;
    %Give matlab access to exiftool
    setenv('PATH', getenv('PATH')+":/usr/local/bin")
    % Go to directory with images
    cd(dirin)
    % Perform exiftool and pipe it to metadata folder in dirout 
    status = ['for i in ','*.JPG; do exiftool -b -n -csv $i > ','Metadata/${i%%.JPG}.csv; done'];
    stat = system(status);
    cd(current_folder)
end

% Extract image files
files = dir(strcat(dirin,'Metadata/','*.csv'));
files_raw = dir(strcat(dirin,'*.JPG'));

if isempty(stop)
    stop = length(files);
end

for i=start:stop
    filename = files(i).name; 
    disp(filename)
    % Next read in this .csv file as a table
    table = readtable(strcat(dirin,'Metadata/',filename));
    image = imread(strcat(dirin,files_raw(i).name)); %save images
    calibration = []; %empty variable for saving nonthermal images
    
    %Do thermal data only for odd (thermal) images
    try table.ThermalData;
        therm = 1; %Marker for therm
        % Pull out relevant data (ThermalData)
        data_encoded = table.ThermalData;
        calibration_encoded = table.ThermalCalibration;
        % Decode this data
        if strcmp(data_encoded{1}(1:7),'base64:')
            % Note: Start from position 8 to remove the 'base64:'
            data = matlab.net.base64decode(data_encoded{1}(8:end));
        else
            %For some reason some of the data is encoded in UTF-8
            data = unicode2native(data_encoded{1},'UTF-8');
        end
        if length(data) ~= 655360
            data = unicode2native(data_encoded{1},'latin1');
        end
        calibration = matlab.net.base64decode(calibration_encoded{1}(8:end));
        % Note: This data should have twice the number of bits,
        %       correct for this by the following formula
        a = double(data(1:2:end)); %odd values
        b = double(data(2:2:end)); %even values
        data_corrected = b*256 + a;
        a = double(calibration(1:2:end));
        b = double(calibration(2:2:end));
        calibration = b*256 + a;
        % Resize this image according to its mxn size
        image = reshape(data_corrected,640,512);
        % Orient properly
        image = rot90(fliplr(image));
        
    catch
        therm = 0; %Marker for therm
        
    end
    
    %Undistort images
    if undistort
        try
            if therm
                [image,~] = undistortImage(image,cameraParams_therm);
            else
                [image,~] = undistortImage(image,cameraParams_vis);
            end
        catch
            disp('image unable to be undistorted')
        end
    end
    
    [~,w_vis,~] = size(image);
    
    %Calculate spatial resolution
    %If there's RTK, use RTK, otherwise, use relative altitude
    try
        rtk = table.RtkFlag;
        altitude = table.AbsoluteAltitude + 33.69; %subtract geoid
    catch
        altitude = table.RelativeAltitude;
    end
    %Modify altitude for slight roll and pitch variations
    roll = deg2rad(table.GimbalRollDegree);
    pitch = deg2rad(table.GimbalPitchDegree+90);
    altitude_corrected = altitude*sqrt(1+sqrt((sin(roll)^2 + sin(pitch)^2)));
    % Next, calculate the resolutions (m/pix) for each image
    if therm
        res = (12.15e-3)*(altitude_corrected/9);
    else
        res = (6.17/w_vis)*(altitude_corrected/4.5);
    end
    
    
    save(strcat(dirout,filename(1:end-4)),'table','image','calibration','res','altitude_corrected');

end