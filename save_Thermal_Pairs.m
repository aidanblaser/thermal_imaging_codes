function save_Thermal_Pairs(dirin,dirout,start,stop,undistort,metadata,geotiff)
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
% start - (optional) the starting index of visible images (defaults to 1)
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
% geotiff (optional) - if 1, will also save the images as geotiffs in a
%                      folder in the outgoing directory Geotiff/. Defaults
%                      to 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Aidan Blaser (ablaser@ucsd.edu)
% Last Edited: 09/13/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<7
    geotiff = 0;
end

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
if geotiff
    mkdir(strcat(dirout,'Geotiff'))
end
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
counter = 0;
for i=start:2:stop-1
    counter = counter + 1;
    filename_vis = files(i).name; 
    filename_therm = files(i+1).name;
    disp(filename_vis)
    % Next read in this .csv file as a table for visible image
    table_vis = readtable(strcat(dirin,'Metadata/',filename_vis));
    image_vis = imread(strcat(dirin,files_raw(i).name)); %save images
    % Do the same for thermal
    table_therm = readtable(strcat(dirin,'Metadata/',filename_therm));
    
    %Do thermal data only for odd (thermal) images
    data_encoded = table_therm.ThermalData;
    calibration_encoded = table_therm.ThermalCalibration;
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
    image_therm = reshape(data_corrected,640,512);
    % Orient properly
    image_therm = rot90(fliplr(image_therm));
    
    %Undistort images
    if undistort
        [image_vis,~] = undistortImage(image_vis,cameraParams_vis);
        [image_therm,~] = undistortImage(image_therm,cameraParams_therm);
    end
    
    % Note their sizes
    [h_vis,w_vis] = size(image_vis(:,:,1));
    [h_therm,w_therm] = size(image_therm);
    
    %Calculate spatial resolution
    %If there's RTK, use RTK, otherwise, use relative altitude
    try
        rtk = table_vis.RtkFlag;
        altitude_vis = table_vis.AbsoluteAltitude + 35.019; %subtract geoid
        altitude_therm = table_therm.AbsoluteAltitude + 35.019;
    catch
        altitude_vis = table_vis.RelativeAltitude;
        altitude_therm = table_therm.RelativeAltitude;
    end
    %Modify altitude for slight roll and pitch variations
    roll_vis = (table_vis.GimbalRollDegree) + table_vis.FlightRollDegree;
    pitch_vis = table_vis.GimbalPitchDegree+90 + table_vis.FlightPitchDegree;
    yaw_vis = table_vis.FlightYawDegree;
    altitude_vis_corrected = altitude_vis*sqrt(1+sqrt((sind(roll_vis)^2 + sind(pitch_vis)^2)));
    % Do the same for thermal
    roll_therm = (table_therm.GimbalRollDegree) + table_therm.FlightRollDegree;
    pitch_therm = (table_therm.GimbalPitchDegree+90) + table_therm.FlightPitchDegree;
    yaw_therm = table_therm.FlightYawDegree;
    altitude_therm_corrected = altitude_therm*sqrt(1+sqrt((sind(roll_therm)^2 + sind(pitch_therm)^2)));
    % Next, calculate the resolutions (m/pix) for each image
    res_therm = (12e-3)*(altitude_therm_corrected/9);
    res_vis = (6.17/w_vis)*(altitude_vis_corrected/4.5);
    
    % Next find Northings-Eastings coordinates for each pixel
    % Pull out center lat-lon data from metadata
    lat_deg_vis = table_vis.GPSLatitude;
    lon_deg_vis = table_vis.GPSLongitude;
    % Pull out center lat-lon data from metadata
    lat_deg_therm = table_therm.GPSLatitude;
    lon_deg_therm = table_therm.GPSLongitude;
    % Convert this to UTM
    [vis_loc_x,vis_loc_y] = ll2utm(lat_deg_vis,lon_deg_vis);
    [therm_loc_x,therm_loc_y] = ll2utm(lat_deg_therm,lon_deg_therm);
    
    % Using resolution, calculate utm coordinate for each pixel
    x_vis = ((0:w_vis-1)-w_vis/2)*(res_vis);
    y_vis = ((0:h_vis-1)-h_vis/2)*(res_vis);
    [X_vis,Y_vis] = meshgrid(x_vis,y_vis);
    % Rotate values and add center
    offset_x_vis = -altitude_vis_corrected*tand(roll_vis);
    offset_y_vis = -altitude_vis_corrected*(tand(pitch_vis));
    Eastings_vis = (X_vis+offset_x_vis).*cosd(yaw_vis) - (Y_vis+offset_y_vis).*sind(yaw_vis) + vis_loc_x;
    Northings_vis = -(Y_vis+offset_y_vis).*cosd(yaw_vis) - (X_vis+offset_x_vis).*sind(yaw_vis)+ vis_loc_y;
    % Do the same with thermal 
    x_therm = ((0:w_therm-1)-w_therm/2)*(res_therm);
    y_therm = ((0:h_therm-1)-h_therm/2)*(res_therm);
    [X_therm,Y_therm] = meshgrid(x_therm,y_therm);
    % Rotate values and add center (150 values determined empirically,
    % distance in visible pixels between center of thermal and visible image
    % respectively based on drone camera geometry)
    offset_x_therm = -altitude_therm_corrected*tand(roll_therm) +150*res_vis;
    offset_y_therm = -altitude_therm_corrected*(tand(pitch_therm)) - 150*res_vis;
    Eastings_therm = (X_therm+offset_x_therm).*cosd(yaw_therm) + (Y_therm+offset_y_therm).*sind(yaw_therm) + therm_loc_x;
    Northings_therm = (Y_therm+offset_y_therm).*cosd(yaw_therm) - (X_therm+offset_x_therm).*sind(yaw_therm)+ therm_loc_y;
    
    % Lastly, get approximate temperature values for image_therm
    image_therm_celsius = sensor_vals_to_temp(image_therm,0.98,altitude_therm_corrected,22);
    
    % Save therm as structure
    therm.filename = filename_therm(1:end-4);
    therm.time = datetime(table_therm.DateTimeOriginal,'InputFormat','yyyy:MM:dd HH:mm:ss');
    therm.image = image_therm;
    therm.image_celsius = image_therm_celsius;
    therm.altitude = altitude_therm_corrected;
    therm.table = table_therm;
    therm.calibration = calibration;
    therm.roll = roll_therm;
    therm.pitch = pitch_therm;
    therm.yaw = yaw_therm;
    therm.eastings = Eastings_therm;
    therm.northings = Northings_therm;
    therm.latitutde = lat_deg_therm;
    therm.longitude = lon_deg_therm;
    therm.resolution = res_therm;
    
    
    % Save vis as structure
    visible.filename = filename_vis(1:end-4);
    visible.time = datetime(table_vis.DateTimeOriginal,'InputFormat','yyyy:MM:dd HH:mm:ss');
    visible.image = image_vis;
    visible.altitude = altitude_vis_corrected;
    visible.table = table_vis;
    visible.roll = roll_vis;
    visible.pitch = pitch_vis;
    visible.yaw = yaw_vis;
    visible.eastings = Eastings_vis;
    visible.northings = Northings_vis;
    visible.latitude = lat_deg_vis;
    visible.longitude = lon_deg_vis;
    visible.resolution = res_vis;
    
    
    % Combine into single structure
    DJI_data.visible = visible;
    DJI_data.thermal = therm;
    DJI_data.source = dirin;
    
    % Name of files
    name = sprintf('%04d',counter);
    % Get folder name of dirin for naming
    slashes = strfind(dirin,'/');
    dirin_save = dirin(slashes(end-1)+1:end-1);
    
    
    if geotiff
        name2 = strcat(name,'_therm');
        % save into geotiff
        im_vis = imrotate(image_vis,-yaw_vis ,'loose');
        [lat_vis_min,lon_vis_min] = utm2ll(min(Eastings_vis(:)),min(Northings_vis(:)),11);
        [lat_vis_max,lon_vis_max] = utm2ll(max(Eastings_vis(:)),max(Northings_vis(:)),11);
        R_vis = georasterref('LatitudeLimits',[lat_vis_min lat_vis_max],...
            'LongitudeLimits',[lon_vis_min lon_vis_max],...
            'RasterSize',size(im_vis),'ColumnsStartFrom','north');
        geotiffwrite(strcat(dirout,'Geotiff/',dirin_save,'_processed_',name),im_vis,R_vis);
        % Same for thermal
        % First, convert thermal to color image
        % Regrid thermal image to [0,1]
        image_therm = image_therm-min(image_therm(:));
        image_therm = image_therm/(max(image_therm(:)));
        % Convert to color image
        image_therm = uint8(255*grs2rgb(image_therm,parula));
        % Then do same for thermal
        im_therm = imrotate(image_therm,-yaw_therm,'loose');
        [lat_therm_min,lon_therm_min] = utm2ll(min(Eastings_therm(:)),min(Northings_therm(:)),11);
        [lat_therm_max,lon_therm_max] = utm2ll(max(Eastings_therm(:)),max(Northings_therm(:)),11);
        R_therm = georasterref('LatitudeLimits',[lat_therm_min lat_therm_max],...
            'LongitudeLimits',[lon_therm_min,lon_therm_max],...
            'RasterSize',size(im_therm),'ColumnsStartFrom','north');
        geotiffwrite(strcat(dirout,'Geotiff/',dirin_save,'_processed_',name2),im_therm,R_therm);
    end
    
    
    save(strcat(dirout,dirin_save,'_processed_',name),'DJI_data');

end