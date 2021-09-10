%Script for superimposing visible and thermal image atop one another
format long

%Load camera distortion parameters
load('../Camera_Calibration/Camera_Calibration.mat');
cameraParams_therm = cameraParams;
load('../Camera_Calibration_Vis/camera_parameters.mat');

% First, load the visible image
vis = imread('../Ground_Control/DJI_0016.JPG');
vis = undistortImage(vis,cameraParams);
vis_undisturbed = vis;
load('/Users/aidanblaser/Desktop/DJI_Mavic2EA_IR_test/Ground_Extracted/DJI_0016.mat')
%load('/Volumes/LaCie/TempTest/DJI_0116.mat')
%vis = image;
%vis = undistortImage(vis,cameraParams);
%vis_undisturbed = vis;

% Parameters for visible image
% Pull out center lat-lon data from metadata
lat = table.GPSLatitude;
string_lat = erase(lat{1},['deg ',"'",'" N']);
lat_dms = sscanf(string_lat,'%f');
lat_deg_vis = lat_dms(1) + (lat_dms(2)/60) + (lat_dms(3)/3600);

lon = table.GPSLongitude;
string_lon = erase(lon{1},['deg ',"'",'" N']);
lon_dms = sscanf(string_lon,'%f');
lon_deg_vis = lon_dms(1) + (lon_dms(2)/60) + (lon_dms(3)/3600);

angle_vis = table.FlightYawDegree;
roll_vis = table.FlightRollDegree;
pitch_vis = table.FlightPitchDegree;

% Load the thermal
load('/Users/aidanblaser/Desktop/DJI_Mavic2EA_IR_test/Ground_Extracted/DJI_0017.mat')
%load('/Volumes/LaCie/TempTest/DJI_0117.mat')
% Parameters
yaw = table.FlightYawDegree;
roll_therm = table.FlightRollDegree;
pitch_therm = table.FlightPitchDegree;
% Pull out center lat-lon data from metadata
lat = table.GPSLatitude;
string_lat = erase(lat{1},['deg ',"'",'" N']);
lat_dms = sscanf(string_lat,'%f');
lat_deg_therm = lat_dms(1) + (lat_dms(2)/60) + (lat_dms(3)/3600);

lon = table.GPSLongitude;
string_lon = erase(lon{1},['deg ',"'",'" N']);
lon_dms = sscanf(string_lon,'%f');
lon_deg_therm = lon_dms(1) + (lon_dms(2)/60) + (lon_dms(3)/3600);

altitude = table.RelativeAltitude;
angle_diff = angle_vis - table.FlightYawDegree;


% Undistort the image
[image,~] = undistortImage(image,cameraParams_therm);


% Note their sizes
[h_vis,w_vis] = size(vis(:,:,1));
[h_therm,w_therm] = size(image);

% Next, calculate the resolutions (m/pix) for each image
res_vis = (6.17/w_vis)*(altitude/4.5);
%res_therm = (12.56e-3)*(altitude/9);
res_therm = (12.15e-3)*(altitude/9);

% Calculate the thermal extents in meters
h_extent = h_therm * res_therm;
w_extent = w_therm * res_therm;

% Convert this to pixels in visible image (half extent)
hpix = round(h_extent / (2*res_vis));
wpix = round(w_extent / (2*res_vis));

% Resize thermal image to match these values
im = imresize(image,[2*hpix,2*wpix]);

% Remap thermal image to [0,1]
im = im-min(im(:));
im = im/(max(im(:)));

% Convert to color image
therm_insert = grs2rgb(im,parula);

% Calculate difference spatially between two images
[vis_loc_x,vis_loc_y] = ll2utm(lat_deg_vis,-lon_deg_vis);
[therm_loc_x,therm_loc_y] = ll2utm(lat_deg_therm,-lon_deg_therm);
dist_m = sqrt((vis_loc_x-therm_loc_x)^2 + (vis_loc_y-therm_loc_y)^2); %maybe divide by 2, maybe not, check later
dist_pix = (dist_m / (1*res_vis)); %maybe divide by 2, maybe not, check later
angle = atan((vis_loc_y-therm_loc_y)/(vis_loc_x-therm_loc_x)); %angle between locations
if vis_loc_x < therm_loc_x %want to find angle from thermal to vis
        angle = angle + pi;
end
if isnan(angle) %in case they're at same location
    angle = 0;
end


%Find how far to modify centerpoint
dist_pix_x = round(dist_pix*cos(angle-yaw+deg2rad(angle_diff)));
dist_pix_y = round(dist_pix*sin(angle-yaw+deg2rad(angle_diff)));

% Roll and Pitch
dist_pix_x = round(dist_pix_x + altitude*(tan(deg2rad(roll_therm))-tan(deg2rad(roll_vis)))/res_vis);
dist_pix_y = round(dist_pix_y + altitude*(tan(deg2rad(pitch_therm))-tan(deg2rad(pitch_vis)))/res_vis);

% Account for discrepency in location between visible and thermal cameras
dist_pix_x = dist_pix_x + 150;
dist_pix_y = dist_pix_y - 150;

%Rotate visible image
vis = imrotate(vis,-angle_diff);

% Insert thermal into visible image
vis(h_vis/2 - hpix+1-dist_pix_y:h_vis/2 + hpix - dist_pix_y ,w_vis/2 - wpix+1+dist_pix_x:w_vis/2 + wpix+dist_pix_x,:) = 255*therm_insert;

% Rotate back
vis = imrotate(vis,angle_diff);

imshow(vis)