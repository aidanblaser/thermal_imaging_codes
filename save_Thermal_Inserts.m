function save_Thermal_Inserts(dirin,dirout)
% This function takes the thermal data and superimposes it on the visual
% data
%
% Parameters:
% dirin - the path to the directory containing the image files
% dirout - the name and location of the directory where images will be
%           saved
%
% Written by: Aidan Blaser (ablaser@ucsd.edu)
% Last Edited: 09/02/2021

% First, create the outgoing directory
mkdir(dirout);

% Extract matlab fil
files_images = dir(strcat(dirin,'*.mat'));

%Load camera distortion parameters
load('../Camera_Calibration/Camera_Calibration.mat','cameraParams');
cameraParams_therm = cameraParams;
load('../Camera_Calibration_Vis/camera_parameters.mat','cameraParams');

for i=101:2:length(files_images)
    disp(files_images(i).name)
    
    %Script for superimposing visible and thermal image atop one another
    format long
    %First, load the visible image
    %vis = imread(strcat(dirin,files_images(i).name));
    load(strcat(dirin,files_images(i).name),'image','table')
    vis = image;
    [vis,~] = undistortImage(vis,cameraParams);
    vis_undisturbed = vis;

    % Parameters for visible image
    % Pull out center lat-lon data from metadata
    lat_deg_vis = table.GPSLatitude;
    lon_deg_vis = table.GPSLongitude;
    % Pull out roll, pitch, yaw
    angle_vis = table.FlightYawDegree;
    roll_vis = table.GimbalRollDegree;
    pitch_vis = table.GimbalPitchDegree+90;

    % Load the thermal
    load(strcat(dirin,files_images(i+1).name),'image','table')
    % Parameters
    roll_therm = table.GimbalRollDegree;
    pitch_therm = table.GimbalPitchDegree+90;
    % Pull out center lat-lon data from metadata
    lat_deg_therm = table.GPSLatitude;
    lon_deg_therm = table.GPSLongitude;
    %altitude = table.RelativeAltitude;
    altitude = table.AbsoluteAltitude + 33.69; %subtract geoid
    yaw = table.FlightYawDegree;
    angle_diff = angle_vis - yaw;

    % Undistort the image
    [image,~] = undistortImage(image,cameraParams_therm);


    % Note their sizes
    [h_vis,w_vis] = size(vis(:,:,1));
    [h_therm,w_therm] = size(image);
    
    % Correct altitude for pitch and roll
    altitude = altitude*sqrt(1+sqrt((cos(roll)^2 + cos(pitch)^2)));

    % Next, calculate the resolutions (m/pix) for each image
    res_vis = (6.17/w_vis)*(altitude/4.5);
    res_therm = (12.15e-3)*(altitude/9);

    % Calculate the thermal extents in meters
    h_extent = h_therm * res_therm;
    w_extent = w_therm * res_therm;

    % Convert this to pixels in visible image (half extent)
    hpix = round(h_extent / (2*res_vis));
    wpix = round(w_extent / (2*res_vis));

    % Resize thermal image to match these values
    im = imresize(image,[2*hpix,2*wpix]);

    % Regrid thermal image to [0,1]
    im = im-min(im(:));
    im = im/(max(im(:)));

    % Convert to color image
    therm_insert = grs2rgb(im,parula);

    % Calculate difference spatially between two images
    [vis_loc_x,vis_loc_y] = ll2utm(lat_deg_vis,-lon_deg_vis);
    [therm_loc_x,therm_loc_y] = ll2utm(lat_deg_therm,-lon_deg_therm);
    dist_m = sqrt((vis_loc_x-therm_loc_x)^2 + (vis_loc_y-therm_loc_y)^2); 
    dist_pix = (dist_m / (res_vis)); 
    angle = atan((vis_loc_y-therm_loc_y)/(vis_loc_x-therm_loc_x)); %angle between locations
    if vis_loc_x < therm_loc_x %want to find angle from thermal to vis
        angle = angle + pi;
    end
    if isnan(angle) %in case they're at same location
        angle = 0;
    end

    %Find how far to modify centerpoint
    dist_pix_x = round(dist_pix*cos(angle+yaw+deg2rad(angle_diff)));
    dist_pix_y = round(dist_pix*sin(angle+yaw+deg2rad(angle_diff)));
    
    % Roll and Pitch Adjustments
    dist_pix_x = round(dist_pix_x - altitude*(tan(deg2rad(roll_therm))-tan(deg2rad(roll_vis)))/res_vis);
    dist_pix_y = round(dist_pix_y - altitude*(tan(deg2rad(pitch_therm))-tan(deg2rad(pitch_vis)))/res_vis);
    
    % Account for discrepency in location between visible and thermal cameras
    dist_pix_x = dist_pix_x + 140;
    dist_pix_y = dist_pix_y - 160;

    %Rotate visible image
    vis = imrotate(vis,-angle_diff);
    %[hnew,wnew] = size(vis); %new size of rotated image
    

    % Insert thermal into visible image
    vis(h_vis/2 - hpix+1-dist_pix_y:h_vis/2 + hpix - dist_pix_y ,w_vis/2 - wpix+1+dist_pix_x:w_vis/2 + wpix+dist_pix_x,:) = 255*therm_insert;

    % Rotate back
    vis = imrotate(vis,angle_diff);
    imshow(vis)

    save(strcat(dirout,files_images(i).name(1:end-3),'mat'),'vis','vis_undisturbed')
end
end