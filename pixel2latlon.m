function coords = pixel2latlon(image,table)
% This function takes in an image with metadata from table and calculates a
% latitude and longitude position for each pixel, saved as coords
%
% Parameters:
% image - the image matrix
% table - the metadata table extracted by save_Thermal_Images
%
% Written by: Aidan Blaser (ablaser@ucsd.edu)
% Last Edited: 08/11/2021

%Constants for now
resolution = 0.0035; %meters per pixel
yaw = deg2rad(table.FlightYawDegree+90); %flight angle from due east
a = 6371000; %radius of earth in meters

%Pull out center lat-lon data from metadata
lat = table.GPSLatitude;
string_lat = erase(lat{1},['deg ',"'",'" N']);
lat_dms = sscanf(string_lat,'%f');
lat_deg = lat_dms(1) + (lat_dms(2)/60) + (lat_dms(3)/3600);

lon = table.GPSLongitude;
string_lon = erase(lon{1},['deg ',"'",'" N']);
lon_dms = sscanf(string_lon,'%f');
lon_deg = lon_dms(1) + (lon_dms(2)/60) + (lon_dms(3)/3600);

%Find size of image in pixels
[length,width,~] = size(image);
coords = zeros(length,width,2);

for i=1:width
    for j=1:length
        %Find vector from center to each point
        displacement = [i-width/2;length/2-j];
        %Convert vector from pixels to meters
        disp_meters = displacement*resolution;
        %Rotate vector
        rotated = [cos(yaw),sin(yaw);-sin(yaw),cos(yaw)]*disp_meters;
        %Convert into lat-lon
        coords(j,i,1) = -lon_deg + rad2deg(rotated(1)/(a*cos(deg2rad(lat_deg))));
        coords(j,i,2) = lat_deg + rad2deg(rotated(2)/a);
    end
end
        
        
        
        