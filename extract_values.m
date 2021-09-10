% Script for pulling out Thermal Data from image metadata
%
% Written by: Aidan Blaser
% Last Edited: 08/09/2021

% First specify .csv filename from EXIFTOOL from shell
filename = '../test.csv'; 
% Next read in this .csv file as a table
A = readtable(filename); 
% Pull out relevant data (ThermalData)
data_encoded = A.ThermalData;
% Decode this data from base64 to numbers
% Note: Start from position 8 to remove the 'base64:'
data = matlab.net.base64decode(data_encoded{1}(8:end));
% Note: This data should have twice the number of bits,
%       correct for this by the following formula
a = double(data(1:2:end)); %odd values
b = double(data(2:2:end)); %even values
data_corrected = b*256 + a;
% Resize this image according to its mxn size
image = reshape(data_corrected,640,512);
% Orient properly
image = rot90(fliplr(image));
% Plot data and add colorbar
imagesc(image)
colorbar
