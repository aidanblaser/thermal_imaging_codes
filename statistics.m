% Script for doing statistical analysis on thermal image

% Load image
close all
name = '0417';
load(['../TempTest/DJI_',name,'.mat'],'image','res','altitude_corrected')

%Rotate image so that waves are vertical
angle = 95; %degrees
image = imrotate(image,angle,'crop');
%Make all 0's from rotate into NaN
t = image==0;
image(t) = NaN;

[height,width] = size(image);

signal = zeros(1,length(image));
% Do some statistics on strips
for i=1:length(image)
    signal(i) = mean(image(:,i),'omitnan');
end

signal_raw = signal;

%Demean signal
signal = signal - mean(signal(:),'omitnan');
%No negative values
signal = signal + abs(min(signal(:)));

%Choosable pixel parameters
fit_start = 450;
fit_end = 500;

x = (0:width-1)*res;
y = 0:height-1 * res;

lny = log(signal(fit_start:fit_end));
p = polyfit(x(fit_start:fit_end),lny,1);
fit = exp(p(2))*exp(p(1)*x(fit_start:fit_end));

figure1 = figure('Position',[100,100,1024,1200]);

sb1 = subplot(2,1,1);
sb1.Position(4) = sb1.Position(4) + 0.05;
imagesc(x,y,image)
ylabel('distance (m)','FontSize',18)

sb2 = subplot(2,1,2); 
sb2.Position(4) = sb2.Position(4) + 0.05;
plot(x,signal)
xlim([0,max(x)])
title(['Image: ',name],'FontSize',18)
hold on
plot(x(fit_start:fit_end),fit)
xlabel('distance from wave (m)','FontSize',18)
ylabel('Y-averaged thermal measurement (arb. units)','FontSize',18)