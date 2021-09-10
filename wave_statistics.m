function wave_statistics(file)
% Computes wave statistics on a thermal image specified by file
%
% Written by: Aidan Blaser (ablaser@ucsd.edu)
% Last Edited: 09/07/2021

close all
%First, load the file
load(file,'image','res')

[~,width] = size(image);

%Vertically average signal
signal = zeros(1,width);
for i=1:width
    signal(i) = mean(image(:,i));
end

%Demean signal
signal = signal - mean(signal(:));
%No negative values
signal = signal + abs(min(signal(:)));

%Generate x
x = (0:width-1)*res;

%Visualize signal
plot(x,signal)
xlabel('distance (m)','FontSize',18)
ylabel('Y-averaged thermal measurement (arb. units)','FontSize',18)
figure(2)
imagesc(image)