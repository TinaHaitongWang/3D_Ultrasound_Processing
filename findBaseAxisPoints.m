function [basePoint,apexPoint, axisNormVec] = findBaseAxisPoints(volume)
%findBaseAxisPoints: asks user for cursor inputs to find base point and
%axis orientation and outputs them in physical coordinates and units
% Reminder: the corner of the volume (pixel(1,1,1))is (0,0,0) in physical space 
% Reminder 2: unit of distance is (cm)
% How the function works: ij

%% Extract slices
numLS = 6; %% number of long slices to generate
stepSize = 0.1; %% step size for displayed long slices
[prelimLS,myMaps] = getPrelimLongSlices(volume,numLS,stepSize);
%% Show slices to user and accept input
x = zeros(2*numLS,1);
y = zeros(2*numLS,1);
for k = 1:numLS
    figure(k)
    imagesc(prelimLS(:,:,k))
    title('Specify base, then apex')
    colormap gray
    [x(2*k-1:2*k),y(2*k-1:2*k)] = ginputY(2);
    hold on
    plot(x(2*k-1:2*k),y(2*k-1:2*k),'-*y')
    pause(1)
end
%% Find base and apex that maximize their Cartesian distance
RCs = [y,x];
maxDist = 0;
Base = [];
Apex = [];
bestg = 1;
for g = 1:numLS
    base = RCs(2*g-1,:);
    apex = RCs(2*g,:);
    if norm(base-apex) > maxDist
       Base = base;
       Apex = apex;
       maxDist = norm(base-apex);
       bestg = g;
    end
end
%% Transform data back to physical domain
myMap = myMaps{bestg};
basePoint = myMap(Base);
apexPoint = myMap(Apex);
axisNormVec = (apexPoint - basePoint)/norm(apexPoint-basePoint);
close all
end

