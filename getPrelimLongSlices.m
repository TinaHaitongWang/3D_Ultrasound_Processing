function [longSlices, toPSmaps] = getPrelimLongSlices(volume,numLS,outStepSize)
%getPrelimLongSlices: Helper function for findBaseAxisPoints.m
% Generates vertical slices rotated around the z axis
% that are used to get user input of base point and axis direction
% ***Arguments***
% volume: dcm volume to be used (get with readDICOM3D.m)
% numLS: number of rotated long slices
% outStepSize: cm/pixel - this argument is specified by
% findBaseAxisPoints.m
% ***Outputs***
% longSlices: n x n x numLS 3D matrix of long slices
% toPSmaps: 1 x numLS cell array consisting of maps that transform user input back to
% physical space

% Note: this is a subfunction - you do not need to use it independently!
%% Prep work
% Extract dimensions of physical space (in cm)
width = volume.widthspan;
height = volume.heightspan;
depth = volume.depthspan;
% Decide size (# of pixels) of output and initialize output slices
outSize = round(max(norm([width,height],2),depth)/outStepSize + 1);
%extract frame and initialize items
inVol= volume.data(:,:,:,1);
longSlices = zeros(outSize,outSize,numLS);
da = 180/numLS;
basePoint = [width/2,height/2,depth/2]';
toPSmaps = cell(1,numLS);
%% Map from output space to physical space 
XYZmap1 = @(x) (x-1)*((outSize-1)*outStepSize)/(outSize-1)- ((outSize-1)*outStepSize/2);

%% Map from physical space to input space
Win = size(inVol,1);
Hin = size(inVol,2);
Din = size(inVol,3);
Xmap2 = @(x) x*(Win-1)/width + 1;
Ymap2 = @(y) y*(Hin-1)/height + 1;
Zmap2 = @(z) z*(Din-1)/depth + 1;
%% Fill in long slices
% Start filling long slices
for u = 1:numLS
     % Generate rotation matrix
            a = (u-1)*da;
            Mrot = [cosd(a),-sind(a),0;
                    sind(a),cosd(a),0;
                    0,       0,     1];
            % Starting point is [1 0 0]
            % First, rotate around z axis
            % Second, rotate to goal coordinate system
            normVec = Mrot*[1 0 0]';
            thisNV = normVec';
            % Find rotation angle
            thisTheta = acosd(thisNV(3)/norm(thisNV,2));
            thisPhi = atan2d(thisNV(2),thisNV(1));
            % rotation of +thisTheta degrees around y (2nd) axis
            thisA = [cosd(thisTheta),0,sind(thisTheta);
                     0,1,0;
                   -sind(thisTheta),0,cosd(thisTheta)];
            % rotation of +thisPhi degrees around z (3rd) axis
            thisB = [cosd(thisPhi),-sind(thisPhi),0;
                     sind(thisPhi),cosd(thisPhi),0;
                     0,0,1];
            %%% ~~~This section stores the maps from output to physical space
            thisMap = @(RC) thisB*thisA*[XYZmap1(RC(1)),XYZmap1(RC(2)),0]' + basePoint;
            toPSmaps{u} = thisMap;           
            %%% In order to let findBaseAxisPoints.m output in the correct
            %%% coordinates
    for r = 1:outSize
        for c = 1:outSize
            % Map to input volume
            pos_rotated = [XYZmap1(r) XYZmap1(c) 0]';
            pos_original = thisB*thisA*pos_rotated  + basePoint; % rotate, then translate
            
            w = Xmap2(pos_original(1));
            h = Ymap2(pos_original(2));
            d = Zmap2(pos_original(3));
        
            if floor(w) > Win || ceil(w) < 1 || ...
               floor(h) > Hin || ceil(h) < 1 || ...
               floor(d) > Din || ceil(d) < 1
                    longSlices(r,c,u) = 0;
                continue
            end
        
            % Interpolate
            % Generate upper and lower bounds
            x2 = min(ceil (w),Win);
            x1 = max(x2-1,1);

            y2 = min(ceil(h),Hin);
            y1 = max(y2-1,1);
            
            z2 = min(ceil(d),Din);
            z1 = max(z2-1,1);
        
            xd = (w-x1);
            yd = (h-y1);
            zd = (d-z1);
            
            %x- interpolation
            V00 = inVol(x1,y1,z1)*(1-xd) + inVol(x2,y1,z1)*xd;
            V01 = inVol(x1,y1,z2)*(1-xd) + inVol(x2,y1,z2)*xd;
            V10 = inVol(x1,y2,z1)*(1-xd) + inVol(x2,y2,z1)*xd;
            v11 = inVol(x1,y2,z2)*(1-xd) + inVol(x2,y2,z2)*xd;

            %y - interpolation
            V0 = V00*(1-yd) + V10*yd;
            V1 = V01*(1-yd) + v11*yd;
            
            %z interpolation
            longSlices(r,c,u) = V0*(1-zd) + V1*zd;
        end
    end
end
end
