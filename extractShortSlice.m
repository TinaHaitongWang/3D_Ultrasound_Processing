function [shortSlice] = extractShortSlice(volume,subdir,frame,basePoint,normVec,...
                                                        dist, outStepSize, outSize)
%extractAnatSlices Extract isotropic long and short axis slices from ultrasound
%                    volume, which is not necessarily isotropic.
% volume: input volume struct obtained from readDicom3D
% basePoint: location of long axis base, as (x,y,z) in PHYSICAL SPACE
% normVect: normal vector along long axis of volume in PHYSICAL SPACE
% LSangle = angle of rotation about long axis for long slice
% dist: distance of translation (in cm) along long axis for short slice
% frame: specific frame to be used 
% outStepSize(in units of cm/pixel): physical distance per pixel in output
% subdir: type of volume data to slice, filtered or non-filtered 

%% Prep work
%extract image and dimensions of physical space
width = volume.widthspan;
height = volume.heightspan;
depth = volume.depthspan;
vol = getfield(volume,subdir); 
inVol= vol(:,:,:,frame);

%initialize slice
shortSlice = zeros(outSize,outSize);

%% Forward rotation maps in physical space
% Find theta and phi from normal vector
x0 = normVec(1);
y0 = normVec(2);
z0 = normVec(3);
theta = acosd(z0/norm(normVec,2));
phi = atan2d(y0,x0);
% rotation of +theta degrees around y (2nd) axis
A = [cosd(theta),0,sind(theta);
      0,1,0;
     -sind(theta),0,cosd(theta)];
% rotation of +phi degrees around z (3rd) axis
B = [cosd(phi),-sind(phi),0;
     sind(phi),cosd(phi),0;
     0,0,1];
    
% Find normal vector AFTER the two rotations, in physical space
nv = B*A*[0 0 1.0]';
  
%% Map from output space to physical space 
XYZmap1 = @(x) (x-1)*(outSize*outStepSize)/(outSize-1)- (outSize*outStepSize/2);
%% Map from physical space to input space
Win = size(inVol,1);
Hin = size(inVol,2);
Din = size(inVol,3);
Xmap2 = @(x) x*(Win-1)/width + 1;
Ymap2 = @(y) y*(Hin-1)/height + 1;
Zmap2 = @(z) z*(Din-1)/depth + 1;
%% Fill in short slices
    for r = 1:outSize
        for c = 1:outSize
            % Map to input volume
            pos_rotated = [XYZmap1(r) XYZmap1(c) 0]';
            pos_original = B*A*pos_rotated + (dist*nv) + basePoint; % rotate, then translate
            w = Xmap2(pos_original(1));
            
            h = Ymap2(pos_original(2));
            d = Zmap2(pos_original(3));
        
            if floor(w) > Win || ceil(w) < 1 || ...
               floor(h) > Hin || ceil(h) < 1 || ...
               floor(d) > Din || ceil(d) < 1
                    shortSlice(r,c) = 0;
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
            shortSlice(r,c) = V0*(1-zd) + V1*zd;
            
        end
    end

end

