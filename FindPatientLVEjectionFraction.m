function [VolumeLV, ejectionFraction] = FindPatientLVEjectionFraction(...
    PatientFileName, thickness, outStepSize,sliceSpace)
% Cardiac Ultrasound Imaging and Function Final project 
% Team: Alex Noboa, Gehua Tong, Haitong Wang 

% This is the main function for this project 
% This function allows user to input the patient file name (.dicom) and
% compute the volume of left ventricle and its ejetion graction with
% specified pixel size (of each slice of 4D data) and thickness of slice. 

% PatientFileName: string, indicate the name of patient data file 
% thickness : slice thickness in cm, and slice are taken at each point of
% the normal vector from the base of the heart to the apex of the heart
% outStepSize: determine the physical size of each pixel in the slice (cm)
%% User Input
% In this section, user needs to input the patient data that they wish to
% compute and define the thickness and output pixel size in cm for the
% short slice image 
volume = readDicom3D(PatientFileName);
volume.data = single(volume.data);
%% Find the base point and axis orientation of the volume 
% This function allows us to define the base of the heart and the normal
% vector that go through the parastrenal long axis 
[base,apex, vec] = findBaseAxisPoints(volume); 

%% Pre work 
basePT = base'; 
nv = vec'; % normal vector 
diagonal = norm([volume.widthspan volume.heightspan volume.depthspan],2);
maxoutSize = round(diagonal/outStepSize);
% estimated ventricle length 
Vlength = sqrt((base(1)-apex(1))^2+(base(2)-apex(2))^2+(base(3)-apex(3))^2);
LAlength = floor(Vlength/thickness)*thickness; %length of long axis
slicenumber = LAlength/thickness; % number of slice 
stepnumber = thickness/outStepSize;% number step in each slice 
width = volume.widthspan/outStepSize; 
height = volume.heightspan/outStepSize;

% filter the volumetric data for later edge detection 
volume = filterImageData(volume); 

% define the ratio of research radius for the edge of left ventricle 
% in each frame 
medianFrame = floor(median(1:volume.NumVolumes));
ratioFrame1 = linspace(0.6,1,medianFrame);
if mod(volume.NumVolumes,2)
    ratioFrame = [fliplr(ratioFrame1) ratioFrame1(2:end)];
else 
    ratioFrame = [fliplr(ratioFrame1) ratioFrame1];
end 
% define the ratio to adjust for each slice 
ratioSlice = fliplr(linspace(0.7,1.2,slicenumber)); 

%% Generate Slice and Calculate Volume 
% uncomment count to track the iteration 
% count = 1
for i = 1:volume.NumVolumes
    frame = i; 
    % for regular slice 
    subdir = 'dataShortAxisFiltered';
    % define the basic research radius for the slice 
    baseSearchRadius = round((width+height)/15)*ratioFrame(i);
    for j = 1:slicenumber
        basedist = (j-1)*thickness;
        outSize = round(sqrt(1-(basedist/LAlength))*maxoutSize);
        sumimage = zeros(outSize,outSize);
        % for each step, extract a slice and add them up 
        for k = 1:stepnumber
            dist = basedist + ((k-1)*sliceSpace);
            shortSlice = extractShortSlice(volume,subdir,frame,basePT, nv, ...
                dist, outStepSize,outSize);
            sumimage = sumimage + shortSlice; 
        end
        % find the average slice over number of steps to generate the thick
        % slice 
        avgimage{i,j} = sumimage/(thickness/sliceSpace);
        % define the Slice of interest, center of the slice and the search
        % radius of the slice 
        Slice = avgimage{i,j}; center = [size(Slice,1)/2,size(Slice,2)/2];
        searchRadius = baseSearchRadius*ratioSlice(j);
        % Find the edge of the left ventricle 
        [PointEdge,Area,DiskVolume] = findVentricleEdgeArea(Slice,center,...
            searchRadius,outStepSize,thickness);
        PointsLV{i,j} = PointEdge;
        SliceArea(i,j) = Area; 
        TotalVolume(i,j) = DiskVolume; 
%       count = count+1
    end
    
    % for the remaining slice, with different thickness, redo previous
    % calculation 
    n = round((Vlength-LAlength)/sliceSpace);
    if n >0
        remthickness = n * sliceSpace;
        outSize = round(sqrt(1-(LAlength/Vlength))*maxoutSize);
        sumimage = zeros(outSize,outSize);
        % for each step, extract a slice and add them up 
        for k = 1:n+1
            dist = LAlength + ((k-1)*sliceSpace);
            shortSlice = extractShortSlice(volume,subdir,frame,basePT, nv, ...
                dist, outStepSize,outSize);
            sumimage = sumimage + shortSlice;
        end
        % find the average slice over number of steps to generate the thick
        % slice 
        avgimage{i,slicenumber+1} = sumimage/(thickness/sliceSpace);
        Slice = avgimage{i,slicenumber+1}; 
        searchRadius = searchRadius* 0.9;
        % Find the edge of the left ventricle 
        [PointEdge,Area,DiskVolume] = findVentricleEdgeArea(Slice,center,...
            searchRadius,outStepSize,remthickness);
        PointsLV{i,slicenumber+1} = PointEdge;
        SliceArea(i,slicenumber+1) = Area; 
        TotalVolume(i,slicenumber+1) = DiskVolume; 
%       count = count+1
    end
    
end 
%% Visualization of the volume change over frames 
for i = 1:volume.NumVolumes
   y(i) = sum(TotalVolume(i,:));
end 
VolumeLV = y; 
plot(y)
xlabel('Frame')
ylabel('Left Ventricle Volume (ml)')
title(PatientFileName)
%% Calculate the ejectionFraction 
ejectionFraction = (max(y)-min(y))/max(y); 

end 