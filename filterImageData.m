function outVolume = filterImageData(inVolume)
%% filterImageData function
% This function allows user to pick type of filter for process the image
% data. 
% inVolme: volume data from readDicom3D (struct)
% outVolume: volume data (struct) with extra field for fillterd data 
% Apply depth averaging filter and remove speckles filte
% along the depth direction 
% Note: this code is adapted from Andy's resource, section 1,2,3: Introduction.
% ,Filtering, and Edge Detection. 
% Reference: Andrew J Wald, Project Tutorial: Finding The Mitral Annulus, 
% BME 543, DUKE BME, http://people.duke.edu/~jcm21/bme265/andy/MitralAnnulus_Tut3.pdf

% Check number of inputs 
if nargin < 1
    error('Not enough input argument');
end 
%% depth averaging filter
% this filter use the blackman windows to filter through the depth (z)
% direction throuhg each slice (width * height)
N = 21; % emprically chosen, can be refined to optimize the filtering 
for i = 1: inVolume.NumVolumes 
    filter = zeros(1,1,N); 
    filter(1,1,:) = blackman(N); % this filter uses a balckman window 
    inVolume.dataDepthfiltered(:,:,:,i) = convn(inVolume.data(:,:,:,i),filter,'same'); 
end 
%% normalized data 
% this function normalize the data x 
normData = @(x) (x-min(min(x)))/(max(max(x))-min(min(x)));
%% filter remove speckle in short axis slices
% this filter us wiener window and second order statistic filter to remove
% speckle in the short axis slice 
inVolume.dataShortAxisFiltered = inVolume.dataDepthfiltered;
for i = 1: inVolume.NumVolumes 
    for j = 1:inVolume.depth
        s = squeeze(inVolume.dataShortAxisFiltered(:,:,j,i)); 
        s = normData(s);
        s = wiener2(s,[7 7]); % apply wiener filter 
        s = normData(s);
        s = imadjust(s, [0.2 0.7],[0 1]);% adjust the intensity to between 0 and 1 
        s = ordfilt2(s, 14, ones(5,5),'symmetric'); % 2d order statistic filter
        inVolume.dataShortAxisFiltered(:,:,j,i)  = normData(s);
    end 
end
%% output volume 
outVolume = inVolume; 
end 
