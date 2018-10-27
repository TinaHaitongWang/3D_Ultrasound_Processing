function [PointEdge,Area,DiskVolume] = findVentricleEdgeArea(Slice,center,...
                                           searchRadius,stepSize,thickness)
%% ventricle edge detection function 
% this function uses Andy's resource for applying different combinations of 
% gaussian blur and sobel filter, then conduct a radial search to detect
% the edge 
% slice : short axis slice 
% center: center of the left ventricle location in the image [x,y] or [row, column]
% searchRadius : maximum distance to look for edge from the center of
% ventricle  eg.round((width+height)/15); 
% stepSize: size of a pixel in physical dimension 

% output: 
% Points: location in x, y in image domain of edge of ventricle wall
% Area : area of left ventricle approximated by those points in physical
% dimension 

% Note: this code is adapted from Andy's resource, section 3: Edge
% Detection. 
% Reference: Andrew J Wald, Project Tutorial: Finding The Mitral Annulus
% Part3: Feature Detection, BME 543, DUKE BME, http://people.duke.edu/~jcm21/bme265/andy/MitralAnnulus_Tut3.pdf

%% Initialize parameters for search 
[width, height] = size(Slice);
numAngles = 60 ; % number of lines drawn at equally spaced angles around a 
% with an origin at the center 
% Number of combinations of gaussian blur and sobel filter for detecting
% the edge 
numBlurVariants = 6; % determine the degree of blurring 
numEdgeDetectionVariants = numBlurVariants; 
% initialize the edge points 
PointEdge = zeros(2,numEdgeDetectionVariants,numAngles); 
% parameters for guassian blurs 
dsigma = 0.003*((width+height)/2);
Slice = squeeze(Slice); 
onEdgeDetectionVariant = 1; 

%% Edge detections by trying different blur filters 
for blurIndex = 1: numBlurVariants
    % apply gaussian blur filter 
    if(blurIndex == 1)
        sigma = 0; 
        sliceFilteredGB = Slice; 
    else 
        sigma = (blurIndex-1)*dsigma; 
        filter = fspecial('gaussian',ceil(8*sigma),sigma);
        sliceFilteredGB = conv2(Slice,filter,'same'); 
    end 
    % apply sobel filter 
    horizontalFilter = conv2(sliceFilteredGB,fspecial('sobel'),'same');
    verticalFilter = conv2(sliceFilteredGB,fspecial('sobel')','same'); 
    sliceFiltered = sqrt(horizontalFilter.^2+ verticalFilter.^2); 
    
    % use the filtered slice to compare with the original slice to conduct
    % a radial search from the center of slice to define the edge of the
    % left ventricle 
    centerColn = round(center(2)); 
    centerRow = round(center(1)); 
    dtheta = 2*pi/numAngles; 
    % initial an auxiliary slice that has the same dimension as the Slice
    % to prevent picking points that co-occur on downhills as we radially
    % search out along each line from the center 
    AuxSlice = sliceFiltered; 
    AuxSliceMin = min(min(sliceFiltered)); AuxSliceMax = max(max(sliceFiltered));
    AuxSliceRange = AuxSliceMax -AuxSliceMin; 
    AuxSliceThreshold = AuxSliceMin + 0.62 * AuxSliceRange;
    M = searchRadius; % default to ensure all pixels are included in the search 
    
    % search for each radial line 
    for thetaIndex = 1:numAngles 
        theta = thetaIndex * dtheta; 
        MaximumLine = 0; % maximum along each line 
        AuxSlice_IntensityAtFirstLastPos = AuxSlice(centerRow,centerColn);
        AuxSlice_IntensityAtSecondLastPos = AuxSlice(centerRow,centerColn);
        
        % search at each radius along the line 
        for radius = linspace(1,searchRadius,M); 
            x1 = round( centerColn + (radius*cos(theta)) );
            y1 = round( centerRow + (radius*sin(theta)) );
            
            if ((x1>height)||(x1<1)||(y1>width)||(y1<1)) 
                % when the location is outside image 
                break; 
            elseif ((Slice(y1,x1)>MaximumLine)&&...
                    (AuxSlice(y1,x1)>AuxSlice_IntensityAtFirstLastPos)&&...
                    (AuxSlice(y1,x1)>AuxSlice_IntensityAtSecondLastPos)&&...
                    (AuxSlice(y1,x1)>AuxSliceMin)&&...
                    (AuxSlice(y1,x1)<AuxSliceThreshold))
                % when the slice has intensity is larger than that in the
                % auxiliary slice 
                % update the comparing parameters 
                MaximumLine = Slice(y1,x1);
                points(1,thetaIndex) = y1; % row
                points(2,thetaIndex) = x1; % column
            end 
            
            % update the auxilary slice 
            AuxSlice_IntensityAtSecondLastPos = AuxSlice_IntensityAtFirstLastPos;
            AuxSlice_IntensityAtFirstLastPos = AuxSlice(y1,x1); 
        end 
        
        % for the case, there is not a maximum for the line, assign the
        % point to be outside the image 
        if (MaximumLine ==0)
            points(:,thetaIndex) = [-1;-1];
        end 
    end 
       % record the points on the maximum points on each radial line 
       PointEdge(1,onEdgeDetectionVariant,:) = points(1,:); % row
       PointEdge(2,onEdgeDetectionVariant,:) = points(2,:); % column
       onEdgeDetectionVariant = onEdgeDetectionVariant+1; 
end 

%% Edge points processing 

% organize points on edge into a single 2 x n matrix 
startPoint = 1; 
endPoint = size(PointEdge,3); 
for i = 1:size(PointEdge,2)
    singleMatrixPoint(1,startPoint:endPoint) = squeeze(PointEdge(1,i,:))';
    singleMatrixPoint(2,startPoint:endPoint) = squeeze(PointEdge(2,i,:))';
    startPoint = endPoint+1;
    endPoint = endPoint+ size(PointEdge,3);
end 

% remove edge points that are outside the image 
for i = 1: size(singleMatrixPoint,2)
    if(singleMatrixPoint(:,i)==[-1;-1])
        singleMatrixPoint(:,i) = [0;0];
    end 
end 

% binarized the area thresholding the area 
% remove repeat points and index the point in ascending order in y
% direction 
pointTran = floor(abs(singleMatrixPoint'));
[C,ia,ic] = unique(pointTran,'rows'); % y 
newPoint = pointTran(ia,:); 
sortPoint = sortrows(newPoint,2); % sorting according each column, x direction 
% check if the first row of sortPoint is zero
if ~any(sortPoint(1,:))
    sortPoint(1,:) = [];
end 

% initialize the binary image 
SliceBW = zeros(size(Slice)); 
% conduct a column research for the maximum and minimum location of the
% edge point in the colunm(y) to estimate the actual edge of left ventricle
startColn = min(sortPoint(:,2));endColn = max(sortPoint(:,2)); 
currentColn = startColn; rangeColn = endColn - startColn+1;
minPreviousRow = 0; maxPreviousRow = 0;
binarizedRange = zeros(rangeColn,3); 

for i = 1: rangeColn
    if(find(sortPoint(:,2)==currentColn))
        indexColn = (find(sortPoint(:,2)==currentColn)); 
        for j = 1:size(indexColn,1)
            currentRow(j) = sortPoint(indexColn(j),1); 
        end 
        minCurrentRow = min(currentRow); 
        maxCurrentRow = max(currentRow); 
         % check if the minX to minY is spans off the radius length 
        if (i == 1)||(i == 2)||(i == 3)
        % for the first 3 column, just record the maximum and minimum found
        % previousely 
        else
            % in other cases, check if the difference between maximum and
            % minimum is smaller than the research radius, if it is smaller
            % than the search radius, take the previous row maximum or
            % minimum 
            if ((maxCurrentRow-minCurrentRow) <= (searchRadius/1.5))
                diffCurrent = maxCurrentRow-minCurrentRow;
                if (abs(minPreviousRow- minCurrentRow)>diffCurrent) 
                    minCurrentRow = minPreviousRow; 
                elseif (abs(maxPreviousRow - maxCurrentRow)>diffCurrent)
                    maxCurrentRow = maxPreviousRow;
                end 
            end 
        end 
        % record the range location in row and column 
        binarizedRange(i,1) = minCurrentRow; % minY
        binarizedRange(i,2) = maxCurrentRow; % maxY
        binarizedRange(i,3) = currentColn; % X
        currentColn = currentColn +1; 
    else
        currentColn = currentColn +1;
        continue; 
    end 
    minPreviousRow = minCurrentRow; maxPreviousRow = maxCurrentRow;
end

zeroRow = any(binarizedRange,2);
% fill the binarized slice pixel with 1 if pixel location is within the
% binarized range for left ventricle 
for i = 1:size(binarizedRange,1)
    if zeroRow(i)
        x = binarizedRange(i,3);
        yMin = binarizedRange(i,1); yMax = binarizedRange(i,2);
        SliceBW(yMin:yMax,x) = 1; 
    else 
        x = x+1;
        SliceBW(yMin:yMax,x) = 1;
    end 
end 
 
%% Calculate the area 
numPixel = nnz(SliceBW); 
areaPixel = stepSize^2; 
Area = numPixel * areaPixel; 
DiskVolume = Area * thickness;% cm 
%% visualization
% uncomment to visualize overlapping area and the left ventricle 
% figure; imshowpair(Slice,SliceBW); 

end 