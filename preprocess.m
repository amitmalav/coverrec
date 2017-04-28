%%%%%%%%%%%% rotation %%%%%%%%%%%%%%%%%%


% image read

[FileName,Path] = uigetfile('*.jpg','Choose fixed image file');
FI = imread(fullfile(Path,FileName));

% grayscale edge detection

% gray
grayImg = rgb2gray(FI);

%edge
edges = edge(grayImg,'Canny');

% hough transform
h = hough(edges);

% get peaks 
P = houghpeaks(h,8);
angles = P(:,2); % get angles of peak

% check if we have enough angles and take the major one
if length(unique(angles)) > 4
	angleRot = angles(1);
else
	angleRot = mode(P(:,2));
end

while(angleRot > 45 || angleRot < -45)
	if angleRot > 0
		angleRot = angleRot - 90;
	else
		angleRot = angleRot + 90;
	end
end

rtImage = imrotate(FI,angleRot);

%%%%%%%%%%%%%%% End of Preprocessing %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% MSER %%%%%%%%%%%%%%%%%%%%%%%

grayImg = rgb2gray(rtImage);
[mser] = detectMSERFeatures(grayImg,'RegionAreaRange',[200 8000],'ThresholdDelta',4);
sz = size(grayImg);

% mserStats = regionprops(mser, 'BoundingBox', 'Eccentricity', ... 'Solidity', 'Extent', 'Euler', 'Image');
% [mser] = detectMSERFeatures(grayImg);
plist = cellfun(@(xy)sub2ind(size(grayImg), xy(:,2), xy(:,1)),mser.PixelList, 'UniformOutput', false);
mser1.Connectivity = 8;
mser1.ImageSize = size(grayImg);
mser1.NumObjects = mser.Count;
mser1.PixelIdxList = plist;

mserstat = regionprops(mser1,'BoundingBox','Eccentricity','Solidity','Extent','Euler','Image');

% Compute aspect ratio w/h

bBox = vertcat(mserstat.BoundingBox);
w = bBox(:,3); 
h = bBox(:,4);
aspectRatio = w./h;
%detect usless mser regions 

filteredRegions = aspectRatio' > 3 | [mserstat.Eccentricity] > .995 | [mserstat.Solidity] < .3 | [mserstat.Extent] < .2 | [mserstat.Extent] > .9 | [mserstat.EulerNumber] < -4;

%remove useless regions

mserstat(filteredRegions) = [];
mser(filteredRegions) = [];

% filter regions based on stroke width
%Text regions tend to have little stroke width variation, whereas non-text regions tend to have larger variations.

[m, n] = size(mserstat);

for i = 1:numel(mserstat)
	%padding the region image for filtering
	msImage = padarray(mserstat(i).Image, [1,1], 0);
	%computing distance of nearest non zero pixel for each pixel
	distanceMatrix = bwdist(~msImage);
	%applies thin morph op
	morphmat  = bwmorph(msImage, 'thin', inf);
	sValue = distanceMatrix(morphmat);
	sMat = std(sValue)/mean(sValue);
	sfilteredRegions(i) = sMat > 0.6;
end

mserstat(sfilteredRegions) = [];
mser(sfilteredRegions) = [];

%create boundry boxes for every word


bBox = vertcat(mserstat.BoundingBox);
w = bBox(:,3); 
h = bBox(:,4);

%get end points of boundry boxes

xlow = bBox(:, 1);
ylow = bBox(:, 2);
xhi = xlow + w - 1;
yhi = ylow + h - 1;

%expand bboxes by little amount

xlow = (1 - 0.02) * xlow;
ylow = (1 - 0.02) * ylow;
xhi = (1 + 0.02) * xhi;
yhi = (1 + 0.02) * yhi;

%clip end points

xlow = max(xlow, 1);
ylow = max(ylow, 1);
xhi = min(xhi, size(rtImage, 2));
yhi = min(yhi, size(rtImage, 1));

bBoxMat = [xlow, ylow, xhi - xlow + 1, yhi - ylow + 1];

%draw bbox

bboxdraw = insertShape(rtImage,'Rectangle',bBoxMat,'LineWidth',3);

%remove bboxes which are not aligned

oratio = bboxOverlapRatio(bBoxMat, bBoxMat);

% Set the overlap ratio between a bounding box and itself to zero to

n = size(oratio, 1);
for i = 1:n
	oratio(i, i) = 0;
end

g = graph(oratio);
componentIndices = conncomp(g);
xlow = accumarray(componentIndices', xlow, [], @min);
ylow = accumarray(componentIndices', ylow, [], @min);
xhi = accumarray(componentIndices', xhi, [], @max);
yhi = accumarray(componentIndices', yhi, [], @max);

bBigBox =  [xlow, ylow, xhi - xlow + 1, yhi - ylow + 1];
numRegionsInGroup = histcounts(componentIndices);
bBigBox(numRegionsInGroup == 1, :) = [];
ITextRegion = insertShape(rtImage, 'Rectangle', bBigBox,'LineWidth',3);


ocrtxt = ocr(rtImage, bBigBox);
rectext = ocrtxt.Text;