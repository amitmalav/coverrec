clear;

%%%%%%%%%%%% rotation %%%%%%%%%%%%%%%%%%

% image read
[FileName,Path] = uigetfile('*.jpg','Choose fixed image file');
FI = imread(fullfile(Path,FileName));

% grayscale edge detection

% gray
grayImg = rgb2gray(FI);

%edge
edges = edge(grayImg,'sobel');

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

% need processing for MSER regions 
plist = cellfun(@(doubMat)sub2ind(size(grayImg), doubMat(:,2), doubMat(:,1)),mser.PixelList, 'UniformOutput', false);
mser1.Connectivity = 8;
mser1.ImageSize = size(grayImg);
mser1.NumObjects = mser.Count;
mser1.PixelIdxList = plist;

mserReg = regionprops(mser1,'BoundingBox','Eccentricity','Solidity','Extent','Euler','Image');

% Compute aspect ratio w/h
containerBox = vertcat(mserReg.BoundingBox);
w = containerBox(:,3); 
h = containerBox(:,4);
aspectRatio = w./h;

%detect usless mser regions 
filteredRegions = aspectRatio' > 3.5 | [mserReg.Eccentricity] > 0.995 | [mserReg.Solidity] < 0.3 | [mserReg.Extent] < 0.2 | [mserReg.Extent] > 0.9 | [mserReg.EulerNumber] < -4;

%remove the detected useless regions
mserReg(filteredRegions) = [];
mser(filteredRegions) = [];

% filter regions based on stroke width
%Text regions tend to have little stroke width variation, whereas non-text regions tend to have larger variations.

[m, n] = size(mserReg);

for i = 1:m

	%padding the region image for filtering
	msImage = padarray(mserReg(i).Image, [1,1]);

	%computing distance of nearest non zero pixel for each pixel
	distanceMatrix = bwdist(~msImage);

	%applies thin morph op
	morphmat  = bwmorph(msImage, 'thin', inf);
	sValue = distanceMatrix(morphmat);
	sMat = std(sValue)/mean(sValue);
	sfilteredRegions(i) = sMat > 0.55;
end

mserReg(sfilteredRegions) = [];
mser(sfilteredRegions) = [];

%create boundry boxes for every word
containerBox = vertcat(mserReg.BoundingBox);
w = containerBox(:,3); 
h = containerBox(:,4);

%get end points of boundry boxes
x_small = containerBox(:, 1);
y_small = containerBox(:, 2);
x_big = x_small + w - 1;
y_big = y_small + h - 1;

%expand bboxes by little amount
x_small = (1 - 0.02) * x_small;
y_small = (1 - 0.02) * y_small;
x_big = (1 + 0.02) * x_big;
y_big = (1 + 0.02) * y_big;

%clip end points
[size1,size2] = size(rtImage);
x_small = max(x_small, 1);
y_small = max(y_small, 1);
x_big = min(x_big, size2);
y_big = min(y_big, size1);

bBoxMat = [x_small - 1, y_small - 1, x_big - x_small + 1, y_big - y_small + 1];

%remove bboxes which are not aligned
oratio = bboxOverlapRatio(bBoxMat, bBoxMat);

% Set the overlap ratio between a bounding box and itself to zero to
n = size(oratio, 1);
for i = 1:n
	oratio(i, i) = 0;
end

g = graph(oratio);
componentIndices = conncomp(g);
x_small = accumarray(componentIndices', x_small, [], @min);
y_small = accumarray(componentIndices', y_small, [], @min);
x_big = accumarray(componentIndices', x_big, [], @max);
y_big = accumarray(componentIndices', y_big, [], @max);

bBigBox =  [x_small, y_small, x_big - x_small + 1, y_big - y_small + 1];
numRegionsInGroup = histcounts(componentIndices);
bBigBox(numRegionsInGroup == 1, :) = [];
ITextRegion = insertShape(rtImage, 'Rectangle', bBigBox,'LineWidth',4);

% get words using ocr
ocrResult = ocr(rtImage, bBigBox);

% check for valid text from ocrout
ocrSelectResult = {};
trueOut = {};
true1 = {};
for i = 1:size(ocrResult,1)
	ocrSelectResult{i} = ocrResult(i).Text;
	select_ocr_text = ((ocrResult(i).CharacterConfidences) < 0.74);
	ocrSelectResult{i}(select_ocr_text) = [];
	c = isletter(ocrSelectResult{i});
	text_ratio = sum(c)/length(c);
	trueOut{end+1} = ocrResult(i).Text; 
	true1{end+1} = ocrSelectResult;
end

for i =1:size(ocrSelectResult,2)
	m = ocrSelectResult(1,i);
	a = m{:};
	if size(a,2) > 2 
		if ~isspace(a(1,1)) fprintf('%c',a(1,1)); end
		for j = 2:size(a,2)
			if ~isspace(a(1,j)) 
				fprintf('%c',a(1,j));
			elseif ~isspace(a(1,j-1))
				fprintf('%c',a(1,j));
			end
		end
	end
end