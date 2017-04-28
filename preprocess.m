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
    [mser] = detectMSERFeatures(grayImg,'RegionAreaRange',[200 8000],'ThresholdDelta',0.1);
sz = size(grayImg);

% mserStats = regionprops(mser, 'BoundingBox', 'Eccentricity', ... 'Solidity', 'Extent', 'Euler', 'Image');

plist = cellfun(@(xy)sub2ind(sz, xy(:,2), xy(:,1)),mser.PixelList, 'UniformOutput', false);
