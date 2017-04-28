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

% orientations = P(:,2);

% if length(unique(orientations)) > 3    
%     temp_orientation = orientations(1);
%     rotate_angle = 0;
% else
%     temp_orientation = mode(P(:,2));
%     rotate_angle = transpose(temp_orientation)-90;
% end

% rotate_angle = transpose(temp_orientation)-90;

% while(abs(rotate_angle) > 45)
%     rotate_angle = rotate_angle - sign(rotate_angle)*90;
% end

% angle = rotate_angle;

% rotate_img = imrotate(FI,rotate_angle,'bilinear');
% %figure;
% imshow(rotate_img);