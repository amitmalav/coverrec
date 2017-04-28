%%%%%%%%%%%% rotation %%%%%%%%%%%%%%%%%%


% image read

[FileName,Path] = uigetfile('*.jpg','Choose fixed image file');
I = imread(fullfile(Path,FileName));

% grayscale edge detection

% gray
grayImg = rgb2gray(I);

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

rtImage = imrotate(I,angleRot);

%%%%%%%%%%%%%%% End of Preprocessing %%%%%%%%%%%%%%%%%

        filter_lower = 5;


        min_area = 50;

        max_area = 4000;

        expansionAmount = 0.02;
    
        strokeWidthThreshold = 0.6;
    
        matlabocr = 1;
        teocr = 0;
    filter_lower = 6;
    min_area = 50;
    max_area = 4000;
    expansionAmount = 0.03;
    strokeWidthThreshold = 0.6;
    matlabocr = 1;
    teocr = 0;
    colorImage = rtImage;
%% Step 1: Detect Candidate Text Regions Using MSER
    %% Step 1: Detect Candidate Text Regions Using MSER
     I = rgb2gray(colorImage);
    sz = size(I);
    % Detect MSER regions.
    [mserRegions] = detectMSERFeatures(I,'RegionAreaRange',[min_area max_area],'ThresholdDelta',0.1);
    pixelIdxList = cellfun(@(xy)sub2ind(sz, xy(:,2), xy(:,1)), ...
        mserRegions.PixelList, 'UniformOutput', false);
    mserConnComp.Connectivity = 8;
    mserConnComp.ImageSize = sz;
    mserConnComp.NumObjects = mserRegions.Count;
    mserConnComp.PixelIdxList = pixelIdxList;
    % Detect MSER regions.
%     figure(3)
%     imshow(I)
%     hold on;
%     plot(mserRegions, 'showPixelList', true,'showEllipses',false)
%     title('MSER regions')
%     hold off

    %% Step 2:Remove Non-Text Regions Based On Basic Geometric Properties
    % Use regionprops to measure MSER properties
    mserStats = regionprops(mserConnComp, 'BoundingBox', 'Eccentricity', ...
        'Solidity', 'Extent', 'Euler', 'Image');

    % Compute the aspect ratio using bounding box data.
    bbox = vertcat(mserStats.BoundingBox);
    w = bbox(:,3);
    h = bbox(:,4);
    aspectRatio = w./h;
    % Threshold the data to determine which regions to remove. These thresholds
    % may need to be tuned for other images.
    filterIdx = aspectRatio' > 3; 
    filterIdx = filterIdx | [mserStats.Eccentricity] > .995 ;
    filterIdx = filterIdx | [mserStats.Solidity] < .3;
    filterIdx = filterIdx | [mserStats.Extent] < 0.1 | [mserStats.Extent] > 1;
    filterIdx = filterIdx | [mserStats.EulerNumber] < -4;

    % Remove regions
    mserStats(filterIdx) = [];
    mserRegions(filterIdx) = [];

    % Show remaining regions
%     figure(4)
%     imshow(I)
%     hold on;
%     plot(mserRegions, 'showPixelList', true,'showEllipses',false)
%     title('After Removing Non-Text Regions Based On Geometric Properties')
%     hold off

    %% Step 3: Remove Non-Text Regions Based On Stroke Width Variation
    for j = 1:numel(mserStats)

        regionImage = mserStats(j).Image;
        regionImage = padarray(regionImage, [1 1], 0);

        distanceImage = bwdist(~regionImage);
        skeletonImage = bwmorph(regionImage, 'thin', inf);

        strokeWidthValues = distanceImage(skeletonImage);

        strokeWidthMetric = std(strokeWidthValues)/mean(strokeWidthValues);

        strokeWidthFilterIdx(j) = strokeWidthMetric > strokeWidthThreshold;

    end

    % Remove regions based on the stroke width variation
    mserRegions(strokeWidthFilterIdx) = [];
    mserStats(strokeWidthFilterIdx) = [];

    % Show remaining regions
%     figure(5)
%     imshow(I)
%     hold on
%     plot(mserRegions, 'showPixelList', true,'showEllipses',false)
%     title('After Removing Non-Text Regions Based On Stroke Width Variation')
%     hold off
    %% Step 4: Merge Text Regions For Final Detection Result
    % At this point, all the detection results are composed of individual text
    % characters. To use these results for recognition tasks, such as OCR, the
    % individual text characters must be merged into words or text lines. This
    % enables recognition of the actual words in an image, which carry more
    % meaningful information than just the individual characters. For example,
    % recognizing the string 'EXIT' vs. the set of individual characters
    % {'X','E','T','I'}, where the meaning of the word is lost without the
    % correct ordering.
    %
    % One approach for merging individual text regions into words or text lines
    % is to first find neighboring text regions and then form a bounding box
    % around these regions. To find neighboring regions, expand the bounding
    % boxes computed earlier with |regionprops|. This makes the bounding boxes
    % of neighboring text regions overlap such that text regions that are part
    % of the same word or text line form a chain of overlapping bounding boxes.

    % Get bounding boxes for all the regions
    bboxes = vertcat(mserStats.BoundingBox);

    % Convert from the [x y width height] bounding box format to the [xmin ymin
    % xmax ymax] format for convenience.
    xmin = bboxes(:,1);
    ymin = bboxes(:,2);
    xmax = xmin + bboxes(:,3) - 1;
    ymax = ymin + bboxes(:,4) - 1;

    % Expand the bounding boxes by a small amount.
    xmin = (1-expansionAmount) * xmin;
    ymin = (1-expansionAmount) * ymin;
    xmax = (1+expansionAmount) * xmax;
    ymax = (1+expansionAmount) * ymax;

    % Clip the bounding boxes to be within the image bounds
    xmin = max(xmin, 1);
    ymin = max(ymin, 1);
    xmax = min(xmax, size(I,2));
    ymax = min(ymax, size(I,1));

    % Show the expanded bounding boxes
    expandedBBoxes = [xmin ymin xmax-xmin+1 ymax-ymin+1];
    IExpandedBBoxes = insertShape(colorImage,'Rectangle',expandedBBoxes,'LineWidth',3);

%     figure(6)
%     imshow(IExpandedBBoxes)
%     title('Expanded Bounding Boxes Text')
    %% Get rid of the text where bounding boxes are not aligned together
    % Filter by the upper bound
    expandedBBoxes_fil = expandedBBoxes;
    y_pos_up = ymin;
    [v_up,pos_up] = hist(y_pos_up,200);
    half_width_up = pos_up(2) - pos_up(1);
    left_bond_up = pos_up-half_width_up;
    right_bond_up = pos_up+half_width_up;
    left_bond_up = left_bond_up(v_up>filter_lower);
    right_bond_up = right_bond_up(v_up>filter_lower);

    % Filter by the lower bound
    y_pos_down = ymax;
    [v_down,pos_down] = hist(y_pos_down,200);
    half_width_down = pos_down(2) - pos_down(1);
    left_bond_down = pos_down-half_width_down;
    right_bond_down = pos_down+half_width_down;
    left_bond_down = left_bond_down(v_down>filter_lower);
    right_bond_down = right_bond_down(v_down>filter_lower);

    % Only keep the bounding boxes that pass both filters
    for i = 1:length(expandedBBoxes)
        flagdown = 0;
        flagup = 0;
        for j = 1:length(left_bond_down)
            if (ymax(i) - left_bond_down(j)) * ...
                    (ymax(i) - right_bond_down(j))<0 
                flagdown = 1;
            end
        end
        for j = 1:length(left_bond_up)
            if (ymin(i) - left_bond_up(j)) * ...
                    (ymin(i) - right_bond_up(j))<0 
                flagup = 1;
            end
        end
        if flagup == 0 || flagdown == 0
                expandedBBoxes_fil(i,:) = 0;
        end
    end
    filter_idx = expandedBBoxes_fil(:,1)==0;
    xmin = xmin(expandedBBoxes_fil(:,1)~=0);
    ymin = ymin(expandedBBoxes_fil(:,1)~=0);
    xmax = xmax(expandedBBoxes_fil(:,1)~=0);
    ymax = ymax(expandedBBoxes_fil(:,1)~=0);
    expandedBBoxes_fil(filter_idx,:) = [];
    expandedBBoxes = expandedBBoxes_fil;
    IExpandedBBoxes = insertShape(colorImage,'Rectangle',expandedBBoxes,'LineWidth',3);
%   figure(7)
%     imshow(IExpandedBBoxes)
%     title('Aligned Bounding Boxes Text')


    %%
    % Now, the overlapping bounding boxes can be merged together to form a
    % single bounding box around individual words or text lines. To do this,
    % compute the overlap ratio between all bounding box pairs. This quantifies
    % the distance between all pairs of text regions so that it is possible to
    % find groups of neighboring text regions by looking for non-zero overlap
    % ratios. Once the pair-wise overlap ratios are computed, use a |graph| to
    % find all the text regions "connected" by a non-zero overlap ratio.
    %
    % Use the |bboxOverlapRatio| function to compute the pair-wise overlap
    % ratios for all the expanded bounding boxes, then use |graph| to find all
    % the connected regions.

%     % Compute the overlap ratio
%     overlapRatio = bboxOverlapRatio(expandedBBoxes, expandedBBoxes);
% 
%     % Set the overlap ratio between a bounding box and itself to zero to
%     % simplify the graph representation.
%     n = size(overlapRatio,1); 
%     overlapRatio(1:n+1:n^2) = 0;
% 
%     % Create the graph
%     g = graph(overlapRatio);
% 
%     % Find the connected text regions within the graph
%     componentIndices = conncomp(g);
[box_hist,region_center] = hist(expandedBBoxes(:,2),100);
idx = 1;
region_label = 1;
all_label = zeros(1,length(box_hist));
while(idx<=length(box_hist))
    if(box_hist(idx)>0)
        while(box_hist(idx)>0)
            all_label(idx) = region_label;
            if idx == length(box_hist)
                break
            else
                idx = idx + 1;
            end
        end
        region_label = region_label + 1;
    end
    idx = idx + 1;
end
% Now get the median position for different group
num_group = max(all_label);
pos_mat = zeros(num_group,2);
for i = 1:num_group
    idx_ = find(all_label==i);
    candi = region_center(idx_);
    pos_mat(i,1) = i;
    pos_mat(i,2) = mean(candi);
end
componentIndices = zeros(1,length(expandedBBoxes(:,2)));
for i = 1:length(componentIndices)
    diff = abs(pos_mat(:,2) - expandedBBoxes(i,2));
    group_idx = find(diff==min(diff));
    componentIndices(i) = group_idx(1);
end
% Assign bounding boxes to different groups
    %%
    % The output of |conncomp| are indices to the connected text regions to
    % which each bounding box belongs. Use these indices to merge multiple
    % neighboring bounding boxes into a single bounding box by computing the
    % minimum and maximum of the individual bounding boxes that make up each
    % connected component.

    % Merge the boxes based on the minimum and maximum dimensions.
    xmin = accumarray(componentIndices', xmin, [], @min);
    ymin = accumarray(componentIndices', ymin, [], @min);
    xmax = accumarray(componentIndices', xmax, [], @max);
    ymax = accumarray(componentIndices', ymax, [], @max);
    fil_idx = zeros(1,length(xmin));
    for i = 1:length(xmin)
        xmin_c = xmin(i);
        ymin_c = ymin(i);
        xmax_c = xmax(i);
        ymax_c = ymax(i);
        for j = 1:length(xmin)
            if xmin_c>xmin(j) && xmax_c<xmax(j) && ymin_c>ymin(j) && ymax_c<ymax(j) 
                fil_idx(i) = 1;
            end
        end
    end
    fil_idx_ = find(fil_idx==1);
    xmin(fil_idx_) = [];
    ymin(fil_idx_) = [];
    xmax(fil_idx_) = [];
    ymax(fil_idx_) = [];
    % Compose the merged bounding boxes using the [x y width height] format.
    textBBoxes = [xmin ymin xmax-xmin+1 ymax-ymin+1];

    %%
    % Finally, before showing the final detection results, suppress false text
    % detections by removing bounding boxes made up of just one text region.
    % This removes isolated regions that are unlikely to be actual text given
    % that text is usually found in groups (words and sentences).
    
    % Remove bounding boxes that only contain one text region
    numRegionsInGroup = histcounts(componentIndices);
    %textBBoxes(numRegionsInGroup == 1, :) = [];

    % Show the final text detection result.
    ITextRegion = insertShape(colorImage, 'Rectangle', textBBoxes,'LineWidth',3);

%     figure(8)
    imshow(ITextRegion)
    title('Detected Text')

    %% Step 5: Recognize Detected Text Using OCR
    % After detecting the text regions, use the |ocr| function to recognize the
    % text within each bounding box. Note that without first finding the text
    % regions, the output of the |ocr| function would be considerably more
    % noisy.
    [sort_r,idx]= sort(textBBoxes(:,2));
    Detection_results = {};
        ocrtxt = ocr(I, textBBoxes);
        true_txt = {};
        for i = 1:numel(ocrtxt)
            fil_idx = (ocrtxt(i).CharacterConfidences<0.7);
            Detection_results{i} = ocrtxt(i).Text;  
            Detection_results{i}(fil_idx) = [];
            character = isletter(Detection_results{i});
            real_ratio = sum(character)/length(character);
            if real_ratio > 0.8
                true_txt{end+1} = ocrtxt(i).Text;
            end
        display('Using matlab ocr');
        % Alternative Step 5: Using tessract to perform better ocr
    end
    % Reformat the result according to the position of the textbox
    Detection_result = '';
    Searching_result = '';
    for i = 1:length(Detection_results)
        current_s = Detection_results{idx(i)};
        current_s = current_s(1:end-2);
        if i == 1
%             Detection_result = word_correction(Detection_results{idx(i)});
            Detection_result = strrep(current_s,sprintf('\n'),' ');
%             Searching_result = word_correction(Detection_results{idx(i)});
            Searching_result = strrep(current_s,sprintf('\n'),' ');
        else
%             Detection_result = strcat(Detection_result,'\n',word_correction(Detection_results{idx(i)}));
            Detection_result = strcat(Detection_result,'\n',strrep(current_s,sprintf('\n'),' '));
%             Searching_result = [Searching_result,' ',word_correction(Detection_results{idx(i)})];
            Searching_result = [Searching_result,' ',strrep(current_s,sprintf('\n'),' ')];
        end
    end