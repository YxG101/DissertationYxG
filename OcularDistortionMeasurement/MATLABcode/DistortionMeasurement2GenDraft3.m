%% Pre-processing_cropping
clc;
clear;
close all;

Name = "292.84mAnaked2";
I = imread(Name+".bmp");
 % figure,
 % imshow(I);

I2 = imresize(I,3);
 % figure,
 % imshow(I2);

[centersBeforeCrop,radiiBeforeCrop] = imfindcircles(I2,[7 50],'ObjectPolarity','bright');
% viscircles(centersBeforeCrop, radiiBeforeCrop,'EdgeColor','b');

centersYcoorBeforeCrop = centersBeforeCrop(:,2);
centersXcoorBeforeCrop = centersBeforeCrop(:,1);

XmaxForCrop = max(centersXcoorBeforeCrop);
XminForCrop = min(centersXcoorBeforeCrop);
YmaxForCrop = max(centersYcoorBeforeCrop);
YminForCrop = min(centersYcoorBeforeCrop);

I4 = imcrop(I2, [XminForCrop-200 YminForCrop-200 XmaxForCrop-XminForCrop+400 YmaxForCrop-YminForCrop+400]);
I3 = I4;
 % figure;
 % imshow(I4);
% figure;
% imshow(I3);
%% Find circles and sort
% 1.1 find the grid center 
[centers,radii] = imfindcircles(I3,[7 20],'ObjectPolarity','bright');
% viscircles(centers, radii,'EdgeColor','b');

radiiMean = mean(radii);

centersYcoor = centers(:,2);
centersXcoor = centers(:,1);

minX = min(centersXcoor); 
maxX = max(centersXcoor);
minY = min(centersYcoor); 
maxY = max(centersYcoor);

XCenter = (minX + maxX)/2;
YCenter = (minY + maxY)/2;

centers1 = [centersXcoor,centersYcoor];
GridCenterIndex1 = dsearchn(centers1,[XCenter YCenter]);
GridCenter = centers1(GridCenterIndex1,:);
% viscircles(GridCenter, 13,'EdgeColor','r');
foundCentroidsSorted = zeros(121,2);

%  1.2 find the 6 row including center points, sort and assign orders
Grid6throw = zeros(11,2);
m = 1;
for n  = 1:121
    if abs(centers1(n,2) - GridCenter(:,2)) < max(radii)+5
       Grid6throw(m,:) = centers1(n,:);
       m = m+1;
    end
end
% for n = 1:11
%     viscircles(Grid6throw(n,:), 20,'EdgeColor','g');
% end

Grid6throwSorted = sort(Grid6throw,1);
Grid6throwSortedCalibrated = [Grid6throwSorted(:,1)-GridCenter(1,1),-Grid6throwSorted(:,2)+GridCenter(1,2)];
m = 1;
for n = 56:66
    foundCentroidsSorted(n,:) = Grid6throwSortedCalibrated(m,:);
    m = m+1;
end

%  1.3 using 6th row to eliminate tilt
tantheta1 = zeros(10,1);
for a = 1:10
    tantheta1(a,1) = (Grid6throwSortedCalibrated(a+1,2)-Grid6throwSortedCalibrated(a,2))./ (Grid6throwSortedCalibrated(a+1,1)-Grid6throwSortedCalibrated(a,1));
end

tanaverage1 = sum(tantheta1)./10;
theta1 = atan(tanaverage1).*180./pi;
I5 = imrotate(I3,theta1,'bicubic');
figure
imshow(I5);

%%  2.1 find row 6 after tilt
[centers,radii] = imfindcircles(I5,[7 20],'ObjectPolarity','bright');
% viscircles(centers, radii,'EdgeColor','b');

radiiMean = mean(radii);

centersYcoor = centers(:,2);
centersXcoor = centers(:,1);

minX = min(centersXcoor); 
maxX = max(centersXcoor);
minY = min(centersYcoor); 
maxY = max(centersYcoor);

XCenter = (minX + maxX)/2;
YCenter = (minY + maxY)/2;

centers1 = [centersXcoor,centersYcoor];
GridCenterIndex1 = dsearchn(centers1,[XCenter YCenter]);
GridCenter = centers1(GridCenterIndex1,:);
viscircles(GridCenter, 13,'EdgeColor','r');
foundCentroidsSorted = zeros(121,2);

%  2.2 find the 6 row including center points, sort and assign orders
Grid6throw = zeros(11,2);
m = 1;
for n  = 1:121
    if abs(centers1(n,2) - GridCenter(:,2)) < max(radii)
       Grid6throw(m,:) = centers1(n,:);
       m = m+1;
    end
end
% for n = 1:11
%     viscircles(Grid6throw(n,:), 20,'EdgeColor','g');
% end

Grid6throwX = Grid6throw(:,1);
Grid6throwY = Grid6throw(:,2);
[Grid6throwXSorted,index] = sort(Grid6throwX);
Grid6throwYSorted = Grid6throwY(index);
Grid6throwSorted = [Grid6throwXSorted,Grid6throwYSorted];
Grid6throwSortedCalibrated = [Grid6throwSorted(:,1)-GridCenter(1,1),-Grid6throwSorted(:,2)+GridCenter(1,2)];
m = 1;
for n = 56:66
    foundCentroidsSorted(n,:) = Grid6throwSortedCalibrated(m,:);
    m = m+1;
end

%  2.2 find the 6 column including center points, sort and assign orders

Grid6thcolumn = zeros(11,2);
m = 1;
for n  = 1:121
    if abs(centers1(n,1) - GridCenter(:,1)) < max(radii)
       Grid6thcolumn(m,:) = centers1(n,:);
       m = m+1;
    end
end
% for n = 1:11
%     viscircles(Grid6thcolumn(n,:), 20,'EdgeColor','g');
% end

Grid6thcolumnX = Grid6thcolumn(:,1);
Grid6thcolumnY = Grid6thcolumn(:,2);
[Grid6thcolumnYSorted,index] = sort(Grid6thcolumnY);
Grid6thcolumnXSorted = Grid6thcolumnX(index);
Grid6thcolumnSorted = [Grid6thcolumnXSorted,Grid6thcolumnYSorted];
Grid6thcolumnSortedCalibrated = [Grid6thcolumnSorted(:,1)-GridCenter(1,1),-Grid6thcolumnSorted(:,2)+GridCenter(1,2)];

for n = 1:11
    foundCentroidsSorted(11*(n-1)+6,:) = Grid6thcolumnSortedCalibrated(n,:);
end

%  2.3 delete the center
centersCalibrated1 = [centersXcoor-GridCenter(1,1),-centersYcoor+GridCenter(1,2)];

[indexForGridCenter1,indexForGridCenter2] = find(centersCalibrated1==[0,0]);
centersCalibrated2 = centersCalibrated1;
centersCalibrated2(indexForGridCenter1,:)=[];

%  2.4 delete the 6 row and column
centers2 = centers1;

for n = 1:11
    [indexForGridCenter1,indexForGridCenter2] = find(centers2 == Grid6thcolumnSorted(n,:));
    centers2(indexForGridCenter1,:) = [];
end
[indexForGridCenter1,indexForGridCenter2] = find(Grid6throwSorted==GridCenter);
Grid6throwSorted(indexForGridCenter1,:)=[];
for n = 1:10
    [indexForGridCenter1,indexForGridCenter2] = find(centers2 == Grid6throwSorted(n,:));
    centers2(indexForGridCenter1,:) = [];
end

%  2.5 find the smallest circle 4 neighbor points around center
centers3 = centers2;
centersCalibrated3 = [centers3(:,1)-GridCenter(1,1),-centers3(:,2)+GridCenter(1,2)];

%  2.6 find and sort 6-9 in the first circle
FirstCircle6to9NeighborCalibrated = zeros(4,2);
for n = 1:4
    NeighborIndex1 = dsearchn(centersCalibrated3,[0 0]);
    FirstCircle6to9NeighborCalibrated(n,:) = centersCalibrated3(NeighborIndex1,:);
    centersCalibrated3(NeighborIndex1,:)=[];
end

% for n = 1:4
%     viscircles(FirstCircle6to9Neighbor(n,:), 20,'EdgeColor','g');
% end


 for a = 1:4
     if FirstCircle6to9NeighborCalibrated(a,1) > 0 
         if FirstCircle6to9NeighborCalibrated(a,2)>0
                      GridCentroid11calibrated = FirstCircle6to9NeighborCalibrated(a,:);
                      GridCentroid11 = [GridCentroid11calibrated(1,1)+GridCenter(1,1),-GridCentroid11calibrated(1,2)+GridCenter(1,2)];

         else 
                      GridCentroid1Minus1calibrated = FirstCircle6to9NeighborCalibrated(a,:);
                      GridCentroid1Minus1 = [GridCentroid1Minus1calibrated(1,1)+GridCenter(1,1),-GridCentroid1Minus1calibrated(1,2)+GridCenter(1,2)];

         end
     else
         if FirstCircle6to9NeighborCalibrated(a,2)>0
                      GridCentroidMinus11calibrated = FirstCircle6to9NeighborCalibrated(a,:);
                      GridCentroidMinus11 = [GridCentroidMinus11calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus11calibrated(1,2)+GridCenter(1,2)];

         else
                      GridCentroidMinus1Minus1calibrated = FirstCircle6to9NeighborCalibrated(a,:);
                      GridCentroidMinus1Minus1 = [GridCentroidMinus1Minus1calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus1Minus1calibrated(1,2)+GridCenter(1,2)];
         end
     end
 end

 % viscircles(GridCentroid11, 20,'EdgeColor','y');
 % viscircles(GridCentroid1Minus1, 20,'EdgeColor','y');
 % viscircles(GridCentroidMinus11, 20,'EdgeColor','y');
 % viscircles(GridCentroidMinus1Minus1, 20,'EdgeColor','y');
foundCentroidsSorted(51,:) = GridCentroid11calibrated;
foundCentroidsSorted(73,:) = GridCentroid1Minus1calibrated;
foundCentroidsSorted(49,:) = GridCentroidMinus11calibrated;
foundCentroidsSorted(71,:) = GridCentroidMinus1Minus1calibrated;


%  2.7 find and sort 14-21 in the second circle
SecondCircle14to21NeighborCalibrated = zeros(8,2);

for n = 1:8
    NeighborIndex2 = dsearchn(centersCalibrated3,[0 0]);
    SecondCircle14to21NeighborCalibrated(n,:) = centersCalibrated3(NeighborIndex2,:);
    centersCalibrated3(NeighborIndex2,:)=[];
end

SecondCircle14to21Neighbor = [SecondCircle14to21NeighborCalibrated(:,1)+GridCenter(1,1),-SecondCircle14to21NeighborCalibrated(:,2)+GridCenter(1,2)];
% for n = 1:8
%     viscircles(SecondCircle14to21Neighbor(n,:), 20,'EdgeColor','g');
% end

c = 1; d = 1; e = 1; f = 1;
GridCentroidSection1calibrated = zeros(2,2);
GridCentroidSection2calibrated = zeros(2,2);
GridCentroidSection3calibrated = zeros(2,2);
GridCentroidSection4calibrated = zeros(2,2);
 for b = 1:8
     if SecondCircle14to21NeighborCalibrated(b,1) > 0 
         if SecondCircle14to21NeighborCalibrated(b,2)>0
                      GridCentroidSection1calibrated(c,:) = SecondCircle14to21NeighborCalibrated(b,:);
                      c = c+1;
         else 
                      GridCentroidSection2calibrated(d,:) = SecondCircle14to21NeighborCalibrated(b,:);
                      d = d+1;
         end
     else
         if SecondCircle14to21NeighborCalibrated(b,2)>0
                      GridCentroidSection4calibrated(e,:) = SecondCircle14to21NeighborCalibrated(b,:);
                      e = e+1;
         else
                      GridCentroidSection3calibrated(f,:) = SecondCircle14to21NeighborCalibrated(b,:);
                      f = f+1;
         end
     end
 end

[GridCentroid21calibratedX,index1] = max(GridCentroidSection1calibrated(:,1));
GridCentroid21calibrated = GridCentroidSection1calibrated(index1,:);

GridCentroid21 = [GridCentroid21calibrated(1,1)+GridCenter(1,1),-GridCentroid21calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroid21, 20,'EdgeColor','y');
foundCentroidsSorted(52,:) = GridCentroid21calibrated;
if index1 == 1
    GridCentroid12calibrated = GridCentroidSection1calibrated(2,:);
else 
    GridCentroid12calibrated = GridCentroidSection1calibrated(1,:);
end
GridCentroid12 = [GridCentroid12calibrated(1,1)+GridCenter(1,1),-GridCentroid12calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroid12, 20,'EdgeColor','y');
foundCentroidsSorted(40,:) = GridCentroid12calibrated;



[GridCentroid2Minus1calibratedX,index1] = max(GridCentroidSection2calibrated(:,1));
GridCentroid2Minus1calibrated = GridCentroidSection2calibrated(index1,:);

GridCentroid2Minus1 = [GridCentroid2Minus1calibrated(1,1)+GridCenter(1,1),-GridCentroid2Minus1calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroid2minus1, 20,'EdgeColor','y');
foundCentroidsSorted(74,:) = GridCentroid2Minus1calibrated;
if index1 == 1
    GridCentroid1Minus2calibrated = GridCentroidSection2calibrated(2,:);
else 
    GridCentroid1Minus2calibrated = GridCentroidSection2calibrated(1,:);
end
GridCentroid1Minus2 = [GridCentroid1Minus2calibrated(1,1)+GridCenter(1,1),-GridCentroid1Minus2calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroid1minus2, 20,'EdgeColor','g');
foundCentroidsSorted(84,:) = GridCentroid1Minus2calibrated;



[GridCentroidMinus1Minus2calibratedX,index1] = max(GridCentroidSection3calibrated(:,1));
GridCentroidMinus1Minus2calibrated = GridCentroidSection3calibrated(index1,:);

GridCentroidMinus1Minus2 = [GridCentroidMinus1Minus2calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus1Minus2calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroidMinus1Minus2, 20,'EdgeColor','y');
foundCentroidsSorted(82,:) = GridCentroidMinus1Minus2calibrated;
if index1 == 1
    GridCentroidMinus2Minus1calibrated = GridCentroidSection3calibrated(2,:);
else 
    GridCentroidMinus2Minus1calibrated = GridCentroidSection3calibrated(1,:);
end
GridCentroidMinus2Minus1 = [GridCentroidMinus2Minus1calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus2Minus1calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroidMinus2Minus1, 20,'EdgeColor','g');
foundCentroidsSorted(70,:) = GridCentroidMinus2Minus1calibrated;



[GridCentroidMinus12calibratedX,index1] = max(GridCentroidSection4calibrated(:,1));
GridCentroidMinus12calibrated = GridCentroidSection4calibrated(index1,:);

GridCentroidMinus12 = [GridCentroidMinus12calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus12calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroidMinus12, 20,'EdgeColor','y');
foundCentroidsSorted(38,:) = GridCentroidMinus12calibrated;
if index1 == 1
    GridCentroidMinus21calibrated = GridCentroidSection4calibrated(2,:);
else 
    GridCentroidMinus21calibrated = GridCentroidSection4calibrated(1,:);
end
GridCentroidMinus21 = [GridCentroidMinus21calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus21calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroidMinus21, 20,'EdgeColor','g');
foundCentroidsSorted(48,:) = GridCentroidMinus21calibrated;

%  2.8 find and sort 22-25 in the second circle
SecondCircle22to25NeighborCalibrated = zeros(4,2);
centersCalibrated4 = centersCalibrated3;
for n = 1:4
    NeighborIndex3 = dsearchn(centersCalibrated4,[0 0]);
    SecondCircle22to25NeighborCalibrated(n,:) = centersCalibrated4(NeighborIndex3,:);
    centersCalibrated4(NeighborIndex3,:)=[];
end

SecondCircle22to25Neighbor = [SecondCircle22to25NeighborCalibrated(:,1)+GridCenter(1,1),-SecondCircle22to25NeighborCalibrated(:,2)+GridCenter(1,2)];

%     viscircles(SecondCircle22to25Neighbor, 20,'EdgeColor','g');

 for a = 1:4
     if SecondCircle22to25NeighborCalibrated(a,1) > 0 
         if SecondCircle22to25NeighborCalibrated(a,2)>0
                      GridCentroid22calibrated = SecondCircle22to25NeighborCalibrated(a,:);
                      GridCentroid22 = [GridCentroid22calibrated(1,1)+GridCenter(1,1),-GridCentroid22calibrated(1,2)+GridCenter(1,2)];

         else 
                      GridCentroid2Minus2calibrated = SecondCircle22to25NeighborCalibrated(a,:);
                      GridCentroid2Minus2 = [GridCentroid2Minus2calibrated(1,1)+GridCenter(1,1),-GridCentroid2Minus2calibrated(1,2)+GridCenter(1,2)];

         end
     else
         if SecondCircle22to25NeighborCalibrated(a,2)>0
                      GridCentroidMinus22calibrated = SecondCircle22to25NeighborCalibrated(a,:);
                      GridCentroidMinus22 = [GridCentroidMinus22calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus22calibrated(1,2)+GridCenter(1,2)];

         else
                      GridCentroidMinus2Minus2calibrated = SecondCircle22to25NeighborCalibrated(a,:);
                      GridCentroidMinus2Minus2 = [GridCentroidMinus2Minus2calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus2Minus2calibrated(1,2)+GridCenter(1,2)];
         end
     end
 end

 % viscircles(GridCentroid22, 20,'EdgeColor','y');
 % viscircles(GridCentroid2Minus2, 20,'EdgeColor','y');
 % viscircles(GridCentroidMinus22, 20,'EdgeColor','y');
 % viscircles(GridCentroidMinus2Minus2, 20,'EdgeColor','y');
foundCentroidsSorted(41,:) = GridCentroid22calibrated;
foundCentroidsSorted(85,:) = GridCentroid2Minus2calibrated;
foundCentroidsSorted(37,:) = GridCentroidMinus22calibrated;
foundCentroidsSorted(81,:) = GridCentroidMinus2Minus2calibrated;

%  2.9 find and sort 30-45 in the third circle
ThirdCircle30to45NeighborCalibrated = zeros(8,2);
centersCalibrated5 = centersCalibrated4;
for n = 1:16
    NeighborIndex4 = dsearchn(centersCalibrated5,[0 0]);
    ThirdCircle30to45NeighborCalibrated(n,:) = centersCalibrated5(NeighborIndex4,:);
    centersCalibrated5(NeighborIndex4,:)=[];
end

ThirdCircle30to45Neighbor = [ThirdCircle30to45NeighborCalibrated(:,1)+GridCenter(1,1),-ThirdCircle30to45NeighborCalibrated(:,2)+GridCenter(1,2)];
% 
% for n = 1:16
%     viscircles(ThirdCircle30to46Neighbor(n,:), 20,'EdgeColor','g');
% end

c = 1; d = 1; e = 1; f = 1;
GridCentroidSection1calibrated = zeros(4,2);
GridCentroidSection2calibrated = zeros(4,2);
GridCentroidSection3calibrated = zeros(4,2);
GridCentroidSection4calibrated = zeros(4,2);
 for b = 1:16
     if ThirdCircle30to45NeighborCalibrated(b,1) > 0 
         if ThirdCircle30to45NeighborCalibrated(b,2)>0
                      GridCentroidSection1calibrated(c,:) = ThirdCircle30to45NeighborCalibrated(b,:);
                      c = c+1;
         else 
                      GridCentroidSection2calibrated(d,:) = ThirdCircle30to45NeighborCalibrated(b,:);
                      d = d+1;
         end
     else
         if ThirdCircle30to45NeighborCalibrated(b,2)>0
                      GridCentroidSection4calibrated(e,:) = ThirdCircle30to45NeighborCalibrated(b,:);
                      e = e+1;
         else
                      GridCentroidSection3calibrated(f,:) = ThirdCircle30to45NeighborCalibrated(b,:);
                      f = f+1;
         end
     end
 end

[GridCentroid13calibratedX,index1] = min(GridCentroidSection1calibrated(:,1));
GridCentroid13calibrated = GridCentroidSection1calibrated(index1,:);
GridCentroid13 = [GridCentroid13calibrated(1,1)+GridCenter(1,1),-GridCentroid13calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroid13, 20,'EdgeColor','g');
foundCentroidsSorted(29,:) = GridCentroid13calibrated;
GridCentroidSection1calibrated(index1,:) = [];

[GridCentroid31calibratedY,index2] = min(GridCentroidSection1calibrated(:,2));
GridCentroid31calibrated = GridCentroidSection1calibrated(index2,:);
GridCentroid31 = [GridCentroid31calibrated(1,1)+GridCenter(1,1),-GridCentroid31calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroid31, 20,'EdgeColor','y');
foundCentroidsSorted(53,:) = GridCentroid31calibrated;
GridCentroidSection1calibrated(index2,:) = [];

if GridCentroidSection1calibrated(1,1) > GridCentroidSection1calibrated(2,1) && ...
    GridCentroidSection1calibrated(1,2) < GridCentroidSection1calibrated(2,2)
    GridCentroid32calibrated = GridCentroidSection1calibrated(1,:);
    GridCentroid32 = [GridCentroid32calibrated(1,1)+GridCenter(1,1),-GridCentroid32calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroid32, 20,'EdgeColor','g');
    foundCentroidsSorted(42,:) = GridCentroid32calibrated;
    GridCentroid23calibrated = GridCentroidSection1calibrated(2,:);
    GridCentroid23 = [GridCentroid23calibrated(1,1)+GridCenter(1,1),-GridCentroid23calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroid23, 20,'EdgeColor','y');
    foundCentroidsSorted(30,:) = GridCentroid23calibrated;
else
    GridCentroid32calibrated = GridCentroidSection1calibrated(2,:);
    GridCentroid32 = [GridCentroid32calibrated(1,1)+GridCenter(1,1),-GridCentroid32calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroid32, 20,'EdgeColor','g');
    foundCentroidsSorted(42,:) = GridCentroid32calibrated;
    GridCentroid23calibrated = GridCentroidSection1calibrated(1,:);
    GridCentroid23 = [GridCentroid23calibrated(1,1)+GridCenter(1,1),-GridCentroid23calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroid23, 20,'EdgeColor','y');
    foundCentroidsSorted(30,:) = GridCentroid23calibrated;
end

[GridCentroid1Minus3calibratedX,index1] = min(GridCentroidSection2calibrated(:,1));
GridCentroid1Minus3calibrated = GridCentroidSection2calibrated(index1,:);
GridCentroid1Minus3 = [GridCentroid1Minus3calibrated(1,1)+GridCenter(1,1),-GridCentroid1Minus3calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroid1Minus3, 20,'EdgeColor','g');
foundCentroidsSorted(95,:) = GridCentroid1Minus3calibrated;
GridCentroidSection2calibrated(index1,:) = [];

[GridCentroid3Minus1calibratedY,index2] = max(GridCentroidSection2calibrated(:,2));
GridCentroid3Minus1calibrated = GridCentroidSection2calibrated(index2,:);
GridCentroid3Minus1 = [GridCentroid3Minus1calibrated(1,1)+GridCenter(1,1),-GridCentroid3Minus1calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroid3Minus1, 20,'EdgeColor','y');
foundCentroidsSorted(75,:) = GridCentroid3Minus1calibrated;
GridCentroidSection2calibrated(index2,:) = [];

if GridCentroidSection2calibrated(1,1) > GridCentroidSection2calibrated(2,1) && ...
    GridCentroidSection2calibrated(1,2) > GridCentroidSection2calibrated(2,2)
    GridCentroid3Minus2calibrated = GridCentroidSection2calibrated(1,:);
    GridCentroid3Minus2 = [GridCentroid3Minus2calibrated(1,1)+GridCenter(1,1),-GridCentroid3Minus2calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroid3Minus2, 20,'EdgeColor','g');
    foundCentroidsSorted(86,:) = GridCentroid3Minus2calibrated;
    GridCentroid2Minus3calibrated = GridCentroidSection2calibrated(2,:);
    GridCentroid2Minus3 = [GridCentroid2Minus3calibrated(1,1)+GridCenter(1,1),-GridCentroid2Minus3calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroid2Minus3, 20,'EdgeColor','y');
    foundCentroidsSorted(96,:) = GridCentroid2Minus3calibrated;
else
    GridCentroid3Minus2calibrated = GridCentroidSection2calibrated(2,:);
    GridCentroid3Minus2 = [GridCentroid3Minus2calibrated(1,1)+GridCenter(1,1),-GridCentroid3Minus2calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroid3Minus2, 20,'EdgeColor','g');
    foundCentroidsSorted(86,:) = GridCentroid3Minus2calibrated;
    GridCentroid2Minus3calibrated = GridCentroidSection2calibrated(1,:);
    GridCentroid2Minus3 = [GridCentroid2Minus3calibrated(1,1)+GridCenter(1,1),-GridCentroid2Minus3calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroid2Minus3, 20,'EdgeColor','y');
    foundCentroidsSorted(96,:) = GridCentroid2Minus3calibrated;
end

[GridCentroidMinus1Minus3calibratedX,index1] = max(GridCentroidSection3calibrated(:,1));
GridCentroidMinus1Minus3calibrated = GridCentroidSection3calibrated(index1,:);
GridCentroidMinus1Minus3 = [GridCentroidMinus1Minus3calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus1Minus3calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroidMinus1Minus3, 20,'EdgeColor','g');
foundCentroidsSorted(93,:) = GridCentroidMinus1Minus3calibrated;
GridCentroidSection3calibrated(index1,:) = [];

[GridCentroidMinus3Minus1calibratedY,index2] = max(GridCentroidSection3calibrated(:,2));
GridCentroidMinus3Minus1calibrated = GridCentroidSection3calibrated(index2,:);
GridCentroidMinus3Minus1 = [GridCentroidMinus3Minus1calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus3Minus1calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroidMinus3Minus1, 20,'EdgeColor','y');
foundCentroidsSorted(69,:) = GridCentroidMinus3Minus1calibrated;
GridCentroidSection3calibrated(index2,:) = [];

if GridCentroidSection3calibrated(1,1) < GridCentroidSection3calibrated(2,1) && ...
    GridCentroidSection3calibrated(1,2) > GridCentroidSection3calibrated(2,2)
    GridCentroidMinus3Minus2calibrated = GridCentroidSection3calibrated(1,:);
    GridCentroidMinus3Minus2 = [GridCentroidMinus3Minus2calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus3Minus2calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroidMinus3Minus2, 20,'EdgeColor','g');
    foundCentroidsSorted(80,:) = GridCentroidMinus3Minus2calibrated;
    GridCentroidMinus2Minus3calibrated = GridCentroidSection3calibrated(2,:);
    GridCentroidMinus2Minus3 = [GridCentroidMinus2Minus3calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus2Minus3calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroidMinus2Minus3, 20,'EdgeColor','y');
    foundCentroidsSorted(92,:) = GridCentroidMinus2Minus3calibrated;
else
    GridCentroidMinus3Minus2calibrated = GridCentroidSection3calibrated(2,:);
    GridCentroidMinus3Minus2 = [GridCentroidMinus3Minus2calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus3Minus2calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroidMinus3Minus2, 20,'EdgeColor','g');
    foundCentroidsSorted(80,:) = GridCentroidMinus3Minus2calibrated;
    GridCentroidMinus2Minus3calibrated = GridCentroidSection3calibrated(1,:);
    GridCentroidMinus2Minus3 = [GridCentroidMinus2Minus3calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus2Minus3calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroidMinus2Minus3, 20,'EdgeColor','y');
    foundCentroidsSorted(92,:) = GridCentroidMinus2Minus3calibrated;
end


[GridCentroidMinus13calibratedX,index1] = max(GridCentroidSection4calibrated(:,1));
GridCentroidMinus13calibrated = GridCentroidSection4calibrated(index1,:);
GridCentroidMinus13 = [GridCentroidMinus13calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus13calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroidMinus13, 20,'EdgeColor','g');
foundCentroidsSorted(27,:) = GridCentroidMinus13calibrated;
GridCentroidSection4calibrated(index1,:) = [];

[GridCentroidMinus31calibratedY,index2] = min(GridCentroidSection4calibrated(:,2));
GridCentroidMinus31calibrated = GridCentroidSection4calibrated(index2,:);
GridCentroidMinus31 = [GridCentroidMinus31calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus31calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroidMinus31, 20,'EdgeColor','y');
foundCentroidsSorted(47,:) = GridCentroidMinus31calibrated;
GridCentroidSection4calibrated(index2,:) = [];

if GridCentroidSection4calibrated(1,1) < GridCentroidSection4calibrated(2,1) && ...
    GridCentroidSection4calibrated(1,2) < GridCentroidSection4calibrated(2,2)
    GridCentroidMinus32calibrated = GridCentroidSection4calibrated(1,:);
    GridCentroidMinus32 = [GridCentroidMinus32calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus32calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroidMinus32, 20,'EdgeColor','g');
    foundCentroidsSorted(36,:) = GridCentroidMinus32calibrated;
    GridCentroidMinus23calibrated = GridCentroidSection4calibrated(2,:);
    GridCentroidMinus23 = [GridCentroidMinus23calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus23calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroidMinus23, 20,'EdgeColor','y');
    foundCentroidsSorted(26,:) = GridCentroidMinus23calibrated;
else
    GridCentroidMinus32calibrated = GridCentroidSection4calibrated(2,:);
    GridCentroidMinus32 = [GridCentroidMinus32calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus32calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroidMinus32, 20,'EdgeColor','g');
    foundCentroidsSorted(36,:) = GridCentroidMinus32calibrated;
    GridCentroidMinus23calibrated = GridCentroidSection4calibrated(1,:);
    GridCentroidMinus23 = [GridCentroidMinus23calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus23calibrated(1,2)+GridCenter(1,2)];
    % viscircles(GridCentroidMinus23, 20,'EdgeColor','y');
    foundCentroidsSorted(26,:) = GridCentroidMinus23calibrated;
end

%  2.10 find and sort 46-49 in the third circle
centersCalibrated6 = centersCalibrated5;
NeighborIndex5 = dsearchn(centersCalibrated6,foundCentroidsSorted(37,:));
GridCentroidMinus33calibrated = centersCalibrated6(NeighborIndex5,:);
GridCentroidMinus33 = [GridCentroidMinus33calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus33calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroidMinus33, 20,'EdgeColor','y');
foundCentroidsSorted(25,:) = GridCentroidMinus33calibrated;
centersCalibrated6(NeighborIndex5,:)=[];

NeighborIndex6 = dsearchn(centersCalibrated6,foundCentroidsSorted(41,:));
GridCentroid33calibrated = centersCalibrated6(NeighborIndex6,:);
GridCentroid33 = [GridCentroid33calibrated(1,1)+GridCenter(1,1),-GridCentroid33calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroid33, 20,'EdgeColor','y');
foundCentroidsSorted(31,:) = GridCentroid33calibrated;
centersCalibrated6(NeighborIndex6,:)=[];

NeighborIndex6 = dsearchn(centersCalibrated6,foundCentroidsSorted(85,:));
GridCentroid3Minus3calibrated = centersCalibrated6(NeighborIndex6,:);
GridCentroid3Minus3 = [GridCentroid3Minus3calibrated(1,1)+GridCenter(1,1),-GridCentroid3Minus3calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroid3Minus3, 20,'EdgeColor','y');
foundCentroidsSorted(97,:) = GridCentroid3Minus3calibrated;
centersCalibrated6(NeighborIndex6,:)=[];

NeighborIndex6 = dsearchn(centersCalibrated6,foundCentroidsSorted(81,:));
GridCentroidMinus3Minus3calibrated = centersCalibrated6(NeighborIndex6,:);
GridCentroidMinus3Minus3 = [GridCentroidMinus3Minus3calibrated(1,1)+GridCenter(1,1),-GridCentroidMinus3Minus3calibrated(1,2)+GridCenter(1,2)];
% viscircles(GridCentroidMinus3Minus3, 20,'EdgeColor','y');
foundCentroidsSorted(91,:) = GridCentroidMinus3Minus3calibrated;
centersCalibrated6(NeighborIndex6,:)=[];

% 2.11 find 45 degree diagonal 50 to 53
Diagonal45Tan = (foundCentroidsSorted(31,2)-foundCentroidsSorted(91,2))/(foundCentroidsSorted(31,1)-foundCentroidsSorted(91,1));
Grid45diagonalCalibrated = zeros(4,2);
m = 1;
for n  = 1:64
    if centersCalibrated6(n,2) < Diagonal45Tan * centersCalibrated6(n,1)+20 ...
            &&  centersCalibrated6(n,2) > Diagonal45Tan * centersCalibrated6(n,1)-20
       Grid45diagonalCalibrated(m,:) = centersCalibrated6(n,:);
       m = m+1;
    end
end
Grid45diagonal = [Grid45diagonalCalibrated(:,1)+GridCenter(1,1),-Grid45diagonalCalibrated(:,2)+GridCenter(1,2)];

%     viscircles(Grid45diagonal, 20,'EdgeColor','g');

Grid45diagonalX = Grid45diagonal(:,1);
Grid45diagonalY = Grid45diagonal(:,2);
[Grid45diagonalXSorted,index] = sort(Grid45diagonalX);
Grid45diagonalYSorted = Grid45diagonalY(index);
Grid45diagonalSorted = [Grid45diagonalXSorted,Grid45diagonalYSorted];
Grid45diagonalSortedCalibrated = [Grid45diagonalSorted(:,1)-GridCenter(1,1),-Grid45diagonalSorted(:,2)+GridCenter(1,2)];

foundCentroidsSorted(111,:) = Grid45diagonalSortedCalibrated(1,:);
foundCentroidsSorted(101,:) = Grid45diagonalSortedCalibrated(2,:);
foundCentroidsSorted(21,:) = Grid45diagonalSortedCalibrated(3,:);
foundCentroidsSorted(11,:) = Grid45diagonalSortedCalibrated(4,:);

centersCalibrated7 = centersCalibrated6;
for n = 1:4
    [indexFor45diagonal1,indexFor45diagonal2] = find(centersCalibrated7 == Grid45diagonalSortedCalibrated(n,:));
    centersCalibrated7(indexFor45diagonal1,:) = [];
end

% 2.12 find 135 degree diagonal 54 to 57
Diagonal135Tan = (foundCentroidsSorted(97,2)-foundCentroidsSorted(25,2))/(foundCentroidsSorted(97,1)-foundCentroidsSorted(25,1));
Grid135diagonalCalibrated = zeros(4,2);
m = 1;
for n  = 1:60
    if centersCalibrated7(n,2) < Diagonal135Tan * centersCalibrated7(n,1)+20 ...
            &&  centersCalibrated7(n,2) > Diagonal135Tan * centersCalibrated7(n,1)-20
       Grid135diagonalCalibrated(m,:) = centersCalibrated7(n,:);
       m = m+1;
    end
end
Grid135diagonal = [Grid135diagonalCalibrated(:,1)+GridCenter(1,1),-Grid135diagonalCalibrated(:,2)+GridCenter(1,2)];

%     viscircles(Grid135diagonal, 20,'EdgeColor','g');

Grid135diagonalX = Grid135diagonal(:,1);
Grid135diagonalY = Grid135diagonal(:,2);
[Grid135diagonalXSorted,index] = sort(Grid135diagonalX);
Grid135diagonalYSorted = Grid135diagonalY(index);
Grid135diagonalSorted = [Grid135diagonalXSorted,Grid135diagonalYSorted];
Grid135diagonalSortedCalibrated = [Grid135diagonalSorted(:,1)-GridCenter(1,1),-Grid135diagonalSorted(:,2)+GridCenter(1,2)];

foundCentroidsSorted(1,:) = Grid135diagonalSortedCalibrated(1,:);
foundCentroidsSorted(13,:) = Grid135diagonalSortedCalibrated(2,:);
foundCentroidsSorted(109,:) = Grid135diagonalSortedCalibrated(3,:);
foundCentroidsSorted(121,:) = Grid135diagonalSortedCalibrated(4,:);

centersCalibrated8 = centersCalibrated7;
for n = 1:4
    [indexFor135diagonal1,indexFor135diagonal2] = find(centersCalibrated8 == Grid135diagonalSortedCalibrated(n,:));
    centersCalibrated8(indexFor135diagonal1,:) = [];
end

% 2.13 find and sort 58-61 in the fourth circle first section
FourthCircle58to61NeighborCalibrated = zeros(4,2);
centersCalibrated9 = centersCalibrated8;
for n = 1:4
    NeighborIndex7 = dsearchn(centersCalibrated9,foundCentroidsSorted(31,:));
    FourthCircle58to61NeighborCalibrated(n,:) = centersCalibrated9(NeighborIndex7,:);
    centersCalibrated9(NeighborIndex7,:)=[];
end

FourthCircle58to61Neighbor = [FourthCircle58to61NeighborCalibrated(:,1)+GridCenter(1,1),-FourthCircle58to61NeighborCalibrated(:,2)+GridCenter(1,2)];

%     viscircles(FourthCircle58to61Neighbor, 20,'EdgeColor','g');

FourthCircle58to59NeighborCalibrated = zeros(2,2);
FourthCircle60to61NeighborCalibrated = zeros(2,2);
a = 1; 
b = 1;
for n = 1:4
    if FourthCircle58to61NeighborCalibrated(n,2)> Diagonal45Tan * FourthCircle58to61NeighborCalibrated(n,1)
        FourthCircle58to59NeighborCalibrated(a,:) = FourthCircle58to61NeighborCalibrated(n,:);
        a = a+1;
    else
        FourthCircle60to61NeighborCalibrated(b,:) = FourthCircle58to61NeighborCalibrated(n,:);
        b = b+1;
    end
end

FourthCircle58to59NeighborCalibratedX = FourthCircle58to59NeighborCalibrated(:,1);
FourthCircle58to59NeighborCalibratedY = FourthCircle58to59NeighborCalibrated(:,2);
[FourthCircle58to59NeighborCalibratedXSorted,index] = sort(FourthCircle58to59NeighborCalibratedX);
FourthCircle58to59NeighborCalibratedYSorted = FourthCircle58to59NeighborCalibratedY(index);
FourthCircle58to59NeighborCalibratedSorted = [FourthCircle58to59NeighborCalibratedXSorted,FourthCircle58to59NeighborCalibratedYSorted];
FourthCircle58to59NeighborSorted = [FourthCircle58to59NeighborCalibratedSorted(:,1)+GridCenter(1,1),-FourthCircle58to59NeighborCalibratedSorted(:,2)+GridCenter(1,2)];
foundCentroidsSorted(19,:) = FourthCircle58to59NeighborCalibratedSorted(1,:);
foundCentroidsSorted(20,:) = FourthCircle58to59NeighborCalibratedSorted(2,:);

 % viscircles(FourthCircle58to59NeighborSorted, 20,'EdgeColor','g');

FourthCircle60to61NeighborCalibratedX = FourthCircle60to61NeighborCalibrated(:,1);
FourthCircle60to61NeighborCalibratedY = FourthCircle60to61NeighborCalibrated(:,2);
[FourthCircle60to61NeighborCalibratedYSorted,index] = sort(FourthCircle60to61NeighborCalibratedY);
FourthCircle60to61NeighborCalibratedXSorted = FourthCircle60to61NeighborCalibratedX(index);
FourthCircle60to61NeighborCalibratedSorted = [FourthCircle60to61NeighborCalibratedXSorted,FourthCircle60to61NeighborCalibratedYSorted];
FourthCircle60to61NeighborSorted = [FourthCircle60to61NeighborCalibratedSorted(:,1)+GridCenter(1,1),-FourthCircle60to61NeighborCalibratedSorted(:,2)+GridCenter(1,2)];
foundCentroidsSorted(43,:) = FourthCircle60to61NeighborCalibratedSorted(1,:);
foundCentroidsSorted(32,:) = FourthCircle60to61NeighborCalibratedSorted(2,:);

 % viscircles(FourthCircle60to61NeighborSorted, 20,'EdgeColor','g');

% 2.14 find and sort 62-65 in the fourth circle second section
FourthCircle62to65NeighborCalibrated = zeros(4,2);
centersCalibrated10 = centersCalibrated9;
for n = 1:4
    NeighborIndex8 = dsearchn(centersCalibrated10,foundCentroidsSorted(97,:));
    FourthCircle62to65NeighborCalibrated(n,:) = centersCalibrated10(NeighborIndex8,:);
    centersCalibrated10(NeighborIndex8,:)=[];
end

FourthCircle62to65Neighbor = [FourthCircle62to65NeighborCalibrated(:,1)+GridCenter(1,1),-FourthCircle62to65NeighborCalibrated(:,2)+GridCenter(1,2)];

%     viscircles(FourthCircle62to65Neighbor, 20,'EdgeColor','g');

FourthCircle62to63NeighborCalibrated = zeros(2,2);
FourthCircle64to65NeighborCalibrated = zeros(2,2);
a = 1; 
b = 1;
for n = 1:4
    if FourthCircle62to65NeighborCalibrated(n,2)< Diagonal135Tan * FourthCircle62to65NeighborCalibrated(n,1)
        FourthCircle62to63NeighborCalibrated(a,:) = FourthCircle62to65NeighborCalibrated(n,:);
        a = a+1;
    else
        FourthCircle64to65NeighborCalibrated(b,:) = FourthCircle62to65NeighborCalibrated(n,:);
        b = b+1;
    end
end

FourthCircle62to63NeighborCalibratedX = FourthCircle62to63NeighborCalibrated(:,1);
FourthCircle62to63NeighborCalibratedY = FourthCircle62to63NeighborCalibrated(:,2);
[FourthCircle62to63NeighborCalibratedXSorted,index] = sort(FourthCircle62to63NeighborCalibratedX);
FourthCircle62to63NeighborCalibratedYSorted = FourthCircle62to63NeighborCalibratedY(index);
FourthCircle62to63NeighborCalibratedSorted = [FourthCircle62to63NeighborCalibratedXSorted,FourthCircle62to63NeighborCalibratedYSorted];
FourthCircle62to63NeighborSorted = [FourthCircle62to63NeighborCalibratedSorted(:,1)+GridCenter(1,1),-FourthCircle62to63NeighborCalibratedSorted(:,2)+GridCenter(1,2)];
foundCentroidsSorted(107,:) = FourthCircle62to63NeighborCalibratedSorted(1,:);
foundCentroidsSorted(108,:) = FourthCircle62to63NeighborCalibratedSorted(2,:);

 % viscircles(FourthCircle62to63NeighborSorted, 20,'EdgeColor','g');

FourthCircle64to65NeighborCalibratedX = FourthCircle64to65NeighborCalibrated(:,1);
FourthCircle64to65NeighborCalibratedY = FourthCircle64to65NeighborCalibrated(:,2);
[FourthCircle64to65NeighborCalibratedYSorted,index] = sort(FourthCircle64to65NeighborCalibratedY);
FourthCircle64to65NeighborCalibratedXSorted = FourthCircle64to65NeighborCalibratedX(index);
FourthCircle64to65NeighborCalibratedSorted = [FourthCircle64to65NeighborCalibratedXSorted,FourthCircle64to65NeighborCalibratedYSorted];
FourthCircle64to65NeighborSorted = [FourthCircle64to65NeighborCalibratedSorted(:,1)+GridCenter(1,1),-FourthCircle64to65NeighborCalibratedSorted(:,2)+GridCenter(1,2)];
foundCentroidsSorted(98,:) = FourthCircle64to65NeighborCalibratedSorted(1,:);
foundCentroidsSorted(87,:) = FourthCircle64to65NeighborCalibratedSorted(2,:);

 % viscircles(FourthCircle64to65NeighborSorted, 20,'EdgeColor','g');

 % 2.15 find and sort 66-69 in the fourth circle third section
FourthCircle66to69NeighborCalibrated = zeros(4,2);
centersCalibrated11 = centersCalibrated10;
for n = 1:4
    NeighborIndex8 = dsearchn(centersCalibrated11,foundCentroidsSorted(91,:));
    FourthCircle66to69NeighborCalibrated(n,:) = centersCalibrated11(NeighborIndex8,:);
    centersCalibrated11(NeighborIndex8,:)=[];
end

FourthCircle66to69Neighbor = [FourthCircle66to69NeighborCalibrated(:,1)+GridCenter(1,1),-FourthCircle66to69NeighborCalibrated(:,2)+GridCenter(1,2)];

%     viscircles(FourthCircle66to69Neighbor, 20,'EdgeColor','g');

FourthCircle66to67NeighborCalibrated = zeros(2,2);
FourthCircle68to69NeighborCalibrated = zeros(2,2);
a = 1; 
b = 1;
for n = 1:4
    if FourthCircle66to69NeighborCalibrated(n,2)< Diagonal45Tan * FourthCircle66to69NeighborCalibrated(n,1)
        FourthCircle66to67NeighborCalibrated(a,:) = FourthCircle66to69NeighborCalibrated(n,:);
        a = a+1;
    else
        FourthCircle68to69NeighborCalibrated(b,:) = FourthCircle66to69NeighborCalibrated(n,:);
        b = b+1;
    end
end

FourthCircle66to67NeighborCalibratedX = FourthCircle66to67NeighborCalibrated(:,1);
FourthCircle66to67NeighborCalibratedY = FourthCircle66to67NeighborCalibrated(:,2);
[FourthCircle66to67NeighborCalibratedXSorted,index] = sort(FourthCircle66to67NeighborCalibratedX);
FourthCircle66to67NeighborCalibratedYSorted = FourthCircle66to67NeighborCalibratedY(index);
FourthCircle66to67NeighborCalibratedSorted = [FourthCircle66to67NeighborCalibratedXSorted,FourthCircle66to67NeighborCalibratedYSorted];
FourthCircle66to67NeighborSorted = [FourthCircle66to67NeighborCalibratedSorted(:,1)+GridCenter(1,1),-FourthCircle66to67NeighborCalibratedSorted(:,2)+GridCenter(1,2)];
foundCentroidsSorted(102,:) = FourthCircle66to67NeighborCalibratedSorted(1,:);
foundCentroidsSorted(103,:) = FourthCircle66to67NeighborCalibratedSorted(2,:);

 % viscircles(FourthCircle66to67NeighborSorted, 20,'EdgeColor','g');
 

FourthCircle68to69NeighborCalibratedX = FourthCircle68to69NeighborCalibrated(:,1);
FourthCircle68to69NeighborCalibratedY = FourthCircle68to69NeighborCalibrated(:,2);
[FourthCircle68to69NeighborCalibratedYSorted,index] = sort(FourthCircle68to69NeighborCalibratedY);
FourthCircle68to69NeighborCalibratedXSorted = FourthCircle68to69NeighborCalibratedX(index);
FourthCircle68to69NeighborCalibratedSorted = [FourthCircle68to69NeighborCalibratedXSorted,FourthCircle68to69NeighborCalibratedYSorted];
FourthCircle68to69NeighborSorted = [FourthCircle68to69NeighborCalibratedSorted(:,1)+GridCenter(1,1),-FourthCircle68to69NeighborCalibratedSorted(:,2)+GridCenter(1,2)];
foundCentroidsSorted(90,:) = FourthCircle68to69NeighborCalibratedSorted(1,:);
foundCentroidsSorted(79,:) = FourthCircle68to69NeighborCalibratedSorted(2,:);

 % viscircles(FourthCircle68to69NeighborSorted, 20,'EdgeColor','g');


  % 2.16 find and sort 70-73 in the fourth circle fourth section
FourthCircle70to73NeighborCalibrated = zeros(4,2);
centersCalibrated12 = centersCalibrated11;
for n = 1:4
    NeighborIndex8 = dsearchn(centersCalibrated12,foundCentroidsSorted(25,:));
    FourthCircle70to73NeighborCalibrated(n,:) = centersCalibrated12(NeighborIndex8,:);
    centersCalibrated12(NeighborIndex8,:)=[];
end

FourthCircle70to73Neighbor = [FourthCircle70to73NeighborCalibrated(:,1)+GridCenter(1,1),-FourthCircle70to73NeighborCalibrated(:,2)+GridCenter(1,2)];

% for n = 1:4
%     viscircles(FourthCircle70to73Neighbor(n,:), 20,'EdgeColor','g');
% end

FourthCircle70to71NeighborCalibrated = zeros(2,2);
FourthCircle72to73NeighborCalibrated = zeros(2,2);
a = 1; 
b = 1;
for n = 1:4
    if FourthCircle70to73NeighborCalibrated(n,2)> Diagonal135Tan * FourthCircle70to73NeighborCalibrated(n,1)
        FourthCircle70to71NeighborCalibrated(a,:) = FourthCircle70to73NeighborCalibrated(n,:);
        a = a+1;
    else
        FourthCircle72to73NeighborCalibrated(b,:) = FourthCircle70to73NeighborCalibrated(n,:);
        b = b+1;
    end
end

FourthCircle70to71NeighborCalibratedX = FourthCircle70to71NeighborCalibrated(:,1);
FourthCircle70to71NeighborCalibratedY = FourthCircle70to71NeighborCalibrated(:,2);
[FourthCircle70to71NeighborCalibratedXSorted,index] = sort(FourthCircle70to71NeighborCalibratedX);
FourthCircle70to71NeighborCalibratedYSorted = FourthCircle70to71NeighborCalibratedY(index);
FourthCircle70to71NeighborCalibratedSorted = [FourthCircle70to71NeighborCalibratedXSorted,FourthCircle70to71NeighborCalibratedYSorted];
FourthCircle70to71NeighborSorted = [FourthCircle70to71NeighborCalibratedSorted(:,1)+GridCenter(1,1),-FourthCircle70to71NeighborCalibratedSorted(:,2)+GridCenter(1,2)];
foundCentroidsSorted(14,:) = FourthCircle70to71NeighborCalibratedSorted(1,:);
foundCentroidsSorted(15,:) = FourthCircle70to71NeighborCalibratedSorted(2,:);

 % viscircles(FourthCircle70to71NeighborSorted(1,:), 20,'EdgeColor','g');
 % viscircles(FourthCircle70to71NeighborSorted(2,:), 20,'EdgeColor','y');

FourthCircle72to73NeighborCalibratedX = FourthCircle72to73NeighborCalibrated(:,1);
FourthCircle72to73NeighborCalibratedY = FourthCircle72to73NeighborCalibrated(:,2);
[FourthCircle72to73NeighborCalibratedYSorted,index] = sort(FourthCircle72to73NeighborCalibratedY);
FourthCircle72to73NeighborCalibratedXSorted = FourthCircle72to73NeighborCalibratedX(index);
FourthCircle72to73NeighborCalibratedSorted = [FourthCircle72to73NeighborCalibratedXSorted,FourthCircle72to73NeighborCalibratedYSorted];
FourthCircle72to73NeighborSorted = [FourthCircle72to73NeighborCalibratedSorted(:,1)+GridCenter(1,1),-FourthCircle72to73NeighborCalibratedSorted(:,2)+GridCenter(1,2)];
foundCentroidsSorted(35,:) = FourthCircle72to73NeighborCalibratedSorted(1,:);
foundCentroidsSorted(24,:) = FourthCircle72to73NeighborCalibratedSorted(2,:);

 % viscircles(FourthCircle72to73NeighborSorted(1,:), 20,'EdgeColor','g');
 % viscircles(FourthCircle72to73NeighborSorted(2,:), 20,'EdgeColor','y');

   % 2.17 find and sort 74-75 in the fourth circle first section
FourthCircle74to75NeighborCalibrated = zeros(2,2);
centersCalibrated13 = centersCalibrated12;
for n = 1:2
    NeighborIndex9 = dsearchn(centersCalibrated13,foundCentroidsSorted(41,:));
    FourthCircle74to75NeighborCalibrated(n,:) = centersCalibrated13(NeighborIndex9,:);
    centersCalibrated13(NeighborIndex9,:)=[];
end
FourthCircle74to75Neighbor = [FourthCircle74to75NeighborCalibrated(:,1)+GridCenter(1,1),-FourthCircle74to75NeighborCalibrated(:,2)+GridCenter(1,2)];
% for n = 1:2
%     viscircles(FourthCircle74to75Neighbor(n,:), 20,'EdgeColor','g');
% end
FourthCircle74to75NeighborCalibratedX = FourthCircle74to75NeighborCalibrated(:,1);
FourthCircle74to75NeighborCalibratedY = FourthCircle74to75NeighborCalibrated(:,2);
[FourthCircle74to75NeighborCalibratedYSorted,index] = sort(FourthCircle74to75NeighborCalibratedY);
FourthCircle74to75NeighborCalibratedXSorted = FourthCircle74to75NeighborCalibratedX(index);
FourthCircle74to75NeighborCalibratedSorted = [FourthCircle74to75NeighborCalibratedXSorted,FourthCircle74to75NeighborCalibratedYSorted];
FourthCircle74to75NeighborSorted = [FourthCircle74to75NeighborCalibratedSorted(:,1)+GridCenter(1,1),-FourthCircle74to75NeighborCalibratedSorted(:,2)+GridCenter(1,2)];
foundCentroidsSorted(54,:) = FourthCircle74to75NeighborCalibratedSorted(1,:);
foundCentroidsSorted(18,:) = FourthCircle74to75NeighborCalibratedSorted(2,:);

   % 2.18 find and sort 76-77 in the fourth circle second section
FourthCircle76to77NeighborCalibrated = zeros(2,2);
centersCalibrated14 = centersCalibrated13;
for n = 1:2
    NeighborIndex9 = dsearchn(centersCalibrated14,foundCentroidsSorted(85,:));
    FourthCircle76to77NeighborCalibrated(n,:) = centersCalibrated14(NeighborIndex9,:);
    centersCalibrated14(NeighborIndex9,:)=[];
end
FourthCircle76to77Neighbor = [FourthCircle76to77NeighborCalibrated(:,1)+GridCenter(1,1),-FourthCircle76to77NeighborCalibrated(:,2)+GridCenter(1,2)];
% for n = 1:2
%     viscircles(FourthCircle76to77Neighbor(n,:), 20,'EdgeColor','g');
% end
FourthCircle76to77NeighborCalibratedX = FourthCircle76to77NeighborCalibrated(:,1);
FourthCircle76to77NeighborCalibratedY = FourthCircle76to77NeighborCalibrated(:,2);
[FourthCircle76to77NeighborCalibratedYSorted,index] = sort(FourthCircle76to77NeighborCalibratedY);
FourthCircle76to77NeighborCalibratedXSorted = FourthCircle76to77NeighborCalibratedX(index);
FourthCircle76to77NeighborCalibratedSorted = [FourthCircle76to77NeighborCalibratedXSorted,FourthCircle76to77NeighborCalibratedYSorted];
FourthCircle76to77NeighborSorted = [FourthCircle76to77NeighborCalibratedSorted(:,1)+GridCenter(1,1),-FourthCircle76to77NeighborCalibratedSorted(:,2)+GridCenter(1,2)];
foundCentroidsSorted(106,:) = FourthCircle76to77NeighborCalibratedSorted(1,:);
foundCentroidsSorted(76,:) = FourthCircle76to77NeighborCalibratedSorted(2,:);

  % 2.19 find and sort 78-79 in the fourth circle third section
FourthCircle78to79NeighborCalibrated = zeros(2,2);
for n = 1:2
    NeighborIndex9 = dsearchn(centersCalibrated14,foundCentroidsSorted(81,:));
    FourthCircle78to79NeighborCalibrated(n,:) = centersCalibrated14(NeighborIndex9,:);
    centersCalibrated14(NeighborIndex9,:)=[];
end
FourthCircle78to79Neighbor = [FourthCircle78to79NeighborCalibrated(:,1)+GridCenter(1,1),-FourthCircle78to79NeighborCalibrated(:,2)+GridCenter(1,2)];
% for n = 1:2
%     viscircles(FourthCircle78to79Neighbor(n,:), 20,'EdgeColor','g');
% end
FourthCircle78to79NeighborCalibratedX = FourthCircle78to79NeighborCalibrated(:,1);
FourthCircle78to79NeighborCalibratedY = FourthCircle78to79NeighborCalibrated(:,2);
[FourthCircle78to79NeighborCalibratedYSorted,index] = sort(FourthCircle78to79NeighborCalibratedY);
FourthCircle78to79NeighborCalibratedXSorted = FourthCircle78to79NeighborCalibratedX(index);
FourthCircle78to79NeighborCalibratedSorted = [FourthCircle78to79NeighborCalibratedXSorted,FourthCircle78to79NeighborCalibratedYSorted];
FourthCircle78to79NeighborSorted = [FourthCircle78to79NeighborCalibratedSorted(:,1)+GridCenter(1,1),-FourthCircle78to79NeighborCalibratedSorted(:,2)+GridCenter(1,2)];
foundCentroidsSorted(104,:) = FourthCircle78to79NeighborCalibratedSorted(1,:);
foundCentroidsSorted(68,:) = FourthCircle78to79NeighborCalibratedSorted(2,:);

  % 2.20 find and sort 80-81 in the fourth circle fourth section
FourthCircle80to81NeighborCalibrated = zeros(2,2);
for n = 1:2
    NeighborIndex9 = dsearchn(centersCalibrated14,foundCentroidsSorted(37,:));
    FourthCircle80to81NeighborCalibrated(n,:) = centersCalibrated14(NeighborIndex9,:);
    centersCalibrated14(NeighborIndex9,:)=[];
end
FourthCircle80to81Neighbor = [FourthCircle80to81NeighborCalibrated(:,1)+GridCenter(1,1),-FourthCircle80to81NeighborCalibrated(:,2)+GridCenter(1,2)];
% for n = 1:2
%     viscircles(FourthCircle80to81Neighbor(n,:), 20,'EdgeColor','g');
% end
FourthCircle80to81NeighborCalibratedX = FourthCircle80to81NeighborCalibrated(:,1);
FourthCircle80to81NeighborCalibratedY = FourthCircle80to81NeighborCalibrated(:,2);
[FourthCircle80to81NeighborCalibratedYSorted,index] = sort(FourthCircle80to81NeighborCalibratedY);
FourthCircle80to81NeighborCalibratedXSorted = FourthCircle80to81NeighborCalibratedX(index);
FourthCircle80to81NeighborCalibratedSorted = [FourthCircle80to81NeighborCalibratedXSorted,FourthCircle80to81NeighborCalibratedYSorted];
FourthCircle80to81NeighborSorted = [FourthCircle80to81NeighborCalibratedSorted(:,1)+GridCenter(1,1),-FourthCircle80to81NeighborCalibratedSorted(:,2)+GridCenter(1,2)];
foundCentroidsSorted(46,:) = FourthCircle80to81NeighborCalibratedSorted(1,:);
foundCentroidsSorted(16,:) = FourthCircle80to81NeighborCalibratedSorted(2,:);

% centers14 = [centersCalibrated14(:,1)+GridCenter(1,1),-centersCalibrated14(:,2)+GridCenter(1,2)];
% viscircles(centers14, 20,'EdgeColor','g');

% 2.21 sort the last 32 in the fifth circle first section
FifthCircleFirstSectionCalibrated = zeros(8,2); a = 1;
FifthCircleSecondSectionCalibrated = zeros(8,2); b = 1;
FifthCircleThirdSectionCalibrated = zeros(8,2); c = 1;
FifthCircleFourthSectionCalibrated = zeros(8,2); d = 1; 
for n = 1:32
    if centersCalibrated14(n,1)>0 
        if centersCalibrated14(n,2)>0
        FifthCircleFirstSectionCalibrated(a,:) = centersCalibrated14(n,:);
        a = a+1;
        else 
            FifthCircleSecondSectionCalibrated(b,:) = centersCalibrated14(n,:);
            b = b+1;
        end
    else 
        if centersCalibrated14(n,2)>0
        FifthCircleFourthSectionCalibrated(d,:) = centersCalibrated14(n,:);
        d = d+1;
        else 
            FifthCircleThirdSectionCalibrated(c,:) = centersCalibrated14(n,:);
            c = c+1;
        end
    end
end
FifthCircleFirstSection = [FifthCircleFirstSectionCalibrated(:,1)+GridCenter(1,1),-FifthCircleFirstSectionCalibrated(:,2)+GridCenter(1,2)];
% viscircles(FifthCircleFirstSection, 20,'EdgeColor','g');

FifthCircleSecondSection = [FifthCircleSecondSectionCalibrated(:,1)+GridCenter(1,1),-FifthCircleSecondSectionCalibrated(:,2)+GridCenter(1,2)];
% viscircles(FifthCircleSecondSection, 20,'EdgeColor','g');

FifthCircleThirdSection = [FifthCircleThirdSectionCalibrated(:,1)+GridCenter(1,1),-FifthCircleThirdSectionCalibrated(:,2)+GridCenter(1,2)];
% viscircles(FifthCircleThirdSection, 20,'EdgeColor','g');

FifthCircleFourthSection = [FifthCircleFourthSectionCalibrated(:,1)+GridCenter(1,1),-FifthCircleFourthSectionCalibrated(:,2)+GridCenter(1,2)];
% viscircles(FifthCircleFourthSection, 20,'EdgeColor','g');

FifthCircleFirstSectionCalibratedHor = zeros(4,2);
FifthCircleFirstSectionCalibratedVerti = zeros(4,2);
a = 1; 
b = 1;
for n = 1:8
    if FifthCircleFirstSectionCalibrated(n,2)> Diagonal45Tan * FifthCircleFirstSectionCalibrated(n,1)
        FifthCircleFirstSectionCalibratedHor(a,:) = FifthCircleFirstSectionCalibrated(n,:);
        a = a+1;
    else
        FifthCircleFirstSectionCalibratedVerti(b,:) = FifthCircleFirstSectionCalibrated(n,:);
        b = b+1;
    end
end

FifthCircleFirstSectionCalibratedHorX = FifthCircleFirstSectionCalibratedHor(:,1);
FifthCircleFirstSectionCalibratedHorY = FifthCircleFirstSectionCalibratedHor(:,2);
[FifthCircleFirstSectionCalibratedHorXSorted,index] = sort(FifthCircleFirstSectionCalibratedHorX);
FifthCircleFirstSectionCalibratedHorYSorted = FifthCircleFirstSectionCalibratedHorY(index);
FifthCircleFirstSectionCalibratedHorSorted = [FifthCircleFirstSectionCalibratedHorXSorted,FifthCircleFirstSectionCalibratedHorYSorted];
FifthCircleFirstSectionHorSorted = [FifthCircleFirstSectionCalibratedHorSorted(:,1)+GridCenter(1,1),-FifthCircleFirstSectionCalibratedHorSorted(:,2)+GridCenter(1,2)];
m = 1;
for n = 7:10
foundCentroidsSorted(n,:) = FifthCircleFirstSectionCalibratedHorSorted(m,:);
m = m+1;
end

FifthCircleFirstSectionCalibratedVertiX = FifthCircleFirstSectionCalibratedVerti(:,1);
FifthCircleFirstSectionCalibratedVertiY = FifthCircleFirstSectionCalibratedVerti(:,2);
[FifthCircleFirstSectionCalibratedVertiYSorted,index] = sort(FifthCircleFirstSectionCalibratedVertiY,'descend');
FifthCircleFirstSectionCalibratedVertiXSorted = FifthCircleFirstSectionCalibratedVertiX(index);
FifthCircleFirstSectionCalibratedVertiSorted = [FifthCircleFirstSectionCalibratedVertiXSorted,FifthCircleFirstSectionCalibratedVertiYSorted];
FifthCircleFirstSectionVertiSorted = [FifthCircleFirstSectionCalibratedVertiSorted(:,1)+GridCenter(1,1),-FifthCircleFirstSectionCalibratedVertiSorted(:,2)+GridCenter(1,2)];
m = 1;
for n = 1:4
foundCentroidsSorted(11*(n+1),:) = FifthCircleFirstSectionCalibratedVertiSorted(m,:);
m = m+1;
end

% viscircles(FifthCircleFirstSectionVertiSorted, 20,'EdgeColor','g');

% 2.22 sort the last 32 in the fifth circle second section
FifthCircleSecondSectionCalibratedHor = zeros(4,2);
FifthCircleSecondSectionCalibratedVerti = zeros(4,2);
a = 1; 
b = 1;
for n = 1:8
    if FifthCircleSecondSectionCalibrated(n,2)< Diagonal135Tan * FifthCircleSecondSectionCalibrated(n,1)
        FifthCircleSecondSectionCalibratedHor(a,:) = FifthCircleSecondSectionCalibrated(n,:);
        a = a+1;
    else
        FifthCircleSecondSectionCalibratedVerti(b,:) = FifthCircleSecondSectionCalibrated(n,:);
        b = b+1;
    end
end

FifthCircleSecondSectionCalibratedHorX = FifthCircleSecondSectionCalibratedHor(:,1);
FifthCircleSecondSectionCalibratedHorY = FifthCircleSecondSectionCalibratedHor(:,2);
[FifthCircleSecondSectionCalibratedHorXSorted,index] = sort(FifthCircleSecondSectionCalibratedHorX);
FifthCircleSecondSectionCalibratedHorYSorted = FifthCircleSecondSectionCalibratedHorY(index);
FifthCircleSecondSectionCalibratedHorSorted = [FifthCircleSecondSectionCalibratedHorXSorted,FifthCircleSecondSectionCalibratedHorYSorted];
FifthCircleSecondSectionHorSorted = [FifthCircleSecondSectionCalibratedHorSorted(:,1)+GridCenter(1,1),-FifthCircleSecondSectionCalibratedHorSorted(:,2)+GridCenter(1,2)];
m = 1;
for n = 117:120
foundCentroidsSorted(n,:) = FifthCircleSecondSectionCalibratedHorSorted(m,:);
m = m+1;
end

FifthCircleSecondSectionCalibratedVertiX = FifthCircleSecondSectionCalibratedVerti(:,1);
FifthCircleSecondSectionCalibratedVertiY = FifthCircleSecondSectionCalibratedVerti(:,2);
[FifthCircleSecondSectionCalibratedVertiYSorted,index] = sort(FifthCircleSecondSectionCalibratedVertiY,'descend');
FifthCircleSecondSectionCalibratedVertiXSorted = FifthCircleSecondSectionCalibratedVertiX(index);
FifthCircleSecondSectionCalibratedVertiSorted = [FifthCircleSecondSectionCalibratedVertiXSorted,FifthCircleSecondSectionCalibratedVertiYSorted];
FifthCircleSecondSectionVertiSorted = [FifthCircleSecondSectionCalibratedVertiSorted(:,1)+GridCenter(1,1),-FifthCircleSecondSectionCalibratedVertiSorted(:,2)+GridCenter(1,2)];
m = 1;
for n = 1:4
foundCentroidsSorted(11*(n+6),:) = FifthCircleSecondSectionCalibratedVertiSorted(m,:);
m = m+1;
end

% viscircles(FifthCircleSecondSectionVertiSorted, 20,'EdgeColor','g');

% 2.23 sort the last 32 in the fifth circle third section
FifthCircleThirdSectionCalibratedHor = zeros(4,2);
FifthCircleThirdSectionCalibratedVerti = zeros(4,2);
a = 1; 
b = 1;
for n = 1:8
    if FifthCircleThirdSectionCalibrated(n,2)< Diagonal45Tan * FifthCircleThirdSectionCalibrated(n,1)
        FifthCircleThirdSectionCalibratedHor(a,:) = FifthCircleThirdSectionCalibrated(n,:);
        a = a+1;
    else
        FifthCircleThirdSectionCalibratedVerti(b,:) = FifthCircleThirdSectionCalibrated(n,:);
        b = b+1;
    end
end

FifthCircleThirdSectionCalibratedHorX = FifthCircleThirdSectionCalibratedHor(:,1);
FifthCircleThirdSectionCalibratedHorY = FifthCircleThirdSectionCalibratedHor(:,2);
[FifthCircleThirdSectionCalibratedHorXSorted,index] = sort(FifthCircleThirdSectionCalibratedHorX);
FifthCircleThirdSectionCalibratedHorYSorted = FifthCircleThirdSectionCalibratedHorY(index);
FifthCircleThirdSectionCalibratedHorSorted = [FifthCircleThirdSectionCalibratedHorXSorted,FifthCircleThirdSectionCalibratedHorYSorted];
FifthCircleThirdSectionHorSorted = [FifthCircleThirdSectionCalibratedHorSorted(:,1)+GridCenter(1,1),-FifthCircleThirdSectionCalibratedHorSorted(:,2)+GridCenter(1,2)];
m = 1;
for n = 112:115
foundCentroidsSorted(n,:) = FifthCircleThirdSectionCalibratedHorSorted(m,:);
m = m+1;
end

FifthCircleThirdSectionCalibratedVertiX = FifthCircleThirdSectionCalibratedVerti(:,1);
FifthCircleThirdSectionCalibratedVertiY = FifthCircleThirdSectionCalibratedVerti(:,2);
[FifthCircleThirdSectionCalibratedVertiYSorted,index] = sort(FifthCircleThirdSectionCalibratedVertiY,'descend');
FifthCircleThirdSectionCalibratedVertiXSorted = FifthCircleThirdSectionCalibratedVertiX(index);
FifthCircleThirdSectionCalibratedVertiSorted = [FifthCircleThirdSectionCalibratedVertiXSorted,FifthCircleThirdSectionCalibratedVertiYSorted];
FifthCircleThirdSectionVertiSorted = [FifthCircleThirdSectionCalibratedVertiSorted(:,1)+GridCenter(1,1),-FifthCircleThirdSectionCalibratedVertiSorted(:,2)+GridCenter(1,2)];
m = 1;
for n = 1:4
foundCentroidsSorted(11*(n+5)+1,:) = FifthCircleThirdSectionCalibratedVertiSorted(m,:);
m = m+1;
end

% viscircles(FifthCircleThirdSectionVertiSorted, 20,'EdgeColor','g');

% 2.24 sort the last 32 in the fifth circle fourth section
FifthCircleFourthSectionCalibratedHor = zeros(4,2);
FifthCircleFourthSectionCalibratedVerti = zeros(4,2);
a = 1; 
b = 1;
for n = 1:8
    if FifthCircleFourthSectionCalibrated(n,2)> Diagonal135Tan * FifthCircleFourthSectionCalibrated(n,1)
        FifthCircleFourthSectionCalibratedHor(a,:) = FifthCircleFourthSectionCalibrated(n,:);
        a = a+1;
    else
        FifthCircleFourthSectionCalibratedVerti(b,:) = FifthCircleFourthSectionCalibrated(n,:);
        b = b+1;
    end
end

FifthCircleFourthSectionCalibratedHorX = FifthCircleFourthSectionCalibratedHor(:,1);
FifthCircleFourthSectionCalibratedHorY = FifthCircleFourthSectionCalibratedHor(:,2);
[FifthCircleFourthSectionCalibratedHorXSorted,index] = sort(FifthCircleFourthSectionCalibratedHorX);
FifthCircleFourthSectionCalibratedHorYSorted = FifthCircleFourthSectionCalibratedHorY(index);
FifthCircleFourthSectionCalibratedHorSorted = [FifthCircleFourthSectionCalibratedHorXSorted,FifthCircleFourthSectionCalibratedHorYSorted];
FifthCircleFourthSectionHorSorted = [FifthCircleFourthSectionCalibratedHorSorted(:,1)+GridCenter(1,1),-FifthCircleFourthSectionCalibratedHorSorted(:,2)+GridCenter(1,2)];
m = 1;
for n = 2:5
foundCentroidsSorted(n,:) = FifthCircleFourthSectionCalibratedHorSorted(m,:);
m = m+1;
end

FifthCircleFourthSectionCalibratedVertiX = FifthCircleFourthSectionCalibratedVerti(:,1);
FifthCircleFourthSectionCalibratedVertiY = FifthCircleFourthSectionCalibratedVerti(:,2);
[FifthCircleFourthSectionCalibratedVertiYSorted,index] = sort(FifthCircleFourthSectionCalibratedVertiY,'descend');
FifthCircleFourthSectionCalibratedVertiXSorted = FifthCircleFourthSectionCalibratedVertiX(index);
FifthCircleFourthSectionCalibratedVertiSorted = [FifthCircleFourthSectionCalibratedVertiXSorted,FifthCircleFourthSectionCalibratedVertiYSorted];
FifthCircleFourthSectionVertiSorted = [FifthCircleFourthSectionCalibratedVertiSorted(:,1)+GridCenter(1,1),-FifthCircleFourthSectionCalibratedVertiSorted(:,2)+GridCenter(1,2)];
m = 1;
for n = 1:4
foundCentroidsSorted(11*n+1,:) = FifthCircleFourthSectionCalibratedVertiSorted(m,:);
m = m+1;
end

% viscircles(FifthCircleFourthSectionVertiSorted, 20,'EdgeColor','g');

for n = 1:121
    text(foundCentroidsSorted(n,1)+GridCenter(1,1)-(radiiMean),-foundCentroidsSorted(n,2)+GridCenter(1,2),num2str(n),'Fontsize',8,'FontWeight','bold');
end

figure;
plot(foundCentroidsSorted(:,1),foundCentroidsSorted(:,2),'*');
xlim([-400 400]);
ylim([-400 400]);
title('Original data');

referenceXspacing = foundCentroidsSorted(62,1)-foundCentroidsSorted(60,1);
referenceYspacing = foundCentroidsSorted(50,2)-foundCentroidsSorted(72,2);
referenceSpacing = (referenceXspacing + referenceYspacing)/4;
hold on;
xRef1 = [-5*referenceSpacing,-4*referenceSpacing,-3*referenceSpacing,-2*referenceSpacing,-1*referenceSpacing,0,referenceSpacing,2*referenceSpacing,3*referenceSpacing,4*referenceSpacing,5*referenceSpacing];
yRef1 = [-5*referenceSpacing,-4*referenceSpacing,-3*referenceSpacing,-2*referenceSpacing,-1*referenceSpacing,0,referenceSpacing,2*referenceSpacing,3*referenceSpacing,4*referenceSpacing,5*referenceSpacing];
[XRef1,YRef1] = meshgrid(xRef1,yRef1);
scatter(XRef1,YRef1,'*','r');
legend('found points','reference points');
%% 3. keystone Elimination
% 3.1 y axis keystone
    % 3.1.1 move the center to upper left corner
y1keystone = foundCentroidsSorted(1,2)-foundCentroidsSorted(111,2);
y2keystone = foundCentroidsSorted(11,2)-foundCentroidsSorted(121,2);
y1keystoneHalf = (foundCentroidsSorted(1,2)-foundCentroidsSorted(111,2))/2;
y2keystoneHalf = (foundCentroidsSorted(11,2)-foundCentroidsSorted(121,2))/2;

if y1keystone>y2keystone
    movingPoints1 = [foundCentroidsSorted(1,:);foundCentroidsSorted(11,:);foundCentroidsSorted(111,:);foundCentroidsSorted(121,:);0,0];
    fixedPoints1 = [-y1keystoneHalf,y1keystoneHalf;y1keystoneHalf,y1keystoneHalf;-y1keystoneHalf,(-y1keystoneHalf);y1keystoneHalf,(-y1keystoneHalf);0,0];
else
    movingPoints1 = [foundCentroidsSorted(1,:);foundCentroidsSorted(11,:);foundCentroidsSorted(111,:);foundCentroidsSorted(121,:);0,0];
    fixedPoints1 = [-y2keystoneHalf,y2keystoneHalf;y2keystoneHalf,y2keystoneHalf;-y2keystoneHalf,(-y2keystoneHalf);y2keystoneHalf,(-y2keystoneHalf);0,0];
end
tform1 = fitgeotform2d(movingPoints1,fixedPoints1,"projective");
[x2,y2] = transformPointsForward(tform1,foundCentroidsSorted(:,1),foundCentroidsSorted(:,2));
figure;
scatter(x2,y2,'*');
xlim([-400 400]);
ylim([-400 400]);

referenceXspacing = x2(62)-x2(60);
referenceYspacing = y2(50)-y2(72);
referenceSpacing = (referenceXspacing + referenceYspacing)/4;
hold on;
xRef2 = [-5*referenceSpacing,-4*referenceSpacing,-3*referenceSpacing,-2*referenceSpacing,-1*referenceSpacing,0,referenceSpacing,2*referenceSpacing,3*referenceSpacing,4*referenceSpacing,5*referenceSpacing];
yRef2 = [-5*referenceSpacing,-4*referenceSpacing,-3*referenceSpacing,-2*referenceSpacing,-1*referenceSpacing,0,referenceSpacing,2*referenceSpacing,3*referenceSpacing,4*referenceSpacing,5*referenceSpacing];
[XRef2,YRef2] = meshgrid(xRef2,yRef2);
scatter(XRef2,YRef2,'*','r');
legend('found points','reference points');



x3 = x2-x2(61);
y3 = y2-y2(61);
figure;
plot(x3,y3,'*');
xlim([-400 400]);
ylim([-400 400]);
title('Data after first Keystone Elimination');

referenceXspacing = x3(62)-x3(60);
referenceYspacing = y3(50)-y3(72);
referenceSpacing = (referenceXspacing + referenceYspacing)/4;
hold on;
xRef2 = [-5*referenceSpacing,-4*referenceSpacing,-3*referenceSpacing,-2*referenceSpacing,-1*referenceSpacing,0,referenceSpacing,2*referenceSpacing,3*referenceSpacing,4*referenceSpacing,5*referenceSpacing];
yRef2 = [-5*referenceSpacing,-4*referenceSpacing,-3*referenceSpacing,-2*referenceSpacing,-1*referenceSpacing,0,referenceSpacing,2*referenceSpacing,3*referenceSpacing,4*referenceSpacing,5*referenceSpacing];
[XRef2,YRef2] = meshgrid(xRef2,yRef2);
scatter(XRef2,YRef2,'*','r');
legend('found points','reference points');

% 3.2 x axis keystone
x1keystone = (x3(11)-x3(1));
x2keystone = (x3(121)-x3(111));
x1keystoneHalf = (x3(11)-x3(1))/2;
x2keystoneHalf = (x3(121)-x3(111))/2;
if x1keystone>x2keystone
    movingPoints2 = [x3(1),y3(1);x3(11),y3(11);x3(111),y3(111);x3(121),y3(121)];
    fixedPoints2 = [-x1keystoneHalf,x1keystoneHalf;x1keystoneHalf,x1keystoneHalf;-x1keystoneHalf,(-x1keystoneHalf);x1keystoneHalf,(-x1keystoneHalf)];
else
    movingPoints2 = [x3(1),y3(1);x3(11),y3(11);x3(111),y3(111);x3(121),y3(121)];
    fixedPoints2 = [-x2keystoneHalf,x2keystoneHalf;x2keystoneHalf,x2keystoneHalf;-x2keystoneHalf,(-x2keystoneHalf);x2keystoneHalf,(-x2keystoneHalf)];
end
tform2 = fitgeotform2d(movingPoints2,fixedPoints2,"projective");
[x4,y4] = transformPointsForward(tform2,x3,y3);
% figure;
% scatter(x2,y2,'*','r');
% xlim([-400 400]);
% ylim([-400 400]);
x5 = x4-x4(61);
y5 = y4-y4(61);
figure;
plot(x5,y5,'*');
xlim([-400 400]);
ylim([-400 400]);
title('Data after second Keystone Elimination');
referenceXspacing = x5(62)-x5(60);
referenceYspacing = y5(50)-y5(72);
referenceSpacing = (referenceXspacing + referenceYspacing)/4;
hold on;
xRef3 = [-5*referenceSpacing,-4*referenceSpacing,-3*referenceSpacing,-2*referenceSpacing,-1*referenceSpacing,0,referenceSpacing,2*referenceSpacing,3*referenceSpacing,4*referenceSpacing,5*referenceSpacing];
yRef3 = [-5*referenceSpacing,-4*referenceSpacing,-3*referenceSpacing,-2*referenceSpacing,-1*referenceSpacing,0,referenceSpacing,2*referenceSpacing,3*referenceSpacing,4*referenceSpacing,5*referenceSpacing];
[XRef3,YRef3] = meshgrid(xRef3,yRef3);
scatter(XRef3,YRef3,'*','r');
legend('found points','reference points');


%% Calculation
% 5.1 
RefDia = sqrt((-5*referenceSpacing)^2+(5*referenceSpacing)^2);
RadDistortion1 = (sqrt(x5(1)^2+y5(1)^2)-RefDia)*100/RefDia;
RadDistortion2 = (sqrt(x5(11)^2+y5(11)^2)-RefDia)*100/RefDia;
RadDistortion3 = (sqrt(x5(111)^2+y5(111)^2)-RefDia)*100/RefDia;
RadDistortion4 = (sqrt(x5(121)^2+y5(121)^2)-RefDia)*100/RefDia;

RadDistortion = (RadDistortion1+RadDistortion2+RadDistortion3+RadDistortion4)/4
