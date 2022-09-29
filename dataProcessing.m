% Clear workspace
clearvars;
close all;
clc;

% Set folders
Folder.toolbox      = 'C:\Users\moissene\OneDrive - unige.ch\2022 - ROBOSHOULDER2\Matlab\';
Folder.data         = 'C:\Users\moissene\OneDrive - unige.ch\2022 - ROBOSHOULDER2\Matlab\data\';
Folder.dependencies = [Folder.toolbox,'dependencies\'];
addpath(Folder.toolbox);
addpath(genpath(Folder.dependencies));

% Load data
cd(Folder.data);
load('test.mat');
pointCloud3D = test(~isnan(test(:,1)),:)'*1e3; % mm

% Fit a plane to the point cloud at the least square sense
fobjPlane = planarFit(pointCloud3D);
figure(1);
[hL,hD]   = plot(fobjPlane);
hold on; axis equal; grid on; box on;
title('3D view');
plot3(pointCloud3D(1,:),pointCloud3D(2,:),pointCloud3D(3,:),'Marker','.','Markersize',20,'Linestyle','none','Color','red');

% Define a plane coordinate system
Oplane = mean(pointCloud3D,2);
Zplane = fobjPlane.normal';
Yplane = cross([1 0 0]',Zplane);
Xplane = cross(Yplane,Zplane);
Rplane = [Xplane Yplane Zplane];
quiver3(Oplane(1),Oplane(2),Oplane(3),Xplane(1),Xplane(2),Xplane(3),10,'red');
quiver3(Oplane(1),Oplane(2),Oplane(3),Yplane(1),Yplane(2),Yplane(3),10,'green');
quiver3(Oplane(1),Oplane(2),Oplane(3),Zplane(1),Zplane(2),Zplane(3),10,'blue');

% Project point cloud on this plane
for ipoint = 1:size(pointCloud3D,2)
    pointCloud3D_proj1(:,ipoint) = pointCloud3D(:,ipoint)-dot(fobjPlane.normal,pointCloud3D(:,ipoint)-Oplane)*fobjPlane.normal';
end
figure(1);
plot3(pointCloud3D_proj1(1,:),pointCloud3D_proj1(2,:),pointCloud3D_proj1(3,:),'Marker','.','Markersize',20,'Linestyle','none','Color','green');
for ipoint = 1:size(pointCloud3D,2)
    pointCloud2D_proj1(:,ipoint) = Rplane'*pointCloud3D_proj1(:,ipoint);
end
pointCloud2D_proj1 = pointCloud2D_proj1-mean(pointCloud2D_proj1,2);
figure(2);
title('2D projection 1');
hold on; axis equal; grid on; box on;
plot3(pointCloud2D_proj1(1,:),pointCloud2D_proj1(2,:),pointCloud2D_proj1(3,:),'Marker','.','Markersize',20,'Linestyle','none','Color','green');

% Ellipse fitting with Least-Squares based on orthogonal distance
ellipseFit = fitWithLSO(pointCloud2D_proj1(1,:),pointCloud2D_proj1(2,:));
ellipse2D = getPointsOnEllipse(ellipseFit(1),ellipseFit(2),ellipseFit(3),ellipseFit(4),rad2deg(ellipseFit(5)),100)';
ellipse2D(3,:) = mean(pointCloud2D_proj1(3,:));
ellipse2D = ellipse2D-mean(pointCloud2D_proj1,2);
figure(2);
plot3(ellipse2D(1,:),ellipse2D(2,:),ellipse2D(3,:),'Color','green');
figure(1);
for ipoint = 1:size(ellipse2D,2)
    ellipse3D(:,ipoint) = Rplane*ellipse2D(:,ipoint)+Oplane;
end
plot3([ellipse3D(1,:) ellipse3D(1,1)],[ellipse3D(2,:) ellipse3D(2,1)],[ellipse3D(3,:) ellipse3D(3,1)],'Color','green');

% Define a unit vector aligned with minor axis of the ellipse
Rellipse2 = [cos(ellipseFit(5)+pi/2) -sin(ellipseFit(5)+pi/2) 0; sin(ellipseFit(5)+pi/2) cos(ellipseFit(5)+pi/2) 0; 0 0 1];
u         = Rellipse2*[1 0 0]';
minorAxis = Rplane*u;
figure(1);
quiver3(Oplane(1),Oplane(2),Oplane(3),minorAxis(1),minorAxis(2),minorAxis(3),10,'black');

% Project point cloud on a plane having minorAxis as normal
for ipoint = 1:size(pointCloud3D,2)
    pointCloud3D_proj2(:,ipoint) = pointCloud3D(:,ipoint)-dot(minorAxis,pointCloud3D(:,ipoint)-mean(pointCloud3D,2))*minorAxis;
end
figure(1);
plot3(pointCloud3D_proj2(1,:),pointCloud3D_proj2(2,:),pointCloud3D_proj2(3,:),'Marker','.','Markersize',20,'Linestyle','none','Color','blue');

% Fit a circle to the projected 2D point cloud
Rplane2                      = [fobjPlane.normal' cross(fobjPlane.normal',minorAxis) minorAxis];
for ipoint                   = 1:size(pointCloud3D,2)
    pointCloud2D_proj2(:,ipoint) = Rplane2'*pointCloud3D_proj2(:,ipoint);
end
[xc,yc,Re,a]                 = circfit(pointCloud2D_proj2(1,:),pointCloud2D_proj2(2,:));
th                           = linspace(0,2*pi,100)';
Xcircle                      = Re*cos(th)+xc;
Ycircle                      = Re*sin(th)+yc;
Zcircle                      = repmat(mean(pointCloud2D_proj2(3,:),2),[size(Xcircle,1) 1]);
figure(3);
title('2D projection 2');
hold on; axis equal; grid on; box on;
circle2D = [Xcircle Ycircle Zcircle]';
plot3([circle2D(1,:) circle2D(1,1)],[circle2D(2,:) circle2D(2,1)],[circle2D(3,:) circle2D(3,1)]);
plot3(pointCloud2D_proj2(1,:),pointCloud2D_proj2(2,:),pointCloud2D_proj2(3,:),'Marker','.','Markersize',20,'Linestyle','none','Color','blue');
figure(1);
circle3D = Rplane2*(circle2D-[mean(Xcircle,1)/2 mean(Ycircle,1) mean(Zcircle,1)]')+Oplane;
plot3([circle3D(1,:) circle3D(1,1)],[circle3D(2,:) circle3D(2,1)],[circle3D(3,:) circle3D(3,1)],'Color','blue');

% Parametres à extraire
% 1) Position
% - Orientation (R) de l'ellipse 3D par rapport au repère inertiel (3 paramètres)
% - Position (T) de l'ellipse 3D par rapport au repère inertiel (3 paramètres)
% 2) Forme
% - Longueur (L1) grand axe de l'ellipse
% - Longueur (L2) petit axe de l'ellipse
% - Rayon (r) du cercle
