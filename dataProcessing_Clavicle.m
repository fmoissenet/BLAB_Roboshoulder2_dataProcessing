% Clear workspace
clearvars;
% close all;
system('taskkill /F /IM EXCEL.EXE'); 
clc;

% Set folders
Folder.toolbox      = 'C:\Users\moissene\OneDrive - unige.ch\2022 - ROBOSHOULDER2\BLAB_Roboshoulder2_dataProcessing\';
Folder.data         = 'C:\Users\moissene\OneDrive - unige.ch\Article ROBOSHOULDER2\Données\';
Folder.export       = 'C:\Users\moissene\OneDrive - unige.ch\Article ROBOSHOULDER2\Données\Stats\input\';
Folder.dependencies = [Folder.toolbox,'dependencies\'];
addpath(Folder.toolbox);
addpath(genpath(Folder.dependencies));

% Load data
clavicleList = {'RS001_L','RS001_R','RS002_L','RS002_R','RS003_L','RS003_R','RS004_L','RS004_R','RS005_L','RS005_R'};
unitList     = {'deg','deg','deg','mm','mm','mm','mm','mm','mm'};
cd(Folder.data);
kparameter  = 0;
for iclavicle = 2%1:10
    % Load file
    cd(Folder.data);
    csvFile   = ['Aruco\',clavicleList{iclavicle},'_Clavicle_composite.csv']; 
    stlFile   = ['Scan\',clavicleList{iclavicle},'_Mesh_Clavicle.stl']; 
    fcsvFile1 = ['Scan\',clavicleList{iclavicle},'_Clavicle_Landmarks.fcsv']; 
    fcsvFile2 = ['Scan\',clavicleList{iclavicle},'_Clavicle_SC_Contour.fcsv']; 
    % Get 3D bony landmark location (Aruco)
    CAJa = mean(xlsread(csvFile,1,'B3:D12'),1)'*1e3; % Range defined manually
    CAAa = mean(xlsread(csvFile,1,'E3:G12'),1)'*1e3; % Range defined manually
    CAEa = mean(xlsread(csvFile,1,'H3:J12'),1)'*1e3; % Range defined manually
    CSJa = mean(xlsread(csvFile,1,'K3:M12'),1)'*1e3; % Range defined manually
    CASa = mean(xlsread(csvFile,1,'N3:P12'),1)'*1e3; % Range defined manually
    % Get 3D bony landmark location (Scan)
    temp       = readmatrix(fcsvFile1,'Filetype','text','NumHeaderLines',3,'Delimiter',',');
    pointCloud = temp(:,2:4);
    CAJs       = pointCloud(1,:)';
    CAAs       = pointCloud(2,:)';
    CAEs       = pointCloud(3,:)';
    CSJs       = pointCloud(4,:)';
    CASs       = pointCloud(5,:)';
    clear temp pointCloud;
    % Get SC contour (Scan)
    temp       = readmatrix(fcsvFile2,'Filetype','text','NumHeaderLines',3,'Delimiter',',');
    pointCloud = temp(:,2:4);
    clear temp;
    % Get bone mesh
    temp     = stlread(stlFile);
    vertices = permute(temp.Points,[2,1,3]);
    faces    = permute(temp.ConnectivityList,[2,1,3]);
    clear temp;
    % Compute rigid transformation from scan (x) to aruco (y)
    x = [CAJs'; CAAs'; CAEs'; CSJs'; CASs'];
    y = [CAJa'; CAAa'; CAEa'; CSJa'; CASa'];
    [R,d,rsm] = soder(x,y);
    clear x y;
    % Apply rigid transformation to scan data
    CAJsa = R*CAJs+d;
    CAAsa = R*CAAs+d;
    CAEsa = R*CAEs+d;
    CSJsa = R*CSJs+d;
    CASsa = R*CASs+d;
    for ipointCloud = 1:size(pointCloud,1)
        pointClouda(:,ipointCloud) = R*pointCloud(ipointCloud,:)'+d;
    end
    for ivertices = 1:size(vertices,2)
        verticesa(:,ivertices) = R*vertices(:,ivertices)+d;
    end
    % Set scapula coordinate system (ISB recommendations)
    if contains(clavicleList{iclavicle},'_R') % Right side convention
        Oc = CSJa;
        Zc = (CAJa-CSJa)./norm(CAJa-CSJa);
        Xc = (CASa-mean(pointClouda,2))./norm(CASa-mean(pointClouda,2));
        Yc = cross(Zc,Xc);
        Rc = [Xc Yc Zc];
    elseif contains(clavicleList{iclavicle},'_L') % Left side convention
        Oc = CSJa;
        Zc = -(CAJa-CSJa)./norm(CAJa-CSJa);
        Xc = (CASa-mean(pointClouda,2))./norm(CASa-mean(pointClouda,2));
        Yc = cross(Zc,Xc);
        Rc = [Xc Yc Zc];
    end
    % Compute cartilage contour parameters
    for irater = 1:3 
        for itrial = 1:3
            cd(Folder.data);
            disp(['Clavicle ',num2str(iclavicle),' - Rater ',num2str(irater),' - Trial ',num2str(itrial)]);
            pointCloud3Dfull = nan(300,3*3*3);
            pointCloud3Dfull(1:length(xlsread(csvFile,1,'Z3:AH300')*1e3),1:9)   = xlsread(csvFile,1,'Z3:AH300')*1e3; % mm % Other columns are related to AC joint
            pointCloud3Dfull(1:length(xlsread(csvFile,1,'AI3:AQ300')*1e3),10:18) = xlsread(csvFile,1,'AI3:AQ300')*1e3; % mm % Other columns are related to AC joint
            pointCloud3Dfull(1:length(xlsread(csvFile,1,'AR3:AZ300')*1e3),19:27) = xlsread(csvFile,1,'AR3:AZ300')*1e3; % mm % Other columns are related to AC joint
            temp             = pointCloud3Dfull(:,(irater-1)*9+itrial*3-2:(irater-1)*9+itrial*3); % X, Y, Z coordinates
            pointCloud3Di    = temp(~isnan(temp(:,1)),:)';

            % Express the 3D point cloud in the clavicle coordinate system
            pointCloud3Di = inv(Rc)*(pointCloud3Di-Oc);
            CAJai         = inv(Rc)*(CAJa-Oc);
            CSJai         = inv(Rc)*(CSJa-Oc);
            CASai         = inv(Rc)*(CASa-Oc);
            Oci           = inv(Rc)*(Oc-Oc);
            Xci           = inv(Rc)*Xc;
            Yci           = inv(Rc)*Yc;
            Zci           = inv(Rc)*Zc;
            for ivertices = 1:size(verticesa,2)
                verticesai(:,ivertices) = inv(Rc)*(verticesa(:,ivertices)-Oc);
            end

            % Transform left clavicles as right clavicles
            if contains(clavicleList{iclavicle},'_L')
                pointCloud3Di(3,:) = -pointCloud3Di(3,:);
                CAJai(3,:)         = -CAJai(3,:);
                CSJai(3,:)         = -CSJai(3,:);
                CASai(3,:)         = -CASai(3,:);
                Oci(3,:)           = -Oci(3,:);
                Xci(3,:)           = -Xci(3,:);
                Yci(3,:)           = -Yci(3,:);
                Zci(3,:)           = Zci(3,:);
                verticesai(3,:)    = -verticesai(3,:);
            end

            % Fit a plane to the point cloud at the least square sense
            fobjPlane = planarFit(pointCloud3Di);
            figure(1);
            [hL,hD]   = plot(fobjPlane);
            hold on; axis equal; grid on; box on;
            title('3D view');
            plot3(pointCloud3Di(1,:),pointCloud3Di(2,:),pointCloud3Di(3,:),'Marker','.','Markersize',20,'Linestyle','none','Color','red');
            plot3(CAJai(1,:),CAJai(2,:),CAJai(3,:),'Marker','o','Markersize',10,'Linestyle','none','Color','black');
            plot3(CSJai(1,:),CSJai(2,:),CSJai(3,:),'Marker','o','Markersize',10,'Linestyle','none','Color','black');
            plot3(CASai(1,:),CASai(2,:),CASai(3,:),'Marker','o','Markersize',10,'Linestyle','none','Color','black');
            quiver3(Oci(1),Oci(2),Oci(3),Xci(1),Xci(2),Xci(3),10,'red');
            quiver3(Oci(1),Oci(2),Oci(3),Yci(1),Yci(2),Yci(3),10,'green');
            quiver3(Oci(1),Oci(2),Oci(3),Zci(1),Zci(2),Zci(3),10,'blue');
            patch_array3(faces,verticesai,[0.8 0.8 0.8],'none','gouraud',1);
            
            % Define a plane coordinate system
            Oplane = mean(pointCloud3Di,2);
            Zplane = fobjPlane.normal';
            if Zplane(1) < 0 % Avoid different plane normal directions
                Zplane(1) = -Zplane(1);
                Zplane(2) = -Zplane(2);
            end
            Yplane = -cross([1 0 0]',Zplane)./norm(cross([1 0 0]',Zplane));
            Xplane = cross(Yplane,Zplane)./norm(cross(Yplane,Zplane));
            Rplane = [Xplane Yplane Zplane];
%             figure(1);
%             quiver3(Oplane(1),Oplane(2),Oplane(3),Xplane(1),Xplane(2),Xplane(3),10,'red');
%             quiver3(Oplane(1),Oplane(2),Oplane(3),Yplane(1),Yplane(2),Yplane(3),10,'green');
%             quiver3(Oplane(1),Oplane(2),Oplane(3),Zplane(1),Zplane(2),Zplane(3),10,'blue');
            
            % Project point cloud on this plane
            for ipoint = 1:size(pointCloud3Di,2)
                pointCloud3Di_proj1(:,ipoint) = pointCloud3Di(:,ipoint)-dot(fobjPlane.normal,pointCloud3Di(:,ipoint)-Oplane)*fobjPlane.normal';
            end
%             figure(1);
%             plot3(pointCloud3Di_proj1(1,:),pointCloud3Di_proj1(2,:),pointCloud3Di_proj1(3,:),'Marker','.','Markersize',20,'Linestyle','none','Color','green');
            for ipoint = 1:size(pointCloud3Di,2)
                pointCloud2D_proj1(:,ipoint) = Rplane'*pointCloud3Di_proj1(:,ipoint);
            end
            pointCloud2D_proj1 = pointCloud2D_proj1-mean(pointCloud2D_proj1,2);
%             figure(2);
%             title('2D projection 1');
%             hold on; axis equal; grid on; box on;
%             plot3(pointCloud2D_proj1(1,:),pointCloud2D_proj1(2,:),pointCloud2D_proj1(3,:),'Marker','.','Markersize',20,'Linestyle','none','Color','green');
            
            % Ellipse fitting with Least-Squares based on orthogonal distance
            ellipseFit     = fitWithLSO(pointCloud2D_proj1(1,:),pointCloud2D_proj1(2,:));
            ellipse2D      = getPointsOnEllipse(ellipseFit(1),ellipseFit(2),ellipseFit(3),ellipseFit(4),rad2deg(ellipseFit(5)),100)';
            ellipse2D(3,:) = mean(pointCloud2D_proj1(3,:));
            ellipse2D      = ellipse2D-mean(pointCloud2D_proj1,2);
%             figure(2);
%             plot3(ellipse2D(1,:),ellipse2D(2,:),ellipse2D(3,:),'Color','green');
            for ipoint = 1:size(ellipse2D,2)
                ellipse3D(:,ipoint) = (Rplane*ellipse2D(:,ipoint)+Oplane);
            end
%             figure(1);
%             plot3([ellipse3D(1,:) ellipse3D(1,1)],[ellipse3D(2,:) ellipse3D(2,1)],[ellipse3D(3,:) ellipse3D(3,1)],'Color','green');
            
            % Define a unit vector aligned with major and minor axes of the ellipse
            Rellipse2 = [cos(ellipseFit(5)) -sin(ellipseFit(5)) 0; sin(ellipseFit(5)) cos(ellipseFit(5)) 0; 0 0 1];
            u         = Rellipse2*[1 0 0]';
            v         = Rellipse2*[0 1 0]';
            majorAxis = Rplane*u;
            minorAxis = Rplane*v;
%             figure(1);
%             quiver3(Oplane(1),Oplane(2),Oplane(3),majorAxis(1),majorAxis(2),majorAxis(3),10,'red');
%             quiver3(Oplane(1),Oplane(2),Oplane(3),minorAxis(1),minorAxis(2),minorAxis(3),10,'green');
            
            % Project point cloud on a plane having minorAxis as normal
            for ipoint = 1:size(pointCloud3Di,2)
                pointCloud3Di_proj2(:,ipoint) = pointCloud3Di(:,ipoint)-dot(minorAxis,pointCloud3Di(:,ipoint)-mean(pointCloud3Di,2))*minorAxis;
            end
%             figure(1);
%             plot3(pointCloud3Di_proj2(1,:),pointCloud3Di_proj2(2,:),pointCloud3Di_proj2(3,:),'Marker','.','Markersize',20,'Linestyle','none','Color','blue');
            
            % Fit a circle to the projected 2D point cloud
            Rplane2                      = [fobjPlane.normal' cross(fobjPlane.normal',minorAxis) minorAxis];
            for ipoint                   = 1:size(pointCloud3Di,2)
                pointCloud2D_proj2(:,ipoint) = Rplane2'*pointCloud3Di_proj2(:,ipoint);
            end
            [xc,yc,Re,a]                 = circfit(pointCloud2D_proj2(1,:),pointCloud2D_proj2(2,:));
            th                           = linspace(0,2*pi,100)';
            Xcircle                      = Re*cos(th)+xc;
            Ycircle                      = Re*sin(th)+yc;
            Zcircle                      = repmat(mean(pointCloud2D_proj2(3,:),2),[size(Xcircle,1) 1]);
            circle2D                     = [Xcircle Ycircle Zcircle]';
%             figure(3);
%             title('2D projection 2');
%             hold on; axis equal; grid on; box on;
%             plot3([circle2D(1,:) circle2D(1,1)],[circle2D(2,:) circle2D(2,1)],[circle2D(3,:) circle2D(3,1)]);
%             plot3(pointCloud2D_proj2(1,:),pointCloud2D_proj2(2,:),pointCloud2D_proj2(3,:),'Marker','.','Markersize',20,'Linestyle','none','Color','blue');
            circle3D = Rplane2*(circle2D-[mean(Xcircle,1) mean(Ycircle,1) mean(Zcircle,1)]'+[-Re 0 0]')+Oplane;
%             figure(1);
%             plot3([circle3D(1,:) circle3D(1,1)],[circle3D(2,:) circle3D(2,1)],[circle3D(3,:) circle3D(3,1)],'Color','blue');

            % Store parameters
            Rtemp                   = Rplane*Rellipse2;
            temp                    = R2mobileXYZ_array3(Rtemp);
            kparameter              = kparameter+1;
            for ipoint = 1:size(pointCloud3Di,2)
                temp2(ipoint) = abs(Re-sqrt((pointCloud2D_proj2(1,ipoint)-xc)^2+(pointCloud2D_proj2(2,ipoint)-yc)^2));
            end
            circleRMS(kparameter) = mean(temp2)*100/Re;
            Parameter(1,kparameter) = -rad2deg(temp(1)); % Orientation of the 3D ellipse in SCS (/Xscapula) (deg)
            Parameter(2,kparameter) = -rad2deg(temp(2)); % Orientation of the 3D ellipse in SCS (/Yscapula) (deg)
            Parameter(3,kparameter) = -rad2deg(temp(3)); % Orientation of the 3D ellipse in SCS (/Zscapula) (deg)
            Parameter(4,kparameter) = Oplane(1); % Centre of the 3D ellipse in SCS (/Xscapula) (mm)
            Parameter(5,kparameter) = Oplane(2); % Centre of the 3D ellipse in SCS (/Yscapula) (mm)
            Parameter(6,kparameter) = Oplane(3); % Centre of the 3D ellipse in SCS (/Zscapula) (mm)
            Parameter(7,kparameter) = ellipseFit(3); % Ellipse major axis length (mm)
            Parameter(8,kparameter) = ellipseFit(4); % Ellipse major axis length (mm)
            Parameter(9,kparameter) = Re; % Circle radius (mm)
            clear Rtemp temp temp2;
            
            % Export parameters
%             cd(Folder.export);
% %             close all;
%             for iparameter = 1:9
%                 if kparameter == 0
%                     xlswrite(['Cartilage_clavicle_acromioclavicular_parameter_',num2str(iparameter),'.csv'],{'CARTILAGE'},1,'A1');
%                     xlswrite(['Cartilage_clavicle_acromioclavicular_parameter_',num2str(iparameter),'.csv'],{'OPERATOR'},1,'B1');
%                     xlswrite(['Cartilage_clavicle_acromioclavicular_parameter_',num2str(iparameter),'.csv'],{'DIGITALISATION'},1,'C1');
%                     xlswrite(['Cartilage_clavicle_acromioclavicular_parameter_',num2str(iparameter),'.csv'],{'VALUE'},1,'D1');
%                 end
%                 xlswrite(['Cartilage_clavicle_acromioclavicular_parameter_',num2str(iparameter),'.csv'],[iclavicle irater itrial Parameter(iparameter,kparameter)],1,['A',num2str(kparameter+1)]);
%             end

            % Clear workspace
            clearvars -except iclavicle irater itrial kparameter Folder pointCloudLabel clavicleList csvFile Oc Xc Yc Zc Rc CAJa CSJa CASa verticesa vertices faces Parameter circleRMS;
            
        end
    end

    % Clear workspace
    clearvars -except Folder kparameter clavicleList Parameter circleRMS;
end
hold on;
plot(sort(circleRMS),'red')