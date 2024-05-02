%% Preprocess
close all;
clear;clc;

set(0,'defaultfigurecolor','w');


%% Parameters
freq = 60;      % [Hz]

IMU = [62.5,112.5,165,222]; % [mm] Location of IMUs
L1 = IMU(2);                % [mm] The length of the first segment
L2 = IMU(4)-L1;             % [mm] The length of the second segment

%% Load (time-consuming)
gt = importdata('GroundTruth.mat');

fileName = 'Ex072401.bag';
filepath=fullfile(fileName);
bag=rosbag(fileName);
t_total=bag.EndTime-bag.StartTime;
imu_msg0  = select(bag,'Time',[bag.StartTime bag.StartTime+1.0*t_total],'Topic', '/imu0/data');
imu_msg1  = select(bag,'Time',[bag.StartTime bag.StartTime+1.0*t_total],'Topic', '/imu1/data');
imu_msg2  = select(bag,'Time',[bag.StartTime bag.StartTime+1.0*t_total],'Topic', '/imu2/data');
imu_msg3  = select(bag,'Time',[bag.StartTime bag.StartTime+1.0*t_total],'Topic', '/imu3/data');

imus0 = readMessages(imu_msg0);
imus1 = readMessages(imu_msg1);
imus2 = readMessages(imu_msg2);
imus3 = readMessages(imu_msg3);

%% Attitude Data
sample_size = min([size(imus0,1),size(imus1,1),size(imus2,1),size(imus3,1)]);
time = (1:sample_size)./freq; 
imus = [imus0(1:sample_size),imus1(1:sample_size),imus2(1:sample_size),imus3(1:sample_size)];

pose = zeros(4,sample_size,4);  

initDcm = zeros(3,3,4);

for n = 1:4
    for i = 1:sample_size
        imu = imus(i,n);

        pose(1,i,n) = imu{1}.Orientation.W;   % q0
        pose(2,i,n) = imu{1}.Orientation.X;   % q1
        pose(3,i,n) = imu{1}.Orientation.Y;   % q2
        pose(4,i,n) = imu{1}.Orientation.Z;   % q3
    end
    initDcm(:,:,n) = quat2dcm(pose(:,1,n)');
end

%% Initialize attitude (Yaw)
for n = 1:4
    for i = 1:sample_size
        pose(:,i,n) = dcm2quat( initDcm(:,:,n)\quat2dcm(pose(:,i,n)') ); 
        pose(:,i,n) = [pose(1,i,n);pose(2,i,n);-pose(3,i,n);-pose(4,i,n)]; 
    end
end

for n = 1:4
    for i = 1:sample_size
        pose(:,i,n) = [pose(1,i,n);pose(2,i,n);pose(3,i,n);0]./sqrt(pose(1,i,n)^2+pose(2,i,n)^2+pose(3,i,n)^2);    
    end
end

% Calculate Pose of Segment 2
for n = 3:4
    for i = 1:sample_size
        pose(:,i,n) = dcm2quat( quat2dcm(pose(:,i,n)')/quat2dcm(pose(:,i,2)') ); 
    end
end

% smoothing quaternion
for n = 1:4
    for i = 1:sample_size
        for j = 1:4
            if (i>1)&&(pose(j,i-1,n)*pose(j,i,n)<0)&&(abs(pose(j,i,n))>1e-1) 
                pose(:,i,n) = -pose(:,i,n);
                break; 
            end
        end
    end
end

%% Plot Quaternion
figure(1), hold on

for n = 1:4
    subplot(4,1,n);hold on;
    plot(time,pose(1,:,n), '-k','lineWidth',2); % w q0
    plot(time,pose(2,:,n), '-r','lineWidth',2); % x q1
    plot(time,pose(3,:,n), '-g','lineWidth',2); % y q2
    plot(time,pose(4,:,n), '-b','lineWidth',2); % z q3

    legend(["q_0";"q_1";"q_2";"q_3"]);
end

%% Plot Euler angle
eul=zeros(3,sample_size,4);

rad2deg=180/pi;

figure(2), hold on
for n = 1:4
    for i=1:sample_size
        eul(:,i,n) = rad2deg.*quat2eul(pose(:,i,n)','XYZ'); 
    end
    subplot(4,1,n);hold on;
    plot(time,eul(1,:,n), '-r','lineWidth',2);
    plot(time,eul(2,:,n), '-g','lineWidth',2);
    plot(time,eul(3,:,n), '-b','lineWidth',2);

    legend(["X";"Y";"Z"]);
end


%% Main Process
% Segment 1
s0_1 = IMU(1)/L1;
s1_1 = IMU(2)/L1;

A = [s0_1  1/2*s0_1^2;
     s1_1  1/2*s1_1^2];

th1 = zeros(sample_size,2);
alpha1 = zeros(sample_size,2);
phi1 = zeros(sample_size,1);
for i = 1:sample_size
    w_0 = abs(pose(1,i,1));
    w_1 = abs(pose(1,i,2));
    alpha1(i,1) = 2*asin(sign(pose(3,i,1))*sqrt(1-w_0^2));
    alpha1(i,2) = 2*asin(sign(pose(3,i,2))*sqrt(1-w_1^2));
    
    temp = A\[alpha1(i,1);alpha1(i,2)];
    th1(i,1) = temp(1);  %theta 0
    th1(i,2) = temp(2);  %theta 1
 
    phi1(i,1) = (-atan(pose(2,i,1)/pose(3,i,1))+...
                 -atan(pose(2,i,2)/pose(3,i,2)))/2;

    phi1(i,1) = phi1(i,1)-9*pi/180; % imu mounting deviation
end

figure(3);hold off;

subplot(4,2,1);hold on;grid on;
plot(time,abs(pose(1,:,1)), '-r','lineWidth',2); 
plot(time,abs(pose(1,:,2)), '-b','lineWidth',2); 
legend(['w_0';'w_1']);
title('The real part.');

subplot(4,2,3);hold on;grid on;
plot(time,rad2deg.*alpha1(:,1), '-r','lineWidth',2); 
plot(time,rad2deg.*alpha1(:,2), '-b','lineWidth',2); 
legend(['alpha_0';'alpha_1']);
title('Deflection angle');

subplot(4,2,5);hold on;grid on;
plot(time,th1(:,1), '-r','lineWidth',2);    %theta 0
plot(time,th1(:,2), '-b','lineWidth',2);    %theta 1
legend(['theta_0';'theta_1']);
title('Configuration');

subplot(4,2,7);hold on;grid on;
plot(time,rad2deg.*phi1(:,1), '-b','lineWidth',2);    %theta 1
title('Direction of deflection');


% Segment 2
s0_2 = (IMU(3)-L1)/L2;
s1_2 = (IMU(4)-L1)/L2;

A = [s0_2  1/2*s0_2^2;
     s1_2  1/2*s1_2^2];


th2 = zeros(sample_size,2);
alpha2 = zeros(sample_size,2);
phi2 = zeros(sample_size,2);
for i = 1:sample_size
    w_0 = abs(pose(1,i,3));
    w_1 = abs(pose(1,i,4));
    alpha2(i,1) = 2*asin(sign(pose(3,i,3))*sqrt(1-w_0^2));
    alpha2(i,2) = 2*asin(sign(pose(3,i,3))*sqrt(1-w_1^2));
    
    temp = A\[alpha2(i,1);alpha2(i,2)];
    th2(i,1) = temp(1);  %theta 0
    th2(i,2) = temp(2);  %theta 1
    phi2(i,1) =  - atan(pose(2,i,3)/pose(3,i,3));
    
    phi2(i,1) =  phi2(i,1) - 10*pi/180; % imu mounting deviation
end

figure(3); hold off;

subplot(4,2,2);hold on;grid on;
plot(time,abs(pose(1,:,1)), '-r','lineWidth',2); 
plot(time,abs(pose(1,:,2)), '-b','lineWidth',2); 
legend(['w_3';'w_4']);
title('The real part.');

subplot(4,2,4);hold on;grid on;
plot(time,rad2deg.*alpha2(:,1), '-r','lineWidth',2); 
plot(time,rad2deg.*alpha2(:,2), '-b','lineWidth',2); 
legend(['alpha_2';'alpha_3']);
title('Deflection angle');

subplot(4,2,6);hold on;grid on;
plot(time,th2(:,1), '-r','lineWidth',2);    %theta 0
plot(time,th2(:,2), '-b','lineWidth',2);    %theta 1
legend(['theta_0';'theta_1']);
title('Configuration');

subplot(4,2,8);hold on;grid on;
plot(time,rad2deg.*phi2(:,1), '-b','lineWidth',2);  
title('Direction of deflection');


%% Display Shape 
write_video = 0;
if write_video
    writerObj = VideoWriter(['demo\',fileName,'.mp4'],'MPEG-4'); 
    writerObj.FrameRate = 15;
    open(writerObj);  
end

gt_start = 1; 
offset = 3*freq;
GT_3d = gt(1:size(gt,1),:); % 60 [Hz]

fig = figure(4);fig.WindowState = 'maximized';

for i = 2:20:sample_size  
    Arm1 = plotArm(th1(i,1),th1(i,2),0,L1);
    Arm2 = plotArm(th2(i,1),th2(i,2),0,L2); 
    
    % @ View 1
    subplot(1,2,1); hold off;
    % Segment 1
    plot3(Arm1(1,:)*cos(phi1(i)),Arm1(1,:)*sin(phi1(i)),Arm1(2,:),'r','LineWidth',1.5);
    hold on;
    % Segment 2
    R = angle2dcm(phi1(i),alpha1(i,2),-phi1(i),'ZYZ')';
    t = [Arm1(1,size(Arm1,2))*cos(phi1(i));Arm1(1,size(Arm1,2))*sin(phi1(i));Arm1(2,size(Arm1,2))];
    plot3(t(1),t(2),t(3),'xb','LineWidth',1.5);
    points = [Arm2(1,:)*cos(phi2(i));Arm2(1,:)*sin(phi2(i));Arm2(2,:)];
    points = R*points + repmat(t,1,size(points,2));
    plot3(points(1,:),points(2,:),points(3,:),'b','LineWidth',1.5);
    axis equal;grid on;
    axis([-100 100 -100 100 0 250])
    xlabel("x")
    ylabel("y")
    zlabel("z")
    
    % Ground Truth
    gt_ = reshape(GT_3d(i+offset,2:22),[3,7]);
    plot3(gt_(1,:)',gt_(2,:)',gt_(3,:)','-xk','MarkerSize',18);
    
    view(45,20)
     
    set(gca,'ZDir','reverse');
    set(gca,'Fontsize',12)

    % @ View 2
    subplot(1,2,2); hold off;
    % Segment 1
    plot3(Arm1(1,:)*cos(phi1(i)),Arm1(1,:)*sin(phi1(i)),Arm1(2,:),'r','LineWidth',1.5);
    hold on;
    % Segment 2
    R = angle2dcm(phi1(i),alpha1(i,2),-phi1(i),'ZYZ')';
    t = [Arm1(1,size(Arm1,2))*cos(phi1(i));Arm1(1,size(Arm1,2))*sin(phi1(i));Arm1(2,size(Arm1,2))];
    plot3(t(1),t(2),t(3),'xb','LineWidth',1.5);
    points = [Arm2(1,:)*cos(phi2(i));Arm2(1,:)*sin(phi2(i));Arm2(2,:)];
    points = R*points + repmat(t,1,size(points,2));
    plot3(points(1,:),points(2,:),points(3,:),'b','LineWidth',1.5);
    axis equal;grid on;
    axis([-100 100 -100 100 0 250])
    xlabel("x")
    ylabel("y")
    zlabel("z")
    
    % Ground Truth
    gt_ = reshape(GT_3d(i+offset,2:22),[3,7]);
    plot3(gt_(1,:)',gt_(2,:)',gt_(3,:)','-xk','MarkerSize',18);
    
    view(0,90)
     
    set(gca,'ZDir','reverse');
    set(gca,'Fontsize',12)
    
    sgtitle( sprintf("Time = %.2f s / %.2f s" ,time(i) ,time(size(time,2)) ) );
    
    if write_video
        frame = getframe(gcf);          
        writeVideo(writerObj,frame);
    else
        pause(0.0001);
    end
end
fprintf('Motion End\n');
if write_video
    close(writerObj);
end


