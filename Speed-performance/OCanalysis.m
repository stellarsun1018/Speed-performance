% copy column contents:
% 1,2: target x and y in wac pixels 
% 3: switch time, NaN means no switch made
% 4,5: the onset and end time of the reach
% 6,7: endpoint x and y in wac pixels
% 8,9: start position in wac pixels
% 10: target distance in mm
% 11,12: target x and y in mm
% 13,14: endpoint x and y in mm
% 15: target size in mm
% 16: actual duration of the reach
% 17: error size in mm
% 18: switch logical array
% 19,20: start position in mm
% 21: actual reach distance
% 22: average speed
% 23: error along the reach direction (vector projection) in mm
% 24,25: relative target position in mm
% 26: score
% 27: trial order
% 28: hit or not
% 29: error orthogonal to the reach direction (vector rejection) in mm
%%
% participant = 1;
% if participant == 1
%     load('JM_practice_traj_S1_06-Apr-2023_tform.mat')
%     load('JM_practice_traj_S1c_06-Apr-2023_rawtotal.mat')
%     load('JM_practice_traj_S1c_06-Apr-2023_traXtotal.mat')
%     load('JM_practice_traj_S1c_06-Apr-2023_traYtotal.mat')
% elseif participant == 0
%     load('pilot_practice_traj_S1_31-Mar-2023_tform.mat')
%     load('rawComplete.mat')
%     load('xTrajComplete.mat')
%     load('yTrajComplete.mat')
% end

%%

load(['data_onlineConf/pilot/pilot_practice_traj_S1_08-Aug-2024_calibration.mat'])
load(['data_onlineConf/pilot/pilot_practice_traj_S1_08-Aug-2024_tform.mat'])
load(['data_onlineConf/pilot/pilot_practice_traj_S1_08-Aug-2024_trialdata.mat'])
load(['data_onlineConf/pilot/pilot_practice_traj_S1c_08-Aug-2024_rawtotal.mat'])
load(['data_onlineConf/pilot/pilot_practice_traj_S1c_08-Aug-2024_traXtotal.mat'])
load(['data_onlineConf/pilot/pilot_practice_traj_S1c_08-Aug-2024_traYtotal.mat'])
%%
index = NaN(size(data,1),1);
for i = 1:size(data,1)
    index(i) = ~isnan(sum(data(i,:))); % 把NaN都去除掉 
end
valid = data(index==true,:);
validTraX = traXtotal(index==true,:);
validTraY = traYtotal(index==true,:);
%%
pixellength = 0.248;
copy = valid;
copy(:,[1,2]) = transformPointsInverse(tform,copy(:,[1,2]));
copy(:,[8,9]) = transformPointsInverse(tform,copy(:,[8,9]));
copy(:,10) = sqrt(sum((copy(:,1:2) - copy(:,8:9)).^2,2)) .* pixellength;
copy(:,[11,12]) = [copy(:,1)*pixellength (1080 - copy(:,2))*pixellength];
copy(:,[13,14]) = [copy(:,6)*pixellength (1080 - copy(:,7))*pixellength]; % 1080 = tablet pixel height
copy(:,15) = valid(:,10) .* pixellength;
copy(:,16) = copy(:,5) - copy(:,4);
copy(:,17) = sqrt( (copy(:,13)-copy(:,11)).^2 + (copy(:,14)-copy(:,12)).^2 );
copy(:,27) = 1:size(copy,1);
copy(:,19:20) = (copy(:,6:7) - copy(:,8:9)) .* pixellength;
copy(:,21) = sqrt(sum((copy(:,6:7) - copy(:,8:9)).^2,2)) .* pixellength;
copy(:,22) = copy(:,21) ./ copy(:,16);
copy(:,[24,25]) = (copy(:,1:2) - copy(:,8:9)) .* pixellength;% relative target coordinate
copy(:,23) = (abs(dot(copy(:,19:20),copy(:,24:25),2) ./ dot(copy(:,24:25),copy(:,24:25),2)) - 1) .*copy(:,10);
copy(:,26) = valid(:,11);
copy(:,28) = copy(:,26) ~= 0;
endPoints = (copy(:,6:7) - copy(:,8:9)) .* pixellength;% relative endpoint coordinate
projScale = (abs(dot(copy(:,19:20),copy(:,24:25),2) ./ dot(copy(:,24:25),copy(:,24:25),2)));
rejections = endPoints - projScale.* copy(:,24:25);
rejLength = sqrt(rejections(:,1).^2 + rejections(:,2).^2);
copy(:,29) = rejLength;
for i = size(copy,1)
    copy(i,29) = copy(i,29) * (rem(1,2));
end
%
reCenteredTrajX = NaN(size(copy,1),size(validTraX,2));
reCenteredTrajY = NaN(size(copy,1),size(validTraX,2));

for i = 1:size(copy,1)
    reCenteredTrajX(i,:) = (validTraX(i,:) - copy(i,8)) .* pixellength;
    reCenteredTrajY(i,:) = (validTraY(i,:) - copy(i,9)) .* pixellength;
end

maxDuration = sum((sum(~isnan(reCenteredTrajX),1))~=0);
traProjection = NaN(size(copy,1),maxDuration);
for i = 1:maxDuration
    traProjection(:,i) = (abs(dot([reCenteredTrajX(:,i),reCenteredTrajY(:,i)],copy(:,24:25),2) ./ dot(copy(:,24:25),copy(:,24:25),2)));
end
speedProjection = traProjection(:,2:end) - traProjection(:,1:end-1);
accProjection = speedProjection(:,2:end) - speedProjection(:,1:end-1);

endTime = sum(~isnan(reCenteredTrajX),2);
b0 = [1,-1,45];
sigmoidFit = NaN(3,size(copy,1));
for i = 1:size(copy,1)
    x = 50:endTime(i);
    y = traProjection(i,x);
    x = x - 50; % 50 can be changed as parameter
    fun = @(b) (b(1)./(1+exp(b(2)*(x-b(3)))))-y; %sigmoidFit's formula 
    b = lsqnonlin(fun,b0); % lsqnonlin is the funtion of sigmoidFit  %改成bads lower bound/higher bound
    sigmoidFit(:,i) = b; % save b(b,a,c,3 parameters) into a table named signoidFit
end
maxSpeed = -sigmoidFit(2,:)'.*0.25.*copy(:,10).*60; % sigmoid parameter a * 0.25 = how much percent of distance/per frame, then * distance, and *60 to get 1 second 
%%
for i = 1:3
    x = 50:endTime(i);
    y = traProjection(i,x);
    x = x - 50;
    b = sigmoidFit(:,i);
    adjustedX = x - b(3);
    figure %open the new figure window in everytime loop
    plot(adjustedX,b(1)./(1+exp(b(2)*(x-b(3)))),'-')
    hold on
    plot(adjustedX,y,'o')
    xline(0);
    yline(1);
    yline(0);
    hold off
    xlim([-45 45])
    ylim([-0.1,1.2])
    pause(0.2)
end


%% max speed vs euclidean error
figure
subplot(1,2,1)
plot(maxSpeed,(copy(:,23)),'o')
xlim([0,2000])
xlabel('Max speed (mm/s)')
ylabel('Gain error (mm)')

subplot(1,2,2)
plot(maxSpeed,copy(:,29),'o')
xlim([0,2000])
xlabel('Max speed (mm/s)')
ylabel('Orthogonal error (mm)')

sgtitle('Endpoint error v.s. max speed')
%%
figure
plot(copy(:,11) - copy(:,13),copy(:,12) - copy(:,14),'o')
xline(0,'--')
yline(0,'--')
xlabel('Horizontal (mm)')
ylabel('Vertical (mm)')
title('Endpoint scatter relative to target center')
xlim([-100,100])
ylim([-100,100])
%%
% copy = ZL;
% tSizes = unique(copy(:,15));
% SizeDistHitRate = NaN(5,3);
% for i = 1:5
%     for j = 1:3
%         SizeDistHitRate(i,j) = mean(copy(copy(:,15)==tSizes(i) & round(copy(:,10),-2) == j*100,28));
%     end
% end
% figure
% for i = 1:3
%     plot(tSizes,SizeDistHitRate(:,i),'-o')
%     hold on
% end
% 
% plot(tSizes,0.3:0.1:0.7,'r--')
% hold off
% xlabel('Target Size (mm)')
% ylabel('Hit Rate')
% legend('~100 mm','~200 mm','~300 mm','Expectation','Location','Northwest')
% title('Hit Rates of Each Distance Range')

%% grouping by order (e.g. first 30, second 30, etc.)
% [~,ind] = sort(copy(:,22));
% SpeedSortedTrials = copy(ind,:);
% group_n = 20;
% stds = NaN(1,group_n);
% hitrate = NaN(1,group_n);
% means = NaN(1,group_n);
% speedStd = NaN(1,group_n);
% speedPoint = NaN(1,group_n);
% groupSize = length(copy)/group_n;
% for i = 1:group_n
%     stds(i) = std(SpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,29));
%     speedPoint(i) = mean(SpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,22));
%     hitrate(i) = sum(SpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,28));
%     means(i) = mean(SpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,23));
%     speedStd(i) = std(SpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,22));
% end
% figure
% errorbar(speedPoint,stds,speedStd,'horizontal','o')
% xlabel('Average Speed (mm/s)')
% ylabel('Standard Error (mm)')
% title('Average Speed vs Error')
% mld = fitlm(speedPoint,stds)
%% max speed order, fixed group size
% [sortedMaxSpeed,ind] = sort(maxSpeed);
% maxSpeedSortedTrials = copy(ind,:);
% group_n = 20;
% stds = NaN(1,group_n);
% hitrate = NaN(1,group_n);
% means = NaN(1,group_n);
% speedStd = NaN(1,group_n);
% speedPoint = NaN(1,group_n);
% groupSize = length(copy)/group_n;
% for i = 1:group_n
%     stds(i) = std(maxSpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,29));
%     speedPoint(i) = mean(sortedMaxSpeed(i*groupSize-groupSize+1:i*groupSize));
%     hitrate(i) = sum(maxSpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,28));
%     means(i) = mean(maxSpeedSortedTrials(i*groupSize-groupSize+1:i*groupSize,23));
%     speedStd(i) = std(sortedMaxSpeed(i*groupSize-groupSize+1:i*groupSize));
% end
% figure
% errorbar(speedPoint,stds,speedStd,'horizontal','o')
% xlabel('Max Speed (mm/s)')
% ylabel('Standard Error (mm)')
% title('Max Speed vs Error')
% mld = fitlm(speedPoint,stds)

%%
% compare max speed and average speed