% copy column contents:
% 1,2: target x and y in wac pixels 
% 3: switch time, NaN means no switch made
% 4,5: the onset and end time of the reach      % Time
% 6,7: endpoint x and y in wac pixels
% 8,9: start position in wac pixels
% 10: target distance in mm
% 11,12: target x and y in mm             % endpoint - target = error(reach/orth)
% 13,14: endpoint x and y in mm
% 15: target size in mm
% 16: actual duration of the reach     %？ 
% 17: error size in mm               %？
% 18: switch logical array
% 19,20: start position in mm
% 21: actual reach distance
% 22: average speed            % Ave Speed
% 23: error along the reach direction (vector projection) in mm  % Reach
% Error (overshoot and undershoot) 
% 24,25: relative target position in mm
% 26: score
% 27: trial order
% 28: hit or not
% 29: error orthogonal to the reach direction (vector rejection) in mm % Orthog Error
% 30: target shrinking duration(second)
% 31: target shrinking speed(mm/s)

%%
% load(['data_onlineConf/pilot/pilot_practice_traj_S1_08-Aug-2024_calibration.mat'])
% load(['data_onlineConf/pilot/pilot_practice_traj_S1_08-Aug-2024_tform.mat'])
% load(['data_onlineConf/pilot/pilot_practice_traj_S1_08-Aug-2024_trialdata.mat'])
% load(['data_onlineConf/pilot/pilot_practice_traj_S1c_08-Aug-2024_rawtotal.mat'])
% load(['data_onlineConf/pilot/pilot_practice_traj_S1c_08-Aug-2024_traXtotal.mat'])
% load(['data_onlineConf/pilot/pilot_practice_traj_S1c_08-Aug-2024_traYtotal.mat'])
%%
index = NaN(size(data,1),1);
for i = 1:size(data,1)
    index(i) = ~isnan(sum(data(i,:)));
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
copy(:,23) = (abs(dot(copy(:,19:20),copy(:,24:25),2) ./ dot(copy(:,24:25),copy(:,24:25),2)) - 1) .*copy(:,10); %neg+positive
copy(:,26) = valid(:,11);
copy(:,28) = copy(:,26) ~= 0;
endPoints = (copy(:,6:7) - copy(:,8:9)) .* pixellength;% relative endpoint coordinate
projScale = (abs(dot(copy(:,19:20),copy(:,24:25),2) ./ dot(copy(:,24:25),copy(:,24:25),2)));
rejections = endPoints - projScale.* copy(:,24:25);
rejLength = sqrt(rejections(:,1).^2 + rejections(:,2).^2);
copy(:,29) = rejLength; % non-negative

% 创建第30列：“target shrinking duration”
copy(:,30) = NaN; % 初始化为NaN

% 根据 trial order 赋值
copy(copy(:,27)>=1 & copy(:,27)<=60, 30)   = 1.0;
copy(copy(:,27)>=61 & copy(:,27)<=120, 30) = 0.6;
copy(copy(:,27)>=121 & copy(:,27)<=180,30) = 0.8;
copy(copy(:,27)>=181 & copy(:,27)<=240,30) = 0.4;

% 在copy中新增第31列：target shrinking speed (mm/s)
copy(:,31) = copy(:,15) ./ copy(:,30);

%%

for i = 1:size(copy,1)
    if  rem(i,2)
        copy(i,29) = copy(i,29);
    else 
        copy(i,29) = -copy(i,29);
    end
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
bUB = [2,-0.1,100];
bLB = [0.1,-1.5,0];

sigmoidFit = NaN(3,size(copy,1));
for i = 1:size(copy,1)
    x = 1:endTime(i);
    y = traProjection(i,x);
    x = 1:sum(y<1.5);
    y = traProjection(i,x);
    fun = @(b) sum(((b(1)./(1+exp(b(2)*(x-b(3)))))-y).^2); %sigmoidFit's formula 
    b = bads(fun,b0,bLB,bUB); % lsqnonlin is the funtion of sigmoidFit
    sigmoidFit(:,i) = b; % save b(b,a,c,3 parameters) into a table named signoidFit
end
figure(1)
    for i = 1:10
        trial_i = randi(length(copy),1);
        x = 1:endTime(trial_i);
        y = traProjection(trial_i,x);
        b = sigmoidFit(:,trial_i);
        adjustedX = x - b(3);
        plot(adjustedX,b(1)./(1+exp(b(2)*(x-b(3)))),'-')
        hold on
        plot(adjustedX,y,'o')
        xline(0);
        yline(1);
        yline(0);
        hold off
        xlim([-45 45])
        ylim([-0.1,1.2])
        xlabel('Recentered Time (s)')
        ylabel('Standardized Trajectory (mm)')
        title(['Trial #' num2str(trial_i)])
        pause(0.2)
    end
maxSpeed = -sigmoidFit(2,:)'.*0.25.*copy(:,10).*60; % sigmoid parameter a * 0.25 = how much percent of distance/per frame, then * distance, and *60 to get 1 second 
copy(:,30) = maxSpeed;
%%
for i = 100
    x = 1:endTime(i);
    y = traProjection(i,x);
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
scatter(maxSpeed,(copy(:,23)),'filled')
hold on 
yline(0,'--','LineWidth',2)
hold off
xlim([0,2000])
ylim([-75,75])
xlabel('Max speed (mm/s)')
ylabel('Gain error (mm)')

subplot(1,2,2)
scatter(maxSpeed,copy(:,29),'filled')
hold on 
yline(0,'--','LineWidth',2)
hold off
xlim([0,2000])
ylim([-75,75])
xlabel('Max speed (mm/s)')
ylabel('Orthogonal error (mm)')

sgtitle('Endpoint error v.s. max speed')

%%
xTotal = copy(:,23);

yTotal = copy(:,29);

theta0 = [0,0,0.0001,5];
UB = [1,30,1,50];
LB = [-1,-30,0.0001,1];

f = @(theta,speed,error) -(log(1/sqrt(2*pi)) * length(speed) + sum(-log(theta(3).*speed + theta(4)) - (((error - (theta(1).*speed + theta(2))).^2)./ (2.*(theta(3).*speed+theta(4)).^2))));

xParams = NaN(3,4);
xTotalNLL = @(theta) f(theta,copy(:,30),xTotal);
xParams(1,:) = bads(xTotalNLL,theta0,LB,UB);

xTotalAIC = 2 * 4 - 2 * -xTotalNLL(xParams(1,:));
xTotalBIC = 4 * log(length(copy(:,30))) - 2 * -xTotalNLL(xParams(1,:));

theta0 = [0,0,0.0001,5];
UB = [1,30,1,50];
LB = [-1,-30,0.0001,1];

yParams = NaN(3,4);
yTotalNLL = @(theta) f(theta,copy(:,30),yTotal);
yParams(1,:) = bads(yTotalNLL,theta0,LB,UB);

yTotalAIC = 2 * 4 - 2 * -yTotalNLL(yParams(1,:));
yTotalBIC = 4 * log(length(copy(:,30))) - 2 * -yTotalNLL(yParams(1,:));



T = table([xTotalAIC;yTotalAIC],...
    [xTotalBIC;yTotalBIC],...
    'VariableNames',{'AIC','BIC'},...
    'RowName',{'X Common','Y Commom'...
    }); 
disp(T)
%% Graphs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphaRange = [0.1,0.9];
theme_color = [1 0 0];
map = theme_color .* linspace(alphaRange(1),alphaRange(2),100)' + [1 1 1] .* (1-linspace(alphaRange(1),alphaRange(2),100)');

figure(1)
alphaRange = [0.1,0.9];
speed2alpha = (copy(:,30)-min(copy(:,30))) ./ (max(copy(:,30))-min(copy(:,30))) .* (alphaRange(2)-alphaRange(1)) + alphaRange(1);
d = scatter(copy(:,29),copy(:,23),'.','SizeData',250);
d.MarkerFaceColor = 'flat';
d.CData = theme_color .* speed2alpha + [1 1 1] .* (1-speed2alpha);
ylims = [-70,70];
xticks(-60:20:60)
yticks(-60:20:60)
xlim(ylims)
ylim(ylims)

xlabel('Orthogonal Error (mm)')
ylabel('Along Error (mm)')
title('All Trials')
xline(0)
yline(0)
nw_label = 'Overshoot                                 ';
sw_label = 'Undershoot                                 ';
text(0.05, max(ylims), nw_label, 'FontWeight', 'bold', 'HorizontalAlignment','right','VerticalAlignment', 'top');
text(0.05, min(ylims), sw_label, 'FontWeight', 'bold', 'HorizontalAlignment','right','VerticalAlignment', 'bottom');
nw_label = 'left                                ';
sw_label = '                                right';
text(0.05, min(ylims)-12, nw_label, 'FontWeight', 'bold', 'HorizontalAlignment','right','VerticalAlignment', 'bottom');
text(0.05, min(ylims)-12, sw_label, 'FontWeight', 'bold', 'HorizontalAlignment','left','VerticalAlignment', 'bottom');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
indices = 1:1800;
theme_color = [0 0.4470 0.7410];
central_linfit = (xParams(1,1) .* indices) + xParams(1,2);
predict_bias = (xParams(1,1) .* copy(:,30)) + xParams(1,2);
predict_std = (xParams(1,3) * copy(:,30)) + xParams(1,4);
z_scores = (xTotal-predict_bias)./predict_std;
ylims = [-100,100];
shades_n = 7;
z_gr = linspace(0,3.5,shades_n+1);
shade_gr = linspace(0.8,0.05,shades_n);
trials_bin = NaN(1,shades_n *2);
hold on
for i = 1:shades_n
    x_UUB = central_linfit + z_gr(i+1) .* (xParams(1,3) * indices + xParams(1,4));
    x_ULB = central_linfit + z_gr(i) .* (xParams(1,3) * indices + xParams(1,4));
    x_LUB = central_linfit - z_gr(i) .* (xParams(1,3) * indices + xParams(1,4));
    x_LLB = central_linfit - z_gr(i+1) .* (xParams(1,3) * indices + xParams(1,4));
    patch([indices,fliplr(indices)],[x_UUB,fliplr(x_ULB)], theme_color, 'FaceAlpha', shade_gr(i), 'EdgeColor','none')
    patch([indices,fliplr(indices)],[x_LUB,fliplr(x_LLB)], theme_color, 'FaceAlpha', shade_gr(i), 'EdgeColor','none')
    trials_bin(shades_n-i+1) = sum(z_scores(:) > z_gr(i) & z_scores(:) < z_gr(i+1));
    trials_bin(shades_n+i) = sum(z_scores(:) < -z_gr(i) & z_scores(:) > -z_gr(i+1));
end
plot(copy(:,30),xTotal,'.b','MarkerSize',20)
plot(indices,central_linfit,'--','Color',theme_color,'LineWidth',2)
hold off

xlabel('Speed (mm/s)')
ylabel('Endpoint Error (mm)')
title('All Trials')
xlim([0,1800])
ylim(ylims)
xticks(0:500:1500)

% Add labels to the northwest and southwest corners of the plot
nw_label = 'Overshoot        ';
sw_label = 'Undershoot       ';
text(0.05, max(ylims), nw_label, 'FontWeight', 'bold', 'HorizontalAlignment','right','VerticalAlignment', 'top');
text(0.05, min(ylims), sw_label, 'FontWeight', 'bold', 'HorizontalAlignment','right','VerticalAlignment', 'bottom');

%%
figure(3)

indices = 1:1800;
theme_color = [0 0.4470 0.7410];
central_linfit = (yParams(1,1) .* indices) + yParams(1,2);
predict_bias = (yParams(1,1) .* copy(:,30)) + yParams(1,2);
predict_std = (yParams(1,3) * copy(:,30)) + yParams(1,4);
z_scores = (yTotal-predict_bias)./predict_std;
ylims = [-100,100];
shades_n = 7;
z_gr = linspace(0,3.5,shades_n+1);
shade_gr = linspace(0.8,0.05,shades_n);
trials_bin = NaN(1,shades_n *2);
hold on
for i = 1:shades_n
    y_UUB = central_linfit + z_gr(i+1) .* (yParams(1,3) * indices + yParams(1,4));
    y_ULB = central_linfit + z_gr(i) .* (yParams(1,3) * indices + yParams(1,4));
    y_LUB = central_linfit - z_gr(i) .* (yParams(1,3) * indices + yParams(1,4));
    y_LLB = central_linfit - z_gr(i+1) .* (yParams(1,3) * indices + yParams(1,4));
    patch([indices,fliplr(indices)],[y_UUB,fliplr(y_ULB)], theme_color, 'FaceAlpha', shade_gr(i), 'EdgeColor','none')
    patch([indices,fliplr(indices)],[y_LUB,fliplr(y_LLB)], theme_color, 'FaceAlpha', shade_gr(i), 'EdgeColor','none')
    trials_bin(shades_n-i+1) = sum(z_scores(:) > z_gr(i) & z_scores(:) < z_gr(i+1));
    trials_bin(shades_n+i) = sum(z_scores(:) < -z_gr(i) & z_scores(:) > -z_gr(i+1));
end
plot(copy(:,30),yTotal,'.b','MarkerSize',20)
plot(indices,central_linfit,'--','Color',theme_color,'LineWidth',2)
hold off

xlabel('Speed (mm/s)')
ylabel('Endpoint Error (mm)')
title('All Trials')
xlim([0,1800])
ylim(ylims)
xticks(0:500:1500)

% Add labels to the northwest and southwest corners of the plot
nw_label = 'Left        ';
sw_label = 'Right       ';
text(0.05, max(ylims), nw_label, 'FontWeight', 'bold', 'HorizontalAlignment','right','VerticalAlignment', 'top');
text(0.05, min(ylims), sw_label, 'FontWeight', 'bold', 'HorizontalAlignment','right','VerticalAlignment', 'bottom');

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
