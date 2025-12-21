% copy column contents:
% 1,2: target x and y in wac pixels 
% 3: life span by block
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
% 30: maxSpeed
% 31: target shrinking duration(second)
% 32: target shrinking speed(mm/s)

%%
lifespan = [0.6,0.6*3^(0.25),0.6*3^(0.5)]; %[1.0,0.6,0.8,0.4]; %[1.1,0.9,1.0,0.8,0.6,0.7] 设定各blocks中target的不同时长  %lifespan控制了受试者实际可用的、逐渐减少的目标"可见e时间窗，这一时间越短，任务难度越高（因为受试者必须更快速地完成任务以取得更高分数）。
for i = 1:3
    copy((1+(i-1)*240):(i*240),3) = lifespan(i);
end
%%
clear all
participant = 'LC';
session = 3;
fname_preamble = sprintf('data_onlineConf/%s/%s_sptatialTemporalCostFunc_S%d*.mat',participant,participant,session);
files = dir(fname_preamble);
for k = 1:numel(files)
    f = fullfile(files(k).folder, files(k).name);
    load(f);
end

%%
index = NaN(size(data,1),1);
for i = 1:size(data,1)
    index(i) = ~isnan(sum(data(i,:)));
end
valid = data(index==true,:);
validTraX = traXtotal(index==true,:);
validTraY = traYtotal(index==true,:);
%%
% addpath('/Users/stellar/Desktop/Landy lab/Speed-Accuracy-2D Shrinking Circle/bads-master'); % 更改为你解压的路径
% savepath; % 保存路径，下次启动 MATLAB 自动加载
which bads

%%
pixellength = 0.248;
Affine2d =tform.T(1:2,1:2);
[~,s,~] = svd(Affine2d);
proj2tablet = 1./mean([s(1,1),s(2,2)]);  
mmPerProjPx = proj2tablet .* pixellength; % mmPerProjPx是从projector-实际tablet上的mm 之倍数  % 每次calibration能得到新的该值，每次每个人得到的target size都不一样

copy = valid;
copy(:,[1,2]) = transformPointsInverse(tform,copy(:,[1,2]));
copy(:,[8,9]) = transformPointsInverse(tform,copy(:,[8,9]));
copy(:,10) = sqrt(sum((copy(:,1:2) - copy(:,8:9)).^2,2)) .* pixellength;
copy(:,[11,12]) = [copy(:,1)*pixellength (1080 - copy(:,2))*pixellength];
copy(:,[13,14]) = [copy(:,6)*pixellength (1080 - copy(:,7))*pixellength]; % 1080 = tablet pixel height
copy(:,15) = valid(:,10) .* mmPerProjPx;
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
% Calculate the projection scale (no abs here to preserve sign)
projScale = dot(copy(:,19:20), copy(:,24:25), 2) ./ dot(copy(:,24:25), copy(:,24:25), 2);

% Calculate the projection vector
projections = projScale .* copy(:,24:25);

% Compute the rejections (subtract projections from original vectors)
rejections = copy(:,19:20) - projections;

% Calculate the signed rejection length using determinant for directionality
rejLength = sign(copy(:,19).*copy(:,25) - copy(:,20).*copy(:,24)) .* sqrt(rejections(:,1).^2 + rejections(:,2).^2);

% Store the result in copy(:,29)
copy(:,29) = rejLength;



% 创建第30列：“target shrinking duration”
copy(:,31) = NaN; % 初始化为NaN

% 根据 trial order 赋值
copy(copy(:,27)>=1 & copy(:,27)<=60, 31)   = 1.0;
copy(copy(:,27)>=61 & copy(:,27)<=120, 31) = 0.6;
copy(copy(:,27)>=121 & copy(:,27)<=180,31) = 0.8;
copy(copy(:,27)>=181 & copy(:,27)<=240,31) = 0.4;

% 在copy中新增第32列：target shrinking speed (mm/s)
copy(:,32) = copy(:,15) ./ copy(:,31);

%%  polar coordination 

angle_error_in_radian = asin(copy(:,29) ./ copy(:,21));

angle_error_in_degree = rad2deg(angle_error_in_radian);

copy(:,33) = angle_error_in_radian;
copy(:,34) = angle_error_in_degree;

%% 3 blocks - conditions graphs

lim_scale = 1.2;
figure()
for i = 1:3
    block_ind = zeros(size(copy,1),1);
    block_ind((1+(i-1)*240):(i*240)) = 1;
    subplot(1,3,i)

    distances = copy(block_ind==1,10);
    avg_speed = copy(block_ind==1,22);
    x = linspace(unique(copy(block_ind==1,15)),lim_scale*max(copy(:,10)),2);
    y = x ./ unique(copy(block_ind==1,3));
    
    plot(distances,avg_speed,'o');
    hold on

% 原始 condition line: speed = distance / lifespan
    plot(x,y,'--','Color',[0 0.4 0.8],'LineWidth',1.5); % 蓝色

% 新增：半大小的 condition line（红色）
    y_half = x ./ (0.5 * unique(copy(block_ind==1,3)));
    plot(x, y_half, '--r','LineWidth',1.5); % 红色虚线

    hold off
    xlabel("Target Distance (mm)");
    ylabel("Average Speed (mm/s)");
    legend('Trial data','Min speed (full target)', 'Min speed (half-size target)','Location','northwest')
    xlim([0,lim_scale * max(copy(:,10))]);
    ylim([0,lim_scale * max(copy(:,22))]);

    % hold off
    % xlabel("Target Distance (mm)");
    % ylabel("Average Speed (mm/s)");
    % legend('Trial data','Minimum speed')
    % xlim([0,lim_scale * max(copy(:,10))]);
    % ylim([0,lim_scale * max(copy(:,22))]);
end


%% Overlay all 3 blocks with regression lines through origin
figure;
colors = lines(3);  % 三种颜色

% 初始化 legend 对象数组
h_scatter = gobjects(3,1);
h_line = gobjects(3,1);

for i = 1:3
    block_inds = (1+(i-1)*240):(i*240);
    distances = copy(block_inds,10);
    avg_speed = copy(block_inds,22);

    % 画散点图，并保存 scatter 句柄
    h_scatter(i) = scatter(distances, avg_speed, 20, colors(i,:), 'filled'); hold on;

    % 原点回归：y = kx
    k = distances \ avg_speed;
    x_fit = linspace(min(distances), max(distances), 100);
    y_fit = k * x_fit;

    % 虚线拟合线，并保存句柄
    h_line(i) = plot(x_fit, y_fit, '--', 'Color', colors(i,:), 'LineWidth', 2);

    % 显示 slope
    fprintf('Block %d: slope = %.2f, implied lifespan = %.2f sec\n', i, k, 1/k);
end

xlabel('Target Distance (mm)');
ylabel('Average Speed (mm/s)');
title('Distance vs Speed (3 Blocks with Origin-Fixed Regression)');

% 用回归线句柄生成 legend（不是散点）
legend(h_line, ...
    {'Block 1 (slope)', 'Block 2 (slope)', 'Block 3 (slope)'}, ...
    'Location', 'northwest');

grid on;

%%
% Add red dashed lines for half-size threshold
lifespans = [0.6, 0.6*3^(0.25), 0.6*3^(0.5)];
for i = 1:3
    half_time = lifespans(i) / 2;  % time to shrink to half
    x_half = linspace(0, max(copy(:,10)), 100);
    y_half = x_half ./ half_time;
    plot(x_half, y_half, 'r--', 'LineWidth', 1.5);
end

%%  Duration histograms by block (overlapped) + disappearance time lines
% copy(:,16) = duration(s); copy(:,3) = lifespan(s)（每个block恒定）

assert(size(copy,1) >= 3*240, '预期每个block约240 trial，请确认数据尺寸。');

trialsPerBlock = 240;
blocks = 1:3;

% 颜色（与要求一致：红/蓝/黄）
blockColors = [1 0 0; 0 0.45 0.95; 1 0.9 0]; % red, blue, yellow

% 收集各 block 的 duration
Dur = cell(1,3);
Tdisappear = nan(1,3);
for i = blocks
    r = (1+(i-1)*trialsPerBlock) : (i*trialsPerBlock);
    r = r(r <= size(copy,1));      % 防越界
    Dur{i} = copy(r,16);           % duration (s)
    ls = unique(copy(r,3));        % lifespan(s) —— 目标消失时间
    ls = ls(~isnan(ls));
    if isempty(ls)
        Tdisappear(i) = NaN;
    else
        % 若存在极少数异常，可用中位数稳健估计
        Tdisappear(i) = median(ls);
    end
end

% 统一 bin（避免直方图形状受bin不同影响）
allDur = copy(:,16);
allDur = allDur(isfinite(allDur));
edges = linspace(prctile(allDur,1), prctile(allDur,99), 30);

figure('Name','Duration Histograms by Block');
hold on;
hH = gobjects(1,3);
hV = gobjects(1,3);
for i = blocks
    hH(i) = histogram(Dur{i}, edges, ...
        'Normalization','pdf', ...
        'FaceColor', blockColors(i,:), ...
        'EdgeColor','none', ...
        'FaceAlpha', 0.35);
    if isfinite(Tdisappear(i))
        hV(i) = xline(Tdisappear(i), '--', ...
            'Color', blockColors(i,:), 'LineWidth', 2);
    end
end
xlabel('Reach Duration T (s)');
ylabel('Density');
title('Overlapped Duration Histograms with Disappearance Times');
grid on;

% 图例：区分直方图与其对应的消失时间
lgd = legend( ...
    [hH(1) hH(2) hH(3) hV(1) hV(2) hV(3)], ...
    {'Block 1 durations','Block 2 durations','Block 3 durations', ...
     sprintf('Block 1 disappear @ %.2fs', Tdisappear(1)), ...
     sprintf('Block 2 disappear @ %.2fs', Tdisappear(2)), ...
     sprintf('Block 3 disappear @ %.2fs', Tdisappear(3))}, ...
    'Location','northwest');
lgd.Box = 'off';
hold off;

%% Compare mean durations across blocks (bar + errorbars) + print stats

meansT = nan(1,3);
stdT   = nan(1,3);
nT     = nan(1,3);
semT   = nan(1,3);

for i = blocks
    di = Dur{i};
    di = di(isfinite(di));
    meansT(i) = mean(di);
    stdT(i)   = std(di);
    nT(i)     = numel(di);
    semT(i)   = stdT(i)/sqrt(max(nT(i),1));
end

% 柱状图 （+ SEM 误差条）
figure('Name','Mean Duration by Block');
bh = bar(1:3, meansT, 'FaceColor','flat'); 
for i = blocks
    bh.CData(i,:) = blockColors(i,:);
end
% hold on;
% errorbar(1:3, meansT, semT, 'k', 'LineStyle','none', 'LineWidth',1.5, 'CapSize',8);

% 在柱子上方写均值
for i = 1:3
    text(i, meansT(i) + 0.005, sprintf('%.3f', meansT(i)), ...
        'HorizontalAlignment','center', 'FontSize', 12, 'FontWeight','bold');
end

xticks(1:3); xticklabels({'Block 1','Block 2','Block 3'});
ylabel('Mean Duration (s)');
title('Average Reach Duration per Block');
grid on; box off;

% 控制台打印统计
fprintf('\n=== Duration Summary by Block ===\n');
for i = blocks
    fprintf('Block %d: N=%d, Mean=%.4f s, SD=%.4f s, SEM=%.4f s, Disappear=%.4f s\n', ...
        i, nT(i), meansT(i), stdT(i), semT(i), Tdisappear(i));
end


%%

%% 3d dot plot + linear reg
distances = copy(:,10);
avg_speed = copy(:,22);
errors = copy(:,17);
plot3(distances,avg_speed,errors,'o');
xlabel("Target Distance (mm)");
ylabel("Average Speed (mm/s)");
zlabel("Endpoint Error (mm)");

%% 3D reg plane - Euclidean Error 
% Extract variables
reach_distances = copy(:,21);
avg_speed = copy(:,22);
errors = copy(:,17);

% Prepare the design matrix (adding a column of ones for intercept)
X = [ones(size(reach_distances)), reach_distances, avg_speed];

% Perform multivariate linear regression
coeffs = regress(errors, X);

% Display regression coefficients
fprintf('Intercept: %.4f\n', coeffs(1));
fprintf('Distance coefficient: %.4f\n', coeffs(2));
fprintf('Average speed coefficient: %.4f\n', coeffs(3));

% 3D scatter plot of the original data
figure;
plot3(reach_distances, avg_speed, errors, 'o');
hold on;

xlim([0 max(reach_distances)]);
ylim([0 max(avg_speed)]);

% Generate grid for regression plane
[dist_grid, speed_grid] = meshgrid(linspace(min(reach_distances), max(reach_distances), 20), ...
                                   linspace(min(avg_speed), max(avg_speed), 20));

% Predict errors using regression model
error_fit = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*speed_grid;

% Plot regression plane
mesh(dist_grid, speed_grid, error_fit);
xlabel('Reach Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Endpoint Euc Error (mm)');
title('Multivariate Linear Regression: Error ~ Reach Distance + Speed');
grid on;
hold off;

%%
% Euclidean Error 
% 3D scatter plot of the original data

figure;
plot3(reach_distances, avg_speed, errors, 'o');
hold on;

xlim([0 max(reach_distances)]);
ylim([0 max(avg_speed)]);

% Generate grid for regression plane
[dist_grid, speed_grid] = meshgrid(linspace(min(reach_distances), max(reach_distances), 50), ...
                                   linspace(min(avg_speed), max(avg_speed), 50));

% Predict errors using regression model
error_fit = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*speed_grid;

% Plot regression plane
mesh(dist_grid, speed_grid, error_fit);

% Add contour projection to bottom (Z = min)
contour3(dist_grid, speed_grid, error_fit, 5, 'k', 'LineWidth', 1);

% Axis labels and formatting
xlabel('Reach Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Reach Error (mm)');
title('Multivariate Linear Regression: Reach Error ~ Reach Distance + Speed');
grid on;
hold off;
%%
% %% 3D reg plane - Endpoint(x) distribution of diff Reaching distance & ave speed
% % Extract variables
% reach_distances = copy(:,21);
% avg_speed = copy(:,22);
% endpoints_x = copy(:,23);
% 
% % Prepare the design matrix (adding a column of ones for intercept)
% X = [ones(size(reach_distances)), reach_distances, avg_speed];
% 
% % Perform multivariate linear regression
% coeffs = regress(endpoints_x, X);
% 
% % Display regression coefficients
% fprintf('Intercept: %.4f\n', coeffs(1));
% fprintf('Distance coefficient: %.4f\n', coeffs(2));
% fprintf('Average speed coefficient: %.4f\n', coeffs(3));
% 
% % 3D scatter plot of the original data
% figure;
% plot3(reach_distances, avg_speed, endpoints_x, 'o');
% hold on;
% 
% % Generate grid for regression plane
% [dist_grid, speed_grid] = meshgrid(linspace(min(reach_distances), max(reach_distances), 20), ...
%                                    linspace(min(avg_speed), max(avg_speed), 20));
% 
% % Predict endpoints_x using regression model
% error_fit = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*speed_grid;
% 
% % Plot regression plane
% mesh(dist_grid, speed_grid, error_fit);
% xlabel('Reach Distance (mm)');
% ylabel('Average Speed (mm/s)');
% zlabel('endpoints x (mm)');
% title('Multivariate Linear Regression: endpoints x ~ Reach Distance + Speed');
% grid on;
% hold off;
% 
% %%
% %% 3D reg plane - Endpoint(y) distribution of diff Reaching distance & ave speed
% % Extract variables
% reach_distances = copy(:,21);
% avg_speed = copy(:,22);
% endpoints_y = copy(:,14);
% 
% % Prepare the design matrix (adding a column of ones for intercept)
% X = [ones(size(reach_distances)), reach_distances, avg_speed];
% 
% % Perform multivariate linear regression
% coeffs = regress(endpoints_y, X);
% 
% % Display regression coefficients
% fprintf('Intercept: %.4f\n', coeffs(1));
% fprintf('Distance coefficient: %.4f\n', coeffs(2));
% fprintf('Average speed coefficient: %.4f\n', coeffs(3));
% 
% % 3D scatter plot of the original data
% figure;
% plot3(reach_distances, avg_speed, endpoints_y, 'o');
% hold on;
% 
% % Generate grid for regression plane
% [dist_grid, speed_grid] = meshgrid(linspace(min(reach_distances), max(reach_distances), 20), ...
%                                    linspace(min(avg_speed), max(avg_speed), 20));
% 
% % Predict endpoints_y using regression model
% error_fit = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*speed_grid;
% 
% % Plot regression plane
% mesh(dist_grid, speed_grid, error_fit);
% xlabel('Reach Distance (mm)');
% ylabel('Average Speed (mm/s)');
% zlabel('endpoints y (mm)');
% title('Multivariate Linear Regression: endpoints y ~ Reach Distance + Speed');
% grid on;
% hold off;

%%  Z axis: Gain error ~ reach distance + speed
% Extract variables
reach_distances = copy(:,21);
avg_speed = copy(:,22);
Reach_errors = copy(:,23);

% Prepare the design matrix (adding a column of ones for intercept)
X = [ones(size(reach_distances)), reach_distances, avg_speed];

% Perform multivariate linear regression
coeffs = regress(Reach_errors, X);

% Display regression coefficients
fprintf('Intercept: %.4f\n', coeffs(1));
fprintf('Distance coefficient: %.4f\n', coeffs(2));
fprintf('Average speed coefficient: %.4f\n', coeffs(3));

% 3D scatter plot of the original data
figure;
plot3(reach_distances, avg_speed, Reach_errors, 'o');
hold on;

xlim([0 max(reach_distances)]);
ylim([0 max(avg_speed)]);

% Generate grid for regression plane
[dist_grid, speed_grid] = meshgrid(linspace(min(reach_distances), max(reach_distances), 20), ...
                                   linspace(min(avg_speed), max(avg_speed), 20));

% Predict errors using regression model
error_fit = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*speed_grid;

% Plot regression plane
mesh(dist_grid, speed_grid, error_fit);
xlabel('Reach Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Gain Error (mm)');
title('Multivariate Linear Regression: Gain Error ~ Reach Distance + Speed');
grid on;
hold off;

%% add Contours 加入等高线
% 3D scatter plot of the original data
figure;
plot3(reach_distances, avg_speed, Reach_errors, 'o');
hold on;

xlim([0 max(reach_distances)]);
ylim([0 max(avg_speed)]);

% Generate grid for regression plane
[dist_grid, speed_grid] = meshgrid(linspace(min(reach_distances), max(reach_distances), 50), ...
                                   linspace(min(avg_speed), max(avg_speed), 50));

% Predict errors using regression model
error_fit = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*speed_grid;

% Plot regression plane
h_mesh = mesh(dist_grid, speed_grid, error_fit);
alpha(h_mesh, 0.7); % 让mesh半透明以便观察等高线

% Add contour projection on bottom (z = min(error_fit))
z_level = min(error_fit(:)); % 投影平面位置
contour3(dist_grid, speed_grid, error_fit, 5, 'k', 'LineWidth', 1); % 黑色等高线

% Axis labels and formatting
xlabel('Reach Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Gain Error (mm)');
title('Multivariate Linear Regression: Gain Error ~ Reach Distance + Speed');
grid on;
view(3); % 确保是3D视角
hold off;

%% step 2: residuals 
errors_predicted = X * coeffs;

% Calculate residuals
residuals1 = Reach_errors - errors_predicted;

% 3D scatter plot of residuals
figure;
plot3(reach_distances, avg_speed, residuals1, 'ro');
xlim([0 max(reach_distances)]);
ylim([0 max(avg_speed)]);

xlabel('Reach Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Residuals (mm)');
title('Residuals of Multivariate Linear Regression: Reach Error ~ Distance + Speed');
grid on;

%% reach error qq plot
figure;
qqplot(residuals1)

%% Z axis: angle_error_in_degree (orthognal) ~ reach distance + speed

% Extract variables
reach_distances = copy(:,21);
avg_speed = copy(:,22);
angle_error_in_degree = copy(:,34);

% Prepare the design matrix (adding a column of ones for intercept)
X = [ones(size(reach_distances)), reach_distances, avg_speed];

% Perform multivariate linear regression
coeffs = regress(angle_error_in_degree, X);

% Display regression coefficients
fprintf('Intercept: %.4f\n', coeffs(1));
fprintf('Distance coefficient: %.4f\n', coeffs(2));
fprintf('Average speed coefficient: %.4f\n', coeffs(3));

% 3D scatter plot of the original data
figure;
plot3(reach_distances, avg_speed, angle_error_in_degree, 'o');
hold on;

xlim([0 max(reach_distances)]);
ylim([0 max(avg_speed)]);

% Generate grid for regression plane
[dist_grid, speed_grid] = meshgrid(linspace(min(reach_distances), max(reach_distances), 20), ...
                                   linspace(min(avg_speed), max(avg_speed), 20));

% Predict errors using regression model
error_fit = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*speed_grid;

% Plot regression plane
mesh(dist_grid, speed_grid, error_fit);
xlabel('Reach Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Angle error (degree)');
title('Multivariate Linear Regression: Angle error ~ Reach Distance + Speed');
grid on;
hold off;

%%
% 3D scatter plot of the original data
figure;
plot3(reach_distances, avg_speed, angle_error_in_degree, 'o');
hold on;

xlim([0 max(reach_distances)]);
ylim([0 max(avg_speed)]);

% Generate grid for regression plane
[dist_grid, speed_grid] = meshgrid(linspace(min(reach_distances), max(reach_distances), 50), ...
                                   linspace(min(avg_speed), max(avg_speed), 50));

% Predict errors using regression model
error_fit = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*speed_grid;

% Plot regression plane
h = mesh(dist_grid, speed_grid, error_fit);
alpha(h, 0.7);  % 半透明 mesh，方便观察 contour

% Add contour lines projected onto the bottom
contour3(dist_grid, speed_grid, error_fit, 5, 'k', 'LineWidth', 1);

% Labels and aesthetics
xlabel('Reach Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Angle Error (degree)');
title('Multivariate Linear Regression: Angle Error ~ Distance + Speed');
grid on;
view(3);
hold off;

%% step 2: orthognal error residuals 
errors_predicted = X * coeffs;

% Calculate residuals
residuals2 = angle_error_in_degree - errors_predicted;

% 3D scatter plot of residuals
figure;
plot3(reach_distances, avg_speed, residuals2, 'ro');
xlim([0 max(reach_distances)]);
ylim([0 max(avg_speed)]);

xlabel('Reach Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Residuals (mm)');
title('Residuals of Multivariate Linear Regression: Angle Error ~ Distance + Speed');
grid on;

%% qq plot
figure;
qqplot(residuals2)

%%



%% Step 1: reaching direction 平行误差的线性拟合 + bias plane
%拟合平行误差 (copy(:,23)) 作为距离和平均速度的函数
% 变量提取
distances = copy(:,10);      % 第10列：目标距离 (mm)
avg_speed = copy(:,22);      % 第22列：平均速度 (mm/s)
parallel_error = copy(:,23); % 第23列：平行方向误差 (reach direction error, mm)

% 构建回归设计矩阵 (加上截距项)
X = [ones(size(distances)), distances, avg_speed];

% 使用多元线性回归拟合平行误差
coeffs = regress(parallel_error, X);  % 拟合系数 [intercept, beta1, beta2]

% 显示拟合系数（可选）
fprintf('Bias model (Parallel error):\n');
fprintf('  Intercept: %.4f\n', coeffs(1));
fprintf('  Distance coefficient: %.4f\n', coeffs(2));
fprintf('  Speed coefficient: %.4f\n', coeffs(3));

% 用拟合模型预测误差值（系统偏差）
predicted_bias = X * coeffs;

% 保存系统偏差面到 copy 的新列中（第33列）
copy(:,33) = predicted_bias;

% 计算残差（真实误差减去bias）
residual_bias = parallel_error - predicted_bias;

% 保存残差（清除bias后的纯粹误差）到 copy 的第34列
copy(:,34) = residual_bias;

%% 可视化：3D scatter + bias plane（拟合面）
figure;
scatter3(distances, avg_speed, parallel_error, 30, 'b', 'filled');
hold on;

% 创建网格用于绘制拟合面
[dist_grid, speed_grid] = meshgrid(linspace(min(distances), max(distances), 30), ...
                                   linspace(min(avg_speed), max(avg_speed), 30));
bias_plane = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*speed_grid;

% 绘制 bias plane
mesh(dist_grid, speed_grid, bias_plane);
alpha(0.5); % 设置半透明
colormap(parula);

xlabel('Target Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Reach gain error (mm)');
title('Linear Fit of Parallel Error (Bias Plane)');
legend('Observed Data','Fitted Bias Plane');
grid on;
view(-45, 30);
hold off;

%%
% 构建回归设计矩阵 (加上截距项)
X = [ones(size(distances)), distances, avg_speed];

% 使用多元线性回归拟合平行误差
coeffs = regress(abs(residual_bias), X);  % 拟合系数 [intercept, beta1, beta2]

figure;
scatter3(distances, avg_speed, abs(residual_bias), 30, 'b', 'filled');
hold on;

% 创建网格用于绘制拟合面
[dist_grid, speed_grid] = meshgrid(linspace(min(distances), max(distances), 30), ...
                                   linspace(min(avg_speed), max(avg_speed), 30));
bias_plane = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*speed_grid;

% 绘制 bias plane
mesh(dist_grid, speed_grid, bias_plane);
alpha(0.5); % 设置半透明
colormap(parula);

xlabel('Target Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Residual (mm)');
title('Linear Fit of Parallel Error Residual');
legend('Observed Data','Fitted Bias Plane');
grid on;
view(-45, 30);
hold off;

%%
% 线性高斯模型：使用 BADS 同时拟合 均值 与 方差（Reach 方向）

%%
distances = copy(:,10);      % 第10列：目标距离 (mm)
avg_speed = copy(:,22);      % 第22列：平均速度 (mm/s)
parallel_error = copy(:,23); % absolute unbiased gain error

b0 = [0,0,0,0.01,0.01,10];
bUB = [0.2,0.2,30,0.5,0.5,30]; % [b(1)上下值，slope, maxspeed x轴值]  % b(1) - 决定 sigmoid 曲线的最大值
bLB = [-0.2,-0.2,-30,0,0,1];
fun = @(b) -sum(log(normpdf(parallel_error,b(1)*distances + b(2)*avg_speed + b(3),b(4)*distances + b(5)*avg_speed + b(6))));
b = bads(fun,b0,bLB,bUB);

%% 画图：Linear Gaussian Model fitted Bias

figure(9)
% 创建网格用于绘制拟合面
[dist_grid, speed_grid] = meshgrid(linspace(min(distances), max(distances), 30), ...
                                   linspace(min(avg_speed), max(avg_speed), 30));
bias_plane = b(3) + b(1)*dist_grid + b(2)*speed_grid;

% 绘制 bias plane
mesh(dist_grid, speed_grid, bias_plane);
alpha(0.5); % 设置半透明
colormap(parula);
xlabel('Target Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Linear Gaussian Model fitted Bias (mm)');
title("Gain error")

figure(10)
[dist_grid, speed_grid] = meshgrid(linspace(min(distances), max(distances), 30), ...
                                   linspace(min(avg_speed), max(avg_speed), 30));
std_plane = b(6) + b(4)*dist_grid + b(5)*speed_grid;

% 绘制 bias plane
mesh(dist_grid, speed_grid, std_plane);
alpha(0.5); % 设置半透明
colormap(parula);
xlabel('Target Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Linear Gaussian Model fitted Standard Deviation (mm)');
title("Gain error")

%%
distances = copy(:,10);      % 第10列：目标距离 (mm)
avg_speed = copy(:,22);      % 第22列：平均速度 (mm/s)
parallel_error = copy(:,29); % absolute unbiased gain error

b0 = [0,0,0,0.01,0.01,10];
bUB = [0.2,0.2,30,0.5,0.5,30]; % [b(1)上下值，slope, maxspeed x轴值]  % b(1) - 决定 sigmoid 曲线的最大值
bLB = [-0.2,-0.2,-30,0,0,1];
fun = @(b) -sum(log(normpdf(parallel_error,b(1)*distances + b(2)*avg_speed + b(3),b(4)*distances + b(5)*avg_speed + b(6))));
b = bads(fun,b0,bLB,bUB);

%% 画图：Linear Gaussian Model fitted Bias

figure(11)
% 创建网格用于绘制拟合面
[dist_grid, speed_grid] = meshgrid(linspace(min(distances), max(distances), 30), ...
                                   linspace(min(avg_speed), max(avg_speed), 30));
bias_plane = b(3) + b(1)*dist_grid + b(2)*speed_grid;

% 绘制 bias plane
mesh(dist_grid, speed_grid, bias_plane);
alpha(0.5); % 设置半透明
colormap(parula);
xlabel('Target Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Linear Gaussian Model fitted Bias (mm)');
title("Orthogonal")

figure(12)
[dist_grid, speed_grid] = meshgrid(linspace(min(distances), max(distances), 30), ...
                                   linspace(min(avg_speed), max(avg_speed), 30));
std_plane = b(6) + b(4)*dist_grid + b(5)*speed_grid;

% 绘制 bias plane
mesh(dist_grid, speed_grid, std_plane);
alpha(0.5); % 设置半透明
colormap(parula);
xlabel('Target Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Linear Gaussian Model fitted Standard Deviation (mm)');
title("Orthogonal")

%%


%% Step 2: reaching direction 平行误差的线性拟合 + bias plane
%拟合平行误差 (copy(:,23)) 作为距离和平均速度的函数
% 变量提取
distances = copy(:,10);      % 第10列：目标距离 (mm)
avg_speed = copy(:,22);      % 第22列：平均速度 (mm/s)
parallel_error = abs(copy(:,23) - mean(copy(:,23))); % absolute unbiased gain error

% 构建回归设计矩阵 (加上截距项)
X = [ones(size(distances)), distances, avg_speed];

% 使用多元线性回归拟合平行误差
coeffs = regress(parallel_error, X);  % 拟合系数 [intercept, beta1, beta2]

% 显示拟合系数（可选）
fprintf('Bias model (Parallel error):\n');
fprintf('  Intercept: %.4f\n', coeffs(1));
fprintf('  Distance coefficient: %.4f\n', coeffs(2));
fprintf('  Speed coefficient: %.4f\n', coeffs(3));

% 用拟合模型预测误差值（系统偏差）
predicted_bias = X * coeffs;

% 保存系统偏差面到 copy 的新列中（第33列）
copy(:,33) = predicted_bias;

% 计算残差（真实误差减去bias）
residual_bias = parallel_error - predicted_bias;

% 保存残差（清除bias后的纯粹误差）到 copy 的第34列
copy(:,34) = residual_bias;

%% 可视化：3D scatter + bias plane（拟合面）
figure;
scatter3(distances, avg_speed, parallel_error, 30, 'b', 'filled');
hold on;

% 创建网格用于绘制拟合面
[dist_grid, speed_grid] = meshgrid(linspace(min(distances), max(distances), 30), ...
                                   linspace(min(avg_speed), max(avg_speed), 30));
bias_plane = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*speed_grid;

% 绘制 bias plane
mesh(dist_grid, speed_grid, bias_plane);
alpha(0.5); % 设置半透明
colormap(parula);

xlabel('Target Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Absolute Unbiased Gain Error (mm)');
title('Linear Fit of Parallel Error (Bias Plane)');
legend('Observed Data','Fitted Bias Plane');
grid on;
view(-45, 30);
hold off;

%% 可视化：残差 QQ 图
figure;
qqplot(residual_bias);
title('QQ Plot of Residuals (Parallel Error - Bias)');


%% Step 2：正交误差的线性拟合 + 残差提取
% 提取正交误差
orth_error = copy(:,29);

% 回归分析
coeffs_orth = regress(orth_error, X);

% 预测值 + 残差
predicted_orth = X * coeffs_orth;
residual_orth = orth_error - predicted_orth;

% 存储残差
copy(:,34) = residual_orth; % 第34列用于存 residual of orthogonal error

% qqplot 检查残差是否近似正态
figure; qqplot(residual_orth);
title('QQ Plot of Residual (Orthogonal Error)');

%%
% bias面
error_fit_parallel = coeffs_parallel(1) + coeffs_parallel(2)*dist_grid + coeffs_parallel(3)*speed_grid;
mesh(dist_grid, speed_grid, error_fit_parallel);


%% change y axis to duration copy(:,16) 
% Extract variables
distances = copy(:,10);
reach_duration = copy(:,16);
errors = copy(:,17);

% Prepare the design matrix (adding a column of ones for intercept)
X = [ones(size(distances)), distances, reach_duration];

% Perform multivariate linear regression
coeffs = regress(errors, X);

% Display regression coefficients
fprintf('Intercept: %.4f\n', coeffs(1));
fprintf('Distance coefficient: %.4f\n', coeffs(2));
fprintf('Reach Duration coefficient: %.4f\n', coeffs(3));

% 3D scatter plot of the original data
figure;
plot3(distances, reach_duration, errors, 'o');
hold on;

% Generate grid for regression plane
[dist_grid, duration_grid] = meshgrid(linspace(min(distances), max(distances), 20), ...
                                   linspace(min(reach_duration), max(reach_duration), 20));

% Predict errors using regression model
error_fit = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*duration_grid;

% Plot regression plane
mesh(dist_grid, duration_grid, error_fit);
xlabel('Target Distance (mm)');
ylabel('Duration (s)');
zlabel('Endpoint Error (mm)');
title('Multivariate Linear Regression: Error ~ Distance + Duration');
grid on;
hold off;

%% step 2: residuals 
errors_predicted = X * coeffs;

% Calculate residuals
residuals = errors - errors_predicted;

% 3D scatter plot of residuals
figure;
plot3(distances, avg_speed, residuals, 'ro');
xlabel('Target Distance (mm)');
ylabel('Duration(s)');
zlabel('Residuals (mm)');
title('Residuals of Multivariate Linear Regression: Error ~ Distance + Duration');
grid on;


%% 3D plot of Distance vs. Shrinking Speed vs. Avg Error

% 明确定义目标距离区间的范围（修改后的区间）
distance_bins = [50,150; 150,250; 250,350]; % 三个距离区间 (mm)
speed_values = [11.6560, 14.5700, 19.4267, 29.1400]; % 四个速度水平 (mm/s)

% 将每个距离区间的中心点用于绘图
distance_labels = {'100mm','200mm','300mm'};
distance_centers = [100, 200, 300]; % 中心点分别对应区间中心点

% 创建绘图的网格坐标点
[X,Y] = meshgrid(speed_values, distance_centers);

% 初始化Z矩阵以存放平均误差
Z = NaN(size(X));

% 计算每个grid内的平均误差值
for i = 1:size(distance_bins,1)
    for j = 1:length(speed_values)
        
        % 提取当前bin的数据索引
        idx = copy(:,10) >= distance_bins(i,1) & copy(:,10) < distance_bins(i,2) & ...
              abs(copy(:,32) - speed_values(j)) < 0.01;
        
        % 计算平均误差大小
        Z(i,j) = mean(copy(idx,17),'omitnan');
    end
end

% 绘制3D曲面图 (平均误差大小)
figure;
surf(X,Y,Z);
% 设置颜色映射和颜色条
colormap('jet');
colorbar;

% 设置轴标签和图标题
xlabel('Target Shrinking Speed (mm/s)','FontSize',12);
ylabel('Target Distance Bin','FontSize',12);
zlabel('Average Error Size (mm)','FontSize',12);
title('3D Surface: Distance Bin vs. Shrinking Speed vs. Avg Error','FontSize',14);

% 调整Y轴刻度标签，便于理解距离区间
yticks(distance_centers);
yticklabels(distance_labels);

% 设置网格线和平滑显示
grid on;
shading interp; 
view(-45,30);

% 在每个数据点上标记具体的误差数值
hold on;
for i = 1:numel(X)
    if ~isnan(Z(i))
        text(X(i), Y(i), Z(i)+0.2, sprintf('%.2f',Z(i)),...
            'FontSize',10,'FontWeight','bold',...
            'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');
    end
end
hold off;


%% 3D plot of Distance vs. Shrinking Speed vs. Sd of Error (reach direction)


figure;

% 明确定义目标距离区间的范围
distance_bins = [50,150; 150,250; 250,300]; % 三个距离区间 (mm)
speed_values = [11.6560, 14.5700, 19.4267, 29.1400]; % 四个速度水平 (mm/s)

% 定义每个距离区间的中心点和标签用于绘图
distance_centers = [100, 200, 275]; % 中心点 (最后一个区间250-300取275)
distance_labels = {'50-150mm','150-250mm','250-300mm'};

% 创建绘图网格坐标
[X,Y] = meshgrid(speed_values, distance_centers);

% 初始化Z矩阵用于存放标准差
Z_sd_reach_dir = NaN(size(X));

% 计算每个grid内reach方向误差(copy(:,23))的标准差
for i = 1:size(distance_bins,1)
    for j = 1:length(speed_values)
        
        % 提取当前bin的数据索引
        idx = copy(:,10) >= distance_bins(i,1) & copy(:,10) < distance_bins(i,2) & ...
              abs(copy(:,32) - speed_values(j)) < 0.01;  % 注意：这里是copy(:,32)
        
        % 计算标准差(Standard Deviation)
        Z_sd_reach_dir(i,j) = std(copy(idx,23),'omitnan');
    end
end

% 绘制3D曲面图 (Standard Deviation of Reach Direction Error)
surf(X,Y,Z_sd_reach_dir);

% 设置颜色映射和颜色条
colormap('jet');
colorbar;

% 设置轴标签和图标题
xlabel('Target Shrinking Speed (mm/s)','FontSize',12);
ylabel('Target Distance Bin','FontSize',12);
zlabel('STD of Error Along Reach Direction (mm)','FontSize',12);
title('3D Surface: Distance Bin vs. Shrinking Speed vs. STD of Reach Direction Error','FontSize',14);

% 调整Y轴刻度标签
yticks(distance_centers);
yticklabels(distance_labels);

% 设置网格线和平滑显示
grid on;
shading interp; 
view(-45,30);

% 在每个数据点上标记具体的标准差数值
hold on;
for i = 1:numel(X)
    if ~isnan(Z_sd_reach_dir(i))
        text(X(i), Y(i), Z_sd_reach_dir(i)+0.1, sprintf('%.2f',Z_sd_reach_dir(i)),...
            'FontSize',10,'FontWeight','bold',...
            'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');
    end
end
hold off;

%% 3D plot of Distance vs. Shrinking Speed vs. Sd of Error (orthognal)

% 创建新的图窗口（避免覆盖之前的图）
figure;

% 定义目标距离区间范围
distance_bins = [50,150; 150,250; 250,300]; % 三个距离区间 (mm)
speed_values = [11.6560, 14.5700, 19.4267, 29.1400]; % 四个速度水平 (mm/s)

% 定义每个区间的中心点和标签（用于绘图）
distance_centers = [100, 200, 275]; % 最后一个区间中心取275
distance_labels = {'50-150mm','150-250mm','250-300mm'};

% 创建绘图网格
[X,Y] = meshgrid(speed_values, distance_centers);

% 初始化Z矩阵用于存放正交误差标准差
Z_sd_orthogonal = NaN(size(X));

% 计算每个grid内正交误差(copy(:,29))的标准差
for i = 1:size(distance_bins,1)
    for j = 1:length(speed_values)
        
        % 提取当前bin的数据索引
        idx = copy(:,10) >= distance_bins(i,1) & copy(:,10) < distance_bins(i,2) & ...
              abs(copy(:,32) - speed_values(j)) < 0.01; % 注意使用的是copy(:,32)

        % 计算标准差(Standard Deviation)
        Z_sd_orthogonal(i,j) = std(copy(idx,29),'omitnan');
    end
end

% 绘制3D曲面图 (正交误差标准差)
surf(X,Y,Z_sd_orthogonal);

% 设置颜色映射和颜色条
colormap('jet');
colorbar;

% 设置轴标签和图标题
xlabel('Target Shrinking Speed (mm/s)','FontSize',12);
ylabel('Target Distance Bin','FontSize',12);
zlabel('STD of Orthogonal Error (mm)','FontSize',12);
title('3D Surface: Distance Bin vs. Shrinking Speed vs. STD of Orthogonal Error','FontSize',14);

% 调整Y轴刻度标签
yticks(distance_centers);
yticklabels(distance_labels);

% 设置网格线和平滑显示
grid on;
shading interp; 
view(-45,30);

% 在每个数据点上标记具体的标准差数值
hold on;
for i = 1:numel(X)
    if ~isnan(Z_sd_orthogonal(i))
        text(X(i), Y(i), Z_sd_orthogonal(i)+0.1, sprintf('%.2f',Z_sd_orthogonal(i)),...
            'FontSize',10,'FontWeight','bold',...
            'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');
    end
end
hold off;


%%
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
    traProjection(:,i) = (abs(dot([reCenteredTrajX(:,i),reCenteredTrajY(:,i)],copy(:,24:25),2) ./ dot(copy(:,24:25),copy(:,24:25),2))); %逐时间点投影计算
end
speedProjection = traProjection(:,2:end) - traProjection(:,1:end-1);
accProjection = speedProjection(:,2:end) - speedProjection(:,1:end-1); %思考：这个有用吗？

endTime = sum(~isnan(reCenteredTrajX),2); 
b0 = [1,-1,45];
bUB = [2,-0.1,100]; % [b(1)上下值，slope, maxspeed x轴值]  % b(1) - 决定 sigmoid 曲线的最大值
bLB = [0.1,-1.5,0];

sigmoidFit = NaN(3,size(copy,1));
for i = 1:size(copy,1)
    x = 1:endTime(i); % x是时间
    y = traProjection(i,x);  % y是走过的比例 
    x = 1:sum(y<1.5); % x给了一个限制（限制y要小于1.5倍的距离
    y = traProjection(i,x); % 在该限制下的y取值
    fun = @(b) sum(((b(1)./(1+exp(b(2)*(x-b(3)))))-y).^2); %sigmoidFit's formula 
    % b = bads(fun,b0,bLB,bUB); % lsqnonlin is the funtion of sigmoidFit
    b = bads(fun,b0,bLB,bUB); % lsqnonlin is the funtion of sigmoidFit
    sigmoidFit(:,i) = b; % save b(b,a,c,3 parameters) into a table named signoidFit
end
%%
 figure(1) %画sigmoid的图出来
    for i = 1:5
        figure(i)
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
% % 图 1: maxSpeed vs copy(:,23)
% figure;
% scatter(maxSpeed, copy(:,23), 'filled'); % 绘制散点图
% xlabel('Max Speed (mm/s)'); % 横轴标签
% ylabel('Error Along the Reach Direction (mm)'); % 纵轴标签
% title('Max Speed vs Reach Direction Error'); % 图标题
% grid on; % 添加网格线

% 图 1: maxSpeed vs copy(:,23)
figure;
scatter(maxSpeed, copy(:,23), 'filled'); % 绘制散点图
xlabel('Max Speed (mm/s)'); % 横轴标签
ylabel('Error Along the Reach Direction (mm)'); % 纵轴标签
title('Max Speed vs Reach Direction Error'); % 图标题
grid on; % 添加网格线
hold on; % 保持当前图形

% 计算线性回归
coeffs_rd = polyfit(maxSpeed, copy(:,23), 1); % 一阶多项式拟合 (线性拟合)
x_fit_rd = linspace(min(maxSpeed), max(maxSpeed), 100); % 生成拟合线的横坐标范围
y_fit_rd = polyval(coeffs_rd, x_fit_rd); % 计算回归线的 y 值

% 绘制线性回归直线
plot(x_fit_rd, y_fit_rd, 'r-', 'LineWidth', 2); % 红色回归线

% 计算相关系数
corr_rd = corrcoef(maxSpeed, copy(:,23)); % 计算相关矩阵
corr_value_rd = corr_rd(1,2); % 取出相关系数

% 在图上显示相关系数
text(mean(maxSpeed), max(copy(:,23)), sprintf('r = %.2f', corr_value_rd), ...
    'FontSize', 12, 'Color', 'black', 'FontWeight', 'bold');

% 添加图例
legend('Data Points', 'Linear Fit', 'Location', 'Best');
hold off;

%%

% 图: Max Speed vs Orthogonal Error
figure;
scatter(maxSpeed, copy(:,29), 'filled'); % 绘制散点图
xlabel('Max Speed (mm/s)'); % 横轴标签
ylabel('Orthogonal Error (mm)'); % 纵轴标签
title('Max Speed vs Orthogonal Error'); % 图标题
grid on; % 添加网格线
hold on; % 保持当前图形

% 计算线性回归
coeffs_oe = polyfit(maxSpeed, copy(:,29), 1); % 一阶多项式拟合 (线性拟合)
x_fit_oe = linspace(min(maxSpeed), max(maxSpeed), 100); % 生成拟合线的横坐标范围
y_fit_oe = polyval(coeffs_oe, x_fit_oe); % 计算回归线的 y 值

% 绘制线性回归直线
plot(x_fit_oe, y_fit_oe, 'r-', 'LineWidth', 2); % 红色回归线

% 计算相关系数
corr_oe = corrcoef(maxSpeed, copy(:,29)); % 计算相关矩阵
corr_value_oe = corr_oe(1,2); % 取出相关系数

% 在图上显示相关系数
text(mean(maxSpeed), max(copy(:,29)), sprintf('r = %.2f', corr_value_oe), ...
    'FontSize', 12, 'Color', 'black', 'FontWeight', 'bold');

% 添加图例
legend('Data Points', 'Linear Fit', 'Location', 'Best');
hold off;

%%

% 图 4: copy(:,23) vs copy(:,29) (颜色由 maxSpeed 决定)
figure; hold on;

% 归一化 maxSpeed，用于颜色映射
normalizedSpeed = (maxSpeed - min(maxSpeed)) / (max(maxSpeed) - min(maxSpeed));

% 创建自定义红色渐变色映射：从粉红到深红
nColors = 100; % 颜色渐变的步数
cmap = [linspace(1, 0.5, nColors)', linspace(0, 0, nColors)', linspace(0.5, 0, nColors)']; 
% R从1到0.5，G为0，B从0.5到0

% 将 normalizedSpeed 映射到自定义颜色表
colorIndex = round(normalizedSpeed * (nColors - 1)) + 1; % 颜色索引
colorIndex(colorIndex > nColors) = nColors; % 确保索引不越界

% 绘制散点图
for i = 1:length(copy(:,23))
    scatter(copy(i,23), copy(i,29), 50, cmap(colorIndex(i), :), 'filled');
end

% 添加图例和标签
xlabel('Error Along the Reach Direction (mm)'); % 横轴标签
ylabel('Orthogonal Error (mm)'); % 纵轴标签
title('Reach Direction Error vs Orthogonal Error (Colored by Max Speed)'); % 图标题
grid on; % 添加网格线

% 添加 color bar，标示速度对应的颜色梯度
colormap(cmap); % 应用自定义的颜色映射
cb = colorbar;
cb.Label.String = 'Max Speed (mm/s)';
cb.Label.FontSize = 12; % 设置字体大小
cb.Ticks = linspace(0, 1, 5); % 颜色条刻度
cb.TickLabels = round(linspace(min(maxSpeed), max(maxSpeed), 5), 0); % 速度范围
hold off;

%% Plot Distance vs. Error

plot(copy(:,10),'o')

plot(copy(:,23),copy(:,10),'o')

plot(copy(:,29),copy(:,10),'o')

%% 图5
% 目标：绘制 copy(:,23) vs copy(:,29)，颜色由 copy(:,10) target distance 决定
% 颜色分类：
% - target distance 50-150 mm：浅蓝色
% - target distance 150-250 mm：湖蓝色
% - target distance 250-350 mm：黑色

figure; hold on;
grid on;

% 定义颜色
color1 = [0, 0, 0]; % 黑
color2 = [0, 0.5, 1];   % 湖蓝色
color3 = [0.6, 0.8, 1];     % 浅蓝

% 获取不同类别的索引
idx1 = copy(:,10) >= 50 & copy(:,10) < 150;
idx2 = copy(:,10) >= 150 & copy(:,10) < 250;
idx3 = copy(:,10) >= 250 & copy(:,10) <= 350;

% 绘制不同类别的散点
scatter(copy(idx1,23), copy(idx1,29), 50, color1, 'filled');
scatter(copy(idx2,23), copy(idx2,29), 50, color2, 'filled');
scatter(copy(idx3,23), copy(idx3,29), 50, color3, 'filled');

% 添加图例和标签
xlabel('Error Along the Reach Direction (mm)'); % 横轴标签
ylabel('Orthogonal Error (mm)'); % 纵轴标签
title('Reach Direction Error vs Orthogonal Error (Colored by Target Distance)'); % 图标题
legend('50-150 mm (Black)', '150-250 mm (Medium Blue)', '250-350 mm (Light Blue)', 'Location', 'Best');

hold off;

%% 图6

% 目标：绘制 copy(:,23) vs copy(:,30)，颜色由 copy(:,10) target distance 决定
% 颜色分类：


figure; hold on;
grid on;

% 定义颜色
color1 = [0, 0, 0];      % 黑色 (50-150 mm)
color2 = [0, 0.5, 0];    % 深绿
color3 = [0.6, 0.98, 0.7];  % 浅绿

% 获取不同类别的索引
idx1 = copy(:,10) >= 50 & copy(:,10) < 150;
idx2 = copy(:,10) >= 150 & copy(:,10) < 250;
idx3 = copy(:,10) >= 250 & copy(:,10) <= 350;

% 绘制不同类别的散点图
scatter(copy(idx1,23), copy(idx1,30), 50, color1, 'filled');
scatter(copy(idx2,23), copy(idx2,30), 50, color2, 'filled');
scatter(copy(idx3,23), copy(idx3,30), 50, color3, 'filled');

% 添加图例和标签
xlabel('Error Along the Reach Direction (mm)'); % X 轴标签
ylabel('Max Speed (mm/s)'); % Y 轴标签 (已修改)
title('Reach Direction Error vs Max Speed (Colored by Target Distance)'); % 图标题
legend('50-150 mm (Black)', '150-250 mm (Medium Blue)', '250-350 mm (Light Blue)', 'Location', 'Best');

hold off;



%%
% for i = 1:size(data,1)
%     x = 50:endTime(i);
%     y = traProjection(i,x);
%     x = x - 50;
%     b = sigmoidFit(:,i);
%     adjustedX = x - b(3);
%     figure %open the new figure window in everytime loop
%     plot(adjustedX,b(1)./(1+exp(b(2)*(x-b(3)))),'-')
%     hold on
%     plot(adjustedX,y,'o')
%     xline(0);
%     yline(1);
%     yline(0);
%     hold off
%     xlim([-45 45])
%     ylim([-0.1,1.2])
%     pause(0.2)
% end
% 

%% grouping by range (e.g. 50-100, 100-150, etc)
% copy = ZL;
% [sortspd,ind] = sort(copy(:,22));
% group_n = 15;
% stds = NaN(1,group_n);
% hitrate = NaN(1,group_n);
% means = NaN(1,group_n);
% speedStd = NaN(1,group_n);
% speedPoint = NaN(1,group_n);
% for i = 1:group_n
%     stds(i) = std(copy((i*50 < copy(:,22) & i*50+50 > copy(:,22)),23));
%     hitrate(i) = sum(copy((i*50 < copy(:,22) & i*50+50 > copy(:,22)),28));
%     means(i) = mean(copy((i*50 < copy(:,22) & i*50+50 > copy(:,22)),23));
%     speedStd(i) = std(copy(i*50 < copy(:,22) & i*50+50 > copy(:,22),22));
%     speedPoint(i) = mean(copy(i*50 < copy(:,22) & i*50+50 > copy(:,22),22));
%     sum(i*50 < copy(:,22) & i*50+50 > copy(:,22))
% end
% 
% mld = fitlm(speedPoint,stds)
% 
% errorstdX = std(copy(:,23));
% bootMax = 1000;
% numTestTrials = 50;
% 
% numConds = 13;
% numCondd = 14;
% bootXs = NaN(numTestTrials,numConds);
% bootXd = NaN(numTestTrials,numCondd);
% stdBoots = NaN(bootMax,numConds);
% meanBoots = NaN(bootMax,numConds);
% 
% for ii = 1:bootMax
%     
%     for jj = 1:numConds
%         i = jj+2;
%         selectionLogic = i*50 < copy(:,22) & i*50+50 > copy(:,22);
%         numTrials = sum(selectionLogic);
%         bootInd = randi(numTrials,numTestTrials,1);
%         subgroup = copy(selectionLogic,23);
%         bootXs(:,jj) = subgroup(bootInd); 
% %         bootY(:,jj) = epMatY(bootInd(:,jj),jj);
%     end
% 
%     bootstdXs = std(bootXs);
%     bootmeanXs = mean(bootXs);
% %     bootstdY = std(bootY);
% %     bootstdXY = sqrt(bootstdX.^2 + bootstdY.^2);
%     stdBoots(ii,:) = bootstdXs;
%     meanBoots(ii,:) = bootmeanXs;
% end
% 
% errorsErrors = NaN(1,numConds);
% for i = 1:numConds
%     errorsErrors(i) = std(stdBoots(:,i));
% end
% 
% figure
% errorbar(speedPoint(3:15),stds(3:15),errorsErrors,'vertical','-bo')
% hold on 
% errorbar(speedPoint(3:15),stds(3:15),speedStd(3:15),'horizontal','-bo')
% hold off
% ylim([0,25])
% ylabel('STD along the reach (mm)')
% xlabel('Average speed of the reach (mm/s)')
% title('Errorbar data from 1000 bootstrap samples of 50 trials each')
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

% plot(sort(maxSpeed),sort(copy(:,22)),'o')
% 
% xlabel('Max Speed mm/s')
% ylabel('Average Speed mm/s')
% title('Max v.s. Average Speed, each point = one trial')
% 
% 
% 
% plot(sort(maxSpeed),sort(copy(:,23)),'o')
% 
% xlabel('Max Speed mm/s')
% ylabel('error along the reach direction (vector projection) in mm')
% title('Max v.s. error along the reach direction, each point = one trial')