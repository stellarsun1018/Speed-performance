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
clear all
participant = 'SX';
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



% 创建第30列："target shrinking duration"
copy(:,31) = NaN; % 初始化为NaN

% 根据 trial order 赋值
copy(copy(:,27)>=1 & copy(:,27)<=60, 31)   = 1.0;
copy(copy(:,27)>=61 & copy(:,27)<=120, 31) = 0.6;
copy(copy(:,27)>=121 & copy(:,27)<=180,31) = 0.8;
copy(copy(:,27)>=181 & copy(:,27)<=240,31) = 0.4;

% 在copy中新增第32列：target shrinking speed (mm/s)
copy(:,32) = copy(:,15) ./ copy(:,31);


%% 3 blocks - conditions graphs SPEED

lim_scale = 1.2;
colors = lines(3); % 蓝/红/黄 metrics

blockNames = {'Block 1 (short distances)', ...
              'Block 2 (medium distances)', ...
              'Block 3 (long distances)'};

figure()
for i = 1:3
    block_ind = zeros(size(copy,1),1);
    block_ind((1+(i-1)*240):(i*240)) = 1;
    
    subplot(1,3,i)

    distances = copy(block_ind==1,10);
    avg_speed = copy(block_ind==1,22);
    x = linspace(0,lim_scale*max(copy(:,10)),2);
    y = x ./ unique(copy(block_ind==1,3));
    
    % plot(distances,avg_speed,'o');
    scatter(distances, avg_speed, 25, colors(i,:), 'o');
    hold on

    mdl = fitlm(distances,avg_speed);
    x_fit = linspace(0, lim_scale * max(copy(:,10)), 2);
    y_fit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) .* x_fit;
    plot(x_fit, y_fit, 'k-', 'LineWidth', 1);


% 原始 condition line: speed = distance / lifespan
    plot(x,y,'--','Color',[0 0.4 0.8],'LineWidth',1.5); % 蓝色

% 新增：半大小的 condition line（红色）
    y_half = x ./ (0.5 * unique(copy(block_ind==1,3)));
    plot(x, y_half, '--r','LineWidth',1.5); % 红色虚线

    hold off
    xlabel("Target Distance (mm)");
    ylabel("Average Speed (mm/s)");
    title(blockNames{i})
    legend('Data','Linear Fit','Disappear deadline', 'Half-target time','Location','southeast')
    xlim([0,lim_scale * max(copy(:,10))]);
    ylim([0,lim_scale * max(copy(:,22))]);

    % hold off
    % xlabel("Target Distance (mm)");
    % ylabel("Average Speed (mm/s)");
    % legend('Trial data','Minimum speed')
    % xlim([0,lim_scale * max(copy(:,10))]);
    % ylim([0,lim_scale * max(copy(:,22))]);
end

% sgtitle("Arriving earlier for shorter distances and later for longer distances")
% 要改颜色
% assert(size(copy,1) >= 3*240, '预期每个block约240 trial，请确认数据尺寸。');
% 
% trialsPerBlock = 240;
% blocks = 1:3;
% 
% % 颜色（与要求一致：红/蓝/黄）
% blockColors = [1 0 0; 0 0.45 0.95; 1 0.9 0]; % red, blue, yellow

%%
%% 3 blocks - conditions graphs DURATION

lim_scale = 1.2;
figure()
for i = 1:3
    block_ind = zeros(size(copy,1),1);
    block_ind((1+(i-1)*240):(i*240)) = 1;
    subplot(1,3,i)

    distances = copy(block_ind==1,10);
    duration = copy(block_ind==1,16);

    plot(distances,duration,'o');
    hold on

    mdl = fitlm(distances,duration);
    x_fit = linspace(0, lim_scale * max(copy(:,10)), 2);
    y_fit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) .* x_fit;
    plot(x_fit, y_fit, '--', 'LineWidth', 2);

    yline(unique(copy(block_ind==1,3)).*0.5,'r--')
    yline(unique(copy(block_ind==1,3)),'b--')
    hold off
    xlabel("Target Distance (mm)");
    ylabel("Duration (s)");
    legend('Trial data','Linear Fit','Disappear Deadline', 'Target Half Life','Location','northwest')
    xlim([0,lim_scale * max(copy(:,10))]);
    ylim([0,lim_scale * max(copy(:,16))]);

    % hold off
    % xlabel("Target Distance (mm)");
    % ylabel("Average Speed (mm/s)");
    % legend('Trial data','Minimum speed')
    % xlim([0,lim_scale * max(copy(:,10))]);
    % ylim([0,lim_scale * max(copy(:,22))]);
end

sgtitle("Arriving earlier for shorter distances and later for longer distances")

%% Overlay all 3 blocks with regression lines through origin
figure;
colors = lines(3);  % 三种颜色
% blockColors = [0 0.45 0.95;   % blue
%                1 0 0;        % red
%                1 0.9 0];     % yellow

% 初始化 legend 对象数组
h_scatter = gobjects(3,1);
h_line = gobjects(3,1);

for i = 1:3
    block_inds = (1+(i-1)*240):(i*240);
    distances = copy(block_inds,10);
    avg_speed = copy(block_inds,22);

    % 画散点图，并保存 scatter 句柄
    h_scatter(i) = scatter(distances, avg_speed, 20, colors(i,:), 'filled'); hold on;

    % 线性回归：
    mdl = fitlm(distances,avg_speed);
    x_fit = linspace(0, max(distances), 100);
    y_fit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) .* x_fit;

    % 虚线拟合线，并保存句柄
    h_line(i) = plot(x_fit, y_fit, 'k-', 'LineWidth', 1.3);

    % 显示 slope
    fprintf('Block %d: slope = %.2f, implied lifespan = %.2f sec\n', i, k, 1/k);
end

xlim([0,lim_scale * max(copy(:,10))]);
ylim([0,lim_scale * max(copy(:,22))]);

xlabel('Target Distance (mm)');
ylabel('Average Speed (mm/s)');
% title('Distance vs Speed (3 Blocks)');

% 用回归线句柄生成 legend（不是散点）
h_fit_legend = plot(nan, nan, 'k-', 'LineWidth', 2); % dummy line for legend

legend([h_scatter(1), h_scatter(2), h_scatter(3), h_fit_legend], ...
       {'Block 1 (fast shrinking)', ...
        'Block 2 (medium shrinking)', ...
        'Block 3 (slow shrinking)', ...
        'Regression fit lines'}, ...
       'Location', 'northwest');

grid on;

%%  Duration histograms by block (overlapped) + disappearance time lines
% copy(:,16) = duration(s); copy(:,3) = lifespan(s)（每个block恒定）

assert(size(copy,1) >= 3*240, '预期每个block约240 trial，请确认数据尺寸。');

trialsPerBlock = 240;
blocks = 1:3;

% 颜色（与要求一致：红/蓝/黄）
blockColors = [0 0.45 0.95; 1 0 0; 1 0.9 0]; % red, blue, yellow

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
% edges = linspace(prctile(allDur,1), prctile(allDur,99), 30);
sample_rate = 1/60;
edges = min(allDur):sample_rate:max(allDur);


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
xlabel('Reach Duration (s)');
ylabel('Density');
% title('Overlapped Duration Histograms with Disappear Times');
grid on;

% 图例：区分直方图与其对应的消失时间
lgd = legend( ...
    [hH(1) hH(2) hH(3) hV(1) hV(2) hV(3)], ...
    {'Block 1 (short distance)','Block 2 (medium distance)','Block 3 (long distance)', ...
     sprintf('Block 1 disappear time %.2fs', Tdisappear(1)), ...
     sprintf('Block 2 disappear time %.2fs', Tdisappear(2)), ...
     sprintf('Block 3 disappear time %.2fs', Tdisappear(3))}, ...
    'Location','northwest');
lgd.Box = 'off';
hold off;

%% Compare mean durations across blocks (bar + errorbars) + print stats
% !!!!!!!add error bar to this, make alpha = 0.7 or something for it's more
% transparent. also do t-tests against each pair (3! tests) and illustrate
% significant pairs with horizontal brackets and star signs for p-value
% error bar use SEM
% meansT = nan(1,3);
% stdT   = nan(1,3);
% nT     = nan(1,3);
% semT   = nan(1,3);
% 
% for i = blocks
%     di = Dur{i};
%     di = di(isfinite(di));
%     meansT(i) = mean(di);
%     stdT(i)   = std(di);
%     nT(i)     = numel(di);
%     semT(i)   = stdT(i)/sqrt(max(nT(i),1));
% end
% 
% % 柱状图 （+ SEM 误差条）
% figure('Name','Mean Duration by Block');
% bh = bar(1:3, meansT, 'FaceColor','flat'); 
% for i = blocks
%     bh.CData(i,:) = blockColors(i,:);
% end
% % hold on;
% % errorbar(1:3, meansT, semT, 'k', 'LineStyle','none', 'LineWidth',1.5, 'CapSize',8);
% 
% % 在柱子上方写均值
% for i = 1:3
%     text(i, meansT(i) + 0.005, sprintf('%.3f', meansT(i)), ...
%         'HorizontalAlignment','center', 'FontSize', 12, 'FontWeight','bold');
% end
% 
% xticks(1:3); xticklabels({'Block 1','Block 2','Block 3'});
% ylabel('Mean Duration (s)');
% title('Average Reach Duration per Block');
% grid on; box off;
% 
% % 控制台打印统计
% fprintf('\n=== Duration Summary by Block ===\n');
% for i = blocks
%     fprintf('Block %d: N=%d, Mean=%.4f s, SD=%.4f s, SEM=%.4f s, Disappear=%.4f s\n', ...
%         i, nT(i), meansT(i), stdT(i), semT(i), Tdisappear(i));
% end
% 
% %% 3D reg plane - Euclidean Error 
% % Extract variables
% reach_distances = copy(:,21);
% avg_speed = copy(:,22);
% errors = copy(:,17);
% 
% % Prepare the design matrix (adding a column of ones for intercept)
% X = [ones(size(reach_distances)), reach_distances, avg_speed];
% 
% % Perform multivariate linear regression
% coeffs = regress(errors, X);
% 
% % Display regression coefficients
% fprintf('Intercept: %.4f\n', coeffs(1));
% fprintf('Distance coefficient: %.4f\n', coeffs(2));
% fprintf('Average speed coefficient: %.4f\n', coeffs(3));
% 
% % 3D scatter plot of the original data
% figure;
% plot3(reach_distances, avg_speed, errors, 'o');
% hold on;
% 
% xlim([0 max(reach_distances)]);
% ylim([0 max(avg_speed)]);
% 
% % Generate grid for regression plane
% [dist_grid, speed_grid] = meshgrid(linspace(min(reach_distances), max(reach_distances), 20), ...
%                                    linspace(min(avg_speed), max(avg_speed), 20));
% 
% % Predict errors using regression model
% error_fit = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*speed_grid;
% 
% % Plot regression plane
% mesh(dist_grid, speed_grid, error_fit);
% xlabel('Reach Distance (mm)');
% ylabel('Average Speed (mm/s)');
% zlabel('Endpoint Euclidean Error (mm)');
% title('Multivariate Linear Regression: Error ~ Reach Distance + Speed');
% grid on;
% hold off;

% Compute mean, std, SEM
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

%% === Bar Plot with SEM Error Bars ===
figure('Name','Mean Duration by Block');
bh = bar(1:3, meansT, 'FaceColor','flat', 'FaceAlpha', 0.7); % 透明度0.7

for i = blocks
    bh.CData(i,:) = blockColors(i,:);
end
hold on;

% 添加误差条（SEM）
er = errorbar(1:3, meansT, semT, 'k', 'LineStyle','none', ...
    'LineWidth',1.5, 'CapSize',8, 'Color',[0 0 0 0.8]); % 半透明黑色误差线

% 在柱子上方标出均值
for i = 1:3
    text(i, meansT(i) + 0.005, sprintf('%.3f', meansT(i)), ...
        'HorizontalAlignment','center', 'FontSize', 12, 'FontWeight','bold');
end

xticks(1:3);
xticklabels({'Block 1','Block 2','Block 3'});
ylabel('Mean Duration (s)');
title('Average Reach Duration per Block');
grid on; box off;

%% === T-tests among the 3 pairs ===
pairs = [1 2; 1 3; 2 3];
y_max = max(meansT + semT) + 0.02; % 起始高度
increment = 0.02; % 每条横线间隔高度

p_values = nan(1,3);

for p = 1:3
    a = pairs(p,1);
    b = pairs(p,2);
    [~, pval] = ttest2(Dur{a}, Dur{b});
    p_values(p) = pval;

    % 画显著性横线
    x1 = a; x2 = b;
    y = y_max + (p-1)*increment;
    plot([x1 x2], [y y], 'k', 'LineWidth', 1.5); % 横线
    plot([x1 x1], [y-0.005 y], 'k', 'LineWidth', 1);
    plot([x2 x2], [y-0.005 y], 'k', 'LineWidth', 1);

    % 根据显著性画星号
    if pval < 0.001
        stars = '***';
    elseif pval < 0.01
        stars = '**';
    elseif pval < 0.05
        stars = '*';
    else
        stars = 'n.s.';
    end
    text(mean([x1 x2]), y + 0.005, stars, 'HorizontalAlignment', 'center', ...
         'FontSize', 14, 'FontWeight', 'bold');
end

% 控制台打印结果
fprintf('\n=== Duration Summary by Block ===\n');
for i = blocks
    fprintf('Block %d: N=%d, Mean=%.4f s, SD=%.4f s, SEM=%.4f s\n', ...
        i, nT(i), meansT(i), stdT(i), semT(i));
end

fprintf('\n=== Pairwise t-tests ===\n');
for p = 1:3
    fprintf('Block %d vs Block %d: p = %.4g\n', pairs(p,1), pairs(p,2), p_values(p));
end

%% Duration

%% Z axis: Gain error ~ reach distance + Duration copy(:,16) 
% Extract variables
reach_distances = copy(:,21);
reach_dur = copy(:,16);
Reach_errors = copy(:,23);

% Prepare the design matrix (adding a column of ones for intercept)
X = [ones(size(reach_distances)), reach_distances, reach_dur];

% Perform multivariate linear regression
coeffs = regress(Reach_errors, X);

% Display regression coefficients
fprintf('Intercept: %.4f\n', coeffs(1));
fprintf('Distance coefficient: %.4f\n', coeffs(2));
fprintf('Reach duration coefficient: %.4f\n', coeffs(3));

% 3D scatter plot of the original data
figure;
plot3(reach_distances, reach_dur, Reach_errors, 'o');
hold on;

xlim([0 max(reach_distances)]);
ylim([0 max(reach_dur)]);

% Generate grid for regression plane
[dist_grid, dur_grid] = meshgrid(linspace(min(reach_distances), max(reach_distances), 20), ...
                                   linspace(min(reach_dur), max(reach_dur), 20));

% Predict errors using regression model
error_fit = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*dur_grid;

% Plot regression plane
mesh(dist_grid, dur_grid, error_fit);
xlabel('Reach Distance (mm)');
ylabel('Reach Duration (s)');
zlabel('Gain Error (mm)');
title('Multivariate Linear Regression: Gain Error ~ Reach Distance + Reach Duration');
grid on;
hold off;

%% step 2: residuals 
errors_predicted = X * coeffs;

% Calculate residuals
residuals1 = Reach_errors - errors_predicted;

% 3D scatter plot of residuals
figure;
plot3(reach_distances, reach_dur, residuals1, 'ro');
xlim([0 max(reach_distances)]);
ylim([0 max(reach_dur)]);

xlabel('Reach Distance (mm)');
ylabel('Reach Duration (s)');
zlabel('Residuals (mm)');
title('Residuals of Multivariate Linear Regression: Reach Error ~ Distance + Duration');
grid on;

%% add Countours 加入等高线的residuals图 （version 1） 
errors_predicted = X * coeffs;

% Calculate residuals
residuals1 = Reach_errors - errors_predicted;

% 3D scatter plot of residuals
figure;
plot3(reach_distances, reach_dur, residuals1, 'ro');
hold on;

xlim([0 max(reach_distances)]);
ylim([0 max(reach_dur)]);

% ---- Build grid for residual contour ----
[dist_grid, dur_grid] = meshgrid( ...
    linspace(min(reach_distances), max(reach_distances), 40), ...
    linspace(min(reach_dur), max(reach_dur), 40));

% Interpolate residuals onto grid
residual_grid = griddata( ...
    reach_distances, reach_dur, residuals1, ...
    dist_grid, dur_grid, 'natural');

% ---- Add contour projection on bottom plane ----
z_level = min(residual_grid(:));  % 投影到最低 residual平面？还是说应该投影到零残差平面？
contour3(dist_grid, dur_grid, residual_grid, ...
         6, 'k', 'LineWidth', 1);

xlabel('Reach Distance (mm)');
ylabel('Reach Duration (s)');
zlabel('Residuals (mm)');
title('Residuals of Multivariate Linear Regression: Reach Error ~ Distance + Duration');
grid on;

view(0,90)
colorbar


%% add Countours 加入等高线的residuals图 （version 2） 热力图
errors_predicted = X * coeffs;

% Calculate residuals
residuals1 = Reach_errors - errors_predicted;

% ===== SD heatmap of residuals (adaptive bins) + contours =====
% Required inputs already exist:
%   reach_distances (mm)  -> x
%   reach_dur (s)         -> y
%   residuals1 (mm)       -> z  (treat as raw motor error residual)

% ---------- parameters you may tweak ----------
targetPerBin = 6;     % 统计上希望“平均每个bin大概有多少trial”(越大 -> bin越少 -> 更稳定)
minBin = 8;           % 最少每个维度多少个bin
maxBin = 18;          % 最多每个维度多少个bin（太大容易空bin -> 不好看也不稳）
minN = 3;             % 每个bin至少多少个trial才计算SD（<minN会先置NaN再插值）
doSmooth = true;      % 平滑让contour更顺
smoothSigma_bins = 1.0; % 平滑强度（单位：bin，1~1.5常用）
fineN = 220;          % 显示用的细网格分辨率（越大越细腻）

nLevels = 12;         % contour层数
overlayPoints = true; % 是否叠加原始点云

% ---------- clean vectors ----------
x = reach_distances(:);
y = reach_dur(:);
z = residuals1(:);

valid = isfinite(x) & isfinite(y) & isfinite(z);
x = x(valid); y = y(valid); z = z(valid);

n = numel(z);

% ---------- choose a sensible bin count automatically ----------
nb = round(sqrt(n / targetPerBin));   % e.g. n=720,targetPerBin=6 -> nb~11
nb = max(minBin, min(maxBin, nb));
nbx = nb; nby = nb;

% ---------- edges ----------
xedges = linspace(min(x), max(x), nbx+1);
yedges = linspace(min(y), max(y), nby+1);

% ---------- bin assignment ----------
[~,~,bx] = histcounts(x, xedges);
[~,~,by] = histcounts(y, yedges);
in = (bx>0) & (by>0);
bx = bx(in); by = by(in);
x  = x(in);  y  = y(in);  z  = z(in);

% ---------- per-bin N, sum, sumsq (stable SD) ----------
countGrid = accumarray([by, bx], 1,    [nby, nbx], @sum, 0);
sumGrid   = accumarray([by, bx], z,    [nby, nbx], @sum, 0);
sumSqGrid = accumarray([by, bx], z.^2, [nby, nbx], @sum, 0);

% unbiased sample variance
varGrid = (sumSqGrid - (sumGrid.^2)./max(countGrid,1)) ./ max(countGrid-1, 1);
varGrid(countGrid < minN) = NaN;
varGrid(varGrid < 0) = 0;     % numerical guard
sdGrid = sqrt(varGrid);

% 如果minN太严格导致可用bin太少，自动放宽到2
if nnz(isfinite(sdGrid)) < 12
    minN = 2;
    varGrid = (sumSqGrid - (sumGrid.^2)./max(countGrid,1)) ./ max(countGrid-1, 1);
    varGrid(countGrid < minN) = NaN;
    varGrid(varGrid < 0) = 0;
    sdGrid = sqrt(varGrid);
end

% ---------- bin centers ----------
xc = (xedges(1:end-1) + xedges(2:end)) / 2;   % 1×nbx
yc = (yedges(1:end-1) + yedges(2:end)) / 2;   % 1×nby
[Xc, Yc] = meshgrid(xc, yc);                  % nby×nbx

% ---------- fill empty bins so the heatmap is FULL (关键：铺满渐变颜色) ----------
mask = isfinite(sdGrid);
F = scatteredInterpolant(Xc(mask), Yc(mask), sdGrid(mask), 'natural', 'nearest');
sdFill = F(Xc, Yc);   % now full matrix (no NaN)

% ---------- smooth in bin-space (optional, makes contours nice) ----------
sdSmooth = sdFill;
if doSmooth
    ksz = max(3, 2*ceil(3*smoothSigma_bins)+1);
    r = (-floor(ksz/2)):(floor(ksz/2));
    [kx, ky] = meshgrid(r, r);
    K = exp(-(kx.^2 + ky.^2) / (2*smoothSigma_bins^2));
    K = K / sum(K(:));
    sdSmooth = conv2(sdSmooth, K, 'same');
end

% ---------- upsample to a fine grid for a smooth-looking heatmap ----------
xq = linspace(min(xc), max(xc), fineN);
yq = linspace(min(yc), max(yc), fineN);
[XX, YY] = meshgrid(xq, yq);

% cubic makes it look like your attached smooth gradients
sdFine = interp2(xc, yc, sdSmooth, XX, YY, 'cubic');

% ---------- plot: imagesc + contour (like your example) ----------
figure('Name','Residual SD heatmap (adaptive bins) + contours');

imagesc(xq, yq, sdFine);
set(gca, 'YDir', 'normal');   % axis xy
axis tight;
hold on;

contour(xq, yq, sdFine, nLevels, 'k', 'LineWidth', 1);

if overlayPoints
    scatter(x, y, 10, 'k', 'filled', ...
        'Marker', 's', 'MarkerFaceAlpha', 0.18, 'MarkerEdgeAlpha', 0.18);
end

cb = colorbar;
cb.Label.String = 'SD of residual (mm)';

% contrast: avoid "looks like no color"
cl = prctile(sdFine(:), [2 98]);
caxis(cl);

xlabel('Reach Distance (mm)');
ylabel('Reach Duration (s)');
title(sprintf('Residual SD heatmap (bins=%dx%d, targetPerBin=%g, minN=%d)', nbx, nby, targetPerBin, minN));
grid on;
hold off;



%% Z axis: Orthognal error ~ reach distance + Duration copy(:,16) 
% Extract variables
reach_distances = copy(:,21);
reach_dur = copy(:,16);
Orth_errors = copy(:,29);

% Prepare the design matrix (adding a column of ones for intercept)
X = [ones(size(reach_distances)), reach_distances, reach_dur];

% Perform multivariate linear regression
coeffs = regress(Orth_errors, X);

% Display regression coefficients
fprintf('Intercept: %.4f\n', coeffs(1));
fprintf('Distance coefficient: %.4f\n', coeffs(2));
fprintf('Reach duration coefficient: %.4f\n', coeffs(3));

% 3D scatter plot of the original data
figure;
plot3(reach_distances, reach_dur, Orth_errors, 'o');
hold on;

xlim([0 max(reach_distances)]);
ylim([0 max(reach_dur)]);

% Generate grid for regression plane
[dist_grid, dur_grid] = meshgrid(linspace(min(reach_distances), max(reach_distances), 20), ...
                                   linspace(min(reach_dur), max(reach_dur), 20));

% Predict errors using regression model
error_fit = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*dur_grid;

% Plot regression plane
mesh(dist_grid, dur_grid, error_fit);
xlabel('Reach Distance (mm)');
ylabel('Reach Duration (s)');
zlabel('Orthognal Error (mm)');
title('Multivariate Linear Regression: Orthognal Error ~ Reach Distance + Reach Duration');
grid on;
hold off;

%% step 2: residuals 
errors_predicted = X * coeffs;

% Calculate residuals
residuals1 = Orth_errors - errors_predicted;

% 3D scatter plot of residuals
figure;
plot3(reach_distances, reach_dur, residuals1, 'ro');
xlim([0 max(reach_distances)]);
ylim([0 max(reach_dur)]);

xlabel('Reach Distance (mm)');
ylabel('Reach Duration (s)');
zlabel('Residuals (mm)');
title('Residuals of Multivariate Linear Regression: Orthognal Error ~ Distance + Duration');
grid on;


%% Speed


%%  Z axis: Gain error ~ reach distance + speed (avg_speed = copy(:,22))
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

%% Z axis: Orthognal error ~ reach distance + speed (avg_speed = copy(:,22))
% Extract variables
reach_distances = copy(:,21);
avg_speed = copy(:,22);
Orth_errors = copy(:,29);

% Prepare the design matrix (adding a column of ones for intercept)
X = [ones(size(reach_distances)), reach_distances, avg_speed];

% Perform multivariate linear regression
coeffs = regress(Orth_errors, X);

% Display regression coefficients
fprintf('Intercept: %.4f\n', coeffs(1));
fprintf('Distance coefficient: %.4f\n', coeffs(2));
fprintf('Ave speed coefficient: %.4f\n', coeffs(3));

% 3D scatter plot of the original data
figure;
plot3(reach_distances, avg_speed, Orth_errors, 'o');
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
ylabel(''Average Speed (mm/s)'');
zlabel('Orthognal Error (mm)');
title('Multivariate Linear Regression: Orthognal Error ~ Reach Distance + Speed');
grid on;
hold off;

%% step 2: residuals 
errors_predicted = X * coeffs;

% Calculate residuals
residuals1 = Orth_errors - errors_predicted;

% 3D scatter plot of residuals
figure;
plot3(reach_distances, avg_speed, residuals1, 'ro');
xlim([0 max(reach_distances)]);
ylim([0 max(avg_speed)]);

xlabel('Reach Distance (mm)');
ylabel('Average Speed (mm/s)');
zlabel('Residuals (mm)');
title('Residuals of Multivariate Linear Regression: Orthognal Error ~ Distance + Speed');
grid on;


%%
%% Z axis: angular_error_in_degree ~ reach distance + speed
% 
% % Extract variables
% reach_distances = copy(:,21);
% avg_speed = copy(:,22);
% angle_error_in_degree = copy(:,34);
% 
% % Prepare the design matrix (adding a column of ones for intercept)
% X = [ones(size(reach_distances)), reach_distances, avg_speed];
% 
% % Perform multivariate linear regression
% coeffs = regress(angle_error_in_degree, X);
% 
% % Display regression coefficients
% fprintf('Intercept: %.4f\n', coeffs(1));
% fprintf('Distance coefficient: %.4f\n', coeffs(2));
% fprintf('Average speed coefficient: %.4f\n', coeffs(3));
% 
% % 3D scatter plot of the original data
% figure;
% plot3(reach_distances, avg_speed, angle_error_in_degree, 'o');
% hold on;
% 
% xlim([0 max(reach_distances)]);
% ylim([0 max(avg_speed)]);
% 
% % Generate grid for regression plane
% [dist_grid, speed_grid] = meshgrid(linspace(min(reach_distances), max(reach_distances), 20), ...
%                                    linspace(min(avg_speed), max(avg_speed), 20));
% 
% % Predict errors using regression model
% error_fit = coeffs(1) + coeffs(2)*dist_grid + coeffs(3)*speed_grid;
% 
% % Plot regression plane
% mesh(dist_grid, speed_grid, error_fit);
% xlabel('Reach Distance (mm)');
% ylabel('Average Speed (mm/s)');
% zlabel('Angular error (degree)');
% title('Multivariate Linear Regression: Angle error ~ Reach Distance + Speed');
% grid on;
% hold off;
% 
% %% step 2: Angular error residuals 
% errors_predicted = X * coeffs;
% 
% % Calculate residuals
% residuals2 = angle_error_in_degree - errors_predicted;
% 
% % 3D scatter plot of residuals
% figure;
% plot3(reach_distances, avg_speed, residuals2, 'ro');
% xlim([0 max(reach_distances)]);
% ylim([0 max(avg_speed)]);
% 
% xlabel('Reach Distance (mm)');
% ylabel('Average Speed (mm/s)');
% zlabel('Residuals (mm)');
% title('Residuals of Multivariate Linear Regression: Angular Error ~ Distance + Speed');
% grid on;
% 
