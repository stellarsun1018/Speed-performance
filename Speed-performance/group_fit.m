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
clear; clc;

participants = {"JHL","JH","LC", "LN", "RC", "SM", "ML", "SX", "SY", "WMZ", "YL", "LL"};

% 如果你想按 session 读（像你上面 SX 那样），就用这一段 pattern：
use_session = false;   % true: 用 S%d 方式；false: 用 usable 方式
session = 3;

nBlocks = 3;
nTrialsPerBlock = 240;
lim_scale = 1.2;

pixellength = 0.248;

colors = lines(numel(participants));

% 存回归参数：block x (intercept,slope) x subject
fits_speed = NaN(nBlocks, 2, numel(participants));

% 用于统一坐标轴范围与 condition line
maxDistBlock  = zeros(nBlocks,1);
maxSpdBlock   = zeros(nBlocks,1);
lifeBlockSubs = NaN(nBlocks, numel(participants));   % 每个block的lifespan(列3)，每个受试者记录一个

%% ====== Loop subjects: load -> preprocess -> fit ======
for ip = 1:numel(participants)
    part = participants{ip};

    % --- find files ---
    if use_session
        fname_preamble = sprintf('data_onlineConf/%s/%s_sptatialTemporalCostFunc_S%d*.mat', part, part, session);
    else
        fname_preamble = sprintf('data_onlineConf/usable/%s_sptatialTemporalCostFunc*.mat', part);
    end
    files = dir(fname_preamble);
    if isempty(files)
        warning('No files found for %s using pattern: %s', part, fname_preamble);
        continue;
    end

    % --- load & (optional) concat if multiple files ---
    data_all = [];
    traX_all = [];
    traY_all = [];
    tform_ref = [];
    for k = 1:numel(files)
        f = fullfile(files(k).folder, files(k).name);
        S = load(f);

        if isfield(S,'data'),      data_all = [data_all; S.data]; end %#ok<AGROW>
        if isfield(S,'traXtotal'), traX_all = [traX_all; S.traXtotal]; end %#ok<AGROW>
        if isfield(S,'traYtotal'), traY_all = [traY_all; S.traYtotal]; end %#ok<AGROW>

        if isfield(S,'tform') && isempty(tform_ref)
            tform_ref = S.tform;
        end
    end
    if isempty(data_all) || isempty(tform_ref)
        warning('Missing data/tform for %s. Skip.', part);
        continue;
    end

    % --- valid trials (等价于你写法的向量化版本) ---
    index = ~isnan(sum(data_all,2));
    valid = data_all(index==true,:);

    % ====== Build "copy" exactly like your code ======
    copy = valid;

    tform = tform_ref; %#ok<NASGU>

    Affine2d = tform_ref.T(1:2,1:2);
    [~,s,~] = svd(Affine2d);
    proj2tablet = 1./mean([s(1,1),s(2,2)]);
    mmPerProjPx = proj2tablet .* pixellength;

    copy(:,[1,2]) = transformPointsInverse(tform_ref, copy(:,[1,2]));
    copy(:,[8,9]) = transformPointsInverse(tform_ref, copy(:,[8,9]));
    copy(:,10) = sqrt(sum((copy(:,1:2) - copy(:,8:9)).^2,2)) .* pixellength;

    copy(:,[11,12]) = [copy(:,1)*pixellength, (1080 - copy(:,2))*pixellength];
    copy(:,[13,14]) = [copy(:,6)*pixellength, (1080 - copy(:,7))*pixellength];

    copy(:,15) = valid(:,10) .* mmPerProjPx;
    copy(:,16) = copy(:,5) - copy(:,4);
    copy(:,17) = sqrt( (copy(:,13)-copy(:,11)).^2 + (copy(:,14)-copy(:,12)).^2 );

    copy(:,27) = (1:size(copy,1))';

    copy(:,19:20) = (copy(:,6:7) - copy(:,8:9)) .* pixellength;
    copy(:,21) = sqrt(sum((copy(:,6:7) - copy(:,8:9)).^2,2)) .* pixellength;
    copy(:,22) = copy(:,21) ./ copy(:,16);

    copy(:,[24,25]) = (copy(:,1:2) - copy(:,8:9)) .* pixellength;
    copy(:,23) = (abs(dot(copy(:,19:20),copy(:,24:25),2) ./ dot(copy(:,24:25),copy(:,24:25),2)) - 1) .*copy(:,10);

    copy(:,26) = valid(:,11);
    copy(:,28) = copy(:,26) ~= 0;

    % signed rejection length -> copy(:,29)
    projScale = dot(copy(:,19:20), copy(:,24:25), 2) ./ dot(copy(:,24:25), copy(:,24:25), 2);
    projections = projScale .* copy(:,24:25);
    rejections = copy(:,19:20) - projections;
    rejLength = sign(copy(:,19).*copy(:,25) - copy(:,20).*copy(:,24)) .* sqrt(rejections(:,1).^2 + rejections(:,2).^2);
    copy(:,29) = rejLength;

    % polar coordination (not required for regression lines, but keep consistent)
    angle_error_in_radian = asin(copy(:,29) ./ copy(:,21));
    copy(:,33) = angle_error_in_radian;
    copy(:,34) = rad2deg(angle_error_in_radian);

    % ====== Fit per block: speed vs distance ======
    N = size(copy,1);
    for b = 1:nBlocks
        i1 = 1 + (b-1)*nTrialsPerBlock;
        i2 = b*nTrialsPerBlock;
        if i1 > N, continue; end
        i2 = min(i2, N);

        block_ind = false(N,1);
        block_ind(i1:i2) = true;

        distances = copy(block_ind, 10);
        avg_speed = copy(block_ind, 22);

        if isempty(distances) || all(isnan(distances)) || all(isnan(avg_speed))
            continue;
        end

        mdl = fitlm(distances, avg_speed);
        fits_speed(b,1,ip) = mdl.Coefficients.Estimate(1); % intercept
        fits_speed(b,2,ip) = mdl.Coefficients.Estimate(2); % slope

        maxDistBlock(b) = max(maxDistBlock(b), max(distances,[],'omitnan'));
        maxSpdBlock(b)  = max(maxSpdBlock(b),  max(avg_speed,[],'omitnan'));

        % 记录 lifespan（copy(:,3)）用于 condition line
        life_u = unique(copy(block_ind,3));
        life_u = life_u(~isnan(life_u));
        if ~isempty(life_u)
            lifeBlockSubs(b,ip) = life_u(1);
        end
    end
end

%% ====== Plot: 3 subplots, overlay 12 regression lines ======
figure('Color','w');
for b = 1:nBlocks
    subplot(1,3,b); hold on;

    x_max = lim_scale * maxDistBlock(b);
    if x_max <= 0
        title(sprintf('Block %d (no data)', b));
        continue;
    end
    x_fit = linspace(0, x_max, 2);

    % --- subject regression lines ---
    y_max_fit = 0;
    for ip = 1:numel(participants)
        a = fits_speed(b,1,ip);
        s = fits_speed(b,2,ip);
        if isnan(a) || isnan(s), continue; end

        y_fit = a + s .* x_fit;
        plot(x_fit, y_fit, '--', 'LineWidth', 2, ...
            'Color', colors(ip,:), 'DisplayName', participants{ip});

        y_max_fit = max(y_max_fit, max(y_fit));
    end

    % --- condition lines (use median lifespan across subjects for this block) ---
    life = median(lifeBlockSubs(b, ~isnan(lifeBlockSubs(b,:))), 'omitnan');
    if ~isnan(life) && life > 0
        y_cond = x_fit ./ life;
        plot(x_fit, y_cond, '--', 'Color', [0 0.4 0.8], 'LineWidth', 1.5, ...
            'DisplayName', 'Min speed (match disappear deadline)');

        y_half = x_fit ./ (0.5*life);
        plot(x_fit, y_half, '--r', 'LineWidth', 1.5, ...
            'DisplayName', 'Min speed (half-size target)');
    end

    xlabel("Target Distance (mm)");
    ylabel("Average Speed (mm/s)");
    title(sprintf('Block %d', b));

    % --- axes ---
    y_max = max([maxSpdBlock(b)*lim_scale, y_max_fit, ...
                 (exist('y_cond','var')*max(y_cond) + ~exist('y_cond','var')*0), ...
                 (exist('y_half','var')*max(y_half) + ~exist('y_half','var')*0)]);
    if y_max <= 0, y_max = 1; end
    xlim([0, x_max]);
    ylim([0, lim_scale * y_max]);

    grid on;
end

sgtitle("12 subjects: regression lines (speed vs distance) overlaid per block");

% legend 太大时可以改成 'bestoutside' 或者注释掉
legend('Location','northwest', 'NumColumns', 2);

%%
%%
clear all
% participant = 'SX';
lifespan = [0.6,0.6*3^(0.25),0.6*3^(0.5)]; %[1.0,0.6,0.8,0.4]; %[1.1,0.9,1.0,0.8,0.6,0.7] 设定各blocks中target的不同时长  %lifespan控制了受试者实际可用的、逐渐减少的目标"可见e时间窗，这一时间越短，任务难度越高（因为受试者必须更快速地完成任务以取得更高分数）。
for i = 1:3
    copy((1+(i-1)*240):(i*240),3) = lifespan(i);
end
participants = {"JHL", "JH", "LC", "LN", "RC", "SM", "ML", "SX", "SY", "WMZ", "YL", "LL"};

linear_fits = NaN(4 * 3, numel(participants)); % three blocks * 4 parameters
% % 存储所有被试的拟合参数：linear_fits_group(block_index, param_index, subject_index)
% % param_index: 1=intercept, 2=slope (for duration vs distance)
% linear_fits_group = NaN(3, 2, numel(participants)); 
% 
% % 用于统一x轴范围
% all_max_distances = NaN(1, numel(participants));

for ip = 1:numel(participants)
    part = participants{ip};
    fname_preamble = sprintf('data_onlineConf/usable/%s_sptatialTemporalCostFunc*.mat',part);
    files = dir(fname_preamble);
    for k = 1:numel(files)
        f = fullfile(files(k).folder, files(k).name);
        load(f);
    end

    
    index = NaN(size(data,1),1);
    for i = 1:size(data,1)
        index(i) = ~isnan(sum(data(i,:)));
    end
    valid = data(index==true,:);
    validTraX = traXtotal(index==true,:);
    validTraY = traYtotal(index==true,:);

    
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
    % copy(copy(:,27)>=1 & copy(:,27)<=60, 31)   = 1.0;
    % copy(copy(:,27)>=61 & copy(:,27)<=120, 31) = 0.6;
    % copy(copy(:,27)>=121 & copy(:,27)<=180,31) = 0.8;
    % copy(copy(:,27)>=181 & copy(:,27)<=240,31) = 0.4;

    % 在copy中新增第32列：target shrinking speed (mm/s)
    % copy(:,32) = copy(:,15) ./ copy(:,31);

    % polar coordination

    angle_error_in_radian = asin(copy(:,29) ./ copy(:,21));

    angle_error_in_degree = rad2deg(angle_error_in_radian);

    copy(:,33) = angle_error_in_radian;
    copy(:,34) = angle_error_in_degree;

    % 3 blocks - conditions graphs

    lim_scale = 1.2;

    for i = 1:3
        cond_ind = i-1;
        block_ind = zeros(size(copy,1),1);
        block_ind((1+(i-1)*240):(i*240)) = 1;

        distances = copy(block_ind==1,10);
        avg_speed = copy(block_ind==1,22);
        duration = copy(block_ind==1,16);

        plot(distances,duration,'o');
        hold on

        mdl = fitlm(distances,avg_speed);
        linear_fits(1 + cond_ind * 4,ip) = mdl.Coefficients.Estimate(1);
        linear_fits(2 + cond_ind * 4,ip) = mdl.Coefficients.Estimate(2);


        mdl = fitlm(distances,duration);
        linear_fits(3 + cond_ind * 4,ip) = mdl.Coefficients.Estimate(1);
        linear_fits(4 + cond_ind * 4,ip) = mdl.Coefficients.Estimate(2);



    end

end
C = [participants; num2cell(linear_fits)];
fname = sprintf('results/group_fits_spd_dur_dist.csv');
writecell(C, fname)
% writematrix(processed, ['polarProcessed/polarProcessed' part 'S' num2str(sess) '.csv'])

%%

%%

% figure
% 
% for i = 1:numel(participants)
%     for i = 1:3
%         block_ind = zeros(size(copy,1),1);
%         block_ind((1+(i-1)*240):(i*240)) = 1;
%         subplot(1,3,i)
% 
%         distances = copy(block_ind==1,10);
%         avg_speed = copy(block_ind==1,22);
%         x = linspace(0,lim_scale*max(copy(:,10)),2);
%         y = x ./ unique(copy(block_ind==1,3));
% 
%         plot(distances,avg_speed,'o');
%         hold on
% 
%         mdl = fitlm(distances,avg_speed);
%         x_fit = linspace(0, lim_scale * max(copy(:,10)), 2);
%         y_fit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) .* x_fit;
%         plot(x_fit, y_fit, '--', 'LineWidth', 2);
% 
% 
% % 原始 condition line: speed = distance / lifespan
%         plot(x,y,'--','Color',[0 0.4 0.8],'LineWidth',1.5); % 蓝色
% 
% % 新增：半大小的 condition line（红色）
%         y_half = x ./ (0.5 * unique(copy(block_ind==1,3)));
%         plot(x, y_half, '--r','LineWidth',1.5); % 红色虚线
% 
%         hold off
%         xlabel("Target Distance (mm)");
%         ylabel("Average Speed (mm/s)");
%         legend('Trial data','Linear Fit','Min speed (match disappear deadline)', 'Min speed (half-size target)','Location','northwest')
%         xlim([0,lim_scale * max(copy(:,10))]);
%         ylim([0,lim_scale * max(copy(:,22))]);
%     end
% end


%% 画图：每个 block 的截距-斜率散点 + 每个被试跨 block 的连线
figure

subplot(1,2,1)
hold on
for i = 1:3
    cond_ind = i - 1;
    scatter(linear_fits(1 + cond_ind * 4,:), linear_fits(2 + cond_ind * 4,:),'filled')
end

for i = 1:numel(participants)
    selection = NaN(3,2);
    for j = 1:3
        cond_ind = j - 1;
        selection(j,:) = [linear_fits(1 + cond_ind * 4,i); linear_fits(2 + cond_ind * 4,i)];
    end
    plot(selection(:,1), selection(:,2),'r-')
end

hold off
xlabel("Intercept")
ylabel("Slope")
title("Average Speed vs Distance")
legend("Block 1", "Block 2", "Block 3")


subplot(1,2,2)
hold on
for i = 1:3
    cond_ind = i - 1;
    scatter(linear_fits(3 + cond_ind * 4,:), linear_fits(4 + cond_ind * 4,:),'filled')
end

for i = 1:numel(participants)
    selection = NaN(3,2);
    for j = 1:3
        cond_ind = j - 1;
        selection(j,:) = [linear_fits(3 + cond_ind * 4,i); linear_fits(4 + cond_ind * 4,i)];
    end
    plot(selection(:,1), selection(:,2),'r-')
end

hold off
xlabel("Intercept")
ylabel("Slope")
title("Duration vs Distance")
legend("Block 1", "Block 2", "Block 3")


% delete 
figure


