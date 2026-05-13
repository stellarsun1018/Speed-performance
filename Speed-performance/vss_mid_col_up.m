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
clear; 
clc;

participants = {"JHL","JH","LC", "LN", "RC", "SM", "ML", "LL", "SY", "WMZ", "YL", "SX"};

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

%% ====== Plot: 3 subplots, overlay 12 regression lines ====== version 1
figure('Color','w');

for i = 1:3
    block_ind = zeros(size(copy,1),1);
    block_ind((1+(i-1)*240):(i*240)) = 1;
    
    subplot(2,3,i)

    distances = copy(block_ind==1,10);
    avg_speed = copy(block_ind==1,22);
    x = linspace(0,380,2);
    y = x ./ unique(copy(block_ind==1,3));
    
    % plot(distances,avg_speed,'o');
    scatter(distances, avg_speed, 25, colors(i,:), 'filled', "MarkerEdgeColor","none",'MarkerFaceAlpha', 0.7);
    hold on

    mdl = fitlm(distances,avg_speed);
    x_fit = linspace(0, 380, 2);
    y_fit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) .* x_fit;
    plt = plot(x_fit, y_fit, 'k-', 'LineWidth', 6);
    % plt.Color(4) = 0.3;


% 原始 condition line: speed = distance / lifespan
    plot(x,y,'--','Color', "#00923a",'LineWidth',6); % 蓝色

% 新增：半大小的 condition line（红色）
    y_half = x ./ (0.5 * unique(copy(block_ind==1,3)));
    plot(x, y_half, '--','Color', "#8fc31f",'LineWidth',6); % 红色虚线

    hold off
    xlabel("Target Distance (mm)");
    ylabel("Average Speed (mm/s)");
    % title(blockNames{i})
    % legend('Data','Linear Fit','Disappear deadline', 'Half-target time','Location','southeast')
    xlim([0,380]);
    xticks(0:100:300)
    ylim([0,600]);
    yticks(0:200:600)

    % hold off
    % xlabel("Target Distance (mm)");
    % ylabel("Average Speed (mm/s)");
    % legend('Trial data','Minimum speed')
    % xlim([0,lim_scale * max(copy(:,10))]);
    % ylim([0,lim_scale * max(copy(:,22))]);
end

for b = 1:nBlocks
    subplot(2,3,b+3); hold on;

    x_max = lim_scale * maxDistBlock(b);
    if x_max <= 0
        title(sprintf('Block %d (no data)', b));
        continue;
    end
    x_fit = linspace(0, 380, 2);

    % --- subject regression lines ---
    y_max_fit = 0;
    for ip = 1:numel(participants)
        a = fits_speed(b,1,ip);
        s = fits_speed(b,2,ip);
        if isnan(a) || isnan(s), continue; end

        y_fit = a + s .* x_fit;
        plt = plot(x_fit, y_fit, '-k', 'LineWidth', 6, 'DisplayName', participants{ip});
        % plt.Color(4) = 0.3;
        y_max_fit = max(y_max_fit, max(y_fit));
    end

    % --- condition lines (use median lifespan across subjects for this block) ---
    life = median(lifeBlockSubs(b, ~isnan(lifeBlockSubs(b,:))), 'omitnan');
    if ~isnan(life) && life > 0
        y_cond = x_fit ./ life;
        plot(x_fit, y_cond, '--', 'Color', "#00923a", 'LineWidth', 6, ...
            'DisplayName', 'Min speed (match disappear deadline)');

        y_half = x_fit ./ (0.5*life);
        plot(x_fit, y_half, '--','Color', "#8fc31f", 'LineWidth', 6, ...
            'DisplayName', 'Min speed (half-size target)');
    end

    xlabel("Target Distance (mm)");
    ylabel("Average Speed (mm/s)");

    % --- axes ---
    y_max = max([maxSpdBlock(b)*lim_scale, y_max_fit, ...
                 (exist('y_cond','var')*max(y_cond) + ~exist('y_cond','var')*0), ...
                 (exist('y_half','var')*max(y_half) + ~exist('y_half','var')*0)]);
    if y_max <= 0, y_max = 1; end
    xlim([0, 380]);
    xticks(0:100:300)
    ylim([0, 1000]);

    % grid on;
end

% sgtitle("12 subjects: regression lines (speed vs distance) overlaid per block");

% legend 太大时可以改成 'bestoutside' 或者注释掉
% legend('Location','northwest', 'NumColumns', 2);

%%


