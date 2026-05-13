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
participants = {"JHL","JH","LC", "LN", "RC", "SM", "ML", "SX", "SY", "WMZ", "YL", "LL"};
% participants = {"JH"};
part_n = numel(participants);
use_session = false;   % true: 用 S%d 方式；false: 用 usable 方式
session = 3;
nBlocks = 3;
nTrialsPerBlock = 240;
lim_scale = 1.2;
pixellength = 0.248;
colors = lines(part_n);

% --- Initialize variables to aggregate all participants' data ---
all_reach_distances = [];
all_avg_speed = [];
all_duration = [];
all_lifespan = [];
all_hits = [];

%% ====== Loop subjects: load -> preprocess ======
for ip = 1:part_n
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
    
    % --- valid trials ---
    index = ~isnan(sum(data_all,2));
    valid = data_all(index==true,:);
    
    % ====== Build "copy" exactly like your code ======
    copy = valid;
    tform = tform_ref;
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
    copy(:,26) = valid(:,11);
    copy(:,28) = copy(:,26) ~= 0; % Hit or not (1 or 0)
    
    % ====== Extract Metrics ======
    reach_distances = copy(:,21);
    avg_speed = copy(:,22);
    duration = copy(:,16);
    hits = copy(:,28);
    lifespan = copy(:,3);
    
    % Aggregate data for the final plot
    all_reach_distances = [all_reach_distances; reach_distances];
    all_avg_speed = [all_avg_speed; avg_speed];
    all_duration = [all_duration; duration];
    all_hits = [all_hits; hits];
    all_lifespan = [all_lifespan; lifespan];
    
   
end

%% ====== Aggregate Plotting (Quantile Binning with Shaded SEM) ======
n_bins = 8; % Using 10 quantiles

% --- Distance Binning ---
edges_dist = quantile(all_reach_distances, linspace(0, 1, n_bins+1));
edges_dist(1) = min(all_reach_distances) - eps; 
edges_dist(end) = max(all_reach_distances) + eps;

binID_dist = discretize(all_reach_distances, edges_dist);
valid_dist = ~isnan(binID_dist);

% Calculate arrays for Distance
mean_dist_centers = accumarray(binID_dist(valid_dist), all_reach_distances(valid_dist), [n_bins 1], @mean, NaN);
mean_hit_dist = accumarray(binID_dist(valid_dist), all_hits(valid_dist), [n_bins 1], @mean, NaN);
std_hit_dist = accumarray(binID_dist(valid_dist), all_hits(valid_dist), [n_bins 1], @std, NaN);
n_dist_counts = accumarray(binID_dist(valid_dist), 1, [n_bins 1]);

% Standard Error of the Mean (SEM) for Hit Rate
sem_hit_dist = std_hit_dist ./ sqrt(n_dist_counts);

% Calculate polygon coordinates for Distance
x_dist_fill = [mean_dist_centers(:)', fliplr(mean_dist_centers(:)')];
y_dist_fill = [(mean_hit_dist(:) + sem_hit_dist(:))', fliplr((mean_hit_dist(:) - sem_hit_dist(:))')];


% --- Speed Binning ---
edges_speed = quantile(all_avg_speed, linspace(0, 1, n_bins+1));
edges_speed(1) = min(all_avg_speed) - eps;
edges_speed(end) = max(all_avg_speed) + eps;

binID_speed = discretize(all_avg_speed, edges_speed);
valid_speed = ~isnan(binID_speed);

% Calculate arrays for Speed
mean_speed_centers = accumarray(binID_speed(valid_speed), all_avg_speed(valid_speed), [n_bins 1], @mean, NaN);
mean_hit_speed = accumarray(binID_speed(valid_speed), all_hits(valid_speed), [n_bins 1], @mean, NaN);
std_hit_speed = accumarray(binID_speed(valid_speed), all_hits(valid_speed), [n_bins 1], @std, NaN);
n_speed_counts = accumarray(binID_speed(valid_speed), 1, [n_bins 1]);

% Standard Error of the Mean (SEM) for Hit Rate
sem_hit_speed = std_hit_speed ./ sqrt(n_speed_counts);

% Calculate polygon coordinates for Speed
x_speed_fill = [mean_speed_centers(:)', fliplr(mean_speed_centers(:)')];
y_speed_fill = [(mean_hit_speed(:) + sem_hit_speed(:))', fliplr((mean_hit_speed(:) - sem_hit_speed(:))')];

% --- Duration Binning ---
edges_duration = quantile(all_duration, linspace(0, 1, n_bins+1));
edges_duration(1) = min(all_duration) - eps;
edges_duration(end) = max(all_duration) + eps;

binID_duration = discretize(all_duration, edges_duration);
valid_duration = ~isnan(binID_duration);

% Calculate arrays for Speed
mean_duration_centers = accumarray(binID_duration(valid_duration), all_duration(valid_duration), [n_bins 1], @mean, NaN);
mean_hit_duration = accumarray(binID_duration(valid_duration), all_hits(valid_duration), [n_bins 1], @mean, NaN);
std_hit_duration = accumarray(binID_duration(valid_duration), all_hits(valid_duration), [n_bins 1], @std, NaN);
n_duration_counts = accumarray(binID_duration(valid_duration), 1, [n_bins 1]);

% Standard Error of the Mean (SEM) for Hit Rate
sem_hit_duration = std_hit_duration ./ sqrt(n_duration_counts);

% Calculate polygon coordinates for Speed
x_duration_fill = [mean_duration_centers(:)', fliplr(mean_duration_centers(:)')];
y_duration_fill = [(mean_hit_duration(:) + sem_hit_duration(:))', fliplr((mean_hit_duration(:) - sem_hit_duration(:))')];

%%
% --- Plotting ---
figure();
plotColor = [0.8500, 0.3250, 0.0980]; % MATLAB orange, distinguishing from previous script

% Plot 1: Hit Rate vs. Reach Distance
subplot(1, 2, 1);
hold on;
% fill(x_dist_fill, y_dist_fill, plotColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(mean_dist_centers, mean_hit_dist, 'o-', 'Color', plotColor, 'LineWidth', 2, 'MarkerFaceColor', plotColor);
hold off;
xlabel('Reach distance (mm)'); 
ylabel('Hit Rate (proportion)');
title('Hit Rate by Distance Quantiles');
ylim([0 1]); % Hit rate is strictly between 0 and 1
grid on;

% Plot 2: Hit Rate vs. Average Speed
subplot(1, 2, 2);
hold on;
% fill(x_speed_fill, y_speed_fill, plotColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(mean_speed_centers, mean_hit_speed, 'o-', 'Color', plotColor, 'LineWidth', 2, 'MarkerFaceColor', plotColor);
hold off;
xlabel('Average speed'); 
ylabel('Hit Rate (proportion)');
title('Hit Rate by Speed Quantiles');
ylim([0 1]); % Hit rate is strictly between 0 and 1
grid on;


% Plot 3: Hit Rate vs. Duration
figure()
% subplot(1, 3, 3);
hold on;
% fill(x_duration_fill, y_duration_fill, plotColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(mean_duration_centers, mean_hit_duration, 'o-', 'Color', plotColor, 'LineWidth', 2, 'MarkerFaceColor', plotColor);
hold off;
xlabel('Duration (s)'); 
ylabel('Hit Rate (proportion)');
title('Hit Rate by Duration Quantiles');
ylim([0.5 0.9]); % Hit rate is strictly between 0 and 1
grid on;


%%
plotColor = [0.8500, 0.3250, 0.0980];
% Extract data based on lifespan conditions
hit_per_block(:,1) = all_hits(round(all_lifespan,2) == 0.60);
hit_per_block(:,2) = all_hits(round(all_lifespan,2) == 0.79);
hit_per_block(:,3) = all_hits(round(all_lifespan,2) == 1.04);

% Calculate Mean and SEM
hit_rate_per_block = mean(hit_per_block, 1);
sem_hit_rate_per_block = std(hit_per_block, [], 1) ./ sqrt(size(hit_per_block, 1));

% --- Plotting ---
figure('Name', 'Hit Rate per Block');
subplot(1,2,1)
hold on;

x_pos = 1:3;
x_labels = {'Fast', 'Medium', 'Slow'};

% Create the bar chart
b = bar(x_pos, hit_rate_per_block, 'FaceColor', 'flat','FaceAlpha',0.7);

% Apply MATLAB's default blue, orange, and yellow colors
b.CData(1, :) = [0, 0.4470, 0.7410];       % Default Blue
b.CData(2, :) = [0.8500, 0.3250, 0.0980];  % Default Orange
b.CData(3, :) = [0.9290, 0.6940, 0.1250];  % Default Yellow

% Add the error bars
% 'k' makes the error bars black, 'none' removes the connecting lines
errorbar(x_pos, hit_rate_per_block, sem_hit_rate_per_block, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Formatting the axes
set(gca, 'XTick', x_pos, 'XTickLabel', x_labels);
ylabel('Hit Rate');
xlabel('Condition')
ylim([0.5 0.9]); % Clamp y-axis to logical bounds for a hit rate
grid off;
hold off;

subplot(1, 2, 2);
hold on;
% fill(x_duration_fill, y_duration_fill, plotColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(mean_duration_centers, mean_hit_duration, 'o-', 'Color', plotColor, 'LineWidth', 2, 'MarkerFaceColor', plotColor);
hold off;
xlim([0.22,0.65])
xlabel('Duration (s)'); 
ylabel('Hit Rate');
ylim([0.5 0.9]); % Hit rate is strictly between 0 and 1
grid off;