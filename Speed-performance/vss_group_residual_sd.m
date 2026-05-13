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
all_residuals = [];

%% ====== Loop subjects: load -> preprocess -> fit ======
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
    
    % polar coordination
    angle_error_in_radian = asin(copy(:,29) ./ copy(:,21));
    copy(:,33) = angle_error_in_radian;
    copy(:,34) = rad2deg(angle_error_in_radian);
    
    % ====== Fit per participant: dir_errors ======
    reach_distances = copy(:,21);
    avg_speed = copy(:,22);
    dir_errors = copy(:,29);
    
    % Prepare the design matrix (adding a column of ones for intercept)
    X = [ones(size(reach_distances)), reach_distances, avg_speed];
    
    % Perform multivariate linear regression for directional errors
    [coeffs, ~] = regress(dir_errors, X);
    errors_predicted = X * coeffs;
    residuals = dir_errors - errors_predicted;
    
    % Aggregate data for the final plot
    all_reach_distances = [all_reach_distances; reach_distances];
    all_avg_speed = [all_avg_speed; avg_speed];
    all_residuals = [all_residuals; residuals];
    
end

%% ====== Aggregate Plotting (Quantile Binning with Error Bars) ======
n_bins = 8; % Using 10 quantiles

% --- Distance Binning ---
edges_dist = quantile(all_reach_distances, linspace(0, 1, n_bins+1));
edges_dist(1) = min(all_reach_distances) - eps; 
edges_dist(end) = max(all_reach_distances) + eps;

binID_dist = discretize(all_reach_distances, edges_dist);
valid_dist = ~isnan(binID_dist);

% Mean distance per bin
mean_dist_centers = accumarray(binID_dist(valid_dist), all_reach_distances(valid_dist), [n_bins 1], @mean, NaN);
% SD of residuals per bin
std_dist = accumarray(binID_dist(valid_dist), all_residuals(valid_dist), [n_bins 1], @std, NaN);
% Count of items per bin
n_dist_counts = accumarray(binID_dist(valid_dist), 1, [n_bins 1]);
% Standard Error of the Standard Deviation
sem_sd_dist = std_dist ./ sqrt(2 .* n_dist_counts);

% --- Speed Binning ---
edges_speed = quantile(all_avg_speed, linspace(0, 1, n_bins+1));
edges_speed(1) = min(all_avg_speed) - eps;
edges_speed(end) = max(all_avg_speed) + eps;

binID_speed = discretize(all_avg_speed, edges_speed);
valid_speed = ~isnan(binID_speed);

% Mean speed per bin
mean_speed_centers = accumarray(binID_speed(valid_speed), all_avg_speed(valid_speed), [n_bins 1], @mean, NaN);
% SD of residuals per bin
std_speed = accumarray(binID_speed(valid_speed), all_residuals(valid_speed), [n_bins 1], @std, NaN);
% Count of items per bin
n_speed_counts = accumarray(binID_speed(valid_speed), 1, [n_bins 1]);
% Standard Error of the Standard Deviation
sem_sd_speed = std_speed ./ sqrt(2 .* n_speed_counts);

%% Polygon
% Calculate polygon coordinates for Distance
x_dist_fill = [mean_dist_centers(:)', fliplr(mean_dist_centers(:)')];
y_dist_fill = [(std_dist(:) + sem_sd_dist(:))', fliplr((std_dist(:) - sem_sd_dist(:))')];

% Calculate polygon coordinates for Speed
x_speed_fill = [mean_speed_centers(:)', fliplr(mean_speed_centers(:)')];
y_speed_fill = [(std_speed(:) + sem_sd_speed(:))', fliplr((std_speed(:) - sem_sd_speed(:))')];
%%
% --- Plotting ---
figure('Name', 'Aggregated Dir Errors Residual S.D. with Error Bars');
plotColor = [143, 195, 31] / 255; % Default MATLAB blue

% Plot 1: Standard Deviation of Residuals vs. Reach Distance
subplot(1, 2, 1);
plot(mean_dist_centers, std_dist, '-', 'LineWidth', 2, 'MarkerFaceColor', 'none');
hold on
fill(x_dist_fill, y_dist_fill, plotColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
scatter(mean_dist_centers, std_dist, 25, plotColor, 'filled', "MarkerEdgeColor","none");
hold off

xlim([75 325])
xticks(100:100:300)
ylim([4.5 10.5])
yticks(4:2:10)

xlabel('Reach distance (mm)'); 
ylabel('S.D. of Residuals');
% title('Residual S.D. by Distance Quantiles');
grid off;

% Plot 2: Standard Deviation of Residuals vs. Average Speed
subplot(1, 2, 2);
plot(mean_speed_centers, std_speed, '-', 'LineWidth', 2, 'MarkerFaceColor', 'none');
hold on
fill(x_speed_fill, y_speed_fill, plotColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
scatter(mean_speed_centers, std_speed, 25, plotColor, 'filled', "MarkerEdgeColor","none");
hold off

xlim([225 675])
xticks(300:100:600)
ylim([4.5 10.5])
yticks(4:2:10)

xlabel('Average speed'); 
ylabel('S.D. of Residuals');
% title('Residual S.D. by Speed Quantiles');
grid off;