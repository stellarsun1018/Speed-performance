% copy column contents:
% 1,2: target x and y in wac pixels
% 3: life span by block / target shrinking longest duration(second)
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
% 31: target shrinking speed(mm/s)


%% Load Data 
clear all
participant = 'RC';
session = 2;
fname_preamble = sprintf('data_onlineConf/%s/%s_sptatialTemporalCostFunc_S%d*.mat',participant,participant,session);
files = dir(fname_preamble);
for k = 1:numel(files)
    f = fullfile(files(k).folder, files(k).name);
    load(f);
end

%% Proprocess Data
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

lifespan = [0.6,0.6*3^(0.25),0.6*3^(0.5)];
for i = 1:3
    copy((1+(i-1)*240):(i*240),3) = lifespan(i);
end

copy(:,31) = copy(:,15) ./ copy(:,3); % 31: target shrinking speed(mm/s)

%% Fit Fitts's Law Parameter
distances = copy(:,10);
gain_error = copy(:,23); % biased/signed
dir_error = copy(:,29); % biased/signed
durations = copy(:,16);
speeds = copy(:,22);
end_size = copy(:,15) .* durations ./ copy(:,3);


thetaUB = [10,2]; % [theta(1) = index of performance IP, theta(2) tuning factor, theta(3) baseline
thetaLB = [-10,eps];
theta0 = rand .* (thetaUB + thetaLB)./2;

fun_std = @(theta,dist,dur) 2.^(log2(dist * 2) - (theta(1) .* dur)) .* theta(2) + 2;
fun_fit_std = @(theta) fun_std(theta,distances,durations);
fun_base = @(theta,error) -sum(log(normpdf(error,0,fun_fit_std(theta))));

fun_gain = @(theta) fun_base(theta,gain_error);
theta_gain = bads(fun_gain,theta0,thetaLB,thetaUB);

fun_dir = @(theta) fun_base(theta,dir_error);
theta_dir = bads(fun_dir,theta0,thetaLB,thetaUB);


%% Calculate Prediction Mapping from Fitted Parameters
step_n = 100;
tar_size = unique(copy(:,15));
smooth_dist = linspace(min(distances),max(distances),step_n);
smooth_dur = linspace(min(durations),max(durations),step_n);
zoom_dist = (max(distances) - min(distances)) ./ step_n;
zoom_dur = (max(durations) - min(durations))./ step_n;

[x,y] = meshgrid(smooth_dist,smooth_dur);
smooth_gain_std = fun_std(theta_gain,x,y);
smooth_dir_std = fun_std(theta_dir,x,y);

target_live_spans = sort(unique(copy(:,3)));
life_span_keys = zeros(length(target_live_spans),size(copy,1));
smooth_tar_size = NaN(length(target_live_spans),length(smooth_dist),length(smooth_dur));

get_tar_size = @(dur,life_span,tar_size,dist) tar_size .* max((1 - dur ./ life_span),0);

for i = 1:length(target_live_spans)
    life_span_keys(i,:) = copy(:,3) == target_live_spans(i);
    tar_size_array = get_tar_size(y,target_live_spans(i),tar_size,x);
    % tar_size_array = tar_size_array';
    smooth_tar_size(i, :, :) = tar_size_array; %repmat(tar_size_array,length(smooth_dist),1);
end

score_per_hit = @(term_size,tar_size,max_score) max_score .* max((term_size ./ tar_size),0);
p_hit = @(term_size,sigma) normcdf(term_size, 0, sigma) - normcdf(-term_size, 0, sigma);
exp_gain = @(term_size,sigma) p_hit(term_size,sigma) .* score_per_hit(term_size,tar_size,10);

%% Plot Specs
column_num = 4;
plot_colorbar = 1;
%% Gain Axis
figure
for i = 1:length(target_live_spans)
    plot_idx = i-1;
    subplot(length(target_live_spans),column_num,plot_idx * column_num + 1)
    current_size_map = squeeze(smooth_tar_size(i,:,:));
    imagesc(current_size_map)
    hold on
    contour(1:step_n, 1:step_n, current_size_map, 'LineColor', 'k', 'LineWidth', 1);
    plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
        (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
    hold off

    if i == length(target_live_spans)
        xlabel('Distance')
    elseif i == 1
        title('Target Size')
    end

    ylabel('Duration')
    xticks(linspace(1, step_n, 5));
    xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 0));
    yticks(linspace(1, step_n, 5));
    yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
    if plot_colorbar; colorbar; end
    axis xy

    subplot(length(target_live_spans),column_num,plot_idx * column_num + 2)
    current_sigma_map = smooth_gain_std;
    imagesc(current_sigma_map)
    hold on
    contour(1:step_n, 1:step_n, current_sigma_map, 'LineColor', 'k', 'LineWidth', 1);
    plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
        (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
    hold off

    if i == length(target_live_spans)
        xlabel('Distance')
    elseif i == 1
        title('Gain Sigma')
    end

    xticks(linspace(1, step_n, 5));
    xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 0));
    yticks(linspace(1, step_n, 5));
    yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
    if plot_colorbar; colorbar; end
    axis xy

    subplot(length(target_live_spans),column_num,plot_idx * column_num + 3)
    current_phit_map = p_hit(current_size_map,current_sigma_map);
    imagesc(current_phit_map)
    hold on
    contour(1:step_n, 1:step_n, current_phit_map, 'LineColor', 'k', 'LineWidth', 1);
    plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
        (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
    hold off

    if i == length(target_live_spans)
        xlabel('Distance')
    elseif i == 1
        title('P(hit)')
    end

    xticks(linspace(1, step_n, 5));
    xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 0));
    yticks(linspace(1, step_n, 5));
    yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
    if plot_colorbar; colorbar; end
    axis xy

    subplot(length(target_live_spans),column_num,plot_idx * column_num + 4)
    current_egain_map = exp_gain(current_size_map,current_sigma_map);
    imagesc(current_egain_map)
    hold on
    contour(1:step_n, 1:step_n, current_egain_map, 'LineColor', 'k', 'LineWidth', 1);
    plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
        (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
    hold off

    if i == length(target_live_spans)
        xlabel('Distance')
    elseif i == 1
        title('Expected Gain')
    end

    xticks(linspace(1, step_n, 5));
    xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 0));
    yticks(linspace(1, step_n, 5));
    yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
    if plot_colorbar; colorbar; end
    axis xy

end
sgtitle(['Participant ' participant ' on Gain Axis'])
%% Directional Axis
figure
for i = 1:length(target_live_spans)
    plot_idx = i-1;
    subplot(length(target_live_spans),column_num,plot_idx * column_num + 1)
    current_size_map = squeeze(smooth_tar_size(i,:,:));
    imagesc(current_size_map)
    hold on
    contour(1:step_n, 1:step_n, current_size_map, 'LineColor', 'k', 'LineWidth', 1);
    plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
        (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
    hold off

    if i == length(target_live_spans)
        xlabel('Distance')
    elseif i == 1
        title('Target Size')
    end

    ylabel('Duration')
    xticks(linspace(1, step_n, 5));
    xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 0));
    yticks(linspace(1, step_n, 5));
    yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
    if plot_colorbar; colorbar; end
    axis xy

    subplot(length(target_live_spans),column_num,plot_idx * column_num + 2)
    current_sigma_map = smooth_dir_std;
    imagesc(current_sigma_map)
    hold on
    contour(1:step_n, 1:step_n, current_sigma_map, 'LineColor', 'k', 'LineWidth', 1);
    plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
        (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
    hold off

    if i == length(target_live_spans)
        xlabel('Distance')
    elseif i == 1
        title('Directional Sigma')
    end

    xticks(linspace(1, step_n, 5));
    xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 0));
    yticks(linspace(1, step_n, 5));
    yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
    if plot_colorbar; colorbar; end
    axis xy

    subplot(length(target_live_spans),column_num,plot_idx * column_num + 3)
    current_phit_map = p_hit(current_size_map,current_sigma_map);
    imagesc(current_phit_map)
    hold on
    contour(1:step_n, 1:step_n, current_phit_map, 'LineColor', 'k', 'LineWidth', 1);
    plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
        (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
    hold off

    if i == length(target_live_spans)
        xlabel('Distance')
    elseif i == 1
        title('P(hit)')
    end

    xticks(linspace(1, step_n, 5));
    xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 0));
    yticks(linspace(1, step_n, 5));
    yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
    if plot_colorbar; colorbar; end
    axis xy

    subplot(length(target_live_spans),column_num,plot_idx * column_num + 4)
    current_egain_map = exp_gain(current_size_map,current_sigma_map);
    imagesc(current_egain_map)
    hold on
    contour(1:step_n, 1:step_n, current_egain_map, 'LineColor', 'k', 'LineWidth', 1);
    plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
        (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
    hold off

    if i == length(target_live_spans)
        xlabel('Distance')
    elseif i == 1
        title('Expected Gain')
    end
    xticks(linspace(1, step_n, 5));
    xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 0));
    yticks(linspace(1, step_n, 5));
    yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
    if plot_colorbar; colorbar; end
    axis xy
end
sgtitle(['Participant ' participant ' on Directional Axis'])
