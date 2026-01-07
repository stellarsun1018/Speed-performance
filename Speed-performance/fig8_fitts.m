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
close all
participant = 'LL';
fname_preamble = sprintf('data_onlineConf/usable/%s_sptatialTemporalCostFunc*.mat',participant);
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
distances = copy(:,10); % using actual or 
gain_error = copy(:,23); % biased/signed
dir_error = copy(:,29); % biased/signed
durations = copy(:,16);
speeds = copy(:,22);
end_size = copy(:,15) .* durations ./ copy(:,3);

%% Recenter by offseting with a lienar bias multivariate prediction

X = [ones(size(distances)), distances, durations];

gain_coeffs = regress(gain_error, X);
gain_bias_offset = gain_coeffs(1) + gain_coeffs(2)*distances + gain_coeffs(3)*durations;
gain_error = gain_error - gain_bias_offset;

dir_coeffs = regress(dir_error, X);
dir_bias_offset = dir_coeffs(1) + dir_coeffs(2)*distances + dir_coeffs(3)*durations;
dir_error = dir_error - dir_bias_offset;


%%
function Sigma = local_cov_pred(stdF, phi, dist, dur)
sig_gain = stdF(phi(1:2), dist, dur);      % σ_g in linear scale
sig_dir = stdF(phi(3:4), dist, dur);      % σ_d in linear scale
rho = phi(5);                        % static correlation ρ

if isscalar(sig_gain) && isscalar(sig_dir)
    Sigma = [sig_gain.^2, rho*sig_gain*sig_dir; rho*sig_gain*sig_dir, sig_dir.^2];
else
    n = numel(sig_gain);
    Sigma = arrayfun(@(i) [sig_gain(i).^2, rho*sig_gain(i)*sig_dir(i); ...
        rho*sig_gain(i)*sig_dir(i), sig_dir(i).^2], ...
        (1:n)', 'UniformOutput', false);
end
end

fun_std = @(theta,dist,dur) 2.^(log2(dist * 2) - (theta(1) .* dur)) .* theta(2) + 2;

stdF = fun_std;  % alias: std in linear scale (Fitts-style)

% Parameters: [IP_g, k_g, IP_d, k_d, rho]
phiUB = [ 10,  2,  10,  2,  0.95];
phiLB = [-10, eps,-10, eps, -0.95];


% Negative log-likelihood for zero-mean 2D Gaussian with sample-wise Sigma
fun_nll2d = @(phi) ...
    (0.5*sum( ...
    log( (stdF(phi(1:2),distances,durations).^2) .* ...
    (stdF(phi(3:4),distances,durations).^2) .* (1 - phi(5).^2) ) + ...
    (gain_error.^2) ./ (stdF(phi(1:2),distances,durations).^2) ...
    - 2*phi(5) .* (gain_error.*dir_error) ./ ...
    (stdF(phi(1:2),distances,durations).*stdF(phi(3:4),distances,durations)) + ...
    (dir_error.^2) ./ (stdF(phi(3:4),distances,durations).^2) + ...
    2*log(2*pi) ));


% Iterative BADS
num_runs = 10;
phi_hat_array = NaN(num_runs,size(phiUB,2));
nll_from_runs = NaN(num_runs,1);
for i = 1:num_runs
    phi0 = [rand*mean([phiUB(1),phiLB(1)]),  ... % IP_g
        rand*mean([phiUB(2),phiLB(2)]),  ... % k_g
        rand*mean([phiUB(3),phiLB(3)]),  ... % IP_d
        rand*mean([phiUB(4),phiLB(4)]),  ... % k_d
        0];

    phi_hat_array(i,:) = bads(fun_nll2d, phi0, phiLB, phiUB);
    nll_from_runs(i) = fun_nll2d(phi_hat_array(i,:));
end

[~,best_ind] = min(nll_from_runs);
phi_hat = phi_hat_array(best_ind,:);

% Covariance predictor:
%  - If dist,dur are scalars: returns a 2x2 matrix.
%  - If dist,dur are vectors (same length): returns a 1x1 cell array per i, each is 2x2.
cov_pred = @(phi, dist, dur) ...
    local_cov_pred(stdF, phi, dist, dur);

%% Calculate Prediction Mapping from Fitted Parameters
step_n = 100;
rho = phi_hat(5);     % static correlation

tar_size = unique(copy(:,15));
smooth_dist = linspace(min(distances),max(distances),step_n);
smooth_dur = linspace(min(durations),max(durations),step_n);
zoom_dist = (max(distances) - min(distances)) ./ step_n;
zoom_dur = (max(durations) - min(durations))./ step_n;

[x,y] = meshgrid(smooth_dist,smooth_dur);
smooth_gain_std = fun_std(phi_hat(1:2),x,y);
smooth_dir_std = fun_std(phi_hat(3:4),x,y);

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
%%
function P = local_p_hit2d_rect(a, sg, sd, rho)
% Compute P(|G|<=a, |D|<=a) where [G;D] ~ N(0, Sigma), Sigma = [sg^2, rho*sg*sd; rho*sg*sd, sd^2]
% Vectorized over matching arrays a, sg, sd. Tries mvncdf; if unavailable, falls back to MC.

% basic input checks / broadcasting guard
if ~isequal(size(a), size(sg)) || ~isequal(size(a), size(sd))
    error('local_p_hit2d_rect: size mismatch among a, sg, sd.');
end

% Flatten to 1D for a simple loop, then reshape back
a_v  = a(:);
sg_v = sg(:);
sd_v = sd(:);
n    = numel(a_v);
P_v  = nan(n,1);

% Try exact bivariate normal CDF over rectangle with mvncdf
can_mvncdf = exist('mvncdf','file') == 2;

if can_mvncdf
    for k = 1:n
        ak  = a_v(k);
        sgk = sg_v(k);
        sdk = sd_v(k);
        % Guard tiny/zero stds
        sgk = max(sgk, realmin('double'));
        sdk = max(sdk, realmin('double'));

        % Lower/Upper bounds and Sigma
        lo = [-ak, -ak];
        hi = [ ak,  ak];
        Sigma = [sgk^2, rho*sgk*sdk; rho*sgk*sdk, sdk^2];

        % mvncdf returns the prob mass in the rectangle
        P_v(k) = mvncdf(lo, hi, [0 0], Sigma);
    end
else
    % Monte Carlo fallback (fast-ish): independent base normals + correlate via Cholesky
    % Note: MC per grid point can be costly; keep N moderate.
    N = 1000;  % adjust if you want smoother estimates vs runtime
    Z = randn(N, 2);  % reused base samples
    for k = 1:n
        ak  = a_v(k);
        sgk = sg_v(k);
        sdk = sd_v(k);
        sgk = max(sgk, realmin('double'));
        sdk = max(sdk, realmin('double'));
        Sigma = [sgk^2, rho*sgk*sdk; rho*sgk*sdk, sdk^2];

        % Cholesky (robustify with jitter if needed)
        [L,p] = chol(Sigma, 'lower');
        if p>0
            % add tiny jitter if Sigma is borderline
            jitter = 1e-12 * max(Sigma(1,1)+Sigma(2,2), 1);
            [L,~] = chol(Sigma + jitter*eye(2), 'lower');
        end
        X = Z * L.';  % N x 2 samples ~ N(0,Sigma)
        P_v(k) = mean( abs(X(:,1)) <= ak & abs(X(:,2)) <= ak );
    end
end

P = reshape(P_v, size(a));
end

p_hit2d = @(a, sg, sd) local_p_hit2d_rect(a, sg, sd, rho);

%%
function EG = local_exp_score2d_mc(term_size, sg, sd, rho, tar_size, max_score, N)
% Monte Carlo expected score per grid cell for the radial scoring rule.
% Shapes: term_size, sg, sd are identical arrays (e.g., step_n x step_n).
% tar_size can be scalar or same-sized array. EG has same shape as inputs.

if ~isequal(size(term_size), size(sg)) || ~isequal(size(term_size), size(sd))
    error('local_exp_score2d_mc: size mismatch among term_size, sg, sd.');
end
if ~isscalar(tar_size) && ~isequal(size(tar_size), size(term_size))
    error('local_exp_score2d_mc: tar_size must be scalar or match term_size size.');
end

% Vectorize over grid by flattening, then reshape back.
a_v  = term_size(:);
sg_v = max(sg(:),  realmin('double'));
sd_v = max(sd(:),  realmin('double'));
if isscalar(tar_size)
    ts_v = repmat(tar_size, numel(a_v), 1);
else
    ts_v = tar_size(:);
end
n  = numel(a_v);
EG_v = zeros(n,1);

% Reuse base standard-normal samples for variance reduction & speed.
if nargin < 7 || isempty(N), N = 2000; end
Z = randn(N, 2); % N x 2, iid N(0, I)

for k = 1:n
    ak  = a_v(k);         % term_size (acceptance radius)
    sgk = sg_v(k);
    sdk = sd_v(k);
    tsk = ts_v(k);

    % Build Sigma and its Cholesky
    Sigma = [sgk^2, rho*sgk*sdk; rho*sgk*sdk, sdk^2];
    [L,p] = chol(Sigma, 'lower');
    if p>0
        jitter = 1e-12 * max(Sigma(1,1)+Sigma(2,2), 1);
        [L,~] = chol(Sigma + jitter*eye(2), 'lower');
    end

    % Draw correlated samples ~ N(0, Sigma)
    X = Z * L.';                 % N x 2
    r = hypot(X(:,1), X(:,2));   % radial distance

    % Radial, linearly decaying score inside acceptance radius
    S = max_score * (1 - r./tsk);
    S(r > ak) = 0;               % zero outside term_size
    S(S < 0)  = 0;               % clamp if term_size < tar_size, etc.

    EG_v(k) = mean(S);
end

EG = reshape(EG_v, size(term_size));
end


% MC expected score maps for each lifespan
max_score = 10;
exp_gain_maps = cell(length(target_live_spans),1);

for i = 1:length(target_live_spans)
    a  = squeeze(smooth_tar_size(i,:,:));   % term_size grid (dist x dur)
    sg = smooth_gain_std;                   % σ_g grid
    sd = smooth_dir_std;                    % σ_d grid

    % tar_size can be a scalar
    exp_gain_maps{i} = local_exp_score2d_mc(a, sg, sd, rho, tar_size, max_score, 2000);
end

% exp_gain2d = @(a) local_exp_score2d_mc(a, smooth_gain_std, smooth_dir_std, rho, tar_size, max_score, 2000);


%% Plot heatmaps
% column_num = 5;
% plot_colorbar = 1;
% figure
% for i = 1:length(target_live_spans)
%     plot_idx = i-1;
%     subplot(length(target_live_spans),column_num,plot_idx * column_num + 1)
%     current_size_map = squeeze(smooth_tar_size(i,:,:));
%     imagesc(current_size_map)
%     hold on
%     contour(1:step_n, 1:step_n, current_size_map, 'LineColor', 'k', 'LineWidth', 1);
%     plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
%         (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
%     hold off
%
%     if i == length(target_live_spans)
%         xlabel('Distance')
%     elseif i == 1
%         title('Target Size')
%     end
%
%     ylabel('Duration')
%     xticks(linspace(1, step_n, 5));
%     xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 0));
%     yticks(linspace(1, step_n, 5));
%     yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
%     if plot_colorbar; colorbar; end
%     axis xy
%
%     subplot(length(target_live_spans),column_num,plot_idx * column_num + 2)
%     imagesc(smooth_gain_std)
%     hold on
%     contour(1:step_n, 1:step_n, smooth_gain_std, 'LineColor', 'k', 'LineWidth', 1);
%     plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
%         (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
%     hold off
%
%     if i == length(target_live_spans)
%         xlabel('Distance')
%     elseif i == 1
%         title('Gain Sigma')
%     end
%
%     xticks(linspace(1, step_n, 5));
%     xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 0));
%     yticks(linspace(1, step_n, 5));
%     yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
%     if plot_colorbar; colorbar; end
%     axis xy
%
%     subplot(length(target_live_spans),column_num,plot_idx * column_num + 3)
%     imagesc(smooth_dir_std)
%     hold on
%     contour(1:step_n, 1:step_n, smooth_dir_std, 'LineColor', 'k', 'LineWidth', 1);
%     plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
%         (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
%     hold off
%
%     if i == length(target_live_spans)
%         xlabel('Distance')
%     elseif i == 1
%         title('Directional Sigma')
%     end
%
%     xticks(linspace(1, step_n, 5));
%     xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 0));
%     yticks(linspace(1, step_n, 5));
%     yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
%     if plot_colorbar; colorbar; end
%     axis xy
%
%     subplot(length(target_live_spans),column_num,plot_idx * column_num + 4)
%     current_phit_map = p_hit2d(current_size_map,smooth_gain_std,smooth_dir_std);
%     imagesc(current_phit_map)
%     hold on
%     contour(1:step_n, 1:step_n, current_phit_map, 'LineColor', 'k', 'LineWidth', 1);
%     plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
%         (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
%     hold off
%
%     if i == length(target_live_spans)
%         xlabel('Distance')
%     elseif i == 1
%         title('P(hit)')
%     end
%
%     xticks(linspace(1, step_n, 5));
%     xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 0));
%     yticks(linspace(1, step_n, 5));
%     yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
%     if plot_colorbar; colorbar; end
%     axis xy
%
%     subplot(length(target_live_spans),column_num,plot_idx * column_num + 5)
%     current_egain_map = exp_gain_maps{i};
%     imagesc(current_egain_map)
%     hold on
%     contour(1:step_n, 1:step_n, current_egain_map, 'LineColor', 'k', 'LineWidth', 1);
%     plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
%         (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
%     hold off
%
%     if i == length(target_live_spans)
%         xlabel('Distance')
%     elseif i == 1
%         title('Expected Gain')
%     end
%
%     xticks(linspace(1, step_n, 5));
%     xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 0));
%     yticks(linspace(1, step_n, 5));
%     yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
%     if plot_colorbar; colorbar; end
%     axis xy
%
% end
% sgtitle(['Participant ' participant ' Ideal Observer'])
%
% saveas(gcf, fullfile('fitts', ['Participant_' participant '_Ideal_Observer.png']));

%%
column_num = 5;
plot_colorbar = 1;

rows = length(target_live_spans);

% Initialize an array to store axes handles
ax_handles = gobjects(rows, column_num);

figure
for i = 1:rows
    plot_idx = i-1;

    %%%%%%%%%% 111111 %%%%%%%%%
    % Assign output of subplot to the handle array (ax_handles)
    ax_handles(i, 1) = subplot(rows, column_num, plot_idx * column_num + 1);
    current_size_map = squeeze(smooth_tar_size(i,:,:)); % if we use final size
    % current_size_map = ones(step_n,step_n) .* copy(1,15); % if we use initial size
    imagesc(current_size_map)
    hold on
    contour(1:step_n, 1:step_n, current_size_map, 'LineColor', 'k', 'LineWidth', 1);
    plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
        (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
    hold off
    if i == rows
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

    %%%%%%%%% 22222 %%%%%%%%%%
    ax_handles(i, 2) = subplot(rows, column_num, plot_idx * column_num + 2);
    current_id_map = log2( 2 .* smooth_dist ./ current_size_map);
    imagesc(current_id_map)
    hold on
    contour(1:step_n, 1:step_n, current_id_map, 'LineColor', 'k', 'LineWidth', 1);
    plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
        (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
    hold off
    if i == rows
        xlabel('Distance')
    elseif i == 1
        title('Index of Difficulty')
    end
    xticks(linspace(1, step_n, 5));
    xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 0));
    yticks(linspace(1, step_n, 5));
    yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
    if plot_colorbar; colorbar; end
    axis xy

    %%%%%%%%% 33333 %%%%%%%%%%
    ax_handles(i, 3) = subplot(rows, column_num, plot_idx * column_num + 3);
    current_id_map = log2( 2 .* smooth_dist ./ current_size_map);
    current_ip_map = current_id_map ./ smooth_dur';
    imagesc(current_ip_map)
    hold on
    contour(1:step_n, 1:step_n, current_ip_map, 'LineColor', 'k', 'LineWidth', 1);
    plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
        (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
    hold off
    if i == rows
        xlabel('Distance')
    elseif i == 1
        title('Index of Performance')
    end
    xticks(linspace(1, step_n, 5));
    xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 0));
    yticks(linspace(1, step_n, 5));
    yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
    if plot_colorbar; colorbar; end
    axis xy

    %%%%%%%%% 44444 %%%%%%%%%%
    ax_handles(i, 4) = subplot(rows, column_num, plot_idx * column_num + 4);
    current_size_map = squeeze(smooth_tar_size(i,:,:));
    current_phit_map = p_hit2d(current_size_map,smooth_gain_std,smooth_dir_std);
    imagesc(current_phit_map)
    hold on
    contour(1:step_n, 1:step_n, current_phit_map, 'LineColor', 'k', 'LineWidth', 1);
    plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
        (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
    hold off
    if i == rows
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

    %%%%%%%%% 55555 %%%%%%%%%%
    ax_handles(i, 5) = subplot(rows, column_num, plot_idx * column_num + 5);
    current_egain_map = exp_gain_maps{i};
    imagesc(current_egain_map)
    hold on
    contour(1:step_n, 1:step_n, current_egain_map, 'LineColor', 'k', 'LineWidth', 1);
    plot((distances(life_span_keys(i,:) == 1)-min(distances))./zoom_dist, ...
        (durations(life_span_keys(i,:) == 1)-min(durations))./zoom_dur,'k.')
    hold off
    if i == rows
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

% Post-processing to synchronize limits per column
for col = 1:column_num
    % Extract all axes handles for the current column
    col_axes = ax_handles(:, col);

    % Get the current limits (automatically set by imagesc)
    % 'CLim' returns a cell array of [min max], convert to matrix
    current_lims = cell2mat(get(col_axes, 'CLim'));

    % Find the global min and max for this entire column
    global_min = min(current_lims(:, 1));
    global_max = max(current_lims(:, 2));

    % Apply the new range to all axes in this column
    set(col_axes, 'CLim', [global_min, global_max]);
end

sgtitle(['Participant ' participant ' Ideal Observer'])
saveas(gcf, fullfile('fitts_plots', ['Participant_' participant '_Ideal_Observer.png']));
close all
%%
%% New Plots: Per-Trial p(hit) and Optimal Duration Analysis
% Compute p(hit) and optimal duration for each trial using the fitted model

p_hits = NaN(size(distances));
optimal_durs = NaN(size(distances));

for blk = 1:length(target_live_spans)
    idx_trials = life_span_keys(blk, :);  % Logical mask for this lifespan block
    lifespan = target_live_spans(blk);
    
    for t = find(idx_trials)
        dist_t = distances(t);
        dur_t = durations(t);
        
        % Compute std devs for this trial
        sg_t = fun_std(phi_hat(1:2), dist_t, dur_t);
        sd_t = fun_std(phi_hat(3:4), dist_t, dur_t);
        
        % Compute terminal target size for this trial
        term_size_t = tar_size * max((1 - dur_t / lifespan), 0);
        
        % Compute p(hit) for this trial
        p_hits(t) = p_hit2d(term_size_t, sg_t, sd_t);
        
        % Find optimal duration: Use precomputed exp_gain map for this block
        % Column index corresponding to dist_t
        col_idx = round((dist_t - min(smooth_dist)) / zoom_dist) + 1;
        col_idx = max(1, min(step_n, col_idx));  % Clamp to valid range
        
        % EG profile over durations for this dist (rows = dur indices)
        eg_profile = exp_gain_maps{blk}(:, col_idx);
        
        % Find duration index that maximizes EG
        [~, max_row_idx] = max(eg_profile);
        
        % Corresponding optimal duration
        optimal_durs(t) = smooth_dur(max_row_idx);
    end
end

% Plot 1: p(hit) vs Distance (colored by lifespan block)
figure('Name', 'Per-Trial p(hit) and Optimal Duration');
colors = lines(length(target_live_spans));  % Color palette for blocks

subplot(1, 2, 1);
hold on;
for blk = 1:length(target_live_spans)
    idx_trials = life_span_keys(blk, :) == 1;
    scatter(distances(idx_trials), p_hits(idx_trials), 50, colors(blk, :), 'filled', 'MarkerFaceAlpha', 0.7);
end
xlabel('Distance (mm)');
ylabel('Predicted P(hit)');
% title('Predicted P(hit) vs Distance');
ylim([0,1])
legend('Fast','Medium','Slow','Location', 'northeast');
grid on;
hold off;

% Plot 2: Actual Duration vs Optimal Duration for Max Expected Gain (colored by lifespan block)
subplot(1, 2, 2);
hold on;
for blk = 1:length(target_live_spans)
    idx_trials = life_span_keys(blk, :) == 1;
    scatter(durations(idx_trials), optimal_durs(idx_trials), 50, colors(blk, :), 'filled', 'MarkerFaceAlpha', 0.7);
end
xlabel('Actual Duration (s)');
ylabel('Optimal Duration (s)');
xlim([0.1,0.9])
ylim([0.1,0.9])
% title('Actual vs Optimal Duration for Max Expected Gain');
% Add reference line (y = x)
refline(1, 0);
% refline('Color', [0.5 0.5 0.5], 'LineStyle', '--');
legend('Fast','Medium','Slow','Location', 'northeast');


grid on;
hold off;


sgtitle(sprintf('Participant %s: Per-Trial Predictions', participant));

% Optional: Save the figure
saveas(gcf, fullfile('fitts_plots', sprintf('Participant_%s_PerTrial_Predictions.png', participant)));
close all