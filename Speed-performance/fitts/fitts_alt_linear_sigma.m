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

%participants = {"JHL", "JH", "LC", "LN", "RC", "SM", "ML", "SX", "SY"};

%% Load Data 
clear all
participant = 'JH';
fname_preamble = sprintf('../data_onlineConf/usable/%s_sptatialTemporalCostFunc*.mat',participant);
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
sig_gain = stdF(phi(1:3), dist, dur);      % σ_g in linear scale
sig_dir = stdF(phi(4:6), dist, dur);      % σ_d in linear scale
rho = phi(7);                        % static correlation ρ

if isscalar(sig_gain) && isscalar(sig_dir)
    Sigma = [sig_gain.^2, rho*sig_gain*sig_dir; rho*sig_gain*sig_dir, sig_dir.^2];
else
    n = numel(sig_gain);
    Sigma = arrayfun(@(i) [sig_gain(i).^2, rho*sig_gain(i)*sig_dir(i); ...
        rho*sig_gain(i)*sig_dir(i), sig_dir(i).^2], ...
        (1:n)', 'UniformOutput', false);
end
end

fun_std = @(theta,dist,dur) max(theta(1) .* dist + theta(2) .* dur + theta(3),eps);

stdF = fun_std;  % alias: std in linear scale (Fitts-style)

% Parameters: [dist_gain, dur_gain, intercept_gain, dist_dir, dur_dir, intercept_dir, rho]
phiUB = [ 2,  20,  10,  2,  20,  10,  0.95];
phiLB = [-2, -20, -10, -2, -20, -10, -0.95];


% Negative log-likelihood for zero-mean 2D Gaussian with sample-wise Sigma
% fun_nll2d = @(phi) -sum(log(mvnpdf([gain_error,dir_error],zeros(1,2),local_cov_pred(stdF, phi, distances, durations))));
fun_nll2d = @(phi) ...
    (0.5*sum( ...
    log( (stdF(phi(1:3),distances,durations).^2) .* ...
    (stdF(phi(4:6),distances,durations).^2) .* (1 - phi(7).^2) ) + ...
    (gain_error.^2) ./ (stdF(phi(1:3),distances,durations).^2) ...
    - 2*phi(7) .* (gain_error.*dir_error) ./ ...
    (stdF(phi(1:3),distances,durations).*stdF(phi(4:6),distances,durations)) + ...
    (dir_error.^2) ./ (stdF(phi(4:6),distances,durations).^2) + ...
    2*log(2*pi) ));

% Iterative BADS
num_runs = 10;
phi_hat_array = NaN(num_runs,size(phiUB,2));
nll_from_runs = NaN(num_runs,1);
for i = 1:num_runs
    phi0 = [rand*mean([phiUB(1),phiLB(1)]),  ... 
        rand*mean([phiUB(2),phiLB(2)]),  ...
        rand*mean([phiUB(3),phiLB(3)]),  ... 
        rand*mean([phiUB(4),phiLB(4)]),  ...
        rand*mean([phiUB(5),phiLB(5)]),  ...
        rand*mean([phiUB(6),phiLB(6)]),  ...
        0];
    phi0 = phi0 + eps;

    phi_hat_array(i,:) = bads(fun_nll2d, phi0, phiLB, phiUB);
    nll_from_runs(i) = fun_nll2d(phi_hat_array(i,:));
end

[~,best_ind] = min(nll_from_runs);
phi_hat = phi_hat_array(best_ind,:)

% Covariance predictor:
%  - If dist,dur are scalars: returns a 2x2 matrix.
%  - If dist,dur are vectors (same length): returns a 1x1 cell array per i, each is 2x2.

cov_pred = @(phi, dist, dur) ...
    local_cov_pred(stdF, phi, dist, dur);


%% Calculate Prediction Mapping from Fitted Parameters
step_n = 100;
rho = phi_hat(7);     % static correlation

tar_size = unique(copy(:,15));
smooth_dist = linspace(min(distances),max(distances),step_n);
smooth_dur = linspace(min(durations),max(durations),step_n);
zoom_dist = (max(distances) - min(distances)) ./ step_n;
zoom_dur = (max(durations) - min(durations))./ step_n;

[x,y] = meshgrid(smooth_dist,smooth_dur);
smooth_gain_std = fun_std(phi_hat(1:3),x,y);
smooth_dir_std = fun_std(phi_hat(4:6),x,y);

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


x = linspace(min(distances),max(distances),100);

k = 14;
tar_size = 20;
y = log2(2 .* x ./ tar_size) ./ k;
tar_size_array = repmat(tar_size(i),size(x));
rho_array = repmat(rho,size(x));
p_of_hit = local_p_hit2d_rect(tar_size_array,fun_std(phi_hat(1:3),x,y),fun_std(phi_hat(4:6),x,y),phi_hat(7));

%%
plot(p_of_hit,'--')
%%
figure();

% --- Subplot 1: Varying k ---
subplot(1,2,1)
x = linspace(min(distances), max(distances), 100);
k_vals = 8:2:16;
tar_size_fixed = 20;
hold on

for i = 1:length(k_vals)
    y = log2(2 .* x ./ tar_size_fixed) ./ k_vals(i);
    
    % Calculate probability array
    p_of_hit = local_p_hit2d_rect(repmat(tar_size_fixed, size(x)), ...
               fun_std(phi_hat(1:3), x, y), ...
               fun_std(phi_hat(4:6), x, y), phi_hat(7));

    % Use patch for gradient lines. We add a NaN to the end of x and y 
    % because patch colors the segments between vertices.
    patch([x NaN], [y NaN], [p_of_hit NaN], 'EdgeColor', 'interp', ...
          'LineWidth', 3, 'LineStyle', '--', 'DisplayName', sprintf('k = %d', k_vals(i)));
end

title('Variation in k');
xlabel('Distance (mm)'); ylabel('Duration (s)');
grid on;
clim([0 1]); % Sets color scale from 0 to 1
legend('Location', 'northeast'); % Automatically uses 'DisplayName'

% --- Subplot 2: Varying target size ---
subplot(1,2,2)
k_fixed = 10;
tar_sizes = 6:4:22; % Updated to match your original legend intent
hold on

for i = 1:length(tar_sizes)
    y = log2(2 .* x ./ tar_sizes(i)) ./ k_fixed;
    
    p_of_hit = local_p_hit2d_rect(repmat(tar_sizes(i), size(x)), ...
               fun_std(phi_hat(1:3), x, y), ...
               fun_std(phi_hat(4:6), x, y), phi_hat(7));

    patch([x NaN], [y NaN], [p_of_hit NaN], 'EdgeColor', 'interp', ...
          'LineWidth', 3, 'LineStyle', '--', 'DisplayName', sprintf('Size = %g', tar_sizes(i)));
end

title('Variation in Target Size');
xlabel('Distance (mm)'); ylabel('Duration (s)');
grid on;
clim([0 1]); 
legend('Location', 'northeast'); % Automatically uses 'DisplayName'

% --- Global Colorbar ---
% This creates one colorbar that applies to the current axes (0 to 1)
cb = colorbar; 
cb.Label.String = 'Probability of Hit (p\_of\_hit)';
% Move colorbar to the right of both plots
cb.Layout.Tile = 'east'; % If using tiledlayout, otherwise default is fine

