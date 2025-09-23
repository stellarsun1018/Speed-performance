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


%%
% participant = 1;
% if participant == 1
%     load('JM_practice_traj_S1_06-Apr-2023_tform.mat')
%     load('JM_practice_traj_S1c_06-Apr-2023_rawtotal.mat')
%     load('JM_practice_traj_S1c_06-Apr-2023_traXtotal.mat')
%     load('JM_practice_traj_S1c_06-Apr-2023_traYtotal.mat')
% elseif participant == 0
%     load('pilot_practice_traj_S1_31-Mar-2023_tform.mat')
%     load('rawComplete.mat')
%     load('xTrajComplete.mat')
%     load('yTrajComplete.mat')
% end
%%
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

%%
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

%%
step_n = 100;
smooth_dist = linspace(min(distances),max(distances),step_n);
smooth_dur = linspace(min(durations),max(durations),step_n);
zoom_dist = (max(distances) - min(distances)) / step_n;
zoom_dur = (max(durations) - min(durations)) / step_n;

[x,y] = meshgrid(smooth_dist,smooth_dur);
smooth_gain_std = fun_std(theta_gain,x,y);
smooth_dir_std = fun_std(theta_dir,x,y);

subplot(1,2,1)
imagesc(smooth_gain_std)
hold on
contour(1:step_n, 1:step_n, smooth_gain_std, 'LineColor', 'k', 'LineWidth', 1);
plot(distances./zoom_dist, durations./zoom_dur,'ko')
hold off
xlabel('Distance')
ylabel('Duration')
title('Gain Sigma')
xticks(linspace(1, step_n, 5));
xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 2));
yticks(linspace(1, step_n, 5));
yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
colorbar
axis xy

subplot(1,2,2)
imagesc(smooth_dir_std)
hold on
contour(1:step_n, 1:step_n, smooth_dir_std, 'LineColor', 'k', 'LineWidth', 1);
plot(distances./zoom_dist, durations./zoom_dur,'ko')
hold off
xlabel('Distance')
ylabel('Duration')
title('Direction Sigma')
xticks(linspace(1, step_n, 5));
xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 2));
yticks(linspace(1, step_n, 5));
yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
colorbar
axis xy

%%
smooth_gain_width = smooth_gain_std ./ theta_gain(2);
smooth_dir_width = smooth_dir_std ./ theta_dir(2);

subplot(1,2,1)
imagesc(smooth_gain_width)
hold on
contour(1:step_n, 1:step_n, smooth_gain_width, 'LineColor', 'k', 'LineWidth', 1);
hold off
xlabel('Distance')
ylabel('Duration')
title('Gain Width')
xticks(linspace(1, step_n, 5));
xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 2));
yticks(linspace(1, step_n, 5));
yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
colorbar
axis xy

subplot(1,2,2)
imagesc(smooth_dir_width)
hold on
contour(1:step_n, 1:step_n, smooth_dir_width, 'LineColor', 'k', 'LineWidth', 1);
hold off
xlabel('Distance')
ylabel('Duration')
title('Direction Width')
xticks(linspace(1, step_n, 5));
xticklabels(round(linspace(min(smooth_dist), max(smooth_dist), 5), 2));
yticks(linspace(1, step_n, 5));
yticklabels(round(linspace(min(smooth_dur), max(smooth_dur), 5), 2));
colorbar
axis xy
%%
gain_mean = zeros(size(distances));
gain_std = fun_fit_std(theta_gain);
dir_mean = zeros(size(distances));
dir_std = fun_fit_std(theta_dir);

figure;
subplot(1,2,1)
plot3(distances,durations,gain_std,'o')
ylabel('Duration')
xlabel('Distances')
zlabel('Gain Std')
title('Gain')
grid on

subplot(1,2,2)
plot3(distances,durations,dir_std,'o')
ylabel('Duration')
xlabel('Distances')
zlabel('Dir Std')
title('Direction')
grid on

%%
example_sizes = linspace(min(end_size),max(end_size),5);
accuracy = linspace(0.5,0.8,5);

figure;
hold on
for i = 1:length(example_sizes)
smooth_mt = log2(2 .* smooth_dist ./ example_sizes(i)) ./ 10;
plot(smooth_dist,smooth_mt,'k--', LineWidth=2)

std_i = example_sizes(1) / (sqrt(2) * erfinv(accuracy(i)));
width_i = ((std_i - theta_dir(3)) ./ theta_dir(2));
dur_i = log2(2 .* smooth_dist ./ width_i) ./ theta_dir(1);
plot(smooth_dist,dur_i,'r--', LineWidth=2)
end 
% plot(distances, durations, 'o')
hold off
grid on
legend('Vary Sizes','Vary Accuracies')

%%
figure;
hold on
for i = 1:length(example_sizes)
smooth_mt = log2(2 .* smooth_dist ./ example_sizes(i)) ./ theta_gain(1);
std_i = example_sizes(1) / (sqrt(2) * erfinv(accuracy(i)));
dur_i = log2(2 .* smooth_dist ./ ((std_i - theta_gain(3)) ./ theta_gain(2)) ) ./ theta_gain(1);
plot(smooth_dist,smooth_mt,'k--', LineWidth=2)
plot(smooth_dist,dur_i,'r--', LineWidth=2)
end 
% plot(distances, durations, 'o')
hold off
grid on
legend('Vary Sizes','Vary Accuracies')
%%
accuracy = 0.7;
half_width = gaussianBracketCenteredOnZero(gain_mean, smooth_gain_std, accuracy);

figure;
plot3(distances,durations,half_width,'o')
hold on
plot3(distances,durations,end_size,'o')
hold off
xlabel('Distance')
ylabel('Duration')
zlabel('Radius')
title(['Aiming for accuracy = ' num2str(accuracy *100) '%'])
legend('Actual radius','Predicted Radius')
grid on


function half_width = gaussianBracketCenteredOnZero(mu, sigma, p)
    [mu, sigma] = deal(mu(:), sigma(:));
    n = numel(mu);
    
    half_width = zeros(n, 1);
    
    for i = 1:n
        targetFun = @(L) normcdf(L, mu(i), sigma(i)) - ...
                         normcdf(-L, mu(i), sigma(i)) - p;
        % initial guess: symmetric case ignoring mean
        L0 = sigma(i) * norminv((1 + p) / 2);
        half_width(i) = fzero(targetFun, L0);
    end
end