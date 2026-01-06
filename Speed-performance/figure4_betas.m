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
part_n = numel(participants);
% 如果你想按 session 读（像你上面 SX 那样），就用这一段 pattern：
use_session = false;   % true: 用 S%d 方式；false: 用 usable 方式
session = 3;

nBlocks = 3;
nTrialsPerBlock = 240;
lim_scale = 1.2;

pixellength = 0.248;

colors = lines(part_n);

beta_gain_0 = NaN(3,part_n);
beta_gain_dist = NaN(3,part_n);
beta_gain_spd = NaN(3,part_n);

beta_dir_0 = NaN(3,part_n);
beta_dir_dist = NaN(3,part_n);
beta_dir_spd = NaN(3,part_n);


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

    % --- valid trials (等价于你写法的向量化版本) ---
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

    % polar coordination (not required for regression lines, but keep consistent)
    angle_error_in_radian = asin(copy(:,29) ./ copy(:,21));
    copy(:,33) = angle_error_in_radian;
    copy(:,34) = rad2deg(angle_error_in_radian);

    % ====== Fit per participant: speed vs distance ======
    reach_distances = copy(:,21);
    avg_speed = copy(:,22);
    gain_errors = copy(:,23);
    dir_errors = copy(:,29);

    
    % Prepare the design matrix (adding a column of ones for intercept)
    X = [ones(size(reach_distances)), reach_distances, avg_speed];

    % Perform multivariate linear regression
    [coeffs, bInt] = regress(gain_errors, X);

    beta_gain_0(:,ip) = [coeffs(1);bInt(1,:)'];
    beta_gain_dist(:,ip) = [coeffs(2);bInt(2,:)'];
    beta_gain_spd(:,ip) = [coeffs(3);bInt(3,:)'];
    
    [coeffs, bInt] = regress(dir_errors, X);

    beta_dir_0(:,ip) = [coeffs(1);bInt(1,:)'];
    beta_dir_dist(:,ip) = [coeffs(2);bInt(2,:)'];
    beta_dir_spd(:,ip) = [coeffs(3);bInt(3,:)'];

end

%% Gain V1
figure 

idx = 1:part_n;

x0    = linspace(-0.5,  0.6, part_n); 
xdist = linspace(-0.5,  0.6, part_n); 
xspd  = linspace(1,  2.1, part_n); 

subplot(1,2,1)
errorbar(x0,beta_gain_0(1,:),beta_gain_0(2,:),beta_gain_0(3,:),'o','MarkerFaceColor','auto','LineWidth',1);

yline(0,'--');
xlim([x0(1)-0.1,x0(end)+0.1])
ylim([-25,25])
xticks(x0)
xticklabels(participants)
legend({'Intercept'}, 'Location','northeast');
xlabel('Participant');
ylabel('Intercept');
box off


subplot(1,2,2)
errorbar(xdist, beta_gain_dist(1,:), beta_gain_dist(2,:), beta_gain_dist(3,:),'o','MarkerFaceColor','auto','LineWidth',1);
hold on
errorbar(xspd,  beta_gain_spd(1,:),  beta_gain_spd(2,:),  beta_gain_spd(3,:),'o','MarkerFaceColor','auto','LineWidth',1); 
hold off

yline(0,'--');
xlim([xdist(1)-0.2,xspd(end)+0.2])
ylim([-0.2,0.2])
xticks([mean(xdist),mean(xspd)])
xticklabels(["Distance","Speed"])
xlabel('Factor');
ylabel('Slope');
legend({'\beta_{distance}','\beta_{speed}'}, 'Location','northeast');
box off
saveas(gcf,fullfile('results', 'plots', 'fig4','beta_gain_v1.png'));
%% Gain V2
figure
xdist = linspace(-0.25,  43.75, part_n); 
xspd  = linspace(0.25,  44.25, part_n); 
errorbar(xdist, beta_gain_dist(1,:), beta_gain_dist(2,:), beta_gain_dist(3,:),'o','MarkerFaceColor','auto','LineWidth',1);
hold on
errorbar(xspd,  beta_gain_spd(1,:),  beta_gain_spd(2,:),  beta_gain_spd(3,:),'o','MarkerFaceColor','auto','LineWidth',1); 
hold off

yline(0,'--');
xlim([xspd(1)-2,xspd(end)+2])
ylim([-0.2,0.2])

xticks(mean([xdist;xspd]))
xticklabels(participants)
xlabel('Participant');
ylabel('Slope');
legend({'\beta_{distance}','\beta_{speed}'}, 'Location','northeast');
box off
hold off

saveas(gcf,fullfile('results', 'plots', 'fig4','beta_gain_v2.png'));

%% Direction V1
figure 

idx = 1:part_n;

x0    = linspace(-0.5,  0.6, part_n); 
xdist = linspace(-0.5,  0.6, part_n); 
xspd  = linspace(1,  2.1, part_n); 

subplot(1,2,1)
errorbar(x0,beta_dir_0(1,:),beta_dir_0(2,:),beta_dir_0(3,:),'o','MarkerFaceColor','auto','LineWidth',1);

yline(0,'--');
xlim([x0(1)-0.1,x0(end)+0.1])
ylim([-25,25])
xticks(x0)
xticklabels(participants)
legend({'Intercept'}, 'Location','northeast');
xlabel('Participant');
ylabel('Intercept');
box off


subplot(1,2,2)
errorbar(xdist, beta_dir_dist(1,:), beta_dir_dist(2,:), beta_dir_dist(3,:),'o','MarkerFaceColor','auto','LineWidth',1);
hold on
errorbar(xspd,  beta_dir_spd(1,:),  beta_dir_spd(2,:),  beta_dir_spd(3,:),'o','MarkerFaceColor','auto','LineWidth',1); 
hold off

yline(0,'--');
xlim([xdist(1)-0.2,xspd(end)+0.2])
ylim([-0.2,0.2])
xticks([mean(xdist),mean(xspd)])
xticklabels(["Distance","Speed"])
xlabel('Factor');
ylabel('Slope');
legend({'\beta_{distance}','\beta_{speed}'}, 'Location','northeast');
box off
saveas(gcf,fullfile('results', 'plots', 'fig4','beta_dir_v1.png'));
%% Direction V2
figure
xdist = linspace(-0.25,  43.75, part_n); 
xspd  = linspace(0.25,  44.25, part_n); 
errorbar(xdist, beta_dir_dist(1,:), beta_dir_dist(2,:), beta_dir_dist(3,:),'o','MarkerFaceColor','auto','LineWidth',1);
hold on
errorbar(xspd,  beta_dir_spd(1,:),  beta_dir_spd(2,:),  beta_dir_spd(3,:),'o','MarkerFaceColor','auto','LineWidth',1); 
hold off

yline(0,'--');
xlim([xspd(1)-2,xspd(end)+2])
ylim([-0.2,0.2])

xticks(mean([xdist;xspd]))
xticklabels(participants)
xlabel('Participant');
ylabel('Slope');
legend({'\beta_{distance}','\beta_{speed}'}, 'Location','northeast');
box off
hold off

saveas(gcf,fullfile('results', 'plots', 'fig4','beta_dir_v2.png'));