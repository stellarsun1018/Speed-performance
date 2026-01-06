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
% participants = {"JHL"};
part_n = numel(participants);

use_session = false;   % true: 用 S%d 方式；false: 用 usable 方式
session = 3;

nBlocks = 3;
nTrialsPerBlock = 240;
lim_scale = 1.2;

pixellength = 0.248;

colors = lines(part_n);

ip_all = NaN(part_n,720);
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

    if sum(copy(:,3) == 0)
        lifespan = [0.6,0.6*3^(0.25),0.6*3^(0.5)]; %[1.0,0.6,0.8,0.4]; %[1.1,0.9,1.0,0.8,0.6,0.7] 设定各blocks中target的不同时长  %lifespan控制了受试者实际可用的、逐渐减少的目标"可见e时间窗，这一时间越短，任务难度越高（因为受试者必须更快速地完成任务以取得更高分数）。
        for i = 1:3
            copy((1+(i-1)*240):(i*240),3) = lifespan(i);
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    distances = copy(:,10);
    durations = copy(:,16);
    end_size = copy(:,15) .* durations ./ copy(:,3);

    id_all = log2( 2 .* distances ./ end_size);
    ip_all(ip,:) = id_all ./ durations;
end

%%

ip_mean = mean(ip_all,2);
ip_sem = std(ip_all,[],2)./ sqrt(size(ip_all,2));


errorbar(ip_mean,ip_sem,'o','MarkerFaceColor','auto')
xticks(1:part_n)
xticklabels(participants)
xlim([0,part_n+1])
xlabel('Participants')
ylabel('Index of Performance')

saveas(gcf,fullfile('results','plots','fig6', 'errobars.png'))

%% Across Distances
for ip = 1:part_n
    part = participants{ip};

    fname_preamble = sprintf('data_onlineConf/usable/%s_sptatialTemporalCostFunc*.mat', part);
    for k = 1:numel(files)
        f = fullfile(files(k).folder, files(k).name);
        load(f);
    end
    subplot(3,4,ip)
    scatter(copy(:,10),ip_all(ip,:),10,'filled')
    title(part)
    if mod(ip,4) == 1
        ylabel('Index pf Performance')
    end

    if ip >= 9
        xlabel('Distances (mm)')
    end

    xlim([65,340])
    ylim([0,37])

    clear copy

end

saveas(gcf,fullfile('results','plots','fig6', 'across_dist.png'))


%% Across Speeds
for ip = 1:part_n
    part = participants{ip};

    fname_preamble = sprintf('data_onlineConf/usable/%s_sptatialTemporalCostFunc*.mat', part);
    for k = 1:numel(files)
        f = fullfile(files(k).folder, files(k).name);
        load(f);
    end
    subplot(3,4,ip)
    scatter(copy(:,22),ip_all(ip,:),10,'filled')
    title(part)
    if mod(ip,4) == 1
        ylabel('Index pf Performance')
    end

    if ip >= 9
        xlabel('Speed (mm/s)')
    end

    xlim([130,630])
    ylim([0,37])

    clear copy

end

saveas(gcf,fullfile('results','plots','fig6', 'across_speed.png'))

