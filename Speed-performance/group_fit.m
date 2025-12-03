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
% participant = 'SX';
lifespan = [0.6,0.6*3^(0.25),0.6*3^(0.5)]; %[1.0,0.6,0.8,0.4]; %[1.1,0.9,1.0,0.8,0.6,0.7] 设定各blocks中target的不同时长  %lifespan控制了受试者实际可用的、逐渐减少的目标"可见e时间窗，这一时间越短，任务难度越高（因为受试者必须更快速地完成任务以取得更高分数）。
for i = 1:3
    copy((1+(i-1)*240):(i*240),3) = lifespan(i);
end
participants = {"JHL", "JH", "LC", "LN", "RC", "SM", "ML", "SX", "SY"};
% participants = {"JH", "JHL", "LN", "LC"};
linear_fits = NaN(4 * 3, numel(participants)); % three blocks * 4 parameters
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
    copy(copy(:,27)>=1 & copy(:,27)<=60, 31)   = 1.0;
    copy(copy(:,27)>=61 & copy(:,27)<=120, 31) = 0.6;
    copy(copy(:,27)>=121 & copy(:,27)<=180,31) = 0.8;
    copy(copy(:,27)>=181 & copy(:,27)<=240,31) = 0.4;

    % 在copy中新增第32列：target shrinking speed (mm/s)
    copy(:,32) = copy(:,15) ./ copy(:,31);

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

        % plot(distances,duration,'o');
        % hold on

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


