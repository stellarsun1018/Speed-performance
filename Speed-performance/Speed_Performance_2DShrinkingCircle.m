% clear all
Screen('Preference', 'SkipSyncTests', 1); 
cd('C:\Users\labadmin\Documents\Qingjie-GitHub\Speed-performance');

subj = 'pilot';  
dateTime = clock;                % get time for seed             
rng(sum(100*dateTime) );      % 也就是给�?组实验/数�?�dataset 编一个编�?，确�?这组实验的�?�track
expName = 'practice_traj';
session = 01;
redoCalib = 0;

[displayInfo] = startExp(subj,datetime,rng);
[displayInfo] = screenVisuals(displayInfo);
if exist(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_tform.mat']) && redoCalib == 0
    load(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_tform.mat']) %load calibration in case of restart
    load(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_calibration.mat'])
    
else
    [tform, calibration,startPhase] = penCalib(displayInfo);
    save(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_tform.mat'],'tform') % tform.mat —— 存储几何�?��?�信�?�（�??标�?��?�）
    save(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_calibration.mat'],'calibration')
end
mode = 0; % lift = 1, slide = 0
%%
start_size = 20;
cursor_size = 5;
pixellength = 0.248; % �?个�?素对应的物�?�长度（�?��? mm）
wait = 0.5;
patience = 0.5;
% story(1) = wait：等待阶段（目标显示但�?能移动）
% story(2) = wait + lifespan(current_block)：目标开始消失
% story(3) = wait + lifespan(current_block) + patience：实验的最�?�时�?

% 例如，第二个block中：
% 
% 等待阶段：wait = 0.5秒（目标出现，无法移动，仅观察）
% 目标�?�?消失阶段：lifespan = 0.6秒（目标尺寸从满值�?�?�?�?到零，得分�?�?从满分�?低到零）
% 最大�?�?等待时长：patience = 0.5秒（�?过0.5秒若�?未完�?动作，试验将强制结�?�）
% 总时间长度 = wait + lifespan + patience = 0.5 + 0.6 + 0.5 = 1.6秒

topBuff = [0 0 displayInfo.screenXpixels displayInfo.screenAdj/2]; %black bar at top of screen 设置顶部和底部的黑色缓冲区域，用于实验界�?�布局
bottomBuff = [0 displayInfo.screenYpixels-displayInfo.screenAdj/2 displayInfo.screenXpixels displayInfo.screenYpixels]; %black bar at bottom of screen

%% Task Parameters

dists_n = 3; % 3 kind of disciances
UniRandRadius = 50; % �?��? pixels，�?机扰动范围
edgesize = 50;
hitrates = [0.3];
 
rep = 3; % repe at 10 times of 3 (kind of distances)* 2(directions) 

scorebar_length = 200;

mmsigma = [30]; % !! needs to be extraced from previous data %目标大�?的标准差
target_sizes = tSizeGen(mmsigma,hitrates,pixellength);
target_sizes = repmat(target_sizes,1,dists_n*rep);
target_sizes = target_sizes';
target_sizes = target_sizes(:)'; 
switch_scale = 1.5;

all_distances = exp(linspace(log(231),log(693),5));
lifespan = [0.6,0.6*3^(0.25),0.6*3^(0.5)]; %[1.0,0.6,0.8,0.4]; %[1.1,0.9,1.0,0.8,0.6,0.7] 设定�?�blocks中target的�?�?�时长  %lifespan控制了�?�试者实际�?�用的�?�?�?�?少的目标"�?��?e时间窗，这一时间越短，任务难度越高（因为�?�试者必须更快速地完�?任务以�?�得更高分数）。
block_n = length(lifespan); % 实验有4个blocks，�?个block有10*3*2个trails
%% Trial
speedthreshold = 10; % pixel per second, equals to 2.48 mm/s
data = [];
traXtotal = NaN(1,round(framerate * (wait+max(lifespan)+patience)));
traYtotal = NaN(size(traXtotal));
testtimes = zeros(1,10000); % 10 seconds
framerate = Screen('NominalFrameRate',displayInfo.window); % 获�?�显示器的帧率
frames = framerate * 5; % start/preparing page time out at 5 seconds % 计算 5 秒钟内的总帧数
instruct = 'Good luck, try hard, and have fun!';
HideCursor;
Screen('FillRect', displayInfo.window, displayInfo.blackVal);
 

while true
    DrawFormattedText(displayInfo.window,instruct,'center','center', displayInfo.whiteVal); 
    Screen('Flip', displayInfo.window);
    [~,~,b] = GetMouse;
    if b(1)
        break
    end
    if KbCheck
        Screen('CloseAll')
        break
    end
end


for current_block = 1:block_n % block_n = which block it is
    switch current_block
        case 1
            distances = all_distances(1:3);
        case 2
            distances = all_distances(2:4);
        case 3
            distances = all_distances(3:5);
    end

    distances = repmat(distances,1,length(hitrates)*rep); % 2:end-1 选�?�去掉第一个和最�?�一个点；然�?�将这三个�?离�?�?10次 distances = [175, 350, 525, 175, 350, 525...
    seeds = [randperm(size(distances,2)), randperm(size(distances,2))];
    randdists = distances(seeds); % 按照 seeds 的顺�?�?新排列 distances，形�? randdists；randdists 代表 �?次实验中的目标�?离（�?机排列�?�）。
    randdists = randdists(:);
    randsizes = target_sizes(seeds);
    randsizes = randsizes(:);
    params = NaN(length(randdists),11); % 会用实验数�?�覆盖 NaN 值。
    trax = NaN(length(randdists),round(framerate * (wait+lifespan(current_block)+patience)));
    tray = NaN(length(randdists),round(framerate * (wait+lifespan(current_block)+patience)));
    trial_n = length(randdists); % trial_n 计算当�? block 内试验次数
    trials = ones(1,trial_n);
    i = 0;
    DrawFormattedText(displayInfo.window,['Next Block: ' num2str(lifespan(current_block)) ' seconds interval'],'center','center',displayInfo.whiteVal); % not sure how to get this centered yet
    Screen('Flip', displayInfo.window);
    pause(2);

    while sum(trials) > 0
        i = i+1;
        stage = 0;
        frame = 0;
        if i == trial_n + 1
            randdists = [randdists ; randdists(trials==true,:)];
            randsizes = [randsizes ; randsizes(trials==true,:)];
            params = [params ; params(trials==true,:)];
            trax = [trax ; trax(trials==true,:)];
            tray = [tray ; tray(trials==true,:)];
            wrong_n = sum(trials);  %计算失败的试验次数
            origin_trial_n = trial_n;
            trial_n = trial_n + wrong_n; %更新试验总数，确�?这些试验在下一轮执行
            trials = zeros(1,trial_n);
            trials(1,origin_trial_n+1:end) = 1;
        end
        if stage == 0
            while true
                [~,~,keyCode] = KbCheck;
                if find(keyCode) == 27 % KbName(27) = 'ESCAPE'
                    Screen('CloseAll');
                    ShowCursor;
                    break
                end
                [x,y,buttons] = GetMouse(displayInfo.window2);
                if rem(i,2)
                    startpos = [displayInfo.windowRect(3)-edgesize,displayInfo.yCenter];
                    theta = -pi/12 + (pi/6) * rand(1); % 这个计算让 theta 在 [-pi/12, pi/12] 之间�?机�?�化，表示目标的角度扰动：角度�??移 (范围：-15° 到 15°)
                    rho = randdists(i)-UniRandRadius + UniRandRadius * 2 * rand(1);
                    [offset(1),offset(2)] = pol2cart(theta,rho);   % pol2cart 是 MATLAB �?�??标转直角�??标的函数:[offset(1),offset(2)]=新目标�??标[x,y]; x=�?⋅cos(θ), y=�?⋅sin(θ)
                    params(i,1:2) = startpos - offset; % params(i,1:2) 存储目标点�??标（包�?�扰动)，代表第 i 个试验的 (x, y) 目标�?置
                else
                    startpos =  [edgesize,displayInfo.yCenter];
                    theta = -pi/12 + (pi/6) * rand(1);
                    rho = randdists(i)-UniRandRadius + UniRandRadius * 2 * rand(1);
                    [offset(1),offset(2)] = pol2cart(theta,rho);
                    params(i,1:2) = startpos + offset;
                end
                params(i,10) = randsizes(i);
                switch_size = switch_scale * params(i,10);
                Screen('DrawDots', displayInfo.window, startpos, start_size, [1 1 1],[],1);
                [xy(1), xy(2)]  = transformPointsForward(tform,x,y);
                Screen('DrawDots', displayInfo.window, xy, 5, [1 0 0],[],1);
                if buttons(1)
                    if sqrt(sum((xy - startpos).^2)) <= start_size % if in start area
                        stage = 1;
                        break
                    end
                end
                Screen('Flip', displayInfo.window);
            end
        end
        % test to get the effective refresh rate ( the rate at which the
        % Screen() flips )
        %     testtimes = nonzeros(testtimes);
        %     ts = [0;testtimes];
        %     t_diff = testtimes - ts(1:end-1);
        %     plot(t_diff)
        %     rfrate = 1/mean(t_diff);
        %     worst_rfrate = 1/max(t_diff);
        if stage == 1
            Screen('FillRect', displayInfo.window, displayInfo.blackVal);
            Screen('Flip', displayInfo.window);
            time = GetSecs;
            story = [wait, wait+lifespan(current_block), wait+lifespan(current_block)+patience]; 
            % wait is the time during which the target stays the same size
            % for motor planning and perception
            % wait + lifespan is the time when the score drops to zero
            % wait + lifespan + patience is the time before the trial force
            % quits because we ran out of patience for waiting for the 
            % participant to move their stylus from the starting point
            t = 1;
            onset_recorded = 0;
            [x,y,~] = GetMouse(displayInfo.window2);
            tSize = params(i,10);
            for frame = 1: framerate * (story(3)+2) 
                cache = [x,y];
                [x, y, buttons] = GetMouse(displayInfo.window2);
                locdiff = sqrt(sum((cache - [x,y]).^2));
                trax(i,frame) = x;
                tray(i,frame) = y;
                [xy(1), xy(2)]  = transformPointsForward(tform,x,y); 
                max_dot_size = 189; % GPU support
                if frame <= framerate * (story(3)+2)
                    Screen('DrawDots', displayInfo.window, startpos, start_size, [1 1 1],[],1); %Starting point
                    if frame >= framerate * story(1)
                        percent_score = max(1-((frame ./ framerate)-story(1)) / lifespan(current_block),0);
                    else 
                        percent_score = 1;
                    end
                        dot_size = max(min(percent_score * max_dot_size, max_dot_size),1); % limit the dot size
                        Screen('DrawDots', displayInfo.window, params(i,1:2), dot_size, [0 0 1], [], 1);

                        %Screen('DrawDots', displayInfo.window, params(i,1:2), 2 * percent_score * scorebar_length, [0 0 1], [], 1);
                        Screen('DrawLine', displayInfo.window, [1 1 1], displayInfo.xCenter - percent_score * scorebar_length,displayInfo.yCenter-200, displayInfo.xCenter + percent_score * scorebar_length,displayInfo.yCenter-200,5);
                        DrawFormattedText(displayInfo.window,['Score = ' num2str(percent_score * 10)],'center',displayInfo.yCenter-220, displayInfo.whiteVal);
                    
                    Screen('Flip', displayInfo.window);
                    
                    if frame > framerate * story(3) 
                        DrawFormattedText(displayInfo.window,'Too Slow!','center','center', displayInfo.whiteVal); % not sure how to get this centered yet
                        Screen('Flip',displayInfo.window);
                        trials(i) = 1;
                        params(i,4) = NaN; 
                        pause(1)
                        break
                    end
                    if norm(xy - startpos) >= start_size
                        if buttons(1)+mode ==0
                            DrawFormattedText(displayInfo.window,'Stylus Lifted!','center','center', displayInfo.whiteVal); % not sure how to get this centered yet
                            Screen('Flip', displayInfo.window);
                            trials(i) = 1;
                            pause(1)
                            break   
                        end
                        if frame < framerate * story(1)
                            DrawFormattedText(displayInfo.window,'Too Early!','center','center', displayInfo.whiteVal); % not sure how to get this centered yet
                            Screen('Flip', displayInfo.window);
                            trials(i) = 1;
                            pause(1)
                            break
                        elseif ~onset_recorded
                            onset_t = frame / framerate;
%                             abs_onset_time = GetSecs;
                            params(i,4) = onset_t;
                            onset_recorded = 1;
                        end
                        [keyIsDown,~,keyCode] = KbCheck;
                        TarUB = params(i,1)- percent_score * scorebar_length;
                        TarLB = params(i,1)+ percent_score * scorebar_length;
                        if (locdiff <= speedthreshold/framerate && ~mode) || (buttons(1) && mode)
                            if norm(xy - startpos) < randdists(i)/2
                                DrawFormattedText(displayInfo.window,'Not Even Close :(','center','center', displayInfo.whiteVal);
                                Screen('Flip', displayInfo.window);
                                trials(i) = 1;
                                pause(1)  
                            else
                                distance_to_target = sqrt((xy(1) - params(i,1))^2 + (xy(2) - params(i,2))^2);

                                if distance_to_target <= dot_size./2
                                    hit = 1; % success hit
                                    bar_color = [1 1 1]; % white
                                else
                                    hit = 0; % no
                                    bar_color = [1 0 0]; % red
                                end                                
                                end_t = frame / framerate; %initialize end_t
                                endpos = [x y];
%                                 hit = 1;
%                                 bar_color = [1 1 1];
%                                 if xy(1) < TarUB || xy(1) > TarLB
%                                     hit = 0;
%                                     bar_color = [1 0 0];  

                                if hit
                                    bar_color = [1 1 0];
%                                         percent_score = max(1-(((frame+k) ./ framerate)-story(1)) / lifespan(current_block),0);
                                        Screen('DrawDots',displayInfo.window, params(i,1:2), dot_size,[0 0 1],[],1);
                                        Screen('DrawDots', displayInfo.window, xy, 5, [1 0 0],[],1);
                                        Screen('DrawLine', displayInfo.window, bar_color, displayInfo.xCenter - percent_score * scorebar_length,displayInfo.yCenter-200, displayInfo.xCenter + percent_score * scorebar_length,displayInfo.yCenter-200,5);
                                        DrawFormattedText(displayInfo.window,['Score = ' num2str(percent_score * 10)],'center',displayInfo.yCenter-220, displayInfo.whiteVal);
                                        DrawFormattedText(displayInfo.window,[num2str(length(trials)-sum(trials)) '/' num2str(length(seeds)) ' finished'],'center','center', displayInfo.whiteVal);
                                        Screen('Flip', displayInfo.window);
                                    score = percent_score * 10 * hit;
                                else
%                                     dot_size = tSize;
                                    Screen('DrawDots', displayInfo.window, params(i,1:2), dot_size, [0 0 1], [], 1);
%                                     Screen('DrawDots',displayInfo.window, params(i,1:2), tSize,[0 1 0],[],1);
%                                     Screen('FrameOval',displayInfo.window, [1 1 0], [params(i,1)-switch_size./2,params(i,2)-switch_size./2,params(i,1)+switch_size./2,params(i,2)+switch_size./2]);
                                    Screen('DrawDots', displayInfo.window, xy, 5, [1 0 0],[],1);
                                    Screen('DrawLine', displayInfo.window, bar_color, displayInfo.xCenter - percent_score * scorebar_length,displayInfo.yCenter-200, displayInfo.xCenter + percent_score * scorebar_length,displayInfo.yCenter-200,5);
                                    DrawFormattedText(displayInfo.window,['Miss :('],'center',displayInfo.yCenter-220, displayInfo.whiteVal);
                                    DrawFormattedText(displayInfo.window,[num2str(length(trials)-sum(trials)+1) '/' num2str(length(seeds)) ' finished'],'center','center', displayInfo.whiteVal);
                                    Screen('Flip', displayInfo.window);
                                    score = percent_score * 10 * hit;
                                end
                                rest_of_trial = story(3) - end_t;
                                pause(rest_of_trial);
                                params(i,3) = 0;
                                params(i,5) = end_t;
                                params(i,6:7) = endpos;
                                params(i,8:9) = startpos;
                                params(i,10) = randsizes(i);
                                params(i,11) = score;
                                trials(i) = 0;
                            end
                            break
                        end                            
                    end
                end
%                 testtimes(t) = GetSecs; % for testing temporal resolution
%                 t = t+1;
                [~,~,keyCode] = KbCheck;
                if find(keyCode) == 27 % KbName(27) = 'ESCAPE'
                    Screen('CloseAll');
                    ShowCursor;
                    break
                end
            end
        end
    end
    data = [data;params];
    save(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) 'c_' date,'_rawtotal.mat'],'data');
    xTrajTelomere = NaN(size(trax,1),size(traXtotal,2));
    xTrajTelomere(1:size(trax,1),1:size(trax,2)) = trax;
    traXtotal = [traXtotal;xTrajTelomere];
    yTrajTelomere = NaN(size(tray,1),size(traYtotal,2));
    yTrajTelomere(1:size(tray,1),1:size(tray,2)) = tray;
    traYtotal = [traYtotal;yTrajTelomere];
    save(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) 'c_' date,'_traXtotal.mat'],'traXtotal')
    save(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) 'c_' date,'_traYtotal.mat'],'traYtotal')
    while true
        DrawFormattedText(displayInfo.window,'Block finished. Press any key to proceed to next block.','center','center', displayInfo.whiteVal); 
        Screen('Flip', displayInfo.window);
        if KbCheck
            break
        end
    end
end
Screen('CloseAll');
ShowCursor;
%%
index = NaN(size(data,1),1);
for i = 1:size(data,1)
    index(i) = ~isnan(sum(data(i,:)));
end
valid = data(index==true,:);
validTraX = traXtotal(index==true,:);
validTraY = traYtotal(index==true,:);
%%
pixellength = 0.248;
copy = valid;
copy(:,[1,2]) = transformPointsInverse(tform,copy(:,[1,2]));
copy(:,[8,9]) = transformPointsInverse(tform,copy(:,[8,9]));
copy(:,10) = sqrt(sum((copy(:,1:2) - copy(:,8:9)).^2,2)) .* pixellength;
copy(:,[11,12]) = [copy(:,1)*pixellength (1080 - copy(:,2))*pixellength];
copy(:,[13,14]) = [copy(:,6)*pixellength (1080 - copy(:,7))*pixellength]; % 1080 = tablet pixel height
copy(:,15) = valid(:,10) .* pixellength;
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
%%
save(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_trialdata.mat'],'copy')
%%
% copy column contents:
% 1,2: target x and y in wac pixels 
% 3: switch time, 0 means no switch made
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