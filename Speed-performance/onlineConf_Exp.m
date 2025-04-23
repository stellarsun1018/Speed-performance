% clear all
Screen('Preference', 'SkipSyncTests', 1); 
cd('C:\Users\labadmin\Documents\onlineConfExperiment');

subj = 'pilot';  
dateTime = clock;                %get  s time for seed  
rng(sum(100*dateTime) );
expName = 'switch';
session = 01;
redoCalib = 0;

[displayInfo] = startExp(subj,datetime,rng);
[displayInfo] = screenVisuals(displayInfo);
if exist(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_tform.mat']) && redoCalib == 0
    load(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_tform.mat']) %load calibration incase of restart
    load(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_calibration.mat'])
    
else
    [tform, calibration,startPhase] = penCalib(displayInfo);
    save(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_tform.mat'],'tform')
    save(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_calibration.mat'],'calibration')
end
mode = 0; % lift = 1, slide = 0
%%
start_size = 20;
cursor_size = 5;
pixellength = 0.248;
wait = 1;
patience = 0.5;
topBuff = [0 0 displayInfo.screenXpixels displayInfo.screenAdj/2]; %black bar at top of screen
bottomBuff = [0 displayInfo.screenYpixels-displayInfo.screenAdj/2 displayInfo.screenXpixels displayInfo.screenYpixels]; %black bar at bottom of screen

%% Task Parameters

dists_n = 3;
UniRandRadius = 50;
edgesize = 50;
hitrates = [0.3;0.4;0.5;0.6;0.7];

rep = 3;
distances = linspace(edgesize,displayInfo.windowRect(3)-edgesize,dists_n+2)-edgesize;
distances = repmat(distances(2:end-1),1,length(hitrates)*rep);
scorebar_length = 200;
penalty = 0.8;

mmsigma = [15]; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! needs to be extraced from previous data
target_sizes = tSizeGen(mmsigma,hitrates,pixellength);
target_sizes = repmat(target_sizes,1,dists_n*rep);
target_sizes = target_sizes';
target_sizes = target_sizes(:)';
switch_scale = 1.5;
lifespan = [3,3,3,3,3,3];
gap_n = length(lifespan);
%% Trial
speedthreshold = 10; % pixel per second, equals to 2.48 mm/s
data = [];
testtimes = zeros(1,10000); % 10 seconds
framerate = Screen('NominalFrameRate',displayInfo.window);
frames = framerate * 5; % start/preparing page time out at 5 seconds
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

for j = 1:gap_n
    earlyrange = wait*framerate:1:(wait + sound_length*2 + lifespan(j)*1 + lifespan(j)*0.8) * framerate;
    seeds = [randperm(size(distances,2)), randperm(size(distances,2))];
    randdists = distances(seeds);
    randdists = randdists(:);
    randsizes = target_sizes(seeds);
    randsizes = randsizes(:);
    params = NaN(length(randdists),18);
%     params(:,3) = lifespan(j);
    trial_n = length(randdists);
    trials = ones(1,trial_n);
    i = 0;
    DrawFormattedText(displayInfo.window,['Next Block: ' num2str(lifespan(j)) ' seconds interval'],'center','center',displayInfo.whiteVal); % not sure how to get this centered yet
    Screen('Flip', displayInfo.window);
    pause(2);
    while sum(trials) > 0
        i = i+1;
        stage = 0;
        frame = 0;
        if i == trial_n + 1
            randdists = [randdists ; randdists(trials==true,:)];
            params = [params ; params(trials==true,:)];
            wrong_n = sum(trials);
            origin_trial_n = trial_n;
            trial_n = trial_n + wrong_n;
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
                    theta = -pi/12 + (pi/6) * rand(1);
                    rho = randdists(i)-UniRandRadius + UniRandRadius * 2 * rand(1);
                    [offset(1),offset(2)] = pol2cart(theta,rho);
                    params(i,1:2) = startpos - offset;
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
            story = [wait, wait+lifespan(j), wait+lifespan(j)+patience];
            t = 1;
            onset_recorded = 0;
            switch_recorded = 0;
            [x,y,~] = GetMouse(displayInfo.window2);
            tSize = params(i,10);
            for frame = 1: framerate * (story(3)+2) % Using the max display frame rate assumes the computer being able to finish the loop within one frame's time
                cache = [x,y];
                [x, y, buttons] = GetMouse(displayInfo.window2);
                locdiff = sqrt(sum((cache - [x,y]).^2));
                [xy(1), xy(2)]  = transformPointsForward(tform,x,y);
                
                if frame <= framerate * (story(3)+2)
                    Screen('DrawDots', displayInfo.window, startpos, start_size, [1 1 1],[],1);
                    if frame >= framerate * story(1)
                        percent_score = max(1-((frame ./ framerate)-story(1)) / lifespan(j),0);
                        Screen('DrawDots',displayInfo.window, params(i,1:2), tSize,[0 1 0],[],1);
                        Screen('FrameOval',displayInfo.window, [1 1 0], [params(i,1)-switch_size./2,params(i,2)-switch_size./2,params(i,1)+switch_size./2,params(i,2)+switch_size./2]);
                        Screen('DrawLine', displayInfo.window, [1 1 1], displayInfo.xCenter - percent_score * scorebar_length,displayInfo.yCenter-200, displayInfo.xCenter + percent_score * scorebar_length,displayInfo.yCenter-200,5);
                        DrawFormattedText(displayInfo.window,['Score = ' num2str(percent_score * 10)],'center',displayInfo.yCenter-220, displayInfo.whiteVal);
                    end
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
                        if keyIsDown && ismember(16,find(keyCode))
                            % switch_time = abs_switch_time - abs_onset_time;
                            switch_time = frame/framerate;
                            Screen('DrawDots',displayInfo.window, params(i,1:2), switch_size,[0 1 0],[],1);
                            Screen('Flip', displayInfo.window);
                            switch_recorded = 1;
                            tSize = switch_size;
                            params(i,18) = 1;
                        elseif ~switch_recorded
                            switch_time = NaN;
                            params(i,18) = 0;
                        end
                        if (locdiff <= speedthreshold/framerate && ~mode) || (buttons(1) && mode)
                            if norm(xy - startpos) < randdists(i)/2
                                DrawFormattedText(displayInfo.window,'Not Even Close :(','center','center', displayInfo.whiteVal);
                                Screen('Flip', displayInfo.window);
                                trials(i) = 1;
                                pause(1)  
                            else
                                end_t = frame / framerate;
                                endpos = [x y];
                                hit = 1;
                                bar_color = [1 1 1];
                                if norm(xy - params(i,1:2)) > tSize/2
                                    hit = 0;
                                    bar_color = [1 0 0];
                                end
                                if switch_recorded && hit
                                    end_t = end_t + penalty;
                                    bar_color = [1 1 0];
                                    for k = 1:framerate * penalty
                                        percent_score = max(1-(((frame+k) ./ framerate)-story(1)) / lifespan(j),0);
                                        Screen('DrawDots',displayInfo.window, params(i,1:2), switch_size,[0 1 0],[],1);
                                        Screen('DrawDots', displayInfo.window, xy, 5, [1 0 0],[],1);
                                        Screen('DrawLine', displayInfo.window, bar_color, displayInfo.xCenter - percent_score * scorebar_length,displayInfo.yCenter-200, displayInfo.xCenter + percent_score * scorebar_length,displayInfo.yCenter-200,5);
                                        DrawFormattedText(displayInfo.window,['Score = ' num2str(percent_score * 10)],'center',displayInfo.yCenter-220, displayInfo.whiteVal);
                                        DrawFormattedText(displayInfo.window,[num2str(length(trials)-sum(trials)) '/' num2str(length(distances)) ' finished'],'center','center', displayInfo.whiteVal);
                                        Screen('Flip', displayInfo.window);
                                    end
                                    score = percent_score * 10 * hit;
                                else
                                    Screen('DrawDots',displayInfo.window, params(i,1:2), tSize,[0 1 0],[],1);
                                    Screen('FrameOval',displayInfo.window, [1 1 0], [params(i,1)-switch_size./2,params(i,2)-switch_size./2,params(i,1)+switch_size./2,params(i,2)+switch_size./2]);
                                    Screen('DrawDots', displayInfo.window, xy, 5, [1 0 0],[],1);
                                    Screen('DrawLine', displayInfo.window, bar_color, displayInfo.xCenter - percent_score * scorebar_length,displayInfo.yCenter-200, displayInfo.xCenter + percent_score * scorebar_length,displayInfo.yCenter-200,5);
                                    if hit
                                        DrawFormattedText(displayInfo.window,['Score = ' num2str(percent_score * 10)],'center',displayInfo.yCenter-220, displayInfo.whiteVal);
                                    else
                                        DrawFormattedText(displayInfo.window,['Miss :('],'center',displayInfo.yCenter-220, displayInfo.whiteVal);
                                    end
                                    DrawFormattedText(displayInfo.window,[num2str(length(trials)-sum(trials)) '/' num2str(length(distances)) ' finished'],'center','center', displayInfo.whiteVal);
                                    Screen('Flip', displayInfo.window);
                                    score = percent_score * 10 * hit;
                                end
                                rest_of_trial = story(3) - end_t;
                                pause(rest_of_trial);
                                params(i,3) = switch_time;
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
    while true
        DrawFormattedText(displayInfo.window,'Block finished. Press any key to proceed to next block.','center','center', displayInfo.whiteVal); % not sure how to get this centered yet
        Screen('Flip', displayInfo.window);
        if KbCheck
            break
        end
    end
end
Screen('CloseAll');
ShowCursor;
%%
index = NaN(length(data),1);
for i = 1:length(data)
    index(i) = ~isnan(sum(data(i,:)));
end
valid = data(index==true,:);
%%
copy = valid;
copy(:,[1,2]) = transformPointsInverse(tform,copy(:,[1,2]));
copy(:,[8,9]) = transformPointsInverse(tform,copy(:,[8,9]));
copy(:,10) = sqrt(sum((copy(:,1:2) - copy(:,8:9)).^2,2)) .* pixellength;
copy(:,[11,12]) = [copy(:,1)*pixellength (displayInfo.windowRect(4) - copy(:,2))*pixellength];
copy(:,[13,14]) = [copy(:,6)*pixellength (displayInfo.windowRect(4) - copy(:,7))*pixellength];
copy(:,15) = valid(:,10) .* pixellength;
copy(:,16) = copy(:,5) - copy(:,4);
copy(:,17) = sqrt( (copy(:,13)-copy(:,11)).^2 + (copy(:,14)-copy(:,12)).^2 );
copy(:,27) = 1:length(copy);
copy(:,19:20) = (copy(:,6:7) - copy(:,8:9)) .* pixellength;
copy(:,21) = sqrt(sum((copy(:,6:7) - copy(:,8:9)).^2,2)) .* pixellength;
copy(:,22) = copy(:,21) ./ copy(:,16);
copy(:,[24,25]) = (copy(:,1:2) - copy(:,8:9)) .* pixellength;% relative target coordinate
copy(:,23) = (abs(dot(copy(:,19:20),copy(:,24:25),2) ./ dot(copy(:,24:25),copy(:,24:25),2)) - 1) .*copy(:,10);
copy(:,26) = valid(:,11);
%%
save(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) 'b_' date,'_trialdata.mat'],'copy')
%%
% copy column contents:
% 1,2: target x and y in wac pixels 
% 3: target duration (the gap between beeps)
% 4,5: the onset and end time of the reach
% 6,7: endpoint x and y in wac pixels
% 8,9: start position in wac pixels
% 10: target distance in mm
% 11,12: target x and y in mm
% 13,14: endpoint x and y in mm
% 15: target size in mm
% 16: actual duration of the reach
% 17: error size in mm
% 18: trial order
% 19,20: start position in mm
% 21: actual reach distance
% 22: average speed
% 23: error along the reach direction (vector projection) in mm
% 24,25: relative target position in mm
% 26: score