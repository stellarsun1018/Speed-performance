% clear all
Screen('Preference', 'SkipSyncTests', 1); 
cd('C:\Users\labadmin\Documents\Qingjie-GitHub\Speed-performance');

subj = 'pilot';  
dateTime = clock;                %get  s time for seed  
rng(sum(100*dateTime) );
expName = 'practice_traj';
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
hitrates = [0.3];

rep = 5;
distances = linspace(edgesize,displayInfo.windowRect(3)-edgesize,dists_n+2)-edgesize;
distances = repmat(distances(2:end-1),1,length(hitrates)*rep);
scorebar_length = 200;
penalty = 0.8;

mmsigma = [30]; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! needs to be extraced from previous data
target_sizes = tSizeGen(mmsigma,hitrates,pixellength);
target_sizes = repmat(target_sizes,1,dists_n*rep);
target_sizes = target_sizes';
target_sizes = target_sizes(:)';
switch_scale = 1.5;
lifespan = [2,2,2];
block_n = length(lifespan);
%% Trial
speedthreshold = 10; % pixel per second, equals to 2.48 mm/s
data = [];
traXtotal = [];
traYtotal = [];
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

for j = 1:block_n
    seeds = [randperm(size(distances,2)), randperm(size(distances,2))];
    randdists = distances(seeds);
    randdists = randdists(:);
    randsizes = target_sizes(seeds);
    randsizes = randsizes(:);
    params = NaN(length(randdists),11);
    trax = NaN(length(randdists),round(framerate * (wait+lifespan(j)+patience)));
    tray = NaN(length(randdists),round(framerate * (wait+lifespan(j)+patience)));
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
            randsizes = [randsizes ; randsizes(trials==true,:)];
            params = [params ; params(trials==true,:)];
            trax = [trax ; trax(trials==true,:)];
            tray = [tray ; tray(trials==true,:)];
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
            for frame = 1: framerate * (story(3)+2) 
                cache = [x,y];
                [x, y, buttons] = GetMouse(displayInfo.window2);
                locdiff = sqrt(sum((cache - [x,y]).^2));
                trax(i,frame) = x;
                tray(i,frame) = y;
                [xy(1), xy(2)]  = transformPointsForward(tform,x,y);
                
                if frame <= framerate * (story(3)+2)
                    Screen('DrawDots', displayInfo.window, startpos, start_size, [1 1 1],[],1);
                    if frame >= framerate * story(1)
                        percent_score = max(1-((frame ./ framerate)-story(1)) / lifespan(j),0);
                        %Screen('DrawDots',displayInfo.window, params(i,1:2), tSize,[0 1 0],[],1);
                        TarUB = params(i,1)- percent_score * scorebar_length;
                        TarLB = params(i,1)+ percent_score * scorebar_length;
                        Screen('FillRect', displayInfo.window, [0 0 1], [TarUB,300,TarLB,600]);
%                         Screen('FrameOval',displayInfo.window, [1 1 0], [params(i,1)-switch_size./2,params(i,2)-switch_size./2,params(i,1)+switch_size./2,params(i,2)+switch_size./2]);
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
                        TarUB = params(i,1)- percent_score * scorebar_length;
                        TarLB = params(i,1)+ percent_score * scorebar_length;
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
                                if xy(1) < TarUB || xy(1) > TarLB
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
%                                     Screen('FrameOval',displayInfo.window, [1 1 0], [params(i,1)-switch_size./2,params(i,2)-switch_size./2,params(i,1)+switch_size./2,params(i,2)+switch_size./2]);
                                    Screen('DrawDots', displayInfo.window, xy, 5, [1 0 0],[],1);
                                    Screen('DrawLine', displayInfo.window, bar_color, displayInfo.xCenter - percent_score * scorebar_length,displayInfo.yCenter-200, displayInfo.xCenter + percent_score * scorebar_length,displayInfo.yCenter-200,5);
                                    if hit
                                        DrawFormattedText(displayInfo.window,['Score = ' num2str(percent_score * 10)],'center',displayInfo.yCenter-220, displayInfo.whiteVal);
                                    else
                                        DrawFormattedText(displayInfo.window,['Miss :('],'center',displayInfo.yCenter-220, displayInfo.whiteVal);
                                    end
                                    DrawFormattedText(displayInfo.window,[num2str(length(trials)-sum(trials)+1) '/' num2str(length(distances)*2) ' finished'],'center','center', displayInfo.whiteVal);
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