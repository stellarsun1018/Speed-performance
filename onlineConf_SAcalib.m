% clear all
Screen('Preference', 'SkipSyncTests', 1); 
cd('C:\Users\labadmin\Documents\onlineConfExperiment');

subj = 'test';  
dateTime = clock;                %get  s time for seed  
rng(sum(100*dateTime) );
expName = 'speedAcc';
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
target_size = 30;
cursor_size = 5;
pixellength = 0.248;
wait = 0.25;
topBuff = [0 0 displayInfo.screenXpixels displayInfo.screenAdj/2]; %black bar at top of screen
bottomBuff = [0 displayInfo.screenYpixels-displayInfo.screenAdj/2 displayInfo.screenXpixels displayInfo.screenYpixels]; %black bar at bottom of screen
% sound_gap = 0.5;
% sound_length = 0.01;
% temp_acc_threshold = 0.1;
%% sound
sampleRate = Snd('DefaultRate');
frequency = 400;
durationSec = 5;
fVolume = 0.4;
sound_length = 0.05;
nSample = sampleRate*durationSec;
soundVec = sin(2*pi*frequency*(1:nSample)/sampleRate);
soundVec = soundVec * fVolume;
soundVec = repmat(soundVec,2,1); % because i have two speakers hence two channels
pahandle = PsychPortAudio('Open');
%% Task Parameters
% dists_n = 2;
% block = yCenter/(dists_n+1); 
% assign = [1,0;0,1;-1,0;0,-1];
% assign = [assign; [sqrt(2),sqrt(2);sqrt(2),-sqrt(2);-sqrt(2),sqrt(2);-sqrt(2),-sqrt(2)]./2];
% assign = [assign; assign.*2];
% target_pos = repmat([xCenter,yCenter],16,1);
% target_pos = target_pos + assign.*block;

dists_n = 1;
edgesize = 50;
repetition = 50;
distances = linspace(edgesize,displayInfo.windowRect(3)-edgesize,dists_n+2)-edgesize;
distances = repmat(distances(2:end-1),1,repetition);

% target_pos = linspace(edgesize,displayInfo.windowRect(3)-edgesize,dists_n+2);
% target_pos = repmat(target_pos(2:end-1),1,5);

% ind = randperm(size(target_pos,1));
% target_pos = target_pos(ind,:);
% [~,index] = sort(ind); % later just use trial(index) to sort each position  block of trials
gap = [0.5,0.4,0.7,0.3,0.9];
gap_n = length(gap);
% randtarpos = NaN(length(gap)*length(target_pos),2);
% for i = 1:length(gap)
%     randtarpos((i-1)*length(target_pos)+1:i*length(target_pos),:) = target_pos(randperm(size(target_pos,1)),:);
% end
% gap = repmat(gap,length(target_pos),1);
% gap = gap(:);

% params = [randtarpos gap];
%% Durations

% gap_n = length(gap);
% gap = repmat(gap,length(target_pos),1);
% target_pos = repmat(target_pos,gap_n,1);
% gap = gap(:);
% params = [target_pos gap];

crosslength = 10;
cross_shape = [0,0,-crosslength,crosslength;-crosslength,crosslength,0,0];
%% visuallize stimuli positions
% plot(target_pos(:,1),target_pos(:,2),'o')
% xlim(displayInfo.windowRect([1 3]))
% ylim(displayInfo.windowRect([2 4]))


%% Trial
speedthreshold = 10; % pixel per second, equals to 2.48 mm/s
data = [];
testtimes = zeros(1,10000); % 10 seconds
framerate = Screen('NominalFrameRate',displayInfo.window);
frames = framerate * 5; % start/preparing page time out at 5 seconds
instruct = 'In each of the following trials, you will be shown a starting point as a greee dot and your cursor as red dot. \n Tap the starting point to start the trial. During the trial, a target as a white cross will appear. \n There will also be three beeps separated by the same time interval. \n Try to start moving your cursor as the third beeps starts and hit the target after the same interval since the start of your movement. \n Each block has its own interval. You will be notified before each block. \nTap the tablet to continue.  ';
HideCursor;
% WinTabMex(0,displayInfo.window);
% WinTabMex(2);
Screen('FillRect', displayInfo.window, displayInfo.blackVal);
 
while true
    DrawFormattedText(displayInfo.window,instruct,'center','center', displayInfo.whiteVal); 
    Screen('Flip', displayInfo.window);
    [~,~,b] = GetMouse;
    if b(1)
        break
    end
%     pkt = WinTabMex(5);
%     if ~isempty(pkt) %ie the queue is empty
%         buttons = pkt(4);
%         if buttons(1) || KbCheck
%             break
%         end
%     end
    if KbCheck
        Screen('CloseAll')
        break 
    end
end
% pktram = [pkt];
% WinTabMex(3);
for j = 1:gap_n
    earlyrange = wait*framerate:1:(wait + sound_length*2 + gap(j)*1 + gap(j)*0.8) * framerate;
    randdists = [distances(randperm(size(distances,2))),distances(randperm(size(distances,2)))];
    randdists = randdists(:);
%     randtarpos = [randtarpos;randtarpos+displayInfo.xCenter];
%     randtarpos = randtarpos(:);
%     randtarpos(:,2) = displayInfo.yCenter;
    params = NaN(length(randdists),9);
%     params(:,1:2) = randtarpos;
    params(:,3) = gap(j);
    trial_n = length(randdists);
    trials = ones(1,trial_n);
    i = 0;
    DrawFormattedText(displayInfo.window,['Next Block: ' num2str(gap(j)) ' seconds interval'],'center','center',displayInfo.whiteVal); % not sure how to get this centered yet
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
%             WinTabMex(2);
            while true
                if KbCheck
                    Screen('CloseAll')
                    ShowCursor
                    break
                end
                [x,y,buttons] = GetMouse(displayInfo.window2);
                if rem(i,2)
                    startpos = [displayInfo.windowRect(3)-edgesize,displayInfo.yCenter];
                    theta = -pi/12 + (pi/6) * rand(1);
                    rho = randdists(i)-75 + 150 * rand(1);
                    [offset(1),offset(2)] = pol2cart(theta,rho);
                    params(i,1:2) = startpos - offset;
                else
                    startpos =  [edgesize,displayInfo.yCenter];
                    theta = -pi/12 + (pi/6) * rand(1);
                    rho = randdists(i)-75 + 150 * rand(1);
                    [offset(1),offset(2)] = pol2cart(theta,rho);
                    params(i,1:2) = startpos + offset;
                end
                
                Screen('DrawDots', displayInfo.window, startpos, start_size, [1 1 0],[],1);
%                 pkt = WinTabMex(5);
%                 pktram = [pktram pkt];
%                 if ~isempty(pkt)
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
%             WinTabMex(3);
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
%             SetMouse(startpos(1),startpos(2),displayInfowindow);
            Screen('FillRect', displayInfo.window, displayInfo.blackVal);
            Screen('Flip', displayInfo.window);
            time = GetSecs;
            audstarts = [wait, wait+sound_length+gap(j), wait+sound_length*2+gap(j)*2];
            %         audstops = audstarts + sound_length;
            t = 1;
            onset_recorded = 0;
            [x,y,~] = GetMouse(displayInfo.window2);
            for frame = 1: framerate * ( audstarts(1) + sound_length*4 + gap(j)*3 + 10) % Using the max display frame rate assumes the computer being able to finish the loop within one frame's time
                cache = [x,y];
                [x, y, buttons] = GetMouse(displayInfo.window2);
                locdiff = sqrt(sum((cache - [x,y]).^2));
                [xy(1), xy(2)]  = transformPointsForward(tform,x,y);
%                     buttons = pkt(4);
%                 end
                if frame <= framerate * (audstarts(1) + sound_length*4 + gap(j)*3 + gap(j)*0.5)
                    Screen('DrawDots', displayInfo.window, startpos, start_size, [1 1 0],[],1);
                    Screen('DrawLines',displayInfo.window, cross_shape, 2, displayInfo.whiteVal,params(i,1:2));
%                     if ~isempty(pkt)
%                         Screen('DrawDots', displayInfo.window, xy, 5, [1 0 0],[],1);
%                     end
                    Screen('Flip', displayInfo.window);
                    if ismember(frame, round(audstarts .* framerate))
                        PsychPortAudio('FillBuffer',pahandle,soundVec);
                        PsychPortAudio('Start',pahandle,0,0,1,GetSecs+sound_length);
                        %                 elseif ismember(frame, round(audstops .* framerate))
                        %                     PsychPortAudio('Stop',pahandle);
                    end


                    if sqrt(sum((xy - startpos).^2)) >= start_size
                        if buttons(1)+mode ==0
                            DrawFormattedText(displayInfo.window,'Stylus Lifted!','center','center', displayInfo.whiteVal); % not sure how to get this centered yet
                            Screen('Flip', displayInfo.window);
                            trials(i) = 1;
                            pause(1)
                            break   
                        end
                        if ismember(frame,earlyrange)
                            DrawFormattedText(displayInfo.window,'Too Early!','center','center', displayInfo.whiteVal); % not sure how to get this centered yet
                            Screen('Flip', displayInfo.window);
                            trials(i) = 1;
                            pause(1)
                            break
                        elseif ~onset_recorded
                            onset_t = frame / framerate;
                            params(i,4) = onset_t;
                            onset_recorded = 1;
                        end
                        if (locdiff <= speedthreshold/framerate && ~mode) || (buttons(1) && mode)
                            if sqrt(sum((xy-startpos).^2)) < randdists(i)/2
                                DrawFormattedText(displayInfo.window,'Not Even Close :(','center','center', displayInfo.whiteVal); % not sure how to get this centered yet
                                Screen('Flip', displayInfo.window);
                                trials(i) = 1;
                                pause(1)  
                            else
                                end_t = frame / framerate;
                                endpos = [x y];
                                Screen('DrawDots', displayInfo.window, xy, 5, [1 0 0],[],1);
                                Screen('DrawLines',displayInfo.window, cross_shape, 2, displayInfo.whiteVal,params(i,1:2));
                                Screen('Flip', displayInfo.window);
                                pause(1.5);
                                params(i,5) = end_t;
                                params(i,6:7) = endpos;
                                params(i,8:9) = startpos;
                                
                                trials(i) = 0;
                            end
                            break
                        end                            
                    end
                elseif frame > framerate * (audstarts(1) + sound_length*4 + gap(j)*3 + gap(j)*0.5) % 0.5 is arbitary, just how much time we expect the time out should be
                    DrawFormattedText(displayInfo.window,'Too Slow!','center','center', displayInfo.whiteVal); % not sure how to get this centered yet
                    Screen('Flip',displayInfo.window);
                    trials(i) = 1;
                    params(i,4) = NaN; 
                    pause(1)
                    break
                end
                testtimes(t) = GetSecs; % for testing temporal resolution
                t = t+1;
                if KbCheck
                    Screen('CloseAll');
                    ShowCursor;
                    break
                end
            end
%             WinTabMex(3);
        end
        
        
    end
    data = [data;params];
%     WinTabMex(3);
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

copy(:,[11,12]) = [copy(:,1)*pixellength (displayInfo.windowRect(4) - copy(:,2))*pixellength];
copy(:,[13,14]) = [copy(:,6)*pixellength (displayInfo.windowRect(4) - copy(:,7))*pixellength];
for i = 1:gap_n
    copy(copy(:,3)==gap(i),15) = i;
end
copy(:,16) = copy(:,5) - copy(:,4);
pointfivemean = mean(copy(copy(:,15)==1,8));
onemean = mean(copy(copy(:,15)==2,8));
copy(:,17) = sqrt( (copy(:,13)-copy(:,11)).^2 + (copy(:,14)-copy(:,12)).^2 );
copy(:,10) = sqrt(sum((copy(:,1:2) - copy(:,8:9)).^2,2)) .* pixellength;
% copy(:,9) = copy(:,9)./copy(:,10);
copy(:,18) = 1:length(copy);
%%
% ttest2(copy(copy(:,15)==1,17),copy(copy(:,15)==2,17));
%%

%%

% %%
% figure(1)
% boxplot([copy(copy(:,15)==1,17),copy(copy(:,15)==2,17)])
% xticklabels(gap);
% xlabel('Interval Duration');
% ylabel('Error Size');
% figure(2)
% plot(copy(copy(:,15)==1,10),copy(copy(:,15)==1,17),'o');
% hold on
% plot(copy(copy(:,15)==2,10),copy(copy(:,15)==2,17),'o');
% 
% xlabel('target distance in mm');
% ylabel('error size in mm');
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
% 15: target duration condition
% 16: actual duration of the reach
% 17: error size in mm
% 18: trial order
%%
% testtimes = nonzeros(testtimes);
% ts = [0;testtimes];
% t_diff = testtimes - ts(1:end-1);
% % plot(t_diff)
% rfrate = 1/median(t_diff(2:end));
%% done
% audio lagging
% temporal accuracy due to refresh rate
% randomize trials
% timeout feedback
% target location in circles
% milimeters for targets and end points
% interval types
% trial chronological order index
%% not done
% save parameters
% too early click
% too close to start
% cut scene
% instruction

%% scratch
% tic
% for i = 0:1000*length(randtarpos)-1
%     [a,b]= transformPointsForward(tform,randtarpos(rem(i,32)+1,1),randtarpos(rem(i,32)+1,2));
% end
% toc
%% sample rate & frame rate
% WinTabMex(2);
% data = NaN(2,600);
% instruct = 'In each of the following trials, you will be shown a screen center as a white dot and your cursor as red dot. \n Click the center to start the trial. During the trial, a target as a white dot will appear. \n There will also be three beeps separated by the same time interval. \n Try to start moving your cursor as the third beeps starts and hit the target after the same interval since the start of your movement. \n Each block has its own interval. You will be notified before each block.';
% HideCursor;
% 
% for i = 1:600
%     pkt = WinTabMex(5);
%     DrawFormattedText(window,instruct,'center','center', white); % not sure how to get this centered yet
%     Screen('Flip', window);
%     data(2,i) = GetSecs;
%     if ~isempty(pkt)
%         data(1,i) = pkt(1);
%     end
% end
% WinTabMex(3);
% Screen('CloseAll') 
% inter = data(:,~isnan(data(1,:)));
% tdiff = inter(2,[8:end]) - inter(2,[7:end-1]);
% 1/mean(tdiff)