function [displayInfo] = startExp(subj,datetime,rng)
%Starts the screen in psychtoolbox. 


displayInfo = struct();             %Setting up structure to save screen stats

displayInfo.subj = subj;            %adding subject info
displayInfo.dateTime = datetime;    %data collection date and time
displayInfo.randSeed = rng;         %saving the randomization seed

PsychDefaultSetup(2);               % Set default settings for psychtoolbox

screens = Screen('Screens');        % Get the screen numbers

displayInfo.screens = screens;


screenNumber = max(screens)-1;       % Draw to the external screen if avaliable

displayInfo.screenNumber = screenNumber;

white = 255;   % Define white color spaces
black = 0;   % Define black color spaces
grey = white/2;                     % Define grey color spaces

displayInfo.whiteVal = white;
displayInfo.blackVal = black;
displayInfo.greyVal = grey;

[window, windowRect] = Screen('OpenWindow', 3, grey); % Open a grey on screen window
[window2, windowRect2] = Screen('OpenWindow', 1, black); %turn on tablet screen

displayInfo.window = window;
displayInfo.window2 = window2;
displayInfo.windowRect = windowRect;
displayInfo.windowRect2 = windowRect2;

[screenXpixels, screenYpixels] = Screen('WindowSize', window); % Get the size of the on screen window
[screenXpixels2, screenYpixels2] = Screen('WindowSize', window2);

screenAdj = screenYpixels - (1080/1920)*screenXpixels;

displayInfo.screenAdj = screenAdj;
displayInfo.screenXpixels = screenXpixels;
displayInfo.screenYpixels = screenYpixels;
displayInfo.screenXpixels2 = screenXpixels2;
displayInfo.screenYpixels2 = screenYpixels2;

ifi = Screen('GetFlipInterval', window); % Query the frame duration

displayInfo.ifi = ifi;
displayInfo.framerate = Screen('NominalFrameRate',displayInfo.window);

[xCenter, yCenter] = RectCenter(windowRect); % Get the centre coordinate of the window
[xCenter2, yCenter2] = RectCenter(windowRect2); % Get the centre coordinate of the window

displayInfo.xCenter = xCenter;
displayInfo.yCenter = yCenter;

displayInfo.xCenter2 = xCenter2;
displayInfo.yCenter2 = yCenter2;

Screen('TextFont', window, 'Arial');    %font for window
Screen('TextSize', window, 15);         %font size for window
 
%HideCursor(screenNumber); %turn off the mouse
end

