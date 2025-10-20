function degPositions = runGazeCalibrationStimulus(simulation_mode, device_num, agc_convergence_wait_s, subjectId, experiment_name, session, heightCm, widthCm, viewingDistCm)
% Displays 26-dot gaze calibration stimulus at fixed visual angles, with a brief beep signaling each dot onset.
% 
% TO DO: figure out why dots are displaying for ~3.267 seconds!!!
% Description: 
%   Displays calibration dots on a screen for the subject to fixate. The
%   dots are at various distances in degrees visual angle from a central
%   point. The regime also starts the light logger recording if the option
%   to record is turned on (i.e., not in simulation mode).
% 
% Inputs:
%   simulation_mode             - String & Enum. Choose whether to run full stimulus or simulate 
%                                 a part for testing. 
%
%   device_num                  - Numeric. The device number to use when recording from a real device.
%
%   agc_convergence_wait_s      - Numeric. The time in seconds to wait for the AGC to converge to appropriate settings.
%
%   experiment_name             - String. The filename the recording will be saved under on the light logger. 
%
%   heightCm                    - Numeric. Display screen height in cm.
%
%   widthCm                     - Numeric. Display screen width in cm.
%
%   viewingDistCm               - Numeric. Distance from participant to the screen in cm.
% 
% Example:
%{
    subjectId = 'FLIC_2002';
    sessionNum = 1;
    runGazeCalibrationStimulus("full", 2, 60, subjectId, 'GazeCalibration',sessionNum, 106.7, 192.4)
%}
                          
    arguments 
        simulation_mode {mustBeMember(simulation_mode, ["full", "visual", "bluetooth"])} = "full";
        device_num {mustBeNumeric} = 1;
        agc_convergence_wait_s {mustBeNumeric} = 60;
        subjectId = "test";
        experiment_name = "GazeCalibration";
        session = 1;
        heightCm = 106.7; % 2nd floor LGTV
        widthCm = 192.4; % 2nd floor LGTV
        viewingDistCm = 100; % 1 meter standard
    end
    % Hard-coded parametersk
   
    outerDotRadiusDeg = 1.5;
    innerDotRadiusDeg = 0.25;
    dotTime = 3; % sec 
    cheatTime = 0.11; % measured cpu delay (Sam's laptop, connected to 2nd floor conference room)
    dotTime = dotTime - cheatTime;
    repetitions = 2;
    bgColor = [0 0 0];
    outerDotColor = [255 255 255];
    innerDotColor = [255, 0, 0];
    redColor  = [255   0   0];
    if(simulation_mode == "full" || simulation_mode == "visual")
        AssertOpenGL;
        screenNum = max(Screen('Screens'));
        % Initialize PsychPortAudio for beep
        InitializePsychSound(1);
        sampleRate = 44100;
        beepFreq = 1000;
        beepLengthSec = 0.1;
        nrchannels = 1; 
        % Generate beep
        beep = MakeBeep(beepFreq, beepLengthSec, sampleRate);
        % Open audio device
        pahandle = PsychPortAudio('Open', [], 1, 1, sampleRate, nrchannels);
        % Fill buffer with beep waveform
        PsychPortAudio('FillBuffer', pahandle, beep);
        % Query physical display size (mm) and convert to cm
        if ~exist('widthCm', 'var')
            [widthCm, heightCm] = Screen('DisplaySize', screenNum);
            warning('Using DisplaySize, which may incorrectly estimate screen height and width.')
        end
        % Open window and get its pixel resolution
        [win, winRect] = Screen('OpenWindow', screenNum, bgColor);
        [xCen, yCen] = RectCenter(winRect);
        winWidthPx = winRect(3) - winRect(1);
        winHeightPx = winRect(4) - winRect(2);
        % Print physical size and pixels to confirm display
        fprintf('→ Screen %d: %.1f×%.1f cm physical → %d×%d px window\n', ...
                screenNum, widthCm, heightCm, winWidthPx, winHeightPx);
        % Compute pixels-per-cm from actual window
        pxPerCmX = winWidthPx  / widthCm;
        pxPerCmY = winHeightPx / heightCm;
        % çree→pixel conversion lambdas
        deg2pxX = @(deg) viewingDistCm * tand(deg) * pxPerCmX;
        deg2pxY = @(deg) viewingDistCm * tand(deg) * pxPerCmY;
        dotRadiusPx = viewingDistCm * tand(outerDotRadiusDeg) * pxPerCmX;
        innerDotRadiusPx = viewingDistCm * tand(innerDotRadiusDeg) * pxPerCmX;
        % Define 13 calibration positions in degrees [xDeg, yDeg]
        % degPositions = [ ...
        %     0,  0;  -20, 20;   -20, -20;   20, 20;   20, -20; ...
        %     0, 20;   0, -20;   -20,   0;   20, 0; ...
        %     -15, 15;  15, 15;   -15, -15;   15, -15; ...
        %     -10,  10;   -10, -10;    10,  10;    10, -10; ...
        %     0,  10;    0, -10;   -10,   0;     10,   0; ...
        %     -5,   5;   5,   5;   -5,  -5;     5,  -5];
        degPositions = [ ...
            0, 0; -15, 15; -15, -15; 15, 15; 15, -15; ...
            0, 15; 0, -15; -15, 0; 15, 0;...
            -7.5, 7.5; -7.5, -7.5; 7.5, 7.5; 7.5, -7.5; ...
            0, 10; 0, -7.5; -7.5, 0; 7.5, 0];
        nDots = size(degPositions, 1);
        % Check CM distance
        xCm = viewingDistCm * tand(degPositions(:,1));
        yCm = viewingDistCm * tand(degPositions(:,2));
        rCm = hypot(xCm, yCm);
        fprintf(' Dot |  x°   |  y°   |   x_cm   |   y_cm   |  r_cm\n');
        fprintf('-----+-------+-------+----------+----------+--------\n');
        for iDot = 1:nDots
        fprintf('%4d | %5.1f° | %5.1f° | %7.2f  | %7.2f  | %6.2f\n', ...
            iDot, degPositions(iDot,1), degPositions(iDot,2), xCm(iDot), yCm(iDot), rCm(iDot));
        end
        fprintf('\n');
        % Convert to pixel coordinates
        positions = zeros(nDots, 2);
        for iDot = 1:nDots
            xOff = deg2pxX(degPositions(iDot,1));
            yOff = deg2pxY(degPositions(iDot,2));
            positions(iDot,:) = [xCen + xOff, yCen - yOff];
        end
        % Instruction text
        text = [
            'You will see dots appear one by one on the screen,' newline ...
            'each signaled by a brief beep. Please fix your gaze on each dot when it appears,' newline ...
            'and do not try to anticipate the location of the next dot.' newline newline ...
            'Press any key to begin and close your eyes until you hear the computer beep.'
        ];
        Screen('TextSize', win, 24);
        DrawFormattedText(win, text, 'center', 'center', outerDotColor);
        
        % Unify key names across platforms
        KbName('UnifyKeyNames');
        KbReleaseWait;
        FlushEvents('keyDown');
        Screen('Flip', win);
        
        % Wait for key press
        KbWait(-1);
    end
    
    % First, import the bluetooth library for light logger communication 
    bluetooth_central = import_pyfile(getpref("lightLoggerAnalysis", "bluetooth_central_path")); 
    % If a key has been pressed and we are not in simulation mode, 
    % start recording on the light logger 
    fileNameLL = [subjectId, experiment_name, num2str(session)];
    if(simulation_mode == "full" || simulation_mode == "bluetooth")
        disp("Main | Starting recording on light logger..."); 
        % We will attempt to start the light logger recording 
        success = start_recording_light_logger(bluetooth_central, fileNameLL, device_num);
        if(~success)
            Screen('CloseAll'); 
            error('Error starting light logger recording'); 
        end 
        % If the light logger was properly started, we need to leave some time for the AGC to converge 
        % to the appropriate settings 
        fprintf("Main | Recording started. Waiting %f seconds for AGC convergence...\n", agc_convergence_wait_s);
        start_rec_time = GetSecs;
        pause(agc_convergence_wait_s);
    else
        start_rec_time = GetSecs;
    end
    
    % Initialize task timing variables regardless of simulation mode
    totalDots = nDots * repetitions;
    % Array to store the *precise* cputime when each dot appears on screen
    dotStartTimes = zeros(totalDots, 1);
    % Array to store the measured duration of each dot
    actualDurations = zeros(totalDots, 1);
    
    if(simulation_mode == "full" || simulation_mode == "visual")
        
        disp('AGC done converging. Task will begin in about 10 seconds!')
        WaitSecs(10);
        HideCursor;
        Priority(MaxPriority(win));
        
        % Get refresh rate
        ifi = Screen('GetFlipInterval', win); % Inter-frame interval
        
        % --- Pre-calculate all rectangle positions before the loop ---
        outerRects = cell(nDots, 1);
        innerRects = cell(nDots, 1);
        for iDot = 1:nDots
            pos = positions(iDot,:);
            outerRects{iDot} = [pos(1)-dotRadiusPx, pos(2)-dotRadiusPx, ...
                pos(1)+dotRadiusPx, pos(2)+dotRadiusPx];
            innerRects{iDot} = [pos(1)-innerDotRadiusPx, pos(2)-innerDotRadiusPx, ...
                pos(1)+innerDotRadiusPx, pos(2)+innerDotRadiusPx];
        end
        
        % Play beep (single repetition)
        PsychPortAudio('Start', pahandle, 1, 0, 0);
        % Flash the screen red so we know when the task starts ont he world
        % camera
        Screen('FillRect', win, redColor, winRect);

        % Flip to display the single-frame marker. The current time is the start of this frame.
        Screen('Flip', win);

        % Now, immediately clear the back buffer (ready for the first dot)
        Screen('FillRect', win, bgColor, winRect);

        % Get initial time to base all future timings on
        % This is the time when the task officially begins (just before the first dot flip)
        start_task_time = GetSecs + 1; %add one sec because we will do some math below and we don't want to miss the first dot
        task_start_delay = start_task_time - start_rec_time;
        fprintf("Task start delay is %f seconds.\n", task_start_delay);
        % --- Timing setup using a CPU-based wait loop as requested ---
        
        % The total duration a dot is on screen before the next one appears
        dotTotalTime = dotTime; 
        
        prevDotStartTime = cputime; % Initialize with cputime before the loop
        
        % Draw loop with beep on each dot onset
        for rep = 1:repetitions
            for iDot = 1:nDots
                linear_index = (rep - 1) * nDots + iDot;
                
                % Draw stimulus
                Screen('FillOval', win, outerDotColor, outerRects{iDot});
                Screen('FillOval', win, innerDotColor, innerRects{iDot});
                
                % Flip to display the dot and get the precise ONSET time
                Screen('Flip', win);
                
                % Store the actual start time
                dotStartTimes(linear_index) = cputime;
                
                % Calculate and store the actual duration of the *previous* dot
                if linear_index > 1
                    % Duration of previous dot = Current Dot Start Time - Previous Dot Start Time
                    actualDurations(linear_index - 1) = dotStartTimes(linear_index) - dotStartTimes(linear_index - 1);
                end
                
                % Wait for the total dot duration before drawing the next one
                dotDoneTime = dotStartTimes(linear_index) + dotTotalTime;
                waitUntil(dotDoneTime);
                
                % Play beep for the *next* dot's onset (if not the last one)
                if linear_index < totalDots
                    PsychPortAudio('Start', pahandle, 1, 0, 0);
                end
            end
        end

        % Flash the screen red so we know when the task ends on the world
        % camera
        Screen('FillRect', win, redColor, winRect);
        
        % Final flip to clear the screen
        Screen('Flip', win);
        
    
        actualDurations(end) = dotStartTimes(end) + dotTotalTime - dotStartTimes(end); % This will be exactly dotTotalTime
        
        % A more practical measure that accounts for any slight overrun in waitUntil 
        % or cleanup:
        timeAfterLastWait = cputime; % Time right after the final waitUntil call
        actualDurations(end) = timeAfterLastWait - dotStartTimes(end);
        
        disp('Actual dot durations:')
        disp(actualDurations)
        
        Priority(0);
        ShowCursor;
        Screen('CloseAll');
    end
    
    % --- SAVING DATA ---
    % 1. Compile all task timing data
    taskData.experiment_name = experiment_name;
    taskData.task_start_delay_s = task_start_delay; % Time camera started to first dot flip
    taskData.dot_durations_s = actualDurations; % Measured duration of each dot (totalDots x 1)
    taskData.dot_nominal_duration_s = dotTime; % The intended duration
    taskData.gaze_target_positions_deg = degPositions; % Target positions (nDots x 2)
    
    % 2. Create unique filename
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    folders = ['/FLIC_data/lightLogger/GazeCalRunFileData/', subjectId];
    subjDir = fullfile(getpref("lightLoggerAnalysis", 'dropboxBaseDir'), folders);
    if ~exist(subjDir, 'dir')
        mkdir(subjDir)
    end
    filename = fullfile(subjDir, sprintf('%s_%s_session%s_%s.mat', subjectId, experiment_name, num2str(session),timestamp));

    % 3. Save the data structure
    try
        save(filename, 'taskData', '-v7.3');
        fprintf('\nData saved to: %s\n', filename);
    catch ME
        warning('Failed to save task data: %s', ME.message);
    end
    
    % If we are not in simulation mode, stop recording from the light logger 
    if(simulation_mode == "full" || simulation_mode == "bluetooth")
        success = stop_recording_light_logger(bluetooth_central, device_num);
        if(~success)
            error("Error stopping light logger recording"); 
        end 
    end
end

% Local function to perform a busy-wait until a specified time
function waitUntil(stopTimeSeconds)
    while cputime() < stopTimeSeconds
        % Busy-wait loop
    end
end
% Local function to start recording on the light logger 
function success = start_recording_light_logger(bluetooth_central,... 
                                                experiment_name,...
                                                device_num...
                                               ) 
% [No change to this local function]
    % Initialize a success variable 
    success = 0; 
    % Generate a message to send to the light logger 
    disp("Main | Generating recording message...")
    light_logger_update_message = generate_light_logger_recording_message(bluetooth_central, experiment_name); 
    % Send a recording message to the light logger
    disp("Main | Sending recording message...") 
    try
        bluetooth_central.message_peripheral_matlab_wrapper(device_num,... 
                                                            light_logger_update_message...
                                                           );
    catch ME 
        print_error_info(ME); 
        return;
    end  
    % If a message was sent successfully, 
    % wait for the light logger to enter the recording state 
    while(true)
        % Read the current state from the light logger 
        disp("Main | Waiting for acknowledgement...")
        try
            lightlogger_state = struct(bluetooth_central.read_peripheral_matlab_wrapper(device_num));
        catch ME
            print_error_info(ME);
            return; 
        end 
        
        state_name = string(char(lightlogger_state.state));     
        % If the state is error, raise an error on this machine 
        % as well 
        if(state_name == "error") 
            % If the light logger has gone into an errored state, 
            % return failure  
            fprintf("ERROR: Unknown error on the light logger\n"); 
            return ; 
        end 
        % If the light logger entered the target state, stop waiting. 
        if(state_name == "science")
            break 
        end
        
        % Pause for some time between reads
        pause(2); 
    end 
    disp("Main | Acknowledged.")
    % If we successfully executed, denote it as such 
    success = 1; 
    return; 
end 
% Local function to generate the recording message for the light logger 
function message = generate_light_logger_recording_message(bluetooth_central, experiment_name)
% [No change to this local function]
    % Then, initialize a struct that will be sent to the light logger 
    message = bluetooth_central.initialize_update_message();
    
    % Import world_util to retrieve the initial settings for the experiment 
    world_util = import_pyfile(getpref("lightLoggerAnalysis", "world_util_path")); 
    % Retrieve the sensor mode (FPS, size)
    sensor_mode = world_util.WORLD_CAMERA_CUSTOM_MODES{1}; 
    % Initialize the sensors that will be needed for this recording
    % Retrieve a low ND filter initial settings for the world camera 
    % so that it starts closer to the convergence target 
    initial_settings = double(world_util.WORLD_NDF_LEVEL_SETTINGS{3});
    
    sensors.W.Again = initial_settings(1); 
    sensors.W.Dgain = initial_settings(2); 
    sensors.W.exposure = py.int(initial_settings(3));  
    sensors.W.agc = true; 
    sensors.W.save_agc_metadata = true; 
    sensors.W.sensor_mode = sensor_mode;
    sensors.W.awb = false; 
    sensors.W.noise_mode = false; 
    sensors.P.agc = false; 
    sensors.P.save_agc_metadata = false; 
    % Next, we will enter a science state to record until the stimulus is done 
    bluetooth_central.generate_science_state(message,...
                                             experiment_name,...
                                             py.int(30),...
                                             false,...  
                                             sensors...
                                            ); 
end 
% Local function to end light logger recording 
function success = stop_recording_light_logger(bluetooth_central, device_num) 
% [No change to this local function]
    % Initialize a success variable 
    success = 0; 
    % Generate a message to send to the light logger 
    disp("Main | Generating stop message...")
    light_logger_update_message = generate_light_logger_stop_message(bluetooth_central); 
    % Send a recording message to the light logger 
    disp("Main | Sending stop message...") 
    try
        bluetooth_central.message_peripheral_matlab_wrapper(device_num,... 
                                                            light_logger_update_message...
                                                           );
    catch ME 
        print_error_info(ME);
        return; 
    end 
    % If a message was sent successfully, 
    % wait for the light logger to enter the recording state 
    while(true)
        disp("Main | Waiting for acknowledgement...")
        % Read the current state from the light logger 
        try
            lightlogger_state = struct(bluetooth_central.read_peripheral_matlab_wrapper(device_num));
        catch ME
            print_error_info(ME);
            return; 
        end 
        
        state_name = string(char(lightlogger_state.state));     
        % If the state is error, raise an error on this machine 
        % as well 
        if(state_name == "error") 
            % If the light logger has gone into an errored state, 
            % return failure  
            fprintf("ERROR: Unknown error on the light logger\n"); 
            return ; 
        end 
        % If the light logger entered the target state, stop waiting 
        if(state_name == "wait")
            break 
        end
        
        % Pause for some time between reads
        pause(2); 
    end 
    disp("Main | Acknowledged.")
    % If we successfully executed, denote it as such 
    success = 1; 
    return; 
end 
% Local function to generate the recording message for the light logger 
function message = generate_light_logger_stop_message(bluetooth_central)
% [No change to this local function]
    % Then, initialize a struct that will be sent to the light logger 
    message = bluetooth_central.initialize_update_message();
    
    % Next, we will enter a wait to end the science state
    bluetooth_central.generate_wait_state(message); 
end 
% Local function to print out error information from a caught error 
function print_error_info(ME)
% [No change to this local function]
    fprintf('Error message: %s\n', ME.message);
    fprintf('Error occurred in:\n');
    for k = 1:length(ME.stack)
        fprintf('  File: %s\n  Function: %s\n  Line: %d\n\n', ...
                ME.stack(k).file, ME.stack(k).name, ME.stack(k).line);
    end
end