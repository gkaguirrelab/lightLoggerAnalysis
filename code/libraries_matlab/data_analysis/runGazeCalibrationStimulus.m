function degPositions = runGazeCalibrationStimulus(simulation_mode, device_num, agc_convergence_wait_s, experiment_name, heightMm, widthMm)
% Displays 26-dot gaze calibration stimulus at fixed visual angles, with a brief beep signaling each dot onset.
% 
% Description: 
%   xxxx
% 
% Inputs:
%   simulation_mode             - Boolean. Whether or not to run the calibration only, wihout
%                                 light logger recording.
%   device_num                  - Numeric. The device number to use when recording from a real device.
%   agc_convergence_wait_s      - Numeric. The time in seconds to wait for the AGC to converge to appropriate settings.
%   experiment_name             - String. The filename the recording will be saved under on the light logger. 
%   heightMm                    - Numeric. Display screen height in mm.
%   widthMm                     - Numeric. Display screen width in mm.
% 
% Example:
%{
    runGazeCalibrationStimulus(true)
    runGazeCalibrationStimulus(false, 2, 60, 'GazeCalib_Run1', 1067, 1924)
%}
                          
    arguments 
        simulation_mode {mustBeNumericOrLogical} = true;
        device_num {mustBeNumeric} = 1;
        agc_convergence_wait_s {mustBeNumeric} = 60;
        experiment_name = "GazeCalibration";
        heightMm = 1067; % 2nd floor LGTV
        widthMm = 1924; % 2nd floor LGTV 
    end

    % Hard-coded parametersk
    viewingDistCm = 100;
    dotRadiusDeg = 0.8;
    dotTime = 1; 
    repetitions = 3;
    bgColor = [0 0 0];
    fgColor = [255 255 255];
    redColor  = [255   0   0];
    innerFrac = 0.3;

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
    if ~exist('widthMm', 'var')
        [widthMm, heightMm] = Screen('DisplaySize', screenNum);
        warning('Using DisplaySize, which may incorrectly estimate screen height and width.')
    end
    widthCm  = widthMm / 10;
    heightCm = heightMm / 10;

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
    dotRadiusPx = viewingDistCm * tand(dotRadiusDeg) * pxPerCmX;

    % Define 13 calibration positions in degrees [xDeg, yDeg]
    degPositions = [ ...
        0,  0;  -20, 20;   -20, -20;   20, 20;   20, -20; ...
        0, 20;   0, -20;   -20,   0;   20, 0; ...
        -15, 15;  15, 15;   -15, -15;   15, -15; ...

        -10,  10;   -10, -10;    10,  10;    10, -10; ...
        0,  10;    0, -10;   -10,   0;     10,   0; ...
        -5,   5;   5,   5;   -5,  -5;     5,  -5;     0,   0];
    nDots = size(degPositions, 1);

    % Check CM distance
    xCm = viewingDistCm * tand(degPositions(:,1));
    yCm = viewingDistCm * tand(degPositions(:,2));
    rCm = hypot(xCm, yCm);

    fprintf(' Dot |  x°   |  y°   |   x_cm   |   y_cm   |  r_cm\n');
    fprintf('-----+-------+-------+----------+----------+--------\n');
    for i = 1:nDots
    fprintf('%4d | %5.1f° | %5.1f° | %7.2f  | %7.2f  | %6.2f\n', ...
        i, degPositions(i,1), degPositions(i,2), xCm(i), yCm(i), rCm(i));
    end
    fprintf('\n');

    % Convert to pixel coordinates
    positions = zeros(nDots, 2);
    for i = 1:nDots
        xOff = deg2pxX(degPositions(i,1));
        yOff = deg2pxY(degPositions(i,2));
        positions(i,:) = [xCen + xOff, yCen - yOff];
    end

    % Instruction text
    text = [
        'You will see 26 calibration dots appear one by one on the screen,' newline ...
        'each signaled by a brief beep. Please fix your gaze on each dot when it appears,' newline ...
        'and do not try to anticipate the location of the next dot.' newline newline ...
        'Press any key to begin.'
    ];
    Screen('TextSize', win, 24);
    DrawFormattedText(win, text, 'center', 'center', fgColor);
    
    % Unify key names across platforms
    KbName('UnifyKeyNames');
    KbReleaseWait;
    FlushEvents('keyDown');
    Screen('Flip', win);
    
    % Wait for key press
    KbWait(-1);

    % First, import the bluetooth library for light logger communication 
    bluetooth_central = import_pyfile(getpref("lightLoggerAnalysis", "bluetooth_central_path")); 

    % If a key has been pressed and we are not in simulation mode, 
    % start recording on the light logger 
    if(~simulation_mode)
        disp("Main | Starting recording on light logger..."); 

        % We will attempt to start the light logger recording 
        success = start_recording_light_logger(bluetooth_central, experiment_name, device_num);
        if(~success)
            Screen('CloseAll'); 
            error('Error starting light logger recording'); 
        end 

        % If the light logger was properly started, we need to leave some time for the AGC to converge 
        % to the appropriate settings 
        fprintf("Main | Recording started. Waiting %f seconds for AGC convergence...\n", agc_convergence_wait_s); 
        pause(agc_convergence_wait_s); 

    end

    disp('AGC done converging. Press any key to continue.')
    KbWait(-1);

    % Draw loop with beep on each dot onset
    HideCursor;
    Priority(MaxPriority(win));
    for rep = 1:repetitions
        for i = 1:nDots
            % Play beep (single repetition)
            PsychPortAudio('Start', pahandle, 1, 0, 0);
            
            % Draw dot
            pos  = positions(i,:);
            % Outer (white) circle
            outerRect = [pos(1)-dotRadiusPx, pos(2)-dotRadiusPx, ...
                         pos(1)+dotRadiusPx, pos(2)+dotRadiusPx];
            Screen('FillOval', win, fgColor, outerRect);
            
            % Inner (red) circle
            innerRadiusPx = dotRadiusPx * innerFrac;
            innerRect = [pos(1)-innerRadiusPx, pos(2)-innerRadiusPx, ...
                         pos(1)+innerRadiusPx, pos(2)+innerRadiusPx];
            Screen('FillOval', win, redColor, innerRect);
            WaitSecs(dotTime);
    
            Screen('Flip', win);
            WaitSecs(dotTime);
        end
        Screen('Flip', win);
        WaitSecs(dotTime);
    end
    Priority(0);
    ShowCursor;
    Screen('CloseAll');

    % If we are not in simulation mode, stop recording from the light logger 
    if(~simulation_mode)
        success = stop_recording_light_logger(bluetooth_central, device_num);
        if(~success)
            error("Error stopping light logger recording"); 
        end 
    end

end

% Local function to start recording on the light logger 
function success = start_recording_light_logger(bluetooth_central,... 
                                                experiment_name,...
                                                device_num...
                                               ) 

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
    % Then, initialize a struct that will be sent to the light logger 
    message = bluetooth_central.initialize_update_message();
    
    % Import world_util to retrieve the initial settings for the experiment 
    world_util = import_pyfile(getpref("lightLoggerAnalysis", "world_util_path")); 

    % Initialize the sensors that will be needed for this recording
    % Retrieve a low ND filter initial settings for the world camera 
    % so that it starts closer to the convergence target 
    initial_settings = double(world_util.WORLD_NDF_LEVEL_SETTINGS{3});
    

    sensors.W.Again = initial_settings(1); 
    sensors.W.Dgain = initial_settings(2); 
    sensors.W.exposure = py.int(initial_settings(3));  
    sensors.W.agc = true; 
    sensors.W.save_agc_metadata = true; 

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
    % Then, initialize a struct that will be sent to the light logger 
    message = bluetooth_central.initialize_update_message();
    
    % Next, we will enter a wait to end the science state
    bluetooth_central.generate_wait_state(message); 

end 

% Local function to print out error information from a caught error 
function print_error_info(ME)
    fprintf('Error message: %s\n', ME.message);
    fprintf('Error occurred in:\n');
    for k = 1:length(ME.stack)
        fprintf('  File: %s\n  Function: %s\n  Line: %d\n\n', ...
                ME.stack(k).file, ME.stack(k).name, ME.stack(k).line);
    end

end 