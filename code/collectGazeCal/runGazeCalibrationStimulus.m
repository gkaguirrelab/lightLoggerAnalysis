function degPositions = runGazeCalibrationStimulus(simulation_mode, device_num, agc_convergence_wait_s, subjectId, experiment_name, session, heightCm, widthCm, viewingDistCm)
% Present the timed gaze-calibration stimulus and optionally record it.
%
% Syntax:
%   degPositions = runGazeCalibrationStimulus(simulation_mode, device_num, agc_convergence_wait_s, subjectId, experiment_name, session, heightCm, widthCm, viewingDistCm)
%
% Description:
%   This function coordinates the visual gaze-calibration task used with
%   the light logger recordings. In visual modes it opens a Psychtoolbox
%   window, displays instructions, converts a predefined set of target
%   locations from degrees to pixels, and then presents the targets one at
%   a time with audio beeps and red synchronization flashes at task onset
%   and offset. In bluetooth-enabled modes it also commands the light
%   logger to enter and leave its science state, waits for AGC convergence
%   before beginning the task, and saves the measured task timing metadata
%   alongside the target positions for downstream synchronization.
%
% Inputs:
%   simulation_mode          - String. Execution mode: "full" runs both
%                              stimulus presentation and bluetooth control,
%                              "visual" runs the on-screen task only, and
%                              "bluetooth" controls the device without
%                              presenting the task.
%   device_num               - Scalar. Index of the bluetooth peripheral to
%                              command through the Python control layer.
%   agc_convergence_wait_s   - Scalar. Number of seconds to wait after the
%                              recording enters science mode so the camera
%                              AGC can settle before the first target.
%   subjectId                - String or char vector. Subject identifier
%                              used in both the recording name and the
%                              saved MATLAB task data.
%   experiment_name          - String or char vector. Label used when
%                              generating the light logger recording state.
%   session                  - Scalar. Session number appended to the light
%                              logger recording name.
%   heightCm                 - Scalar. Physical display height in
%                              centimeters.
%   widthCm                  - Scalar. Physical display width in
%                              centimeters.
%   viewingDistCm            - Scalar. Viewing distance from the observer
%                              to the display in centimeters.
%
% Outputs:
%   degPositions             - Numeric matrix. Ordered [xDeg, yDeg]
%                              calibration target positions used during the
%                              stimulus presentation.
%
% Examples:
%{
    degPositions = runGazeCalibrationStimulus("visual", 2, 60, ...
        "FLIC_2004", "GazeCalibration", 1, 106.7, 192.4, 100);
%}
%{
    subjectId = 'FLIC_2004';
    sessionNum = 1;
    runGazeCalibrationStimulus("visual", 2, 60, subjectId, ...
        'GazeCalibration', sessionNum, 106.7, 192.4, 100);
%}
                          
    arguments 
        simulation_mode {mustBeMember(simulation_mode, ["full", "visual", "bluetooth"])} = "full";
        device_num {mustBeNumeric} = 1;
        agc_convergence_wait_s {mustBeNumeric} = 60;
        subjectId = "test";
        experiment_name = "gazeCalibration";
        session = 1;
        heightCm = 106.7; % 2nd floor LGTV
        widthCm = 192.4; % 2nd floor LGTV
        viewingDistCm = 100; % 1 meter standard
    end
    % Hard-coded parameters
   
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
            0, 7.5; 0, -7.5; -7.5, 0; 7.5, 0];
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
    if(simulation_mode == "full" || simulation_mode == "bluetooth")
        bluetooth_central = import_pyfile(getpref("lightLoggerAnalysis", "bluetooth_central_path"));
    end
    % If a key has been pressed and we are not in simulation mode,
    % start recording on the light logger
    fileNameLL = [subjectId,'_', experiment_name,'_tf_session' num2str(session)];
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
    folders = ['/FLIC_data/lightLogger/scriptedIndoorOutdoor/', subjectId, 'gazeCalibration/temporalFrequency'];
    subjDir = fullfile(getpref("lightLoggerAnalysis", 'dropboxBaseDir'), folders);
    if ~exist(subjDir, 'dir')
        mkdir(subjDir)
    end
    filename = fullfile(subjDir, sprintf('%s_gazeCal_runData.mat', subjectId));

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
% Busy-wait until a target CPU time is reached.
%
% Syntax:
%   waitUntil(stopTimeSeconds)
%
% Description:
%   This local helper spins until MATLAB's CPU timer reaches the supplied
%   stop time. It is used inside the stimulus loop to keep the nominal
%   dwell time of each calibration target close to the desired duration
%   without inserting additional scheduling logic into the presentation
%   loop itself.
%
% Inputs:
%   stopTimeSeconds          - Scalar. Target value of `cputime` at which
%                              the function should return.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See runGazeCalibrationStimulus.m for usage context.
%}

    while cputime() < stopTimeSeconds
        % Busy-wait loop
    end
end
% Local function to start recording on the light logger 
function success = start_recording_light_logger(bluetooth_central,... 
                                                experiment_name,...
                                                device_num...
                                               ) 
% Request that the light logger enter its science recording state.
%
% Syntax:
%   success = start_recording_light_logger(bluetooth_central, experiment_name, device_num)
%
% Description:
%   This helper assembles the appropriate update message, sends it to the
%   selected peripheral, and then repeatedly polls the device state until
%   the light logger reports that it has entered the "science" state. It
%   returns false immediately if message transmission fails or if the
%   peripheral reports an error state while the acknowledgement loop is
%   running.
%
% Inputs:
%   bluetooth_central        - Python module or object. Interface that
%                              exposes the bluetooth messaging wrappers.
%   experiment_name          - String or char vector. Recording name to
%                              embed in the state transition request.
%   device_num               - Scalar. Identifier of the peripheral to
%                              message.
%
% Outputs:
%   success                  - Logical scalar. True if the device confirms
%                              entry into science mode; false otherwise.
%
% Examples:
%{
    % See runGazeCalibrationStimulus.m for usage context.
%}

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
% Build the light logger state-update message for task recording.
%
% Syntax:
%   message = generate_light_logger_recording_message(bluetooth_central, experiment_name)
%
% Description:
%   This helper creates the Python-side update message that places the
%   device into science mode for the gaze-calibration run. It imports the
%   world camera utility module to pull a default sensor mode and a
%   preselected set of world-camera gain and exposure settings so the
%   logger begins recording near a reasonable operating point before AGC
%   convergence.
%
% Inputs:
%   bluetooth_central        - Python module or object used to create and
%                              populate the update message.
%   experiment_name          - String or char vector. Recording label to
%                              attach to the generated science state.
%
% Outputs:
%   message                  - Python message object ready to be sent to
%                              the peripheral.
%
% Examples:
%{
    % See runGazeCalibrationStimulus.m for usage context.
%}

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
% Request that the light logger leave its science recording state.
%
% Syntax:
%   success = stop_recording_light_logger(bluetooth_central, device_num)
%
% Description:
%   This helper sends a wait-state transition message to the peripheral
%   and then polls until the light logger reports that it has returned to
%   the "wait" state. The logic mirrors the start-recording handshake so
%   the calling code can verify that the recording was shut down cleanly.
%
% Inputs:
%   bluetooth_central        - Python module or object. Interface that
%                              exposes the bluetooth messaging wrappers.
%   device_num               - Scalar. Identifier of the peripheral to
%                              message.
%
% Outputs:
%   success                  - Logical scalar. True if the device confirms
%                              re-entry into wait mode; false otherwise.
%
% Examples:
%{
    % See runGazeCalibrationStimulus.m for usage context.
%}

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
% Build the light logger update message that ends recording.
%
% Syntax:
%   message = generate_light_logger_stop_message(bluetooth_central)
%
% Description:
%   This helper creates a Python-side update message and populates it with
%   the state transition that returns the light logger to its idle wait
%   state after the calibration task has completed.
%
% Inputs:
%   bluetooth_central        - Python module or object used to create and
%                              populate the update message.
%
% Outputs:
%   message                  - Python message object ready to be sent to
%                              the peripheral.
%
% Examples:
%{
    % See runGazeCalibrationStimulus.m for usage context.
%}

    message = bluetooth_central.initialize_update_message();
    
    % Next, we will enter a wait to end the science state
    bluetooth_central.generate_wait_state(message); 
end 
% Local function to print out error information from a caught error 
function print_error_info(ME)
% Print a readable summary of a caught MATLAB exception.
%
% Syntax:
%   print_error_info(ME)
%
% Description:
%   This helper prints the exception message and each frame in the stack so
%   bluetooth communication failures can be inspected from the MATLAB
%   console without rethrowing the original error.
%
% Inputs:
%   ME                       - `MException` object caught by a surrounding
%                              `try/catch` block.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See runGazeCalibrationStimulus.m for usage context.
%}

    fprintf('Error message: %s\n', ME.message);
    fprintf('Error occurred in:\n');
    for k = 1:length(ME.stack)
        fprintf('  File: %s\n  Function: %s\n  Line: %d\n\n', ...
                ME.stack(k).file, ME.stack(k).name, ME.stack(k).line);
    end
end
