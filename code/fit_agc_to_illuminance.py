def fit_agc_to_illuminance(video_paths, recording_paths):
    """
    Fit camera AGC settings to environmental illuminance.

    GOAL
    ----
    Derive a calibration function relating world-camera AGC settings
    (analog gain, digital gain, and exposure) to environmental illuminance
    measured simultaneously by the minispectrometer.

    For each supplied environmental recording:

    1. Load world-camera AGC metadata.
    2. Calculate a camera light score:

        score = 1 / (Again * Dgain * exposure)

    3. Load the corresponding minispectrometer data.
    4. Convert minispectrometer measurements to illuminance.
    5. Align camera and minispectrometer timebases.
    6. Estimate the temporal lag of the AGC response by finding the
       time shift that maximizes correlation between camera score
       and illuminance.
    7. Correct for this lag.
    8. Pool valid observations across recordings.

    Finally, fit:

        illuminance = f(camera_score)

    and return the fitted calibration parameters.

    The function accepts an iterable of video paths so individual
    recordings can easily be included or excluded from calibration.

    %% MAPPING
    video_paths = [
        "/path/to/outdoorDemo.avi"
    ]
    recording_paths = [
        "/Users/sophiamirabal/Library/CloudStorage/"
        "Dropbox-Aguirre-BrainardLab/Team Documents/temp_for_sharing/"
        "outdoorDemoSophia"
    ]
    """

    all_camera_scores = []
    all_illuminance = []

    for video_path, recording_path in zip(video_paths, recording_paths):

        # Load world-camera AGC metadata  
        metadata = extract_metadata(video_path)

        camera_time = metadata["time"]
        analog_gain = metadata["Again"]
        digital_gain = metadata["Dgain"]
        exposure = metadata["exposure"]

        # Convert AGC state into one camera "light score"
        camera_score = 1 / (analog_gain * digital_gain * exposure)

        # Load ONLY minispect data from corresponding recording
        ms_chunks = parse_chunks(recording_path, convert_time_units=True, convert_to_float=True, chunk_ranges={"M": (0, None)})

        # TODO: Combine minispect chunks
        # 
        # The parser returns:
        #
        # chunk["M"]["t"]
        # chunk["M"]["v"]
        #
        # Need to determine which minispect sensor/channel is used
        # to calculate illuminance.
        
        ms_time = ...
        ms_values = ...

        # TODO: Obtain minispect illuminance
        illuminance = ...

        # TODO: Put camera score and illuminance on comparable timebases
        camera_time, camera_score, ms_time, illuminance = \
            synchronize_signals(...)

        # TODO: Estimate AGC lag
        best_lag = lag_that_maximizes_correlation(
            camera_score,
            illuminance
        )

        # TODO: Shift signals according to best lag
        aligned_camera_score = shift(camera_score, best_lag)

        # TODO: Save valid paired observations from this recording
        all_camera_scores.append(aligned_camera_score)
        all_illuminance.append(illuminance)

    # Combine recordings
    camera_score_all = np.concatenate(all_camera_scores)
    illuminance_all = np.concatenate(all_illuminance)

    # Fit the final calibration equation
    #
    # illuminance ~= f(camera_score)
    #
    calibration_parameters = fit_model(
        camera_score_all,
        illuminance_all
    )

    # Return everything useful for validation / later application
    return {
        "parameters": calibration_parameters,
        "camera_score": camera_score_all,
        "illuminance": illuminance_all,
        "best_lag": best_lag,
        ...
    }