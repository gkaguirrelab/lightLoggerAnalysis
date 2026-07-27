def fit_agc_to_illuminance(video_paths):
    """
    Derive a calibration relating camera AGC settings
    to environmental illuminance measured by the minispect.
    """

    all_camera_scores = []
    all_illuminance = []

    for video_path in video_paths:

        # 1. Identify the recording associated with this video
        recording_path = find_recording_from_video(video_path)

        # 2. Load world-camera metadata
        # Need:
        #   timestamps
        #   analog gain
        #   digital gain
        #   exposure
        metadata = extract_metadata(video_path)

        camera_time = metadata["time"]
        analog_gain = metadata["Again"]
        digital_gain = metadata["Dgain"]
        exposure = metadata["exposure"]

        # 3. Convert AGC state into one camera "light score"
        camera_score = 1 / (
            analog_gain *
            digital_gain *
            exposure
        )

        # 4. Load ONLY minispect data from corresponding recording
        ms = parse_chunks(
            recording_path,
            sensors="M",
            ...
        )

        # 5. Extract minispect time and convert its measurements
        # into environmental illuminance
        ms_time = ...
        illuminance = ...

        # 6. Put camera score and illuminance on comparable timebases
        camera_time, camera_score, ms_time, illuminance = \
            synchronize_signals(...)

        # 7. Estimate AGC lag
        # Shift camera score relative to minispect illuminance
        # over a range of candidate temporal offsets.
        lags = ...
        correlations = ...

        best_lag = lag_that_maximizes_correlation(
            camera_score,
            illuminance
        )

        # 8. Shift signals according to best lag
        aligned_camera_score = shift(camera_score, best_lag)

        # 9. Save valid paired observations from this recording
        all_camera_scores.append(aligned_camera_score)
        all_illuminance.append(illuminance)

    # 10. Combine recordings
    camera_score_all = concatenate(all_camera_scores)
    illuminance_all = concatenate(all_illuminance)

    # 11. Fit the final calibration equation
    #
    # illuminance ~= f(camera_score)
    #
    calibration_parameters = fit_model(
        camera_score_all,
        illuminance_all
    )

    # 12. Return everything useful for validation / later application
    return {
        "parameters": calibration_parameters,
        "camera_score": camera_score_all,
        "illuminance": illuminance_all,
        "best_lag": best_lag,
        ...
    }