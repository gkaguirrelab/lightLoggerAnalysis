import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Pi_util / parse_chunks
sys.path.append("/Users/sophiamirabal/Documents/MATLAB/projects/" "lightLogger/raspberry_pi_firmware/utility")

# world_util
sys.path.append("/Users/sophiamirabal/Documents/MATLAB/projects/" "lightLoggerAnalysis/code/library/sensor_utility")

from Pi_util import parse_chunks
import world_util

### HELPER FUNCTION 
def zscore_signal(x):
    """
    Normalize a signal to mean 0 and standard deviation 1. This allows us to compare temporal shape even though the camera
    score and minispect channel have different units / scales.
    """
    x = np.asarray(x, dtype=float)

    mean_x = np.nanmean(x)
    std_x = np.nanstd(x)

    if std_x == 0:
        raise ValueError("Cannot z-score a signal with zero variance.")

    return (x - mean_x) / std_x

### HELPER FUNCTION 
def lowpass_first_order(time, signal, tau):
    """
    Apply a first-order low-pass filter to an irregularly sampled signal.
    """
    time = np.asarray(time, dtype=float)
    signal = np.asarray(signal, dtype=float)

    filtered = np.full(signal.shape, np.nan, dtype=float)

    valid = np.isfinite(time) & np.isfinite(signal)

    if np.sum(valid) < 2:
        return filtered

    t = time[valid]
    x = signal[valid]

    y = np.empty_like(x)

    # Initialize filter at first observed value
    y[0] = x[0]

    for ii in range(1, len(x)):
        dt = t[ii] - t[ii - 1]

        if dt <= 0:
            y[ii] = y[ii - 1]
            continue

        # Exact update for first-order exponential low-pass
        alpha = 1.0 - np.exp(-dt / tau)
        y[ii] = (y[ii - 1] + alpha * (x[ii] - y[ii - 1]))

    filtered[valid] = y
    return filtered

### HELPER FUNCTION 
def characterize_lowpass_filter(camera_time, camera_score, ms_time, ms_signal, 
                                lag_range=(-10, 10), tau_range=(0.1, 20),
                                n_lags=201,n_taus=100):
    """
    Find the first-order low-pass time constant and temporal lag that
    best relate a minispect channel to the camera light score.

    The model is: camera_score ~= delayed lowpass(minispect_signal)

    The comparison is performed after z-scoring both signals so their
    different physical units do not affect the temporal comparison.
    """

    camera_time = np.asarray(camera_time, dtype=float)
    camera_score = np.asarray(camera_score, dtype=float)

    ms_time = np.asarray(ms_time, dtype=float)
    ms_signal = np.asarray(ms_signal, dtype=float)


    # REMOVE INVALID VALUES
    # -------------------------------------------------------------
    valid_camera = (np.isfinite(camera_time) & np.isfinite(camera_score))
    valid_ms = (np.isfinite(ms_time) & np.isfinite(ms_signal))

    camera_time = camera_time[valid_camera]
    camera_score = camera_score[valid_camera]

    ms_time = ms_time[valid_ms]
    ms_signal = ms_signal[valid_ms]

    # SORT BY TIME
    # -------------------------------------------------------------
    camera_order = np.argsort(camera_time)
    camera_time = camera_time[camera_order]
    camera_score = camera_score[camera_order]

    ms_order = np.argsort(ms_time)
    ms_time = ms_time[ms_order]
    ms_signal = ms_signal[ms_order]

    # Normalize because the units / magnitudes are different.
    camera_z = zscore_signal(camera_score)
    ms_z = zscore_signal(ms_signal)

    # PUT CAMERA SCORE ON THE MINISPECT TIMEBASE
    # -------------------------------------------------------------
    camera_on_ms_time = np.interp(ms_time, camera_time, camera_z, left=np.nan, right=np.nan)

    # CANDIDATE FILTER PARAMETERS
    # -------------------------------------------------------------
    candidate_lags = np.linspace(lag_range[0], lag_range[1], n_lags)
    candidate_taus = np.linspace(tau_range[0], tau_range[1], n_taus)

    correlations = np.full((len(candidate_taus), len(candidate_lags)), np.nan)

    # TEST EACH LOW-PASS TIME CONSTANT + LAG
    # -------------------------------------------------------------
    for tau_idx, tau in enumerate(candidate_taus):
        filtered_ms = lowpass_first_order(ms_time, ms_z, tau)

        for lag_idx, lag in enumerate(candidate_lags):

            # Positive lag means camera response occurs AFTER the minispect change.
            # Therefore, camera(t) should match minispect(t - lag).
            predicted_camera = np.interp(ms_time - lag, ms_time, filtered_ms, left=np.nan, right=np.nan)

            valid = (np.isfinite(camera_on_ms_time) & np.isfinite(predicted_camera))

            if np.sum(valid) < 3:
                continue

            correlations[tau_idx, lag_idx] = np.corrcoef(camera_on_ms_time[valid], predicted_camera[valid])[0, 1]

    # FIND BEST PARAMETERS
    # -------------------------------------------------------------
    best_flat_idx = np.nanargmax(correlations)

    best_tau_idx, best_lag_idx = np.unravel_index(best_flat_idx, correlations.shape)

    best_tau = candidate_taus[best_tau_idx]
    best_lag = candidate_lags[best_lag_idx]
    best_correlation = correlations[best_tau_idx, best_lag_idx]

    # GENERATE BEST-FIT PREDICTION
    # -------------------------------------------------------------
    best_filtered_ms = lowpass_first_order(ms_time, ms_z, best_tau)
    best_prediction = np.interp(ms_time - best_lag, ms_time, best_filtered_ms, left=np.nan, right=np.nan)

    # Correlation without any filtering or lag correction
    valid_raw = (np.isfinite(camera_on_ms_time) & np.isfinite(ms_z))
    raw_correlation = np.corrcoef(camera_on_ms_time[valid_raw], ms_z[valid_raw])[0, 1]

    return {
        "best_tau": best_tau,
        "best_lag": best_lag,
        "best_correlation": best_correlation,
        "raw_correlation": raw_correlation,

        "candidate_lags": candidate_lags,
        "candidate_taus": candidate_taus,
        "correlations": correlations,

        "ms_time": ms_time,
        "camera_on_ms_time": camera_on_ms_time,
        "ms_z": ms_z,
        "best_filtered_ms": best_filtered_ms,
        "best_prediction": best_prediction
    }

### MAIN
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
    4. Flatten the AS measurements across minispect chunks.
    5. Convert minispectrometer measurements to illuminance.
    6. Align camera and minispectrometer timebases.
    7. Estimate the temporal lag of the AGC response by finding the
       time shift that maximizes correlation between camera score
       and illuminance.
    8. Correct for this lag.
    9. Pool valid observations across recordings.

    Finally, fit:

        illuminance = f(camera_score)

    and return the fitted calibration parameters.

    The function accepts iterables of video paths and corresponding
    recording paths so individual recordings can easily be included
    or excluded from calibration.

    """

    if len(video_paths) != len(recording_paths):
        raise ValueError(
            "video_paths and recording_paths must have the same length."
        )

    all_camera_scores = []
    all_illuminance = []

    results_by_recording = []

    # LOOP THROUGH RECORDINGS
    # =================================================================

    for video_path, recording_path in zip(video_paths, recording_paths):

        print("\n======================================")
        print("PROCESSING RECORDING")
        print("======================================")
        print("Video path:", video_path)
        print("Raw chunk path:", recording_path)

        if not os.path.isdir(recording_path):
            raise FileNotFoundError(
                f"Recording directory does not exist:\n{recording_path}"
            )

        # CAMERA / AGC
        # =============================================================

        original_metadata_cols = world_util.WORLD_AGC_METADATA_COLS

        try:
            world_util.WORLD_AGC_METADATA_COLS = (
                "Again",
                "Dgain",
                "exposure"
            )

            metadata = world_util.world_metadata_from_chunks(
                recording_path,
                convert_to_seconds=True,
                verbose=False
            )

        finally:
            # Restore the current definition so we do not alter world_util
            # for anything else that uses it later.
            world_util.WORLD_AGC_METADATA_COLS = original_metadata_cols


        print("\nWORLD METADATA COLUMNS")
        print("----------------------")
        print(metadata.columns.tolist())


        camera_time = metadata["timestamp"].to_numpy(dtype=float)

        analog_gain = metadata["Again"].to_numpy(dtype=float)
        digital_gain = metadata["Dgain"].to_numpy(dtype=float)
        exposure = metadata["exposure"].to_numpy(dtype=float)

        # Put camera time onto a relative timebase for plotting.
        camera_time_relative = camera_time - camera_time[0]


        # INSPECT CAMERA METADATA
        # =============================================================

        print("\nCAMERA METADATA")
        print("----------------")
        print("Number of samples:", len(camera_time))

        print("Analog gain range:", np.nanmin(analog_gain), np.nanmax(analog_gain))
        print("Digital gain range:",  np.nanmin(digital_gain), np.nanmax(digital_gain))
        print("Exposure range:", np.nanmin(exposure), np.nanmax(exposure))


        # CALCULATE CAMERA LIGHT SCORE
        # =============================================================

        # Proposed camera light score:
        #
        #     1 / (Dgain * Again * exposure)
        #
        # Invalid / missing AGC values remain NaN.

        denominator = (analog_gain * digital_gain * exposure)
        camera_score = np.full(denominator.shape, np.nan, dtype=float)
        valid_camera = (np.isfinite(denominator) & (denominator > 0))
        camera_score[valid_camera] = (1.0 / denominator[valid_camera])

        print("Camera score range:", np.nanmin(camera_score), np.nanmax(camera_score))


        # PLOT AGC PARAMETERS
        # =============================================================

        plt.figure(figsize=(10, 5))

        plt.plot(camera_time_relative, analog_gain, label="Analog gain")
        plt.plot(camera_time_relative, digital_gain, label="Digital gain")

        plt.xlabel("Time from recording start (s)")
        plt.ylabel("Gain")
        plt.title("Camera gain settings over time")
        plt.legend()
        plt.grid()

        plt.tight_layout()
        plt.show()

        # Exposure separately because its numerical scale is
        # very different from the gain values.

        plt.figure(figsize=(10, 5))

        plt.plot(camera_time_relative, exposure)

        plt.xlabel("Time from recording start (s)")
        plt.ylabel("Exposure")
        plt.title("Camera exposure over time")
        plt.grid()

        plt.tight_layout()
        plt.show()

  
        # PLOT CAMERA LIGHT SCORE
        # =============================================================

        plt.figure(figsize=(10, 5))

        plt.plot(camera_time_relative, camera_score)

        plt.xlabel("Time from recording start (s)")
        plt.ylabel("Camera light score")

        plt.title(
            "Camera light score: "
            "1 / (Again × Dgain × exposure)"
        )

        plt.grid()
        plt.tight_layout()
        plt.show()


        # MINISPECT
        # =============================================================

        # Parse ONLY the minispect data.
        #
        # This avoids unnecessarily loading the large world-camera
        # frame buffers. The world metadata were obtained separately
        # above using world_metadata_from_chunks.

        ms_chunks = parse_chunks(
            recording_path,
            convert_time_units=True,
            convert_to_float=True,
            chunk_ranges={
                "M": (0, None)
            }
        )

        print("\nMINISPECT DATA")
        print("--------------")
        print("Number of parsed chunks:", len(ms_chunks))


        # INSPECT FIRST MINISPECT CHUNK
        # =============================================================

        if len(ms_chunks) == 0:
            raise RuntimeError(
                "No minispect chunks were parsed."
            )

        first_chunk = ms_chunks[0]

        print("\nM keys:")
        print(first_chunk["M"].keys())

        print("\nMinispect time keys:")
        print(first_chunk["M"]["t"].keys())

        print("\nMinispect value keys:")
        print(first_chunk["M"]["v"].keys())

        print("\nMinispect array shapes:")

        for sensor in first_chunk["M"]["v"]:

            print(
                sensor,
                "time:",
                np.shape(
                    first_chunk["M"]["t"][sensor]
                ),
                "values:",
                np.shape(
                    first_chunk["M"]["v"][sensor]
                )
            )


        # COMBINE / FLATTEN AS MINISPECT CHUNKS
        # =============================================================

        as_value_chunks = []
        as_time_chunks = []

        for chunk in ms_chunks:

            as_values_chunk = chunk["M"]["v"]["AS"]
            as_time_chunk = chunk["M"]["t"]["AS"]

            # Skip empty chunks
            if len(as_values_chunk) == 0:
                continue

            if len(as_values_chunk) != len(as_time_chunk):
                raise ValueError(
                    "AS value and timestamp lengths do not match "
                    "within a minispect chunk."
                )

            as_value_chunks.append(as_values_chunk)
            as_time_chunks.append(as_time_chunk)


        if len(as_value_chunks) == 0:
            raise RuntimeError(
                "No AS data were found in the minispect chunks."
            )


        # Stack all chunks into one continuous matrix.
        as_values = np.vstack(as_value_chunks)

        # AS timestamps are 1-D, so concatenate instead of vstack.
        as_time = np.concatenate(as_time_chunks).astype(float)

        # Put AS time onto a relative timebase for plotting.
        as_time_relative = as_time - as_time[0]

        print("\nFLATTENED AS DATA")
        print("-----------------")
        print("AS values shape:", as_values.shape)
        print("AS time shape:", as_time.shape)

        print(
            "AS time range:",
            np.nanmin(as_time),
            np.nanmax(as_time)
        )

        print(
            "AS relative duration:",
            np.nanmin(as_time_relative),
            np.nanmax(as_time_relative)
        )


        # OPTIONAL DATAFRAME FOR EASIER INSPECTION
        # =============================================================

        as_column_names = [f"AS_{ii}" for ii in range(as_values.shape[1])]
        as_df = pd.DataFrame(as_values, columns=as_column_names)
        as_df.insert(0, "timestamp", as_time)

        print("\nAS DATAFRAME")
        print("------------")
        print(as_df.head())

        print("\nAS SUMMARY")
        print("----------")
        print(as_df.describe())


        # PLOT RAW AS CHANNELS
        # =============================================================

        plt.figure(figsize=(11, 6))

        for column_number in range(as_values.shape[1]):
            plt.plot(as_time_relative, as_values[:, column_number], label=f"AS {column_number}")

        plt.xlabel("Time from AS recording start (s)")
        plt.ylabel("AS sensor value")
        plt.title("Flattened minispect AS channels")
        plt.legend(ncol=2, fontsize=8)
        plt.grid()
        plt.tight_layout()
        plt.show()

        # =============================================================
        # EXPERIMENT: CHARACTERIZE AGC TEMPORAL RESPONSE
        # =============================================================

        # Before converting the minispect data to illuminance, pick one minispect channel and 
        # characterize how its temporal changes relate to the camera light score.
        #
        # For this initial experiment, I use AS channel 0.
        #
        # This does NOT assume that AS 0 is illuminance. I use it
        # only to characterize the temporal filtering / lag of the AGC.

        as_channel_index = 0
        ms_signal = as_values[:, as_channel_index]

        filter_result = characterize_lowpass_filter(camera_time, camera_score, as_time, ms_signal,
                                                    lag_range=(-10, 10), tau_range=(0.1, 20),
                                                    n_lags=201, n_taus=100)

        print("\n======================================")
        print("AGC TEMPORAL FILTER RESULT")
        print("======================================")

        print("Minispect channel:", as_channel_index)
        print("Raw correlation:", filter_result["raw_correlation"])
        print("Best low-pass time constant (tau):", filter_result["best_tau"], "seconds")
        print("Best temporal lag:", filter_result["best_lag"], "seconds")
        print("Best correlation after filtering / lag:", filter_result["best_correlation"])

        # PLOT RAW NORMALIZED SIGNALS
        plt.figure(figsize=(11, 5))
        plt.plot(filter_result["ms_time"] - filter_result["ms_time"][0], filter_result["ms_z"], label=f"Minispect AS {as_channel_index}")
        plt.plot(filter_result["ms_time"] - filter_result["ms_time"][0], filter_result["camera_on_ms_time"], label="Camera light score")

        plt.xlabel("Time (s)")
        plt.ylabel("Normalized value (z-score)")
        plt.title("Raw minispect channel vs camera light score")

        plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.show()

        # PLOT BEST FILTERED / DELAYED RESULT
        plt.figure(figsize=(11, 5))
        plt.plot(filter_result["ms_time"] - filter_result["ms_time"][0], filter_result["camera_on_ms_time"], label="Observed camera light score")
        plt.plot(filter_result["ms_time"] - filter_result["ms_time"][0], filter_result["best_prediction"],
            label=(f"Filtered AS {as_channel_index} "f"(tau={filter_result['best_tau']:.2f} s, "f"lag={filter_result['best_lag']:.2f} s)"))

        plt.xlabel("Time (s)")
        plt.ylabel("Normalized value (z-score)")
        plt.title("Camera score vs best low-pass minispect prediction")

        plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.show()

        # PLOT CORRELATION VS LAG AT BEST TAU
        best_tau_idx = np.argmin(np.abs(filter_result["candidate_taus"] - filter_result["best_tau"]))
        correlations_at_best_tau = (filter_result["correlations"][best_tau_idx, :])

        plt.figure(figsize=(9, 5))
        plt.plot(filter_result["candidate_lags"], correlations_at_best_tau)
        plt.axvline(filter_result["best_lag"], linestyle="--",label=(f"Best lag = " f"{filter_result['best_lag']:.2f} s"))

        plt.xlabel("Camera lag relative to minispect (s)")
        plt.ylabel("Correlation")
        plt.title("Correlation vs temporal lag " "at best low-pass time constant")

        plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.show()

        # PLOT TAU / LAG SEARCH
        plt.figure(figsize=(9, 6))
        plt.imshow(filter_result["correlations"], aspect="auto", origin="lower",
            extent=[
                filter_result["candidate_lags"][0],
                filter_result["candidate_lags"][-1],
                filter_result["candidate_taus"][0],
                filter_result["candidate_taus"][-1]
            ])

        plt.colorbar(label="Correlation")
        plt.xlabel("Camera lag relative to minispect (s)")
        plt.ylabel("Low-pass time constant tau (s)")
        plt.title("Search for AGC temporal filter parameters")

        plt.tight_layout()
        plt.show()

        # =============================================================
        # TODO: CONVERT AS COUNTS TO ABSOLUTE ILLUMINANCE
        # =============================================================

        # Existing lightLoggerAnalysis code performs:
        #
        #     illuminance = msCounts2Illuminance(ms_values)
        #
        # We still need to confirm / port the Python equivalent of that
        # calibration before calculating Geoff's illuminance-vs-camera-score
        # relationship.
        #
        # illuminance = ...


        # =============================================================
        # TODO: ALIGN CAMERA SCORE AND ILLUMINANCE
        # =============================================================

        # Once illuminance is available:
        #
        # 1. Put camera score and illuminance on a shared time grid.
        # 2. Test a range of temporal offsets.
        # 3. Compute the correlation at every offset.
        # 4. Select the lag that maximizes the correlation.


        # =============================================================
        # TODO: CROSS-CORRELATE ACROSS TEMPORAL OFFSETS
        # =============================================================

        # Example future structure:
        #
        # candidate_lags = np.arange(-10, 10.01, 0.05)
        # correlations = np.full(candidate_lags.shape, np.nan)
        #
        # for ii, lag in enumerate(candidate_lags):
        #
        #     shifted_camera_time = camera_time_relative + lag
        #
        #     interpolated_camera_score = np.interp(
        #         as_time_relative,
        #         shifted_camera_time,
        #         camera_score,
        #         left=np.nan,
        #         right=np.nan
        #     )
        #
        #     valid = (
        #         np.isfinite(interpolated_camera_score) &
        #         np.isfinite(illuminance)
        #     )
        #
        #     correlations[ii] = np.corrcoef(
        #         interpolated_camera_score[valid],
        #         illuminance[valid]
        #     )[0, 1]
        #
        #
        # best_idx = np.nanargmax(correlations)
        # best_lag = candidate_lags[best_idx]


        # =============================================================
        # TODO: PLOT CORRELATION VS LAG
        # =============================================================

        # plt.figure()
        # plt.plot(candidate_lags, correlations)
        # plt.axvline(best_lag, linestyle="--")
        # plt.xlabel("Temporal offset (s)")
        # plt.ylabel("Correlation")
        # plt.title("Camera score vs illuminance: lag search")
        # plt.grid()
        # plt.show()


        # =============================================================
        # TODO: PLOT ALIGNED CAMERA SCORE AND ILLUMINANCE
        # =============================================================


        # =============================================================
        # TODO: FIT ILLUMINANCE ~= f(CAMERA SCORE)
        # =============================================================


        # =============================================================
        # SAVE CURRENT RESULTS FOR THIS RECORDING
        # =============================================================

        recording_result = {
            "video_path": video_path,
            "recording_path": recording_path,

            "camera_metadata": metadata,

            "camera_time": camera_time,
            "camera_time_relative": camera_time_relative,

            "analog_gain": analog_gain,
            "digital_gain": digital_gain,
            "exposure": exposure,

            "camera_score": camera_score,

            "as_time": as_time,
            "as_time_relative": as_time_relative,
            "as_values": as_values,
            "as_dataframe": as_df
        }

        results_by_recording.append(
            recording_result
        )


    # =================================================================
    # TODO: COMBINE RECORDINGS
    # =================================================================

    # Once illuminance and temporal alignment are implemented:
    #
    # camera_score_all = np.concatenate(all_camera_scores)
    # illuminance_all = np.concatenate(all_illuminance)


    # =================================================================
    # TODO: FIT FINAL CALIBRATION EQUATION
    # =================================================================

    #
    # illuminance ~= f(camera_score)
    #


    # =================================================================
    # RETURN EVERYTHING CURRENTLY AVAILABLE
    # =================================================================

    return results_by_recording


# =====================================================================
# RUN / TEST FUNCTION
# =====================================================================

if __name__ == "__main__":

    video_paths = [
        "/Users/sophiamirabal/Library/CloudStorage/"
        "Dropbox-Aguirre-BrainardLab/Team Documents/temp_for_sharing/"
        "walkOutdoorDemo.avi"
    ]

    recording_paths = [
        "/Users/sophiamirabal/Library/CloudStorage/"
        "Dropbox-Aguirre-BrainardLab/Team Documents/temp_for_sharing/"
        "outdoorDemoSophia"
    ]

    results = fit_agc_to_illuminance(
        video_paths,
        recording_paths
    )