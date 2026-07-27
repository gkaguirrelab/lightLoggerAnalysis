import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Pi_util / parse_chunks
sys.path.append(
    "/Users/sophiamirabal/Documents/MATLAB/projects/"
    "lightLogger/raspberry_pi_firmware/utility"
)

# world_util
sys.path.append(
    "/Users/sophiamirabal/Documents/MATLAB/projects/"
    "lightLoggerAnalysis/code/library/sensor_utility"
)

from Pi_util import parse_chunks
import world_util


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

    NOTE
    ----
    The video path is retained as part of the interface because the
    final procedure is intended to operate on a selectable set of
    environmental video recordings. The AGC metadata itself is loaded
    from the associated raw recording chunks.
    """

    if len(video_paths) != len(recording_paths):
        raise ValueError(
            "video_paths and recording_paths must have the same length."
        )

    all_camera_scores = []
    all_illuminance = []

    results_by_recording = []


    # =================================================================
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

        # =============================================================
        # CAMERA / AGC
        # =============================================================

        # This demo recording uses the older world-metadata format:
        #
        # timestamp, Again, Dgain, exposure
        #
        # The current world_util contains a newer five-setting metadata
        # definition, so temporarily use the legacy column definition while
        # loading this recording.

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


        # =============================================================
        # INSPECT CAMERA METADATA
        # =============================================================

        print("\nCAMERA METADATA")
        print("----------------")
        print("Number of samples:", len(camera_time))

        print("Analog gain range:", np.nanmin(analog_gain), np.nanmax(analog_gain))
        print("Digital gain range:",  np.nanmin(digital_gain), np.nanmax(digital_gain))
        print("Exposure range:", np.nanmin(exposure), np.nanmax(exposure))

        # =============================================================
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

        # =============================================================
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

        # =============================================================
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

        # =============================================================
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

        # =============================================================
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

        # =============================================================
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

        # =============================================================
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


        # =============================================================
        # PLOT RAW AS CHANNELS
        # =============================================================

        plt.figure(figsize=(11, 6))

        for column_number in range(as_values.shape[1]):

            plt.plot(
                as_time_relative,
                as_values[:, column_number],
                label=f"AS {column_number}"
            )

        plt.xlabel("Time from AS recording start (s)")
        plt.ylabel("AS sensor value")
        plt.title("Flattened minispect AS channels")
        plt.legend(ncol=2, fontsize=8)
        plt.grid()
        plt.tight_layout()
        plt.show()

        # =============================================================
        # PLOT CAMERA LIGHT SCORE AGAIN FOR COMPARISON
        # =============================================================

        plt.figure(figsize=(11, 5))

        plt.plot(camera_time_relative, camera_score)
        
        plt.xlabel("Time from recording start (s)")
        plt.ylabel("Camera light score")
        plt.title("Camera light score")
        plt.grid()

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