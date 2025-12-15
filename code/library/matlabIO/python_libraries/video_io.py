import numpy as np 
import cv2
import warnings
import time
import queue
from natsort import natsorted
from typing import Literal, Iterable
import os
import hdf5storage
import tqdm
import h5py
import dill
import pathlib
import sys

# Import utility libraries
light_logger_analysis_dir_path: str = os.path.expanduser("~/Documents/MATLAB/projects/lightLoggerAnalysis")
world_util: str = os.path.join(light_logger_analysis_dir_path, "code", "library", "sensor_utility")
for path in (light_logger_analysis_dir_path, world_util):
    assert os.path.exists(path), f"Expected path: {path} does not exist"
    sys.path.append(path)


"""Given a directory of frames, read them in and convert to video"""
def dir_to_video(dir_path: str, output_path: str, fps: float=30) -> None:
    # Gather all of the frame paths 
    frame_paths: list[str] = [os.path.join(dir_path, filename)
                              for filename in natsorted(os.listdir(dir_path))
                              if not filename.startswith(".")
                             ]
    
    assert len(frame_paths) > 0, "Must have non-zero amount of frames"

    # Read in the first file 
    sample_frame: np.ndarray = cv2.imread(frame_paths[0])

    # Initialize a video writer based on info from this frame 
    video_writer: cv2.VideoWriter = cv2.VideoWriter(output_path, cv2.VideoWriter_fourcc(*"FFV1"), 
                                                    float(fps), sample_frame.shape[:2][::-1], isColor=True
                                                   )
    if(not video_writer.isOpened()):
        raise RuntimeError(f"Could not open video writer to path: {output_path}")

    # Iterate over the frames and add them to the video writer 
    for frame_num, frame_path in enumerate(frame_paths):
        frame: np.ndarray | None = cv2.imread(frame_path)
        
        if(frame is None):
            video_writer.release() 
            raise Exception(f"Could not read frame: {frame_path}")
        
        if(frame.shape != sample_frame.shape):
            video_writer.release()
            raise Exception(f"Video frames are of inhomogenous shape")

        video_writer.write(frame)

    # Release the video writer 
    video_writer.release() 

    return 


"""Given the frames to a video, generate an .avi video at the desired location 
   of those frames at the specified FPS
"""
def frames_to_video(frames: np.ndarray | queue.Queue, output_path: str, fps: float, 
                    stop_event: object=None, timeout: float | None=None 
                   ) -> None:
   
    # Retreive a sample frame from the video so we can get some information 
    # about how it should look
    sample_frame: np.ndarray | None = None
    start_time: float = time.time() 
    if(isinstance(frames, np.ndarray)):
        sample_frame = frames[0]
    else:
        # Wait for a sample frame to write
        while(sample_frame is None):
            current_time: float = time.time() 
            elapsed_time: float = current_time - start_time

            # If we have passed the timeout, simply quit
            if(elapsed_time > timeout):
                break
            
            # If we have received a stop signal, exit the function 
            if(stop_event is not None and stop_event.is_set()):
                return

            # Try to read a frame 
            try:
                sample_frame: np.ndarray = frames.get()
            except queue.Empty:
                continue 
    
    if(sample_frame is None):
        raise Exception("Did not receive a frame to write")
    
    # Retrieve the dimensions and if it is color or not
    frame_height, frame_width = sample_frame.shape[:2]
    is_color = sample_frame.ndim == 3
    num_frames: int | float = len(frames) if isinstance(frames, np.ndarray) else float("inf")

    # Initialize the video wrtier
    video_writer: cv2.VideoWriter = cv2.VideoWriter(output_path, 0, fps, (frame_width, frame_height), isColor=is_color ) 

    # Next, write frames to the video 
    frame_num: int = 0
    while(True):
        # Define frame variable 
        frame: np.ndarray | None | int = -1

        # If we have finished writing a defined number of frames, quit 
        if(frame_num >= num_frames):
            break

        # If we have triggered a stop event, stop writing 
        if(stop_event is not None and stop_event.is_set()):
            break 

        # Retrieve the frame either from the frames 
        if(isinstance(frames, np.ndarray)):
            frame = frames[frame_num]
        else:
            # Attempt to get a frame from the frame stream 
            try:
                frame = frames.get(timeout=timeout)
            except queue.Empty:
                # Close the video writer on error 
                video_writer.release()
                raise Exception("Did not receive a frame in time")

        # If we have finished writing an unknown number of frames, quit 
        if(frame is None):
            break

        # Otherwise, write the frame to the video 
        video_writer.write(frame) 

        # Increment the frame number 
        frame_num += 1 

    # Release the video writer 
    video_writer.release()

    return  

"""Given a path to a video, return the size of each frame of the video"""
def inspect_video_framesize(video_path: str) -> tuple:
    assert os.path.exists(video_path), f"Video path: {video_path} does not exist"

    # Open the video stream
    video_stream: cv2.VideoCapture = cv2.VideoCapture(video_path)
    if(not video_stream.isOpened()):
        raise Exception(f"Failed to open video: {video_path}")

    # Read the first frame to get video size 
    success, frame = video_stream.read() 
    if(not success):
        raise Exception("Failed to read sample frame")

    # Close the video stream 
    video_stream.release() 

    return frame.shape[:2]

"""Given a path to a video, return the FPS of that video"""
def inspect_video_FPS(video_path: str) -> float:
    assert os.path.exists(video_path), f"Video path: {video_path} does not exist"

    # Open the video via cv2 
    video_stream: cv2.VideoCapture = cv2.VideoCapture(video_path)

    # Calculate the number of frames 
    fps: float = int(video_stream.get(cv2.CAP_PROP_FPS)) 
    assert fps != 0, "Could not estimate FPS. cv2 returned 0 FPS."

    # Close the capture stream 
    video_stream.release()  

    # Return the number of frames 
    return fps


"""Given a path to a video, return the number of frames in that video"""
def inspect_video_frame_count(video_path: str) -> int:
    assert os.path.exists(video_path), f"Video path: {video_path} does not exist"

    # Open the video via cv2 
    video_stream: cv2.VideoCapture = cv2.VideoCapture(video_path)

    # Immediately see if we can get a result. Sometimes this is erroneously 0 
    frame_count = int(video_stream.get(cv2.CAP_PROP_FRAME_COUNT))
    if(frame_count != 0):
        video_stream.release() 
        return frame_count

    # Next, write frames to the video 
    frame_num: int = 0
    while(True):
        # Attempt to read in a video 
        # from the stream 
        ret, frame = video_stream.read() 

        # If we could not read in the frame 
        # we reached the end of the video
        if(not ret): 
            break

        # Increment the frame number 
        frame_num += 1 

    # Close the capture stream 
    video_stream.release()  

    # Return the number of frames 
    return frame_num

"""Given a path to a video, extract the provided frames idx from it as a numpy array"""
def extract_frames_from_video(video_path: str, frames_idx: Iterable, is_grayscale: bool=False) -> np.ndarray: 
    assert os.path.exists(video_path), f"Video path: {video_path} does not exist"
    
    # First, sort the frame numbers we want to extract 
    frames_idx: list[int] = [int(idx) for idx in frames_idx] 
    frames_idx.sort() 

    # Initialize array to hold frames 
    # post extraction
    extracted_frames: list[int] = []

    # Open the video via cv2 
    video_stream: cv2.VideoCapture = cv2.VideoCapture(video_path)

    # Stream the frames in from the vidoe 
    extracted_frame_idx: int = 0 
    while(extracted_frame_idx < len(frames_idx)):
        # Jump the the target frame in the stream
        video_stream.set(cv2.CAP_PROP_POS_FRAMES, frames_idx[extracted_frame_idx])

        # Attempt to read in a video 
        # from the stream 
        ret, frame = video_stream.read() 

        # If we could not read in the frame 
        # or we finished capturing the desired frames
        # end the loop with error 
        if(not ret): 
            video_stream.release() 
            raise RuntimeError(f"Could not extract frame idx: {frames_idx[extracted_frame_idx]}")

        # if we want a grayscale video, 
        # convert to grayscale
        if(is_grayscale is True): 
            frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

        # Save the extracted frame 
        extracted_frames.append(frame)

        # Increment the number of extracted frames 
        extracted_frame_idx += 1 

    # Release the video stream 
    video_stream.release() 

    return np.array(extracted_frames, dtype=np.uint8)

"""Given a video file, destruct this video into an array of 
   frames within a given bounds of frames in the video (exclusive).
   Return as a 2D or 3D array based on if the original video was grayscale.
   Optionally stream the frames to a queue, not loading in all of them into 
   memory at once
"""
def destruct_video(video_path: str, start_frame: int=0, end_frame: int=float("inf"), 
                   is_grayscale: bool=False, q: queue.Queue=None, stop_event: object=None
                  ) -> None | np.ndarray:
    # Open the video via cv2 
    video_stream: cv2.VideoCapture = cv2.VideoCapture(video_path)

    # Initialize a container to hold the frames 
    frames: list[np.ndarray] = []
    frame_num: int = start_frame

    # Stream the frames in from the vidoe 
    while(frame_num < end_frame):
        # Jump the the target frame in the stream
        video_stream.set(cv2.CAP_PROP_POS_FRAMES, frame_num)

        # Attempt to read in a video 
        # from the stream 
        ret, frame = video_stream.read() 

        # If we could not read in the frame 
        # or we finished capturing the desired frames
        # end the loop
        if(not ret): 
            # If streaming these frames elsewhere, apply 
            # None to the queue 
            if(q is not None): q.put(None)
            break 

        # If this process has been killed, exit gracefully 
        if(stop_event is not None and stop_event.is_set()):
            break

        # if we want a grayscale video, 
        # convert to grayscale
        if(is_grayscale is True): 
            frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

        # Append this frame to the frame list 
        # Or add to the queue if provided for streaming
        if(q is not None):
            q.put(frame)
        else:
            frames.append(frame)

        # Increment the frame number 
        frame_num += 1 

    # Close the capture stream 
    video_stream.release()  

    # If streaming from a queue, simply return None
    if(q is not None): 
        q.put(None)
        return None

    # Convert to standardized np.ndarray 
    frames: np.ndarray = np.stack(frames, axis=0)

    return frames

"""Convert a video to a matlab 7.3 .mat file saved in HDF5 format"""
def video_to_hdf5(video_path: str, output_path: str, 
                 color_mode: Literal["GRAY", "RGB", "BGR"]="RGB",
                 start_frame: int=0, 
                 end_frame: int | float = float("inf"),
                 zeros_as_nans: bool=False,
                 visualize_results: bool=False
                ) -> None:
    assert os.path.exists(video_path), f"Video path: {video_path} does not exist"
    assert output_path.endswith(".hdf5"), f"Output path: {output_path} must end in .hdf5"

    # Open the video via cv2 
    video_stream: cv2.VideoCapture = cv2.VideoCapture(video_path)
    if(not video_stream.isOpened()):
        raise RuntimeError(f"Could not open video: {video_path}")
    
    # Move the video pointer to the start of the target interval 
    video_stream.set(cv2.CAP_PROP_POS_FRAMES, start_frame)
    end_frame = end_frame if end_frame != float("inf") else inspect_video_frame_count(video_path)

    # Open the output file with an empty dataset of the appropriate type 
    # and shape
    frame_height, frame_width = inspect_video_framesize(video_path)
    frame_array_shape_list: list[int] = [0, frame_height, frame_width]
    if(color_mode != "GRAY"):
        frame_array_shape_list += [3]
    frame_array_shape_tuple: tuple = tuple(frame_array_shape_list)

    max_frame_array_shape_list: list[int | None] = frame_array_shape_list.copy() 
    max_frame_array_shape_list[0] = None
    max_frame_array_shape_tuple: tuple[int | None] = tuple(max_frame_array_shape_list)

    chunks_shape_list: list[int] = frame_array_shape_list.copy()
    chunks_shape_list[0] = 1 
    chunks_shape_tuple: tuple[int] = tuple(chunks_shape_list)

    with h5py.File(output_path, "w") as f:
        dset: object = f.create_dataset("video",
                                        shape=frame_array_shape_tuple,
                                        maxshape=max_frame_array_shape_tuple,
                                        chunks=chunks_shape_tuple,
                                        dtype=np.float64
                                      )


        # Read frames from the interval 
        written_frame_idx: int = 0
        iterator: Iterable = tqdm.tqdm(range(start_frame, end_frame)) if visualize_results is True else range(start_frame, end_frame)
        for frame_num in iterator:
            # Attempt to read in a video 
            # from the stream 
            ret, frame = video_stream.read() 

            # If we could not read in the frame 
            # or we finished capturing the desired frames
            # end the loop
            if(not ret):
                # If we could not read the frame, end frame was specified, and we have not reached it
                # throw an error 
                if(end_frame != float("inf")):
                    video_stream.release()
                    raise Exception(f"Did not read all frames from interval [{start_frame}, {end_frame}). Frame: {frame_num} could not be read")

                break 

            # if we want a grayscale video, 
            # convert to grayscale
            if(color_mode == "GRAY"): 
                frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
            elif(color_mode == "RGB"):
                frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
            elif(color_mode == "BGR"):
                pass 
            else:
                video_stream.release() 
                raise Exception(f"Color mode: {color_mode} is unsupported")


            # Convert to floats so we can use NaN 
            if(zeros_as_nans is True):
                frame = frame.astype(np.float64)
                frame[frame == 0] = np.nan

            # Append to the dataset 
            if color_mode.upper() == "GRAY":
                    dset.resize((written_frame_idx + 1, frame_height, frame_width))
                    dset[written_frame_idx, :, :] = frame
            else:
                dset.resize((written_frame_idx + 1, frame_height, frame_width, 3))
                dset[written_frame_idx, :, :, :] = frame

            # Increment the written frame number
            written_frame_idx += 1

    # Close the capture stream 
    video_stream.release()  

    return 


"""Convert the world chunks of a recording into a playable video"""
def world_chunks_to_video(recording_path: str,
                          output_path: str,
                          debayer_images: bool=False, 
                          convert_to_lms: bool=False, 
                          apply_color_correction: bool=False, 
                          apply_fielding_function: bool=False, 
                          apply_digital_gain: bool=False,
                          fill_missing_frames: bool=False,
                          convert_to_seconds: bool=True,
                          embed_timestamps: bool=True
                         ) -> None:
    
    # First let's read the config file of this recording 
    # to gather some information about how the images look 
    config_filepath: str = os.path.join(recording_path, "config.pkl")
    config_data: dict | None = None 
    with open(config_filepath, 'rb') as f:
        config_data = dill.load(f)
    assert(config_data is not None)

    # Retrieve the FPS and the dimensions of the recording 
    frame_height, frame_width = config_data["sensors"]['W']["sensor_mode"]['size'][::-1]
    FPS = config_data["sensors"]['W']["sensor_mode"]['fps']

    # Initialize the video writer to construct our video. Note: here we denote an FPS, 
    # assuming that it is constant, but we will need to handle dropped frames accordingly 
    # in terms of when we change frames. 
    fourcc = cv2.VideoWriter_fourcc(*"FFV1")
    video_writer = cv2.VideoWriter(output_path, 
                                   fourcc, 
                                   FPS, 
                                   (frame_width, frame_height), 
                                   isColor= (debayer_images is True or convert_to_lms is True) 
                                  ) 

    # Initialize a dummy frame to fill in missing frames if desired 
    dummy_frame: np.ndarray =  np.squeeze(np.full((frame_height, frame_width, 3 if debayer_images is True else 1), 0, dtype=np.uint8))

    # Retrive the sorted chunks for this sensor

    # Iterate over the chunks for this sensor 
    previous_chunk_end_time: float = 0 
    current_chunk_start_time: float = 0 

    for chunk_num, (metadata_matrix_path, frame_buffer_path) in enumerate(chunks_paths):    
        # First, we will retrieve the metadata for this buffer 
        metadata: np.ndarray = np.load(metadata_matrix_path)
        
        # Then, we will read in the timestamps and frame buffer for this chunk 
        # World timestamps are in nanoseconds and thus must be converted to seconds 
        t_vector: np.ndarray = np.ascontiguousarray(metadata[:, 0], dtype=np.float64) / ( (10 ** 9) if convert_to_seconds is True else 1)
        frame_buffer: np.ndarray = np.load(frame_buffer_path)

        # Assert the t vector and frame vector are the same size 
        assert(len(t_vector) == len(frame_buffer))
        
        # Assert the frames are of appropriate dtype 
        assert(frame_buffer.dtype == np.uint8)

        # If we somehow got an empty chunk, skip it 
        if(len(t_vector) == 0): 
            continue

        # Now, we will apply the requested transformations to the frames 
        frame_buffer = frame_buffer.astype(np.float64) # First, convert to float for float multiplications 

        # Apply Dgain if requested
        if(apply_digital_gain is True):
            buffer_dgains: np.ndarray = metadata[:, 1]
            buffer_dgains = buffer_dgains[:, np.newaxis, np.newaxis]
            frame_buffer *= buffer_dgains

        # If we want to apply the color weight correction 
        if(apply_color_correction is True):
            frame_buffer[:] = [ world_util.apply_color_weights(frame)
                                for frame in frame_buffer
                                ] 

        # Apply the fielding function if desired 
        if(apply_fielding_function is True):
            frame_buffer[:] = [ world_util.apply_fielding_function(frame)
                               for frame in frame_buffer 
                              ]

        # If we want to convert to LMS, do so 
        if(convert_to_lms is True):
            pass

        # If we want to debayer the image
        if(debayer_images is True):
            # Note: When writing color images to cv2.VideoWriter, it expects a BGR instead of RGB 
            # frame, so also convert now 
            frame_buffer[:] = [cv2.cvtColor(world_util.debayer_image(frame), cv2.COLOR_RGB2BGR) 
                               for frame in frame_buffer 
                              ]

        # Convert the frame buffer back into uint8 format 
        frame_buffer = np.clip(np.round(frame_buffer), 0, 255).astype(np.uint8)

        # If we want to embed the timestamps in the frame buffer, do that now
        if(embed_timestamps is True):
            frame_buffer = [world_util.embed_timestamp(frame, timestamp)
                            for frame, timestamp in zip(frame_buffer, t_vector)
                           ]

        # Retrieve the current chunk start time 
        current_chunk_start_time = t_vector[0]

        # Before we write the frames for this chunk, we must write the dummy frame (if desired)
        # for those frames between the previous chunk and the current chunk 
        
        # Find the total elapsed time in seconds between the two chunks 
        time_between_chunks: int = current_chunk_start_time - previous_chunk_end_time

        # Calculate the number of missed frames as the elapsed time divided 
        # by the seconds per frame, minus one frame as we have to count the current frame 
        # as captured during this interval 
        missed_frames: int = 0 if chunk_num == 0 or fill_missing_frames is False else int( (time_between_chunks / (1/FPS ) ) -  1) 
        missing_timestamps: np.ndarray = np.linspace(previous_chunk_end_time + (1/FPS), current_chunk_start_time, missed_frames, endpoint=False)

        # Write the number of missing frames in between as the previous frame 
        for missing_timestamp in missing_timestamps:
            missing_frame: np.ndarray = world_util.embed_timestamp(dummy_frame, missing_timestamp) if embed_timestamps is True else dummy_frame
            video_writer.write(missing_frame)

        # Initialize variables to track the previous timestamp 
        # We will use this delta with the current timestamp 
        # to track dropped frames and add dummy frames in the meantime
        previous_timestamp: float = t_vector[0]

        # Write frames to the video 
        for frame_num, (timestamp, frame) in enumerate(zip(t_vector, frame_buffer)):
            # Assert the frame is of the appropriate type after transformation 
            assert(frame.dtype == np.uint8 and dummy_frame.dtype == np.uint8)

            # Otherwise, let's calculate the time delta between the current timestamp and the previous timestamp 
            time_between_frames: float = timestamp - previous_timestamp

            # Calculate the number of missed frames as the elapsed time divided 
            # by the seconds per frame, minus one frame as we have to count the current frame 
            # as captured during this interval 
            missed_frames: int = 0 if frame_num == 0 or fill_missing_frames is False else int( (time_between_frames / (1/FPS ) ) -  1) 
            missing_timestamps: np.ndarray = np.linspace(previous_timestamp + (1/FPS), timestamp, missed_frames, endpoint=False)

            # Write the number of missing frames in between as the previous frame 
            for missing_timestamp in missing_timestamps:
                missing_frame: np.ndarray = world_util.embed_timestamp(dummy_frame, missing_timestamp) if embed_timestamps is True else dummy_frame
                video_writer.write(missing_frame)
            
            # Write the current frame 
            video_writer.write(frame)

            # Save the current timestamp as the previous timestamp
            # so we can reference it for the next frame 
            previous_timestamp = timestamp 
        
        # Save the end time of this chunk 
        previous_chunk_end_time = t_vector[-1]

    return 


def main():
    return 

if(__name__ == "__main__"):
    main()