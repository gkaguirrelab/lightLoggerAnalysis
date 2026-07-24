"""Video I/O utilities for light logger recordings.

Provides functions for reading, writing, converting, and transforming
video data from the light logger's world camera sensor, including
chunk-based video assembly, HDF5 export, frame extraction, and
OCR-based event detection.
"""

import numpy as np
import cv2
import warnings
import time
import queue
from natsort import natsorted
from typing import Literal, Iterable
import os
import hdf5storage
from tqdm.auto import tqdm
import h5py
import dill
import pathlib
import sys
from scipy.signal import find_peaks
import pytesseract
import matplotlib.pyplot as plt
import re
import subprocess
import multiprocessing as mp
from multiprocessing import shared_memory
from numba import njit, prange

# Import utility libraries
light_logger_analysis_dir_path: str = os.path.expanduser("~/Documents/MATLAB/projects/lightLoggerAnalysis")
world_util_path: str = os.path.join(light_logger_analysis_dir_path, "code", "library", "sensor_utility")
for path in (light_logger_analysis_dir_path, world_util_path):
    assert os.path.exists(path), f"Expected path: {path} does not exist"
    sys.path.append(path)

# Import world camera utility library
import world_util

def group_sensors_files(recording_path: str) -> dict[str, list[tuple]]:
    """Pair per-sensor chunk metadata files with their data payloads.

    Each sensor chunk is stored as two sibling files: a metadata array whose
    filename contains ``"_metadata"`` and a value file with the same stem
    minus that marker. This helper performs that filename bookkeeping once
    and returns naturally sorted pairs for the world, pupil, and MS
    sensors.

    Args:
        recording_path: Directory containing the chunk files for one
            recording.

    Returns:
        Dictionary keyed by ``"W"``, ``"P"``, and ``"M"`` whose values are
        ordered ``(metadata_path, value_path)`` tuples.
    """
    def group_sensor_files(sensor_name: str) -> list:
        """Collect the ordered chunk-file pairs for one sensor namespace.

        Args:
            sensor_name: Filename fragment identifying the sensor, such as
                ``"world"``, ``"pupil"``, or ``"ms"``.

        Returns:
            Naturally sorted list of ``(metadata_path, value_path)`` tuples
            for that sensor.
        """
        return [ ( os.path.join(recording_path, file), os.path.join(recording_path, file.replace("_metadata", "") ) ) 
                   for file in natsorted(os.listdir(recording_path)) 
                   if sensor_name in file and "metadata" in file
               ]

    return {sensor[0].upper(): group_sensor_files(sensor) # n chunks = [ (metadata_path, frame_buffer_path), ...  ]
            for sensor in ("world", "pupil", "ms")
           }

def dir_to_video(dir_path: str, output_path: str, fps: float=30) -> None:
    """Encode a directory of image files into a single FFV1 video.

    The function natural-sorts the visible files in ``dir_path``, uses the
    first image to establish the frame size, and then writes every frame to
    an OpenCV ``VideoWriter`` after verifying that all images share the same
    shape.

    Args:
        dir_path: Directory containing one image file per frame.
        output_path: Destination video path.
        fps: Frame rate to encode into the output file.
    """
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


def frames_to_video(frames: np.ndarray | queue.Queue, output_path: str, fps: float,
                    stop_event: object=None, timeout: float | None=None
                   ) -> None:
    """Write frames to a video from either a stack or a live queue.

    Two usage patterns are supported. If ``frames`` is a NumPy array, the
    function writes that fixed stack directly. If it is a queue-like object,
    the function waits for the first frame to determine the video geometry
    and then keeps consuming frames until it receives a ``None`` sentinel,
    times out, or sees an optional stop event.

    Args:
        frames: Either a NumPy frame stack or a queue that yields frames
            followed by ``None``.
        output_path: Destination video path.
        fps: Frame rate written into the output container.
        stop_event: Optional event that aborts waiting or writing early.
        timeout: Maximum number of seconds to wait for a queued frame before
            raising an exception.
    """
   
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

def inspect_video_framesize(video_path: str) -> tuple:
    """Return the (height, width) frame dimensions of a video.

    Args:
        video_path: Path to the video file.

    Returns:
        Tuple of (frame_height, frame_width).
    """
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

def inspect_video_FPS(video_path: str) -> float:
    """Return the frames per second of a video.

    Falls back to a warning if cv2 reports 0 FPS.

    Args:
        video_path: Path to the video file.

    Returns:
        Frames per second as a float.
    """
    assert os.path.exists(video_path), f"Video path: {video_path} does not exist"

    # Open the video via cv2 
    video_stream: cv2.VideoCapture = cv2.VideoCapture(video_path)

    # Calculate the number of frames 
    fps: float = int(video_stream.get(cv2.CAP_PROP_FPS)) 
    if(fps == 0):
        warnings.warn("cv2 returned 0 FPS. Falling back to FFMPEG...")

    # Close the capture stream 
    video_stream.release()  

    # Return the number of frames 
    return fps

def inspect_video_frame_count(video_path: str) -> int:
    """Return the total number of frames in a video.

    First attempts to read the frame count from cv2 metadata. If that
    returns zero, falls back to counting frames by iterating through
    the entire video.

    Args:
        video_path: Path to the video file.

    Returns:
        Total number of frames in the video.
    """
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

def extract_frames_from_video(video_path: str, 
                              frames_idx: Iterable, 
                              is_grayscale: bool=False
                            ) -> np.ndarray:
    """Extract specific frames from a video by their indices.

    Seeks to each requested frame index and reads it, converting to
    grayscale or RGB as specified.

    Args:
        video_path: Path to the video file.
        frames_idx: Iterable of frame indices to extract.
        is_grayscale: If True, converts frames to grayscale; otherwise
            converts from BGR to RGB.

    Returns:
        Numpy array of extracted frames with dtype uint8.
    """
    assert os.path.exists(video_path), f"Video path: {video_path} does not exist"
    
    # First, sort the frame numbers we want to extract 
    frames_idx: list[int] = frames_idx if isinstance(frames_idx, list) else [int(idx) for idx in frames_idx] 
    frames_idx.sort() 

    # Initialize array to hold frames 
    # post extraction
    video_frame_shape: tuple[int] = inspect_video_framesize(video_path)
    extracted_frame_shape: tuple[int] = video_frame_shape if is_grayscale is True else tuple(list(video_frame_shape) + [3])
    extracted_frames: np.ndarray = np.empty((len(frames_idx), *extracted_frame_shape), dtype=np.uint8)


    print(f"Expected frames shape: {extracted_frames.shape}")

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
        else:
            frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        
        print(f"Frame shape")

        # Save the extracted frame 
        extracted_frames[extracted_frame_idx] = frame

        # Increment the number of extracted frames 
        extracted_frame_idx += 1 
    
    # Assert that we extracted all of the correct frames 
    assert extracted_frame_idx == len(extracted_frames), f"Error: Incorrect number of frames extracted: {extracted_frame_idx} compared to target: {len(extracted_frames)}"

    # Release the video stream 
    video_stream.release() 

    return extracted_frames

def destruct_video(video_path: str, start_frame: int=0, end_frame: int=float("inf"),
                   is_grayscale: bool=False, q: queue.Queue=None, stop_event: object=None,
                   verbose: bool=False,
                  ) -> None | np.ndarray:
    """Decompose a video into individual frames within a given range.

    Reads frames from start_frame to end_frame (exclusive), optionally
    converting to grayscale. Frames can be returned as a numpy array or
    streamed to a multiprocessing Queue with a None sentinel at the end.

    Args:
        video_path: Path to the video file.
        start_frame: First frame index to read (supports negative indexing).
        end_frame: Last frame index (exclusive). Use float("inf") for all
            remaining frames. Supports negative indexing.
        is_grayscale: If True, converts frames to grayscale.
        q: Optional multiprocessing Queue for streaming frames. When
            provided, frames are put into the queue instead of returned.
        stop_event: Optional multiprocessing Event to signal early stop.
        verbose: If True, displays a progress bar.

    Returns:
        Numpy array of frames if q is None, otherwise None (frames are
        streamed to the queue).
    """
    # Let's get some information about the video 
    frame_height, frame_width = inspect_video_framesize(video_path)
    num_frames: int = inspect_video_frame_count(video_path)

    # Support inf indexing
    end_frame = end_frame if end_frame != float("inf") else num_frames  

    # Support negative indexing
    if(start_frame < 0):
        start_frame = int(num_frames - abs(start_frame))

    if(end_frame < 0):
        end_frame = int(num_frames - abs(end_frame))


    output_shape: tuple = (end_frame - start_frame, frame_height, frame_width) if is_grayscale is True else (end_frame - start_frame, frame_height, frame_width, 3)

    # Open the video via cv2 
    video_stream: cv2.VideoCapture = cv2.VideoCapture(video_path)

    # Initialize a container to hold the frames 
    frames: np.ndarray | None = None if q is not None else np.empty(shape=output_shape, dtype=np.uint8)

    # Stream the frames in from the video 
    iterator: Iterable = tqdm(range(start_frame, end_frame)) if verbose is True else range(start_frame, end_frame)
    insertion_idx: int = 0 

    # Set the cursor at the start of the range we want to explore 
    video_stream.set(cv2.CAP_PROP_POS_FRAMES, start_frame)
    for frame_num in iterator:
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
        else:
            frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)

        # Append this frame to the frame list 
        # Or add to the queue if provided for streaming

        if(q is not None):
            q.put(frame)
        else:
            frames[insertion_idx] = frame

        insertion_idx += 1 

    # Close the capture stream 
    video_stream.release()  

    # If streaming from a queue, simply return None
    if(q is not None): 
        q.put(None)
        return None

    return frames

def video_to_hdf5(video_path: str, output_path: str,
                 color_mode: Literal["GRAY", "RGB", "BGR"]="RGB",
                 start_frame: int=0, 
                 end_frame: int | float = float("inf"),
                 apply_floor_celing: bool=False, 
                 zeros_as_nans: bool=False,
                 ceiling_as_nans: bool=False, 
                 visualize_results: bool=False,
                 ceiling: int=255, 
                 floor: int=0
                ) -> None:
    """Convert a video to an HDF5 file with optional pixel transformations.

    Reads frames from the specified interval, applies optional floor/ceiling
    clipping and NaN substitution, and writes each frame as float64 into a
    resizable HDF5 dataset.

    Args:
        video_path: Path to the input video file.
        output_path: Output path for the HDF5 file (must end in .hdf5).
        color_mode: Color space for output frames: "GRAY", "RGB", or "BGR".
        start_frame: First frame index to process.
        end_frame: Last frame index (exclusive). Use float("inf") for all.
        apply_floor_celing: If True, clips pixels where any channel hits
            the floor or ceiling to fully black or white.
        zeros_as_nans: If True, replaces zero-valued pixels with NaN.
        ceiling_as_nans: If True, replaces ceiling-valued pixels with NaN.
        visualize_results: If True, displays a progress bar.
        ceiling: Upper pixel value threshold for clipping.
        floor: Lower pixel value threshold for clipping.
    """
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
        iterator: Iterable = tqdm(range(start_frame, end_frame)) if visualize_results is True else range(start_frame, end_frame)
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

            # Covert this frame to float right away 
            # for our transformations and for the fact 
            # that it will be saved as a float
            frame = frame.astype(np.float64)


            # If we want to apply floor ceiling, 
            # we clip the values of each pixel (and all channels in that
            # pixel to the DARKNOISE < x <= 255 RANGE)
            if(apply_floor_celing is True):
                # Clip to the ceiling
                frame = np.where((frame >= ceiling).any(axis=-1, keepdims=True),
                                     255,
                                     frame
                                    )
                # Clip to the floor 
                frame = np.where((frame <= floor).any(axis=-1, keepdims=True),
                                    0,
                                    frame
                                )

            # We need to handle zeros to NaNs 
            # and ceilings to NaN differently 
            # based on the shape if the image 
            if(frame.ndim == 2):
                
                # If we are in a 2D image and want 
                # to make zeros nan, that is easy
                if(zeros_as_nans is True):
                    frame[frame == 0] = np.nan

                # Same for ceiling 
                if(ceiling_as_nans is True):
                    frame[frame >= 255] = np.nan

            elif(frame.ndim == 3):
                # If we are 3D image, we want to make 
                # pixels that ALL channels are 0 
                # into NaNs
                if(zeros_as_nans is True):
                    # Create a mask that is all of the pixels whose 
                    # all 3 channels are == 0
                    zero_mask_three_d: np.ndarray = np.all(frame == floor, axis=2)

                    # Set these pixels equal to fully NaN
                    frame[zero_mask_three_d] = np.nan

                # Same for ceiling, but >= 255 
                if(ceiling_as_nans is True):
                    # Create a mask that is all of the pixels whose 
                    # all 3 channels are >= 255
                    ceiling_mask_three_d: np.ndarray = np.all(frame >= ceiling, axis=2)

                    # Set these pixels equal to fully NaN
                    frame[ceiling_mask_three_d] = np.nan


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

def _load_world_chunk_shared_memory(chunk_paths: list[tuple[str, str]],
                                    chunk_indices: list[int],
                                    metadata_shape: tuple[int, ...],
                                    metadata_dtype: np.dtype,
                                    frame_shape: tuple[int, ...],
                                    frame_dtype: np.dtype,
                                    metadata_shm_names: list[str],
                                    frame_shm_names: list[str],
                                    free_buffer_queue: mp.Queue,
                                   ready_buffer_queue: mp.Queue,
                                   stop_event: object
                                   ) -> None:
    # First, let's access the shared memory that we allocated by the main processs 
    # This is cheap: we are only attaching this process to named shared-memory blocks,
    # not copying any chunk data yet. These SharedMemory objects are just handles that
    # let this process see the same byte buffers as the main process.
    """Background loader that stages world chunks into shared memory.

    ``world_chunks_to_video`` uses this helper in a separate process so that
    disk reads for the next chunk can overlap with transformation and video
    writing of the current chunk. The loader waits for the main process to
    announce a free shared-memory slot, memory-maps the next chunk from
    disk, copies its metadata and frame data into that slot, and reports the
    ready slot number and valid lengths back through a small queue message.

    Args:
        chunk_paths: Ordered ``(metadata_path, frame_path)`` tuples for the
            world sensor.
        chunk_indices: Chunk indices that this loader should process.
        metadata_shape: Shape used to wrap each shared metadata buffer.
        metadata_dtype: Dtype used to interpret the metadata bytes.
        frame_shape: Shape used to wrap each shared frame buffer.
        frame_dtype: Dtype used to interpret the frame bytes.
        metadata_shm_names: Names of the shared-memory blocks holding
            metadata.
        frame_shm_names: Names of the shared-memory blocks holding frames.
        free_buffer_queue: Queue of buffer-slot indices that can be reused.
        ready_buffer_queue: Queue used to notify the main process that a
            buffer has been populated.
        stop_event: Event used to request early termination.
    """
    metadata_shms: list[shared_memory.SharedMemory] = [shared_memory.SharedMemory(name=name) for name in metadata_shm_names]
    frame_shms: list[shared_memory.SharedMemory] = [shared_memory.SharedMemory(name=name) for name in frame_shm_names]

    try:
        # Next, load in the metadata and frame buffers
        # These NumPy arrays are views onto the shared-memory byte regions.
        # Constructing them is cheap: NumPy is not allocating another full chunk
        # here, it is only wrapping the existing shared-memory allocation with shape
        # and dtype information so we can index and assign into it normally.
        metadata_buffers: list[np.ndarray] = [np.ndarray(metadata_shape, dtype=metadata_dtype, buffer=shm.buf) for shm in metadata_shms]
        frame_buffers: list[np.ndarray] = [np.ndarray(frame_shape, dtype=frame_dtype, buffer=shm.buf) for shm in frame_shms]

        # Stop processing chunks if desired 
        for chunk_num in chunk_indices:
            if(stop_event.is_set()):
                break
            
            # Wait until the main process tells us that one of the 2 shared-memory
            # slots is free to be overwritten. This blocking wait bounds memory usage
            # to exactly the 2 preallocated chunk buffers instead of allowing an
            # unbounded number of chunks to accumulate in RAM.
            buffer_num: int = free_buffer_queue.get()
            metadata_matrix_path, frame_buffer_path = chunk_paths[chunk_num]

            # Open the chunk files as read-only memmaps instead of immediately reading
            # the whole arrays into ordinary process RAM. This is fast to create and
            # cheap on memory because NumPy only reads file headers up front, while
            # the OS pages chunk data in lazily as the copy below touches it.
            # The tradeoff is that these memmaps are file-backed sources only; the
            # consumer still cannot use them directly across processes, so we still
            # need one explicit copy into shared memory.
            metadata: np.ndarray = np.load(metadata_matrix_path, mmap_mode='r')
            frame_buffer: np.ndarray = np.load(frame_buffer_path, mmap_mode='r')

            # Copy the finalized chunk data into the chosen shared-memory slot.
            # This is the main unavoidable full-data copy in the pipeline:
            #     disk-backed memmap -> shared-memory buffer
            # Once this copy finishes, the main process can read the chunk without
            # doing another disk load and without receiving a large pickled object
            # through a multiprocessing queue.
            metadata_buffers[buffer_num][:len(metadata)] = metadata
            frame_buffers[buffer_num][:len(frame_buffer)] = frame_buffer

            # Tell the main process that this buffer slot is ready immediately after
            # the copy completes. Importantly, the queue message contains only small
            # bookkeeping data: which buffer is ready and how much of it is valid.
            # We do not send the arrays themselves through the queue, which avoids a
            # very expensive serialization/copy step.
            ready_buffer_queue.put((buffer_num, chunk_num, len(metadata), len(frame_buffer)))

        # Signal stop of reading 
        # `None` is a sentinel meaning there will be no further ready buffers.
        ready_buffer_queue.put(None)

    # Release the memory
    finally:
        # Close this process's handles to the shared-memory blocks.
        # This does not free the underlying shared memory itself; the parent process
        # still owns that lifecycle and later unlinks the blocks after processing.
        for shm in metadata_shms + frame_shms:
            shm.close()


"""Convert the world chunks of a recording into a playable video"""
def world_chunks_to_video(recording_path: str,
                          output_path: str,
                          linearize_camera_responsivity: bool=False, 
                          apply_color_weights: bool=False, 
                          apply_floor_ceiling: bool=False,
                          debayer_images: bool=False, 
                          apply_digital_gain: bool=False,
                          apply_fielding_function: bool=False, 
                          fill_missing_frames: bool=False,
                          convert_to_seconds: bool=True,
                          embed_timestamps: bool=False,
                          verbose: bool=False,
                          start_end: tuple[int | float] = (0, float("inf")),
                          floor_ceiling: tuple[int] = (0, 255), 
                          raw_inf_threshold: int | float=250
                        ) -> None: 
    
    ###########################

    # INITIALIZATION AND MEMORY ALLOCATION

    ##########################

    """Reconstruct a playable video from chunked world-camera recordings.

    The function overlaps chunk loading and frame processing with a
    two-buffer shared-memory pipeline. For each chunk it reads timestamps
    and raw Bayer frames, optionally applies digital gain, fielding, and RGB
    scaling, optionally debayers the frames, optionally propagates raw
    floor/ceiling violations into the RGB result, fills timestamp gaps with
    dummy frames when requested, and finally writes the output to an FFV1
    ``.avi`` file.

    Args:
        recording_path: Recording directory containing ``config.pkl`` and
            world chunk files.
        output_path: Destination ``.avi`` path.
        apply_color_weights: Whether to apply
            ``world_util.WORLD_RGB_SCALARS`` to raw Bayer pixels.
        apply_floor_ceiling: Whether to set a debayered pixel to ``floor``
            or ``ceiling`` when any contributing raw Bayer sample crosses
            those thresholds.
        remove_dark_noise: Whether to zero pixels that remain at or below
            ``WORLD_DARK_NOISE`` after clipping.
        debayer_images: Whether to convert Bayer frames to RGB before
            writing.
        apply_digital_gain: Whether to multiply each frame by the stored
            per-frame digital gain.
        apply_fielding_function: Whether to apply the registered fielding
            correction map for the frame size.
        fill_missing_frames: Whether to synthesize black frames when
            timestamp gaps imply dropped frames.
        convert_to_seconds: Whether to convert world timestamps from
            nanoseconds to seconds before timing calculations.
        embed_timestamps: Reserved for byte-level timestamp embedding. The
            current implementation still raises ``NotImplementedError`` if
            this is requested.
        verbose: Whether to show progress bars.
        start_end: Inclusive/exclusive chunk-index range to render.
        ceiling: Upper clipping threshold used with
            ``apply_floor_ceiling``.
        floor: Lower clipping threshold used with
            ``apply_floor_ceiling``.
    """
    if(fill_missing_frames and not convert_to_seconds):
        raise Exception("Haven't handled this case. This will generate millions of extra frames")

    assert output_path.endswith(".avi"), f"Output path: {output_path} must end in .avi"

    # Make the directories to the output path if they do not exist already 
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

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
                                   isColor= debayer_images is True 
                                  ) 
    if(not video_writer.isOpened()):
        video_writer.release()
        raise RuntimeError(f"Failed to open VideoWriter for {output_path}")

    # Initialize a dummy frame to fill in missing frames if desired 
    dummy_frame: np.ndarray = np.squeeze( np.full((frame_height, frame_width, 3 if debayer_images is True else 1), 0, dtype=np.uint8))

    # Retrive the sorted chunks for this sensor
    chunks_paths: dict[str, list[tuple]] = group_sensors_files(recording_path)['W']
    try:
        assert len(chunks_paths) > 0, f"No chunks found in: {recording_path}"
    except Exception as e:
        print(e)
        video_writer.release()
        return 

    # Iterate over the chunks for this sensor 
    previous_chunk_end_time: float = 0 
    current_chunk_start_time: float = 0 

    # Define the iterator we will use to iterate over the chunks 
    start_chunk: int = start_end[0]
    end_chunk: int = start_end[1] if start_end[1] != float("inf") else len(chunks_paths)
    chunk_iterator: Iterable = range(start_chunk, end_chunk) if verbose is False else tqdm(range(start_chunk, end_chunk), desc="Processing chunks", leave=True)

    # We should gather some information about the chunk sizes here. 
    # We know that all chunks before the last chunk will be the same size, 
    # with the last chunk potentially being smaller than the others. 
    # Therefore, we can allocate arrays based on the size of the first chunk, 
    # then just access slices of this if we need to for other things 
    # For instance, we can do this with the final 
    # Use read-only memmaps here because we only need lightweight metadata about the
    # first chunk: shape, dtype, and total byte size. This avoids eagerly loading an
    # entire first chunk into normal process RAM just to size the shared-memory
    # buffers. The memmap object is cheap to create and does not read the full array
    # contents unless they are actually accessed.
    first_metadata_matrix_path, first_frame_buffer_path = chunks_paths[start_chunk]
    first_metadata: np.ndarray = np.load(first_metadata_matrix_path, mmap_mode='r')
    first_frame_buffer: np.ndarray = np.load(first_frame_buffer_path, mmap_mode='r')
    

    # Next, we should allocate some shared memory buffers. 
    # We will use this to read in chunks while another chunk is being processed 
    # We should have 2 buffers, one to be currently processed and another to hold 
    # the data queued to be processed next. 
    # This is the main fixed memory cost of the overlap strategy: enough shared memory
    # for 2 full metadata chunks and 2 full frame chunks. In exchange, the loader can
    # read the next chunk while the main process transforms/writes the current one.
    metadata_shms: list[shared_memory.SharedMemory] = [shared_memory.SharedMemory(create=True, size=first_metadata.nbytes) for _ in range(2)]
    frame_shms: list[shared_memory.SharedMemory] = [shared_memory.SharedMemory(create=True, size=first_frame_buffer.nbytes) for _ in range(2)]
    free_buffer_queue: mp.Queue = mp.Queue()
    ready_buffer_queue: mp.Queue = mp.Queue()
    stop_event: object = mp.Event()
    for buffer_num in range(2):
        # Initially both buffer slots are free, so seed the loader with both ids.
        # This lets it start filling the first available slot immediately.
        free_buffer_queue.put(buffer_num)

    # Launch the background loader process. After this starts, disk I/O for later
    # chunks can overlap with CPU work in the main process. The extra process has
    # some overhead, but it is much smaller than repeatedly blocking the main loop on
    # chunk reads or shipping large arrays through multiprocessing queues.
    loader_process: mp.Process = mp.Process(target=_load_world_chunk_shared_memory,
                                            args=(chunks_paths,
                                                  list(range(start_chunk, end_chunk)),
                                                  first_metadata.shape,
                                                  first_metadata.dtype,
                                                  first_frame_buffer.shape,
                                                  first_frame_buffer.dtype,
                                                  [shm.name for shm in metadata_shms],
                                                  [shm.name for shm in frame_shms],
                                                  free_buffer_queue,
                                                  ready_buffer_queue,
                                                  stop_event
                                                 )
                                           )
    loader_process.start()


    # Declare some global video pieces of knowledge that we 
    # will fill in from the first chunk, 
    # such as the provenance map of debayered pixels to their contributing pixels in the RAW image 
    # and the RGB bayer pixel locations
    bayer_RGB_mask: np.ndarray | None = None
    bayer_RGB_pixel_locations: list[np.ndarray] | None = None
    debayered_provenance_map: np.ndarray | None = None
    try:
        for _ in chunk_iterator:        
            ###########################

            # LOADING AND TRANSFORMATION

            ##########################

            # Here, we should instead receive a processed chunk from the shared memory 
            # that is populated from the other process reading in the chunks. 
            #       
            # Wait for the next finalized shared-memory slot from the loader process.
            # This queue receive is cheap because only a tiny tuple is transferred,
            # not the chunk arrays themselves.
            ready_chunk: None = ready_buffer_queue.get()
            if(ready_chunk is None):
                break
            buffer_num, chunk_num, metadata_len, frame_buffer_len = ready_chunk

            # First, we will retrieve the metadata for this buffer 
            # Reconstruct a NumPy array view directly over the shared-memory bytes and
            # slice down to the valid region for this chunk. This is effectively
            # zero-copy at retrieval time: we are not duplicating the shared-memory
            # contents, only creating a local ndarray wrapper around them.
            metadata: np.ndarray = np.ndarray(first_metadata.shape, dtype=first_metadata.dtype, buffer=metadata_shms[buffer_num].buf)[:metadata_len]
            
            # Then, we will read in the timestamps and frame buffer for this chunk 
            # World timestamps are in nanoseconds and thus must be converted to seconds 
            # This line does allocate a fresh contiguous float64 timestamp vector,
            # because we are converting units and standardizing dtype for downstream
            # math. That cost is usually small compared with the frame-buffer copy and
            # frame processing work.
            t_vector: np.ndarray = np.ascontiguousarray(metadata[:, 0], dtype=np.float64) if convert_to_seconds is False else np.ascontiguousarray(metadata[:, 0], dtype=np.float64) / (10 ** 9) 
            # Like `metadata`, this is a direct view into the ready shared-memory slot.
            # No disk I/O occurs here and no additional full chunk copy is made yet.
            frame_buffer: np.ndarray = np.ndarray(first_frame_buffer.shape, dtype=first_frame_buffer.dtype, buffer=frame_shms[buffer_num].buf)[:frame_buffer_len]
            
            # Let's keep track of the current colorspace the video is in 
            current_color_space: Literal['BAYER', 'RGB', 'BGR'] = 'BAYER'

            # Assert the t vector and frame vector are the same size 
            try:
                assert(len(t_vector) == len(frame_buffer))
            except Exception as e:
                print(e)
                video_writer.release()
                return 

            # Ensure the frame size is consistent 
            buffer_frame_shape: tuple[int] = tuple(frame_buffer.shape[1:3])
            try:
                assert buffer_frame_shape == (frame_height, frame_width), f"Frame shape unequal after initialization: {(frame_height, frame_width)}"
            except Exception as e:
                print(e)
                video_writer.release() 
                return

            # Assert the frames are of appropriate dtype 
            try:
                assert(frame_buffer.dtype == np.uint8)
            except Exception as e:
                print(e)
                video_writer.release()
                return 

            # If we somehow got an empty chunk, skip it 
            if(len(t_vector) == 0): 
                free_buffer_queue.put(buffer_num)
                continue

            # Now, we will apply the requested transformations to the frames 
            # First, convert to float for float multiplications
            # if not in float already
            frame_buffer = frame_buffer.astype(np.float64, copy=False) # copy = False means reuse the same memory if possible (e.g. if it is already float, which it should be)  
            raw_saturated_pixel_mask: np.ndarray | None = None

            # If we want to apply the floor / ceiling procedure, 
            # we need to mark saturated raw pixels as INF. Raw zero pixels
            # intentionally remain zero until debayering so they can propagate
            # to neighboring debayered pixels via the provenance map.
            if(apply_floor_ceiling is True):
                floor, ceiling = floor_ceiling

                raw_saturated_pixel_mask = frame_buffer >= raw_inf_threshold
                frame_buffer[raw_saturated_pixel_mask] = np.inf

            # Linearize the recording data to account for the full well effect of the camera
            if(linearize_camera_responsivity is True):
                assert frame_buffer.dtype == np.float64 and frame_buffer.ndim == 3, f"Frame buffer dtype must be float 64 and bayer format to linearize the camera responsivity. Current dtype is: {frame_buffer.dtype} and ndim is {frame_buffer.ndim}"
                world_util.linearize_camera_responsivity(frame_buffer, int(config_data["sensors"]["W"]["sensor_mode"]["bit_depth"]), dark_noise=world_util.WORLD_DARK_NOISE, dst=frame_buffer)
                
                # We need to reapply the saturated pixel mask because the above operation performed clipping thereby losing our INF values
                if(raw_saturated_pixel_mask is not None):
                    frame_buffer[raw_saturated_pixel_mask] = np.inf

            # Apply the fielding function of the camera that we measured in the planetarium
            if(apply_fielding_function is True):
                assert frame_buffer.dtype == np.float64 and frame_buffer.ndim == 3, f"Frame buffer dtype must be float64 and in bayer format to apply fielding function. Current dtype is {frame_buffer.dtype}, current ndim: {frame_buffer.ndim}"
                
                # This is an in-place operation, no new memory is allocated here
                world_util.apply_fielding_function(frame_buffer)

            # Apply the radiometric correction weights we calculated outdoors with respect to the PR670
            if(apply_color_weights is True):
                assert frame_buffer.ndim == 3 and frame_buffer.dtype == np.float64, f"To apply color weights, buffer must be in bayer form np.float64."

                # Apply the radiometric correction
                # pre-compute the RGB bayer mask once to save time 
                if(bayer_RGB_mask is None):
                    bayer_RGB_mask = world_util.generate_RGB_mask(frame_buffer[0])
                    bayer_RGB_pixel_locations = [ np.argwhere(bayer_RGB_mask == color) for color in "RGB" ]
                assert bayer_RGB_mask.shape == frame_buffer.shape[1:], f"Bayer RGB mask shape: {bayer_RGB_mask.shape} not equal to frame shape: {frame_buffer.shape[1:]}" 

                # Apply the color corrections IN-PLACE 
                assert len(bayer_RGB_pixel_locations) == len(world_util.WORLD_RGB_SCALARS)
                world_util.apply_color_correction(frame_buffer, bayer_RGB_pixel_locations)

            # Apply Dgain if requested
            if(apply_digital_gain is True):
                assert frame_buffer.dtype == np.float64, f"Frame buffer dtype must be float64 to apply dgain. Current dtype is {frame_buffer.dtype}"

                # Extract the dgains of the buffers 
                buffer_dgains: np.ndarray = metadata[:, 2]
                
                # Apply the digital gain IN PLACE
                world_util.apply_digital_gain(frame_buffer, buffer_dgains)

               
            # If we want to debayer the image
            raw_frame_buffer: np.ndarray | None = None
            if(debayer_images is True):
                # Ensure the pixels that were saturated at the start remain saturated
                if(raw_saturated_pixel_mask is not None):
                    frame_buffer[raw_saturated_pixel_mask] = np.inf

                # Preserve the float Bayer image (with Inf ceiling / 0 floor
                # markers) before de-Bayering so the floor/ceiling propagation
                # below can consult the true raw contributors.
                raw_frame_buffer: np.ndarray = frame_buffer

                # OpenCV's debayer needs finite uint16 input. Scale the float
                # Bayer values up first so their fractional precision survives the
                # round-trip through uint16, map Inf to the uint16 ceiling
                # sentinel, then clip / round / convert. Raw zeros remain zero.
                debayer_scale: int = 100
                debayer_ceiling_marker: int = 2 ** 16 - 1
                debayer_saturation_threshold: float = debayer_ceiling_marker / 3

                debayer_input: np.ndarray = frame_buffer * debayer_scale
                
                # We need to replace INF with a sentinel value as debayering only accepts uint16
                debayer_input[np.isposinf(debayer_input)] = debayer_ceiling_marker
                debayer_input = np.round(np.clip(debayer_input, 0, debayer_ceiling_marker)).astype(np.uint16)
                assert current_color_space == "BAYER" and debayer_input.dtype == np.uint16, f"Frame buffer must be in BAYER space and uint16 to debayer. Current format is: {current_color_space} | {debayer_input.dtype}"

                # Debayer, then re-mark any output that a saturated contributor
                # pulled up as Inf, and scale back down to the original value
                # range. Mirrors the notebook's de-Bayer stage.
                frame_buffer = world_util.debayer(debayer_input).astype(np.float64)
                
                # If we plan on doing more processing (that is, the floor ceiling saturated operation)
                # we need to mark all the saturated pixels as INF
                if(apply_floor_ceiling is True):
                    frame_buffer[frame_buffer >= debayer_saturation_threshold] = debayer_ceiling_marker
                    frame_buffer[frame_buffer == debayer_ceiling_marker] = np.inf
                
                frame_buffer = frame_buffer / debayer_scale
                current_color_space = "RGB"
            
            # If we want to apply floor ceiling [LOW, HIGH] that implies 
            # that any output pixel in the debayered image will be either blacked or whited out 
            # if any of its contributing pixels were floor or ceiling
            if(apply_floor_ceiling is True):
                assert frame_buffer.ndim == 4, f"Frames must be debayered before applying floor/ceiling"

                if(debayered_provenance_map is None):
                    debayered_provenance_map = world_util.generate_debayered_provenance_map(frame_buffer[0])

                # Apply the floor / ceiling algorithm IN PLACE. raw_frame_buffer
                # is the float Bayer image, so Inf is the ceiling sentinel.
                world_util.apply_floor_ceiling(
                    frame_buffer,
                    raw_frame_buffer,
                    debayered_provenance_map,
                    (floor_ceiling[0], np.inf),
                    floor_ceiling,
                )

                # Any pixel still marked Inf was flagged as saturated by the
                # debayer saturation threshold but had no floor/ceiling raw
                # contributor to resolve it; treat it as saturated -> ceiling.
                frame_buffer[np.isposinf(frame_buffer)] = floor_ceiling[1]

            # Convert the frame buffer back into uint8 format
            if(frame_buffer.dtype != np.uint8): 
                if(np.issubdtype(frame_buffer.dtype, np.floating)):
                    assert np.all(np.isfinite(frame_buffer)), (
                        "Frame buffer contains NaN or Inf values immediately "
                        "before uint8 conversion. These should have already "
                        "been filtered out."
                    )
                frame_buffer = np.clip(np.round(frame_buffer), 0, 255).astype(np.uint8)

            # ALL TRANSFORMATIONS ARE DONE AT THIS POINT. 
            # WE NEED TO CONVERT TO BGR TO WRITE FRAMES WITH CV2
            if(frame_buffer.ndim > 3 and current_color_space != 'BGR'):
                if(current_color_space == 'RGB'):
                    # Allocate a new frame buffer to hold the BGR frames 
                    new_buffer: np.ndarray = np.empty_like(frame_buffer)
                    for frame_idx in range(len(frame_buffer)):
                        cv2.cvtColor(frame_buffer[frame_idx], cv2.COLOR_RGB2BGR, dst=new_buffer[frame_idx])
                    frame_buffer = new_buffer
                
                else:
                    raise RuntimeError(f"Unsupported colorspace conversion: {current_color_space} -> RBG")

            assert frame_buffer.dtype == np.uint8, f"Frame buffer of type: {frame_buffer.dtype} must be of dtype np.uint8 before writing."


            ###############
            #    WRITING SECTION    
            ##############
            # If we want to embed the timestamps in the frame buffer, do that now
            if(embed_timestamps is True):
                raise NotImplementedError()

            # Retrieve the current chunk start time 
            current_chunk_start_time = t_vector[0]

            # Before we write the frames for this chunk, we must write the dummy frame (if desired)
            # for those frames between the previous chunk and the current chunk 
            
            # Find the total elapsed time in seconds between the two chunks 
            time_between_chunks: float = current_chunk_start_time - previous_chunk_end_time
            assert time_between_chunks >= 0, f"Time btween chunks is somehow negative, Chunk 1: {previous_chunk_end_time}, Chunk 2: {current_chunk_start_time}"

            # Calculate the number of missed frames as the elapsed time divided 
            # by the seconds per frame, minus one frame as we have to count the current frame 
            # as captured during this interval 
            missed_frames: int = 0 if chunk_num == start_end[0] or fill_missing_frames is False else int( (time_between_chunks / (1/FPS ) ) -  1) 
            if(missed_frames > 5 * FPS):
                warnings.warn(f"Missed {missed_frames} between chunks. This may be unnaturally large")
            missing_timestamps: np.ndarray = np.linspace(previous_chunk_end_time + (1/FPS), current_chunk_start_time, missed_frames, endpoint=False)

            # Write the number of missing frames in between as the previous frame 
            for missing_timestamp in missing_timestamps:
                missing_frame: np.ndarray = world_util.embed_timestamp(dummy_frame, missing_timestamp) if embed_timestamps is True else dummy_frame
                assert missing_frame.dtype == np.uint8 and ( len(missing_frame.shape) == 3 if debayer_images is True else len(missing_frame.shape) == 2 ) 
                video_writer.write(missing_frame)

            # Initialize variables to track the previous timestamp 
            # We will use this delta with the current timestamp 
            # to track dropped frames and add dummy frames in the meantime
            previous_timestamp: float = t_vector[0]

            # Write frames to the video 
            frame_iterator: Iterable = range(len(frame_buffer)) if verbose is False else tqdm(range(len(frame_buffer)), desc="Processing frames", leave=False)
            for frame_num in frame_iterator:
                timestamp, frame = t_vector[frame_num], frame_buffer[frame_num]

                # Assert the frame is of the appropriate type after transformation 
                try:
                    assert(frame.dtype == np.uint8 and dummy_frame.dtype == np.uint8)
                except Exception as e:
                    print(e)
                    video_writer.release()
                    return 

                # Otherwise, let's calculate the time delta between the current timestamp and the previous timestamp 
                time_between_frames: float = timestamp - previous_timestamp
                assert time_between_frames >= 0, f"Time between frames is somehow negative. Frame 1: {previous_timestamp}, Frame 2: {timestamp}"

                # Calculate the number of missed frames as the elapsed time divided 
                # by the seconds per frame, minus one frame as we have to count the current frame 
                # as captured during this interval 
                missed_frames: int = 0 if frame_num == 0 or fill_missing_frames is False else int( (time_between_frames / (1/FPS ) ) -  1) 
                if(missed_frames > 5 * FPS):
                    warnings.warn(f"Missed {missed_frames} frames between frames. This may be unaturally large")
                missing_timestamps: np.ndarray = np.linspace(previous_timestamp + (1/FPS), timestamp, missed_frames, endpoint=False)

                # Write the number of missing frames in between as the previous frame 
                for missing_timestamp in missing_timestamps:
                    missing_frame: np.ndarray = world_util.embed_timestamp(dummy_frame, missing_timestamp) if embed_timestamps is True else dummy_frame
                    
                    assert missing_frame.dtype == np.uint8 and ( len(missing_frame.shape) == 3 if debayer_images is True else len(missing_frame.shape) == 2 )
                    video_writer.write(missing_frame)
                
                # Write the current frame 
                assert frame.dtype == np.uint8 and ( len(frame.shape) == 3 if debayer_images is True else len(frame.shape) == 2 )
                video_writer.write(frame)

                # Save the current timestamp as the previous timestamp
                # so we can reference it for the next frame 
                previous_timestamp = timestamp 
            
            # Save the end time of this chunk 
            previous_chunk_end_time = t_vector[-1]

            # Now that we are done processing this buffer, we can 
            # send it back to load another chunk 
            # This handoff point is what prevents the loader from overwriting a buffer
            # while the main process is still reading from it.
            free_buffer_queue.put(buffer_num)

    # In any case of the function, 
    # we want to close the subprocess 
    # and clear the memory we allocated 
    finally:
        stop_event.set()
        
        if(loader_process.is_alive()):
            loader_process.join(timeout=1)
        if(loader_process.is_alive()):
            loader_process.terminate()
            loader_process.join()
        
        for shm in metadata_shms + frame_shms:
            shm.close()
            shm.unlink()

        # Close the video writer 
        video_writer.release()

    return 


def find_events(video_path: str,
                search_text: str,
                min_conf: float = 60.0,
                case_sensitive: bool = False,
                frame_step: int = 5,
                roi: tuple[int, int, int, int] | None = None,
                refine_hits: bool = True,
                refine_radius: int = 5,
                verbose: bool = False,
                show_images: bool=False
               ) -> int | None:
    """Search a video for frames whose OCR output contains target text.

    The search runs in two passes. A coarse pass samples every
    ``frame_step`` frames and uses a fast substring test on Tesseract's full
    text output. When ``refine_hits`` is enabled, the function then revisits
    a small neighborhood around each coarse hit and performs a more precise
    confidence-filtered OCR pass to localize the matching frame more
    accurately.

    Args:
        video_path: Video file to inspect.
        search_text: Text snippet to search for in OCR output.
        min_conf: Minimum Tesseract word confidence used during the
            refinement pass.
        case_sensitive: Whether matching should preserve case.
        frame_step: Stride used for the coarse scan.
        roi: Optional crop window ``(y0, y1, x0, x1)`` applied before OCR.
        refine_hits: Whether to run the dense second pass around coarse
            detections.
        refine_radius: Number of frames on either side of a coarse hit to
            revisit during refinement.
        verbose: Whether to show progress bars and a detection-score plot.
        show_images: When ``verbose`` is also ``True``, display each matched
            frame.

    Returns:
        The last matching frame index, or ``None`` if no frame contains the
        requested text. Returning only the last hit preserves the historical
        behavior of this helper.
    """

    # -------------------------
    # Helper: normalize text
    # -------------------------
    def normalize_text(txt: str) -> str:
        """Normalize OCR output before matching.

        Args:
            txt: Raw text string produced by OCR.

        Returns:
            Text with repeated whitespace collapsed and, unless
            ``case_sensitive`` is ``True``, lowercased.
        """
        txt = re.sub(r"\s+", " ", txt).strip()
        return txt if case_sensitive else txt.lower()

    # -------------------------
    # Helper: preprocess frame
    # -------------------------
    def preprocess_frame(frame: np.ndarray) -> np.ndarray:
        """Crop and binarize a frame before OCR.

        Args:
            frame: Raw BGR frame from OpenCV.

        Returns:
            Preprocessed grayscale image after optional ROI cropping,
            Gaussian blur, and Otsu thresholding.
        """
        if (roi is not None):
            y0, y1, x0, x1 = roi
            frame = frame[y0:y1, x0:x1]

        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        gray = cv2.GaussianBlur(gray, (3, 3), 0)
        _, gray = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

        return gray

    # -------------------------
    # Fast check (binary match)
    # -------------------------
    def frame_contains_text_fast(frame: np.ndarray, target_text: str) -> float:
        """Run a fast binary OCR check on one frame.

        Args:
            frame: Raw BGR frame.
            target_text: String to search for after normalization.

        Returns:
            ``1.0`` if the normalized OCR text contains ``target_text``,
            otherwise ``0.0``.
        """
        gray = preprocess_frame(frame)

        text = pytesseract.image_to_string(gray, config="--psm 6")
        text = normalize_text(text)
        target = normalize_text(target_text)

        return 1.0 if (target in text) else 0.0

    # -------------------------
    # Accurate check (with confidence)
    # -------------------------
    def frame_contains_text_conf(frame: np.ndarray, target_text: str) -> float:
        """Run confidence-aware OCR on one frame.

        Args:
            frame: Raw BGR frame.
            target_text: String to search for after normalization.

        Returns:
            Mean OCR confidence for matched words when the normalized OCR
            text contains ``target_text``, otherwise ``0.0``.
        """
        gray = preprocess_frame(frame)

        ocr_data = pytesseract.image_to_data(
            gray,
            output_type=pytesseract.Output.DICT,
            config="--psm 6"
        )

        words: list[str] = []
        confs: list[float] = []

        for txt, conf in zip(ocr_data["text"], ocr_data["conf"]):
            txt = str(txt).strip()

            try:
                conf = float(conf)
            except ValueError:
                conf = -1.0

            if ((txt != "") and (conf >= min_conf)):
                words.append(txt)
                confs.append(conf)

        if (len(words) == 0):
            return 0.0

        joined_text = normalize_text(" ".join(words))
        target = normalize_text(target_text)

        if (target in joined_text):
            return float(np.mean(confs))

        return 0.0

    # -------------------------
    # Open video
    # -------------------------
    video_stream = cv2.VideoCapture(video_path)

    if (not video_stream.isOpened()):
        raise RuntimeError(f"Error: Could not open video @: {video_path}")

    num_frames: int = int(video_stream.get(cv2.CAP_PROP_FRAME_COUNT))

    # Array storing coarse detection results
    coarse_scores: np.ndarray = np.zeros(num_frames, dtype=np.float64)

    # Determine which frames to sample
    coarse_frame_nums = range(0, num_frames, frame_step)

    if (verbose is True):
        coarse_frame_nums = tqdm(coarse_frame_nums, desc="Coarse OCR")

    # -------------------------
    # Coarse pass
    # -------------------------
    for frame_num in coarse_frame_nums:
        video_stream.set(cv2.CAP_PROP_POS_FRAMES, frame_num)

        ret, frame = video_stream.read()

        if (not ret):
            continue

        coarse_scores[frame_num] = frame_contains_text_fast(frame, search_text)

    video_stream.release()

    # Frames where coarse detection found text
    coarse_hits: np.ndarray = np.flatnonzero(coarse_scores > 0)

    # -------------------------
    # Refinement pass
    # -------------------------
    final_scores: np.ndarray = np.zeros(num_frames, dtype=np.float64)

    if ((refine_hits is True) and (len(coarse_hits) > 0)):

        # Build set of frames around coarse hits
        candidate_frames: set[int] = set()

        for hit in coarse_hits:
            start = max(0, hit - refine_radius)
            end = min(num_frames, hit + refine_radius + 1)
            candidate_frames.update(range(start, end))

        candidate_frames = sorted(candidate_frames)

        video_stream = cv2.VideoCapture(video_path)

        if (not video_stream.isOpened()):
            raise RuntimeError(f"Error: Could not reopen video @: {video_path}")

        frame_iter = candidate_frames

        if (verbose is True):
            frame_iter = tqdm(frame_iter, desc="Refining OCR")

        for frame_num in frame_iter:
            video_stream.set(cv2.CAP_PROP_POS_FRAMES, frame_num)

            ret, frame = video_stream.read()

            if (not ret):
                continue

            final_scores[frame_num] = frame_contains_text_conf(frame, search_text)

        video_stream.release()

    else:
        # If no refinement, just use coarse results
        final_scores = coarse_scores.copy()

    # -------------------------
    # Extract matched frames
    # -------------------------
    matched_frame_nums: np.ndarray = np.flatnonzero(final_scores > 0)

    # -------------------------
    # Visualization (optional)
    # -------------------------
    if (verbose is True):

        plt.figure(figsize=(12, 4))
        plt.plot(np.arange(len(final_scores)), final_scores, label=f'"{search_text}" score')

        plt.scatter(
            matched_frame_nums,
            final_scores[matched_frame_nums],
            marker="x",
            s=80,
            color="red",
            label="Detected text"
        )

        plt.xlabel("Frame Number")
        plt.ylabel("OCR Match Score")
        plt.title("Detected text plotted")
        plt.grid(True)
        plt.legend()
        plt.show()

        if(show_images is True):
            # Re-open video to display detected frames
            video_stream = cv2.VideoCapture(video_path)

            if (not video_stream.isOpened()):
                raise RuntimeError(f"Error: Could not reopen video @: {video_path}")

            for frame_num in matched_frame_nums:
                video_stream.set(cv2.CAP_PROP_POS_FRAMES, int(frame_num))

                ret, frame = video_stream.read()

                if (not ret):
                    continue

                plt.figure(figsize=(8, 6))
                plt.imshow(cv2.cvtColor(frame, cv2.COLOR_BGR2RGB))
                plt.title(
                    f'Frame: {frame_num} | Text: "{search_text}" | Score: {final_scores[frame_num]:.2f}'
                )
                plt.axis("off")
                plt.show()

            video_stream.release()

    # IMPORTANT: preserve original behavior (return last match)
    return matched_frame_nums[-1] if len(matched_frame_nums) > 0 else None


def convert_avi_to_mp4(input_path: str, output_path: str) -> None:
    """Transcode an AVI file to H.264/AAC MP4 using ``ffmpeg``.

    Args:
        input_path: Source AVI path.
        output_path: Destination MP4 path.
    """
    cmd = [
        "ffmpeg",
        "-y",                     # force overwrite
        "-i", input_path,
        "-c:v", "libx264",
        "-crf", "18",
        "-preset", "slow",
        "-c:a", "aac",
        "-b:a", "192k",
        output_path
    ]

    subprocess.run(cmd, check=True)


def main():
    """Run the command-line entry point."""
    return 

if(__name__ == "__main__"):
    main()
