import numpy as np 
import cv2
import warnings
import time
import queue
from natsort import natsorted
from typing import Iterable
import os

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
    assert os.path.exists(video_path), "Video path does not exist"

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
    assert os.path.exists(video_path), "Video path does not exist"

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
    assert os.path.exists(video_path), "Video path does not exist"

    # Open the video via cv2 
    video_stream: cv2.VideoCapture = cv2.VideoCapture(video_path)

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
    assert os.path.exists(video_path), "Video path does not exist"
    
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

def main():
    return 

if(__name__ == "__main__"):
    main()