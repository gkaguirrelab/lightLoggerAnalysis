import numpy as np
import os 
import sys 

# Import relevant custom libraries with helper functions and constants 
light_logger_dir_path: str = os.path.expanduser("~/Documents/MATLAB/projects/lightLogger")
Pi_util_dir_path: str = os.path.join(light_logger_dir_path, "raspberry_pi_firmware", "utility")
sys.path.append(Pi_util_dir_path)

"""Generate the geometric transform used for virtual foveation"""
def calculate_perspective_transform_w2e() -> object:

    return 

"""Perspective transform an individual frame"""
def perspective_transform_w2e(frame: np.ndarray, transformation: object | None=None) -> np.ndarray:
    # If not transform passed in, calculate it 
    if(transformation is None):
        transformation = calculate_perspective_transform_w2e() 



    return 

def main():
    pass 

if(__name__ == "__main__"):
    pass 