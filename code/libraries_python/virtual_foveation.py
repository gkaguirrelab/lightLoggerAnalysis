import numpy as np
import os 
import sys 
import scipy.io 
import multiprocessing as mp
from scipy.optimize import brentq

# Import relevant custom libraries with helper functions and constants 
light_logger_dir_path: str = os.path.expanduser("~/Documents/MATLAB/projects/lightLogger")
Pi_util_dir_path: str = os.path.join(light_logger_dir_path, "raspberry_pi_firmware", "utility")
sys.path.append(Pi_util_dir_path)

import Pi_util 
import extract_eye_features
from scipy.interpolate import griddata
import matlab.engine

"""Load the transformation matrix from the saved MATLAB file"""
def load_transformation_matrix(transformation_filepath: str, 
                               matlab_engine: object | None=None
                               ) -> np.ndarray:

    # initialize the MATLAB engine if not done so already 
    matlab_pre_initialized: bool = matlab_engine is not None
    if(not matlab_pre_initialized):
        matlab_engine = matlab.engine.start_matlab()

    # Ask MATLAB to return the matrix
    matlab_engine.eval(f"data = load('{transformation_filepath}');", nargout=0)
    transformation_matrix: np.ndarray = np.array(matlab_engine.eval("data.perspective_transform.fit.geometric_transform.A", nargout=1), dtype=np.float64)

    # If we started it for just this funciton, close it 
    if(not matlab_pre_initialized):
        matlab_engine.close() 

    return transformation_matrix

"""Load the intrinsics struct from MATLAB"""
def load_intrinsics(intrinsics_filepath: str, 
                    matlab_engine: object | None=None
                    ) -> dict:
    # initialize the MATLAB engine if not done so already 
    matlab_pre_initialized: bool = matlab_engine is not None
    if(not matlab_pre_initialized):
        matlab_engine = matlab.engine.start_matlab()

    # Initialize our return dictionary 
    intrinsics_dict: dict[str, np.ndarray] = {"mapping_coefficients": None, 
                                              "image_size": None,
                                              "distortion_center": None,
                                              "stretch_matrix": None
                                             }   

    # Ask MATLAB to load in the file so we can extract some data from it
    matlab_engine.eval(f"data = load('{intrinsics_filepath}');", nargout=0)     
    matlab_engine.eval(f"intrinsics_struct = data.camera_intrinsics_calibration.results.Intrinsics;", nargout=0)

    # Exact the data from the MATLAB engine 
    for field in intrinsics_dict:
        # Convert the field name to MATLAB standard 
        matlab_field_name: str = "".join([ word.capitalize() for word in field.split("_") ])

        # Extract the value from matlab and convert to numpy array 
        value: np.ndarray = np.array(matlab_engine.eval(f"intrinsics_struct.{matlab_field_name};", nargout=1), dtype=np.float64)

        # Save the value for this field 
        intrinsics_dict[field] = value.flatten() if field != "stretch_matrix" else value

    # If we started it for just this funciton, close it 
    if(not matlab_pre_initialized):
        matlab_engine.close() 

    return intrinsics_dict

"""Given a path to a video, extract the perimeter dict from it
   That is, the perimeter points of the pupil for each frame 
   massaged into the MATLAB struct format required to use Geoff's 
   gaze angle analysis code. Default keypoint threshold is -1 
   so that filtering is done in Geoff's code 
"""
def generate_perimeter_dict(path_to_video: str, 
                            output_path: str | None=None,
                            is_grayscale: bool=False, 
                            visualize_results: bool=False, 
                            safe_execution: bool=True, 
                            keypoint_threshold: float=-1
                            ) -> dict:
    # Analysis the video with pylids to calculate the perimeter dict 
    _, perimeter_dict = extract_eye_features.extract_pupil_features(path_to_video, 
                                                                    is_grayscale=is_grayscale, 
                                                                    visualize_results=visualize_results, 
                                                                    method='pylids',
                                                                    safe_execution=safe_execution,
                                                                    keypoint_threshold=keypoint_threshold
                                                                  )


    # If we choose to output to a file, simply return 
    if(output_path is None):
        return perimeter_dict

    # Otherwise, output to a file 
    scipy.io.savemat(output_path, {"perimeter": perimeter_dict})

    return perimeter_dict

"""Convet a frame from pixel coordinates to visual angle coordinates"""
def convert_sensor_to_angles(sensor_coordinates: np.ndarray, 
                             center_offset: tuple[float],
                             camera_intrinsics: dict
                            ) -> np.ndarray:
    
    # First, extract the distortion center of the image 
    dist_cx, dist_cy = camera_intrinsics["distortion_center"]

    # Shift the image coordinates according to the distortion center 
    xC: np.ndarray = sensor_coordinates[:, 0] - dist_cx
    yC: np.ndarray = sensor_coordinates[:, 1] - dist_cy

    # Extract the tradius of distortion and mapping coefficents of 
    # the intrincis 
    r: np.ndarray = np.sqrt( xC ** 2 + yC ** 2 )
    coeffs: np.ndarray = camera_intrinsics["mapping_coefficients"]
    a1, a2, a3, a4 = coeffs


    # Create an array to store the solved polar angles (theta) for each pixel.
    # Initialize with NaNs so you can detect where the solver fails.
    theta: np.ndarray = np.full_like(r, np.nan, dtype=float)

    # Loop through each radius value r[k] and solve for theta
    for k in range(len(r)):
        # Define the fisheye mapping polynomial minus the observed radius.
        # This function equals zero when the polynomial maps theta -> r[k].
        func = lambda t: a1*t + a2*t**3 + a3*t**4 + a4*t**5 - r[k]
        try:
            # Numerically find the root of func(t) within [0, pi].
            # This finds the polar angle theta that maps to radius r[k].
            theta[k] = brentq(func, 0, np.pi)
        except ValueError:
            # If no valid root was found (no sign change in the interval),
            # leave this entry as NaN.
            theta[k] = np.nan

    # Compute the azimuthal angle phi for each pixel (direction around optical axis)
    phi: np.ndarray = np.arctan2(yC, xC)

    # Assume unit sphere radius (we only care about direction)
    R: float = 1.0

    # Convert from spherical coordinates (R, theta, phi) to Cartesian (X, Y, Z)
    # These represent 3D directions of rays through each sensor pixel.
    X: np.ndarray = R * np.sin(theta) * np.cos(phi)
    Y: np.ndarray = R * np.sin(theta) * np.sin(phi)
    Z: np.ndarray = R * np.cos(theta)

    # Compute azimuth (horizontal angle) in degrees and apply the horizontal offset
    azi: np.ndarray = np.degrees(np.arctan2(Y, Z)) + center_offset[0]

    # Compute elevation (vertical angle) in degrees and apply the vertical offset
    ele: np.ndarray = np.degrees(np.arctan2(X, Z)) + center_offset[1]

    return np.column_stack([azi, ele])

"""Given a single frame, virtually foveate it"""
def perspective_transform_frame(frame: np.ndarray, 
                                center_offset: tuple[float],
                                camera_intrinsics: dict[str, np.ndarray], 
                                transformation: np.ndarray,
                                degrees_eccentricity: float=60 
                                ) -> np.ndarray:

    # First, make a meshgrid 
    n_rows, n_cols = frame.shape[:2]   

    # Create meshgrid of pixel coordinates 
    xg, yg = np.meshgrid(np.arange(1, n_cols + 1),  # 1..nCols #(added this here to keep it consistent from matlab)
                        np.arange(1, n_rows + 1))  # 1..nRows

    # Flatten and stack into N×2 array of [x, y] coordinates
    sensor_coordinates = np.column_stack((xg.ravel(), yg.ravel()))

    # Get the camera visual field positions corresponding to positions of all
    # locations on the camera sensor
    gaze_angle_coordinates: np.ndarray = convert_sensor_to_angles(sensor_coordinates, center_offset, camera_intrinsics)
    gaze_angle_coordinates = np.hstack([gaze_angle_coordinates, np.ones((len(gaze_angle_coordinates), 1))]) 

    # Transform the camera visual field points to eye rotation coordinate space
    # using the previously calculated tform


    # Apply the projective transformation
    transformed_h: np.ndarray = gaze_angle_coordinates @ transformation.T    # matrix multiply by A transpose (same as MATLAB forward transform)

    # Normalize by the third (homogeneous) coordinate
    # The result is equivalent to MATLAB's transformPointsForward output
    eye_rotation_coordinates: np.ndarray = transformed_h[:, :2] / transformed_h[:, [2]]

    # --- 1. Compute vector norm of each 2D coordinate and mask beyond 60 degrees ---
    # Equivalent to MATLAB's vecnorm(eyeRotationCoordinates, 2, 2)
    norms: np.ndarray = np.linalg.norm(eye_rotation_coordinates, axis=1)  # Euclidean norm per point
    idx = norms > degrees_eccentricity  # Boolean mask for points about 60 degrees of eccentricity

    # --- 2. Copy the image and set those pixels to NaN ---
    sub_image: np.ndarray = frame.copy().astype(np.float64)  # ensure we can assign NaN
    sub_image = sub_image.ravel()            # flatten to align with eyeRotationCoordinates
    sub_image[idx] = np.nan
    sub_image = sub_image.reshape(sub_image.shape)   # reshape back to 2D

    # --- 3. Reshape transformed coordinates to image grid ---
    virtually_foveated_X: np.ndarray = eye_rotation_coordinates[:, 0].reshape(n_rows, n_cols)
    virtually_foveated_Y: np.ndarray = eye_rotation_coordinates[:, 1].reshape(n_rows, n_cols)

    # Create evenly spaced query grid along X and Y
    xq: np.ndarray = np.linspace(np.nanmin(virtually_foveated_X), np.nanmax(virtually_foveated_X), n_cols)
    yq: np.ndarray = np.linspace(np.nanmin(virtually_foveated_Y), np.nanmax(virtually_foveated_Y), n_rows)
    Xq, Yq = np.meshgrid(xq, yq)

    # Flatten the inputs like MATLAB's (:)
    points: np.ndarray = np.column_stack((virtually_foveated_X.ravel(), virtually_foveated_Y.ravel()))
    values: np.ndarray = sub_image.ravel()

    # Interpolate onto the new grid
    Vq: np.ndarray = griddata(points, values, (Xq, Yq), method='linear')

    # Finalize the transformation 
    transformed_frame: np.ndarray = np.flipud(np.nan_to_num(Vq, 0))

    return transformed_frame

def main():
    pass 

if(__name__ == "__main__"):
    pass 