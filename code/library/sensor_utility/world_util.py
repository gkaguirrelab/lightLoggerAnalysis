import numpy as np
import matplotlib.pyplot as plt
import cv2
import os 
import sys
from typing import Literal
import pandas as pd
import mat73
import matlab

# Store the scalar multipliers for all of the different colors of pixel's in an image 
# required to equalize the color to the R pixel values in that frame.
# We are currently not using this, but calculated it off of a 10 minute 
# video in a variety of settings (inside (office) outside (street, park))
# Just something to note that these are inherently unequal
WORLD_RGB_SCALARS: np.ndarray = np.array([1, 0.8223775, 0.95367937], dtype=np.float64)

# Define a mapping between frame sizes and fielding functions of the camera
WORLD_FIELDING_FUNCTIONS: dict[tuple[int], np.ndarray] = {(480, 640): np.ones((480, 640), dtype=np.float64)}

"""Debayer a world camera image into an RGB image"""
def debayer_image(image: np.ndarray, visualize_results: bool=False) -> np.ndarray | tuple[np.ndarray, object]:
    # Initialize a variable for the figure handle that 
    # will be used to visualize (if desired)
    fig: object | None = None

    # Ensure the image is correctly rounded 
    guarded_image: np.ndarray = np.clip(np.round(image), 0, 255).astype(np.uint8)

    # Debayer the image according to the sensor's bayer pattern
    debayered_image: np.ndarray = cv2.cvtColor(guarded_image, cv2.COLOR_BayerRG2RGB)

    # Visualize the results if desired
    if(visualize_results is True):
        # Initialize a figure with two axes 
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Debayering (Before / After)", fontweight='bold', fontsize=18)
        
        # Left axis will be the "before"
        # image 
        axes[0].imshow(guarded_image, cmap="gray")
        axes[0].set_title('Before')

        # Middle axis will be the after but without 
        # color correction 
        axes[1].imshow(debayered_image)
        axes[1].set_title("After")

        # Show the plot 
        plt.show() 

        return debayered_image, fig
        
    return debayered_image 


"""Apply the fielding function to a frame"""
def apply_fielding_function(original_frame: np.ndarray, visualize_results: bool=False) -> tuple[np.ndarray, object] | np.ndarray:
    # Initialize a variable for the figure handle that 
    # will be used to visualize (if desired)
    fig: object | None = None

    # Retrieve the fielding function for the current camera configuration
    fielding_function: np.ndarray = WORLD_FIELDING_FUNCTIONS[original_frame.shape]
    assert original_frame.shape == fielding_function.shape, f"Fielding function shape: {fielding_function.shape} is not equal to frame shape: {original_frame.shape}"

    # Apply the fielding function 
    modified_image: np.ndarray = np.clip(np.round(fielding_function.astype(np.float64) * original_frame.astype(np.float64)), 0, 255).astype(np.uint8)

    # Visualize the results if desired 
    if(visualize_results is True):
        # Initialize a figure with two axes 
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Fielding function (Before / After)", fontweight='bold', fontsize=18)
        
        # Left axis will be the "before"
        # image 
        axes[0].imshow(original_frame, cmap="gray")
        axes[0].set_title('Before')

        # Middle axis will be the after but without 
        # color correction 
        axes[1].imshow(modified_image)
        axes[1].set_title("After")

        # Show the plot 
        plt.show() 

        return modified_image, fig

    return modified_image

"""Generate an RGB mask of the bayer pattern 
   for an arbitrary frame size
"""
def generate_RGB_mask(original_frame: np.ndarray, visualize_results: bool=False) -> tuple[np.ndarray, object] | np.ndarray:
    # Initialize a variable for the figure handle that 
    # will be used to visualize (if desired)
    fig: object | None = None

    # Initialize an array of characters. RGB will represent 
    # the indices where there are the respective colors in 
    # bayer pattern 
    mask: np.ndarray = np.full(original_frame.shape[:2], 'x', dtype='<U1')
    
    # Extract the dimensions of the frame 
    height, width = original_frame.shape[:2]

    # Create a list for only the red pixels in a frame
    world_r_pixels: np.ndarray = np.array( [ (r, c) 
                                            for r in range(height) 
                                            for c in range(width)
                                            if(r % 2 != 0 and c % 2 != 0)
                                        ], 
                                        dtype=np.uint64
                                        )

    # Create a list for only the green pixels in a frame 
    world_g_pixels: np.ndarray = np.array( [ (r, c) 
                                            for r in range(height)
                                            for c in range(width)
                                            if(r % 2 == 0 and c % 2 != 0)
                                            or(r % 2 != 0 and c % 2 == 0) 
                                        ],  
                                        dtype=np.uint64
                                        )

    # Create a list for only the blue pixels in a frame
    world_b_pixels: np.ndarray = np.array([ (r, c)
                                            for r in range(height)
                                            for c in range(width)
                                            if(r % 2 == 0 and c % 2 == 0)
                                          ], 
                                          dtype=np.uint64
                                        )
    
    # Set the values in the mask 
    for color, pixel_coords in zip("RGB", (world_r_pixels, world_g_pixels, world_b_pixels)):
        rows: np.ndarray = pixel_coords[:, 0]
        cols: np.ndarray = pixel_coords[:, 1]

        mask[rows, cols] = color

    # Visualize the results if desired
    if(visualize_results is True):
        # Convert the mask to have color 
        colored_image: np.ndarray = np.empty(original_frame.shape, dtype=np.uint8)
        for idx, (color, pixel_coords) in enumerate(zip("RGB", (world_r_pixels, world_g_pixels, world_b_pixels))):
            rows: np.ndarray = pixel_coords[:, 0]
            cols: np.ndarray = pixel_coords[:, 1]
            color_vector: np.ndarray = np.zeros((3,), dtype=np.uint8)
            color_vector[idx] = 255
            colored_image[rows, cols] = color_vector

        # Initialize a figure with two axes 
        fig, axes = plt.subplots(1, 1)

        # Middle axis will be the after but without 
        # color correction 
        axes[0].imshow(colored_image)
        axes[0].set_title(f"{original_frame.shape} Bayer RGB Pattern")

        # Show the plot 
        plt.show() 

        return mask, fig
    
    return mask

"""Apply the per color weights to the color 
   pixels of a frame 
"""
def apply_color_correction(original_frame: np.ndarray, visualize_results: bool=False) -> tuple[np.ndarray, object] | np.ndarray:
    # Initialize a variable for the figure handle that 
    # will be used to visualize (if desired)
    fig: object | None = None
    
    # First, we must cast the original frame to a float array 
    # to apply float scalars 
    frame_as_float: np.ndarray = original_frame.astype(np.float64)

    # Next, we need to generate a bayer pattern for this size of image 
    mask: np.ndarray = generate_RGB_mask(original_frame)
    assert mask.shape[:2] == original_frame.shape[:2], f"Mask: {mask.shape[:2]} and original frame shape {original_frame.shape[:2]} are unequal"

    # Next, we will apply the weights 
    for color, weight in zip("RGB", WORLD_RGB_SCALARS):
        # Find the pixels that match this color 
        pixels: np.ndarray = np.argwhere(mask == color)
        rows: np.ndarray = pixels[:, 0]
        cols: np.ndarray = pixels[:, 1]

        # Apply the weight to the specified pixels 
        frame_as_float[rows, cols] *= weight
    
    # Round and clip values in the 255 range and cast back to uint8
    modified_frame: np.ndarray = np.clip(np.round(frame_as_float), 0, 255).astype(np.uint8)

    # Visualize the results if desired 
    if(visualize_results is True):
        # Initialize a figure with two axes 
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Color correction (Before / After)", fontweight='bold', fontsize=18)
        
        # Left axis will be the "before"
        # image 
        axes[0].imshow(original_frame, cmap="gray")
        axes[0].set_title('Before')

        # Middle axis will be the after but without 
        # color correction 
        axes[1].imshow(modified_frame)
        axes[1].set_title("After")

        # Show the plot 
        plt.show() 

        return modified_frame, fig

    # Return the modified frame
    return modified_frame

"""Embed a world frame's timestamp into the 8 bit image itself"""
def embed_timestamp(original_frame: np.ndarray, timestamp: np.float64, visualize_results: bool=False) -> tuple[np.ndarray, object] | np.ndarray:
    # Initialize a variable for the figure handle that 
    # will be used to visualize (if desired)
    fig: object | None = None
    assert original_frame.dtype == np.uint8, f"To do proper embedding, the frame must be uint8"
    assert type(timestamp) == np.float64, f"To do proper embedding, the timestamp must be float64"

    # Allocate a copy of the image we will edit and flatten it so we can directly 
    # edit the pixels 
    embedded_frame: np.ndarray = original_frame.copy().flatten() 

    # Let's convert the timestamp to its bytes representation 
    timestamp_as_u8: np.ndarray = np.frombuffer(np.array(timestamp, dtype='<f8').tobytes(), dtype=np.uint8) # little endian bytes
    embedded_frame[:len(timestamp_as_u8)] = timestamp_as_u8
    embedded_frame = embedded_frame.reshape(original_frame.shape)

    # Visualize the results if desired 
    if(visualize_results is True):
        # Initialize a figure with two axes 
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Timestamp Embedding (Before / After)", fontweight='bold', fontsize=18)
        
        # Left axis will be the "before"
        # image 
        axes[0].imshow(original_frame, cmap="gray")
        axes[0].set_title('Before')

        # Middle axis will be the after but without 
        # color correction 
        axes[1].imshow(embedded_frame.reshape(original_frame.shape))
        axes[1].set_title("After")

        # Show the plot 
        plt.show() 

        return embedded_frame, fig

    return embedded_frame 


"""Extract the world frame timestamp's from an embedded 8 bit image"""
def extract_timestamp(embedded_frame: np.ndarray) -> np.float64:
    assert embedded_frame.dtype == np.uint8, f"To extract a timestamp, the frame must be uint8"
    
    # First, let's flatten the image to get the bytes as a sequence 
    flattened_image: np.ndarray = embedded_frame.flatten() 

    # Next, we will take the first 8 bytes and view them as np.float64, 
    # interpreting them as little endian 
    timestamp_bytes: np.ndarray = flattened_image[:8]

    # Convert the timestamp back to np.float64, knowing they were assigned as little endian 
    timestamp: np.float64 = timestamp_bytes.view('<f8')[0]

    return timestamp

"""Convert a bayer image to LMS color space"""
def bayer_to_lms(original_image: np.ndarray, 
                 camera: Literal["IMX219", "standard"]="IMX219",
                 path_to_spectral_sensitivities: str="", 
                 visualize_results: bool=False,
                 matlab_engine: object | None=None
                ) -> tuple[object, np.ndarray] | np.ndarray:  
    # We need a Psychtoobox function to complete this (sadge)
    # namely  WlsToS 
    if(matlab_engine is None):
        import matlab.engine
        matlab_engine = matlab.engine.start_matlab()
        matlab_engine.tbUseProject('lightLoggerAnalysis', nargout=0)

    # First, let's make a copy of the original image
    modified_image: np.ndarray = original_image.copy() 

    # Load in the dataframe and then convert to numpy and not evil pandas (because i am less fluent in it ;-;)
    table: np.ndarray | None = None
    if(path_to_spectral_sensitivities.endswith(".mat")):
        table = pd.DataFrame(mat73.loadmat(path_to_spectral_sensitivities)["T"]).to_numpy()
    else:
        table = pd.read_excel(path_to_spectral_sensitivities).to_numpy() 

    # Extract the wavelengths from the table 
    wavelengths: np.ndarray = table[:, 0]
    rgb_values: np.ndarray = table[:, (-3, -2, -1)]

    # Convert wavelengths to sampling format
    samples: np.ndarray = np.array(matlab_engine.WlsToS(matlab.double(wavelengths), nargout=1), dtype=np.float64)
    
    # Get the sensitivities for the foveal cone classes
    field_size_degrees: int = 30
    observer_age_in_years: int = 30
    pupil_diameter_mm: int = 2
    t_receptors: np.ndarray = np.array(matlab_engine.GetHumanPhotoreceptorSS(matlab.double(samples),
                                                                    {'LConeTabulatedAbsorbance2Deg', 'MConeTabulatedAbsorbance2Deg','SConeTabulatedAbsorbance2Deg'},
                                                                    *[matlab.double(item) for item in (field_size_degrees, observer_age_in_years, pupil_diameter_mm)], [], [], [], [], 
                                                                    nargout=1
                                                                   ), dtype=np.float64)

    #Create the spectrum implied by the rgb camera weights, and then project
    # % that on the receptors

    # Splice the RGB pixels from the image into another matrix  [nPixels, [R, G, B] ]

    # Multiply the result below 
   #  cone_vec: np.ndarray = np.transpose(( t_receptors * np.transpose((rgbVec @ rgb_values)) ))

    # Then put them back in their place

    """
    % Convert RGB --> LMS contrast relative to background
    background_RGB = mean(rgbSignal, 1);
    modulation_RGB = rgbSignal - background_RGB;
    modulation_LMS = cameraToCones(modulation_RGB, options.camera);
    background_LMS = cameraToCones(background_RGB, options.camera);
    
    lmsSignal = modulation_LMS ./ background_LMS;
    
    % Select a post-receptoral channel
    switch options.postreceptoralChannel
        case {'LM'}
            signal = (lmsSignal(:,1)+lmsSignal(:,2))/2;
        case {'L-M'}
            signal = (lmsSignal(:,1)-lmsSignal(:,2));
        case {'S'}
            signal = ((lmsSignal(:,3)-lmsSignal(:,1))+lmsSignal(:,2))/2;
    end
    """


    return 

"""Convert an RGB image to LMS color space"""
def rgb_to_lms(original_image: np.ndarray, visualize_results: bool=False) -> tuple[object, np.ndarray] | np.ndarray:


    return 


"""Calculate the per color weights to apply to each color 
   in order to equalize them to the R channel 
"""
def calculate_color_weights(chunks: list[dict] | list[np.ndarray], guard_saturation: bool=False) -> np.ndarray:
    # Initialize a dict to keep track of mean pixel values by color per frame 
    spacial_averages: dict[str, float] = {let: [] 
                                          for let in 'RGB'
                                         }

    # Iterate over each chunk
    for chunk_idx, chunk in enumerate(chunks):    
        # Extract the world frames for this chunk 
        world_frames: np.ndarray = chunk['W']['v'] if isinstance(chunk, dict) else chunk

        # Find the bayer pattern matching this size 
        # of the frames in this chunk 
        RGB_mask: np.ndarray = generate_RGB_mask(world_frames[0])
        r_pixels, g_pixels, b_pixels = [ np.argwhere(RGB_mask == color) for color in "RGB" ]

        # Iterate over the frames in this chunk 
        for frame in world_frames:
            # Iterate over the colors and their pixels in each frame
            for color, pixels in zip(spacial_averages.keys(), (r_pixels, g_pixels, b_pixels)):
                # Find the pixel colors 
                color_pixels: np.ndarray = frame[pixels[:, 0], pixels[:, 1]]
                
                if(guard_saturation is True and np.any(color_pixels >= 255)):
                    raise Exception("Found saturated pixels")
                
                # Calculate the mean of the pixels of this color in this frame 
                color_frame_mean: float = np.mean(color_pixels)

                # Save this value 
                spacial_averages[color].append(color_frame_mean)


    # Construct the average pixel value over time in addition to over space (per frame)
    temporal_averages: dict[str, float] = {let: np.mean(means)
                                           for let, means in spacial_averages.items()
                                          }

    # Construct the scalars in order to equalize all colors to the R pixels 
    scalars: np.ndarray = np.array([1, temporal_averages['R'] / temporal_averages['G'],  temporal_averages['B'] /  temporal_averages['R'] ], dtype=np.float64)

    return scalars

