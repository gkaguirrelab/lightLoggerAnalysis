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
def calculate_color_weights(sorted_calibration_measurements: dict, visualize_results: bool=False) -> np.ndarray:
    # Let's generate the bayer pattern for a 480, 640 frame 
    bayer_pattern: np.ndarray = generate_RGB_mask(np.zeros((480, 640), dtype=np.uint8))
    R_pixel_locations: set[tuple] = set([(y, x) for (y, x) in zip(*np.where(bayer_pattern == 'R')) ])
    G_pixel_locations: set[tuple] = set([(y, x) for (y, x) in zip(*np.where(bayer_pattern == 'G')) ])
    B_pixel_locations: set[tuple] = set([(y, x) for (y, x) in zip(*np.where(bayer_pattern == 'B')) ])

    # Let's splice out only contrast gamma NDF 0 
    contrast_gamma_NDF0: list = sorted_calibration_measurements["contrast_gamma"][0]
    
    # Then, we will iterate over the contrast levels
    RGB_weights: np.ndarray | None = None
        
    num_contrast_levels: int = len(contrast_gamma_NDF0)
    for contrast_idx in range(num_contrast_levels):
        # Then, we will extract just the first frequency (there is only 1)
        num_frequencies: int = len(contrast_gamma_NDF0[contrast_idx])
        for frequency_idx in range(num_frequencies):
            frequency_measurement = contrast_gamma_NDF0[contrast_idx][frequency_idx]

            # Then, we will go over the 3 measurements 
            avg_world_camera_v: None | np.ndarray = None
            
            # Find the min length of the measurements, because they may not all be the same length
            min_world_v_length: float = float("inf")
            num_measurements: int = len(frequency_measurement)
            for measurement_idx in range(num_measurements):
                measurement = frequency_measurement[measurement_idx][0]
                
                # Next, let's take just the world camera V values 
                world_camera_v = measurement['W']['v']
                min_world_v_length = min(min_world_v_length, len(world_camera_v))

            # TODO: In a previous measurement, a measurement was 
            #       weird, so I am just using a single out of the 3 measurements 
            #       as this weird measurement messed with the averages 
            measurements_to_consider: tuple[int] = (1, 2) # exclusive range
            start, end = measurements_to_consider
            for measurement_idx in range(start, end):
                measurement = frequency_measurement[measurement_idx][0]
                
                # Next, let's take just the world camera V values 
                world_camera_v = measurement['W']['v']

                if(avg_world_camera_v is None):
                    avg_world_camera_v = world_camera_v[:min_world_v_length, :, :].astype(np.float64)
                else:
                    avg_world_camera_v += world_camera_v[:min_world_v_length, :, :].astype(np.float64)

            avg_world_camera_v /= (end - start)

            # Next, let's splice out the target region 
            [midpt_y, midpt_x] = np.array(avg_world_camera_v.shape[1:]) // 2

            # Let's record the ROI coords wrt to entire frame 
            roi_frame_coords = [ (midpt_y + dy, midpt_x + dx)
                                for dy in range(-20, 20)
                                for dx in range(-20, 20)
                                ]
            
            roi_R_pixels = np.array([coord for coord in roi_frame_coords 
                            if coord in R_pixel_locations
                        ])


            roi_B_pixels = np.array([coord for coord in roi_frame_coords 
                            if coord in B_pixel_locations
                        ])

            roi_G_pixels = np.array([coord for coord in roi_frame_coords 
                            if coord in G_pixel_locations
                        ])

            # Now, we calculate the temporal means 
            # for each of RGB 
            temporal_means: list = []
            for means_idx, (colorname, pixels) in enumerate(zip("RGB", (roi_R_pixels, roi_G_pixels, roi_B_pixels))):
                rows: np.ndarray = pixels[:, 0]
                cols: np.ndarray = pixels[:, 1]

                roi_color_intensity = np.mean(avg_world_camera_v[:, rows, cols], axis=1) 

                roi_color_temporal_mean = float(np.mean(roi_color_intensity))
                temporal_means.append((colorname, roi_color_intensity, roi_color_temporal_mean))

            # Save the weights for this combination of contrast + frequency 
            if(RGB_weights is None):
                RGB_weights = np.zeros((num_contrast_levels, num_frequencies, 3), dtype=np.float64)
            
            RGB_weights[contrast_idx][frequency_idx][:] = [m for (c, v, m) in temporal_means]


            # If we do not want to visualize, just continue
            if(visualize_results is False):
                continue

            # Now, plot them over time for this measurement
            fig, axes = plt.subplots(1, 2, figsize=(12, 6))
            ax_ts, ax_tm = axes  # time-series, temporal-mean

            frame_numbers: np.ndarray = np.arange(min_world_v_length)
            for (colorname, spatial_mean, temporal_mean) in temporal_means:
                ax_ts.plot(
                    frame_numbers,
                    spatial_mean,
                    color=colorname.lower(),
                    linewidth=2,
                    label=f"{colorname} ROI"
                )

            # ---- Left: time-series plot ----
            ax_ts.set_title(f"ContrastIdx: {contrast_idx} | ROI Pixel Intensities by Frame Number", fontsize=14)
            ax_ts.set_xlabel("Frame Number", fontsize=12)
            ax_ts.set_ylabel("Avg Intensity", fontsize=12)
            ax_ts.set_ylim(0, 255)
            ax_ts.set_xlim(0, min(100, min_world_v_length - 1))
            ax_ts.legend()
            ax_ts.grid(alpha=0.3)

            # ---- Right: temporal mean per channel ----
            x = np.arange(len(temporal_means))

            # make bars match RGB colors
            heights = [m for (c, v, m) in temporal_means]
            labels  = [c for (c, v, m) in temporal_means]
            bar_colors = [c.lower() for c in labels]

            bars = ax_tm.bar(x,
                            heights,
                            color=bar_colors,
                            edgecolor="black",
                            linewidth=1.2
                            )

            ax_tm.set_xticks(x)
            ax_tm.set_xticklabels(labels)

            ax_tm.set_title(f"ContrastIdx: {contrast_idx} | Temporal Mean (ROI)", fontsize=14)
            ax_tm.set_xlabel("Channel", fontsize=12)
            ax_tm.set_ylabel("Temporal Mean Intensity", fontsize=12)
            ax_tm.set_xticks(x)
            ax_tm.set_xticklabels([c.upper() for c in bar_colors])
            ax_tm.set_ylim(0, 255)
            ax_tm.grid(alpha=0.3, axis="y")

            # Optional: annotate bar values
            for xi, val in zip(x, heights):
                ax_tm.text(xi, val + 3, f"{val:.1f}", ha="center", va="bottom", fontsize=10)


            fig.tight_layout()
            plt.show()
            plt.close(fig)

    return RGB_weights



