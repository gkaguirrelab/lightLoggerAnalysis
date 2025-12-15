import numpy as np
import matplotlib.pyplot as plt
import cv2
import os 
import sys

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
        axes[0].imshow(modified_image, cmap="gray")
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
    mask: np.ndarray = np.full(original_frame.size, 'x', dtype='<U1')
    
    # Extract the dimensions of the frame 
    height, width = original_frame.shape[:2]

    # Create a list for only the red pixels in a frame
    world_r_pixels: np.ndarray = np.array( [ (r, c) 
                                            for r in range(height) 
                                            for c in range(width[1])
                                            if(r % 2 != 0 and c % 2 != 0)
                                        ], 
                                        dtype=np.uint64
                                        )

    # Create a list for only the green pixels in a frame 
    world_g_pixels: np.ndarray = np.array( [ (r, c) 
                                            for r in range(height[0])
                                            for c in range(width[1])
                                            if(r % 2 == 0 and c % 2 != 0)
                                            or(r % 2 != 0 and c % 2 == 0) 
                                        ],  
                                        dtype=np.uint64
                                        )

    # Create a list for only the blue pixels in a frame
    world_b_pixels: np.ndarray = np.array([ (r, c)
                                            for r in range(height[0])
                                            for c in range(width[1])
                                            if(r % 2 == 0 and c % 2 == 0)
                                          ], 
                                          dtype=np.uint64
                                        )
    
    # Set the values in the mask 
    for color, pixel_coords in zip("RGB", (world_r_pixels, world_g_pixels, world_b_pixels)):
        rows: np.ndarray = pixel_coords[: 0]
        cols: np.ndarray = pixel_coords[:, 1]

        mask[rows, cols] = color

    # Visualize the results if desired
    if(visualize_results is True):
        # Convert the mask to have color 
        colored_image: np.ndarray = np.empty(original_frame.shape, dtype=np.uint8)
        for idx, (color, pixel_coords) in enumerate(zip("RGB", (world_r_pixels, world_g_pixels, world_b_pixels))):
            rows: np.ndarray = pixel_coords[: 0]
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
    # First, we must cast the original frame to a float array 
    # to apply float scalars 
    frame_as_float: np.ndarray = original_frame.astype(np.float64)

    # Next, we need to generate a bayer pattern for this size of image 
    mask: np.ndarray = generate_RGB_mask(original_frame)

    # Next, we will apply the weights 
    for color, weight in zip("RGB", WORLD_RGB_SCALARS):
        # Find the pixels that match this color 
        pixels: np.ndarray = np.argwhere(mask == color)

        # Apply the weight to the specified pixels 
        frame_as_float[pixels] *= weight
    
    # Round and clip values in the 255 range and cast back to uint8
    modified_frame: np.ndarray = np.clip(np.round(frame_as_float), 0, 255).astype(np.uint8)

    # Return the modified frame
    return modified_frame

"""Embed a world frame's timestamp into the 8 bit image itself"""
def embed_timestamp(original_frame: np.ndarray, timestamp: np.float64, visualize_results: bool) -> tuple[np.ndarray, object] | np.ndarray:
    assert original_frame.dtype == np.uint8, f"To do proper embedding, the frame must be uint8"
    assert type(timestamp) == np.float64, f"To do proper embedding, the timestamp must be float64"

    # Allocate a copy of the image we will edit and flatten it so we can directly 
    # edit the pixels 
    embedded_frame: np.ndarray = original_frame.copy().flatten() 

    # Let's convert the timestamp to its bytes representation 
    timestamp_as_u8: np.ndarray = timestamp.view('<u1')  # view as little endian uint8
    embedded_frame[:len(timestamp_as_u8)] = timestamp_as_u8

    return embedded_frame.reshape(original_frame.shape)


"""Extract the world frame timestamp's from an embedded 8 bit image"""
def extract_timestamp(embedded_frame: np.ndarray) -> np.float64:
    assert embedded_frame.dtype == np.uint8, f"To extract a timestamp, the frame must be uint8"
    
    # First, let's flatten the image to get the bytes as a sequence 
    flattened_image: np.ndarray = embedded_frame.flatten() 

    # Next, we will take the first 8 bytes and view them as np.float64, 
    # interpreting them as little endian 
    timestamp_bytes: np.ndarray = flattened_image[:8]

    # Convert the timestamp back to np.float64 
    timestamp: np.float64 = timestamp_bytes.view('<f8')[0]

    return timestamp

