import numpy as np 

# Store the scalar multipliers for all of the different colors of pixel's in an image 
# required to equalize the color to the R pixel values in that frame.
# We are currently not using this, but calculated it off of a 10 minute 
# video in a variety of settings (inside (office) outside (street, park))
# Just something to note that these are inherently unequal
WORLD_RGB_SCALARS: np.ndarray = np.array([1, 0.8223775, 0.95367937], dtype=np.float64)

def RGB_mask(world_frame_shape: tuple[int]) -> np.ndarray:
    # Create an empty mask for the desired output frame 
    # that denotes which pixels are which color in the bayer 
    # image (R=0, G=1, B=2)
    world_RGB_mask: np.ndarray = np.zeros(world_frame_shape, dtype="U1")

    # Create a list for only the red pixels in a frame
    world_R_pixels: np.ndarray = np.array( [ (r, c) 
                                            for r in range(world_frame_shape[0]) 
                                            for c in range(world_frame_shape[1])
                                            if(r % 2 != 0 and c % 2 != 0)
                                        ], 
                                        dtype=np.uint64
                                        )

    # Create a list for only the green pixels in a frame 
    world_G_pixels: np.ndarray = np.array( [ (r, c) 
                                            for r in range(world_frame_shape[0])
                                            for c in range(world_frame_shape[1])
                                            if(r % 2 == 0 and c % 2 != 0)
                                            or(r % 2 != 0 and c % 2 == 0) 
                                        ],  
                                        dtype=np.uint64
                                        )

    # Create a list for only the blue pixels in a frame
    world_B_pixels: np.ndarray = np.array([ (r, c)
                                            for r in range(world_frame_shape[0])
                                            for c in range(world_frame_shape[1])
                                            if(r % 2 == 0 and c % 2 == 0)
                                        ], 
                                        dtype=np.uint64
                                        )

    # Assemble the RGB mask 
    for pixel_color, pixel_indices in zip("RGB", (world_R_pixels, world_G_pixels, world_B_pixels)):
        # Extract the rows and cols of these pixels 
        rows: np.ndarray = pixel_indices[:, 0]
        cols: np.ndarray = pixel_indices[:, 1]

        # Set the rows and cols to be this pixel color's number
        world_RGB_mask[rows, cols] = pixel_color

    return world_RGB_mask

def main():
    pass 

if(__name__ == "__main__"):
    main() 