import os
import cv2
import numpy as np
from scipy.io import savemat

video_path = "/Users/sophiamirabal/Library/CloudStorage/Dropbox-Aguirre-BrainardLab/Sophia Mirabal/FLIC_admin/Equipment/ArduCam B0392 IMX219 Wide Angle M12/fielding_function/planetarium_fielding_function.avi"

assert os.path.exists(video_path), f"Video not found: {video_path}"

cap = cv2.VideoCapture(video_path)
assert cap.isOpened(), "Could not open video"

fps = cap.get(cv2.CAP_PROP_FPS)
print("FPS:", fps)

# THROUGHOUT VIDEO
# times_sec = [
#     366, 372,   # 6:06, 6:12
#     393, 401,   # 6:33, 6:41
#     427, 435,   # 7:07, 7:15
#     456, 466,   # 7:36, 7:46
#     488, 498,   # 8:08, 8:18
#     518, 528,   # 8:38, 8:48
#     550, 560,   # 9:10, 9:20
#     582, 592,   # 9:42, 9:52
#     614, 624,   # 10:14, 10:24
#     646, 656,   # 10:46, 10:56
#     678, 686,   # 11:18, 11:26
#     708, 718,   # 11:48, 11:58
#     742, 752,   # 12:22, 12:32
#     774, 784,   # 12:54, 13:04
#     806, 814,   # 13:26, 13:34
#     836, 846,   # 13:56, 14:06
#     868, 878,   # 14:28, 14:38
#     912, 920    # 15:12, 15:20
# ]

# SUSPECTED ROTATIONAL PERIOD NO. 1 (6-7 MIN)
# times_sec = [
#     360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370,
#     371, 372, 373, 374, 375, 376, 378, 384, 385, 386, 387, 
#     388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 
#     399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409
# ]

# SUSPECTED ROTATIONAL PERIOD NO. 2 (8-9 MIN)
# times_sec = [
#     480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490,
#     491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 
#     502, 503, 504, 505, 506, 510, 511, 512, 513, 514, 515, 
#     516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526,
#     527, 528, 529, 530, 531, 532, 533, 534, 535, 536
# ]

# SUSPECTED ROTATIONAL PERIOD NO. 3 (10-11 MIN)
# times_sec = [
#     606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616,
#     617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627,
#     628, 629, 630, 631, 632, 638, 639, 640, 641, 642, 643,
#     644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 
#     655, 656, 657, 658, 659, 660, 661, 662, 663, 664
# ]

# SUSPECTED ROTATIONAL PERIOD NO. 4 (12-13 MIN)
times_sec = [
    734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744,
    745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755,
    756, 757, 758, 759, 760, 766, 767, 768, 769, 770, 771, 
    772, 773, 774, 775, 776, 777, 778, 779, 780, 781, 782,
    783, 784, 785, 786, 787, 788, 789, 790, 791, 792
]

frames_idx = [round(t * fps) for t in times_sec]
print("Frame indices:", frames_idx)

frames = []

for idx in frames_idx:
    cap.set(cv2.CAP_PROP_POS_FRAMES, idx)
    ret, frame = cap.read()

    if not ret:
        raise RuntimeError(f"Could not read frame {idx}")

    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    frames.append(gray)

cap.release()

frames = np.array(frames, dtype=np.uint8)

savemat("selected_twilight_frames.mat", {
    "frames": frames,
    "times_sec": np.array(times_sec),
    "frames_idx": np.array(frames_idx),
    "fps": fps
})

print("Saved selected_twilight_frames.mat")
print("frames shape:", frames.shape)