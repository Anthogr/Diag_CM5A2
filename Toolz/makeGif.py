import os
import glob
from PIL import Image

def makeGif(folder, file_name,timestep_frame):
    listFiles = glob.glob(f"{folder}/*.png")
    listFiles.sort()
    frames = [Image.open(image) for image in listFiles]
    frame_one = frames[0]
    frame_one.save(folder + file_name, format="GIF", append_images=frames,
               save_all=True, duration=timestep_frame, loop=0)
    os.system(f"rm -rf {folder}*.png")