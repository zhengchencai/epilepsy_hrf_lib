import os
import pandas as pd
import subprocess
import glob


def is_oblique(img):
    """Return True if the image is oblique, otherwise False."""
    cmd = f"3dinfo -is_oblique {img}"
    result = subprocess.run(
        cmd, universal_newlines=True, check=True, shell=True, stdout=subprocess.PIPE
    )
    return result.stdout.strip() == "1"


def deoblique(img):
    """Deoblique if image is oblique."""
    cmd = f"adjunct_deob_around_origin -input {img} -prefix {img} -overwrite"
    result = subprocess.run(
        cmd, universal_newlines=True, check=True, shell=True, stdout=subprocess.PIPE
    )
    return result.returncode == 0


bids_path = "/Volumes/django/EP_EEGfMRI_deconvolution"

subj_df = pd.read_csv(os.path.join(bids_path, "participants.tsv"), sep="\t")
subj = subj_df["participant_id"]

for i_sub in subj:
    t1_images = glob.glob(f"{bids_path}/{i_sub}/anat/*.nii*")
    for img in t1_images:
        if is_oblique(img):
            print("Oblique image:", img)
            deoblique(img)
        else:
            print("Not oblique image:", img)
