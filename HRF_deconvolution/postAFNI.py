# this script generate files for results summary. Run it after afni_proc.py
# 1) save FS aseg on BOLD space brik to nii
# 2) save cavity mask on BOLD space brik to nii and save a new version fill gaps and holes -dilate_input 1
# 3) save a brik and nii of the raw HRF parameter
# 4) make data table used later for 3dMSS smoothing
import argparse
import os
import pandas as pd
import numpy as np
import nibabel as nib
import subprocess
import itertools
import glob
import re
from decimal import Decimal


def run_command(command):
    p = subprocess.run(command, universal_newlines=True, check=True, shell=True)
    print(p.stdout)


parser = argparse.ArgumentParser(
    description="Usage", formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument("-sub", help="sub ID", default="NA")
parser.add_argument(
    "-root_dir",
    help="dataset root director",
    default="/home/djangoc/projects/def-jgotman/djangoc/EP_EEGfMRI_deconvolution",
)
parser.add_argument(
    "-afni_dir",
    help="afni results director",
    default="/scratch/djangoc/EP_EEGfMRI_deconvolution/derivatives/afni",
)
parser.add_argument(
    "-regress_anaticor", help="use anaticor or not", action="store_true"
)

args = parser.parse_args()

print(args)

root_dir = args.root_dir
scratch_path = args.afni_dir
i_sub = args.sub
regress_anaticor = args.regress_anaticor
# root_dir = "/home/djangoc/projects/def-jgotman/djangoc/EP_EEGfMRI_deconvolution"

afni_dir = os.path.join(root_dir, "derivatives/afni")
# scratch_path = "/scratch/djangoc/EP_EEGfMRI_deconvolution/derivatives/afni"

# set basis function, space, modulation
basis = "TENT"  # CSPLIN TENT dmUBLOCK
space = "orig"  # tlrc orig
modulation = "dm"

if regress_anaticor:
    anaticor_radius = str(20)  # regress_anaticor FWHM defined the local region
    suffix = f"{basis}_{space}_{modulation}_anaticor{anaticor_radius}"
else:
    suffix = f"{basis}_{space}_{modulation}"

# set smoothed HRF grid
start = 0
end = 25
step = 0.25
TR_smoothed = np.arange(start, end + 0.25, step)

i_sub_dir = os.path.join(scratch_path, i_sub)
result_dir = os.path.join(i_sub_dir, suffix)

# 1) save follower ROI to nii
follower = glob.glob(f"{result_dir}/follow_ROI_*.BRIK")
mask_epi_anat = f"{result_dir}/mask_epi_anat.{i_sub}+{space}.BRIK"

for i_follower in follower:
    if "resection" in i_follower:
        if not os.path.isfile(i_follower.replace(".BRIK", ".nii")):
            cmd_3dAFNItoNIFTI = (
                "3dAFNItoNIFTI -prefix "
                + i_follower.replace(".BRIK", ".nii")
                + " "
                + i_follower
            )
            run_command(cmd_3dAFNItoNIFTI)
            # if resection mask fill gaps and holes and dilate 1 -dilate_input 1
            cmd_3dmask_tool = (
                "3dmask_tool -dilate_input 1 -fill_holes -input "
                + i_follower
                + " -prefix "
                + os.path.dirname(i_follower)
                + "/temp.nii"
            )
            run_command(cmd_3dmask_tool)
            # also interect with brain mask used in AFNI
            cmd_3dmask_tool = (
                "3dmask_tool -inter -input "
                + os.path.dirname(i_follower)
                + "/temp.nii "
                + mask_epi_anat
                + " -prefix "
                + os.path.dirname(i_follower)
                + "/temp1.nii"
            )
            run_command(cmd_3dmask_tool)

            cmd_3dmask_tool = (
                "3dmask_tool -inter -input "
                + os.path.dirname(i_follower)
                + "/temp.nii "
                + mask_epi_anat
                + " -prefix "
                + i_follower.replace("resection", "resection_dilated")
            )
            run_command(cmd_3dmask_tool)
            if os.path.isfile(
                i_follower.replace("resection", "resection_dilated") + ".gz"
            ):
                run_command(
                    f'gzip -d {i_follower.replace("resection", "resection_dilated") + ".gz"}'
                )

            os.remove(os.path.dirname(i_follower) + "/temp.nii")
            os.rename(
                os.path.dirname(i_follower) + "/temp1.nii",
                i_follower.replace("resection", "resection_dilated").replace(
                    ".BRIK", ".nii"
                ),
            )
    else:
        if not os.path.isfile(i_follower.replace(".BRIK", ".nii")):
            cmd_3dAFNItoNIFTI = (
                "3dAFNItoNIFTI -prefix "
                + i_follower.replace(".BRIK", ".nii")
                + " "
                + i_follower
            )
            run_command(cmd_3dAFNItoNIFTI)

# 2) save HRF brik and 3) make data table for 3dMSS
df_dataTable = pd.DataFrame(columns=["subject", "IED", "TR", "InputFile"])
IED_all = []

script_dir = os.path.join(i_sub_dir, "script")
filename = f"{script_dir}/proc.{i_sub}_{suffix}"
pattern = rf"{basis}\((.*?),(.*?),(.*?)\)"
match_list = []

with open(filename, "r") as file:
    contents = file.read()
    matches = re.findall(pattern, contents)
    for match in matches:
        nums = [Decimal(num) for num in match]
        match_list.append(nums)
basis_opt = match_list[0]
basis_step = (basis_opt[1] - basis_opt[0]) / (basis_opt[2] - 1)

hrf_stats = f"{result_dir}/stats.{i_sub}_REML+{space}.BRIK"
afni_img = nib.load(hrf_stats.replace(".BRIK", ".HEAD"))
brick_labels = afni_img.header.info["BRICK_LABS"].split("~")
HRF_idx = []
for index, i_brick in enumerate(brick_labels):
    if "_Coef" in i_brick:
        coef_name = os.path.basename(i_brick)
        IED = coef_name.split("#")[0]
        IED_all.append(IED)
        TR = (
            basis_opt[0]
            + Decimal(str(coef_name.split("#")[1].split("_")[0])) * basis_step
        )
        df_dataTable.loc[len(df_dataTable)] = [
            i_sub,
            IED,
            TR,
            hrf_stats + f"[{index}]",
        ]
        HRF_idx.append(index)

# save HRF sub-brik for visualization in AFNI
if os.path.isfile(hrf_stats.replace("stats", "stats.HRF")):
    os.remove(hrf_stats.replace("stats", "stats.HRF"))
    os.remove(hrf_stats.replace("stats", "stats.HRF").replace("BRIK", "HEAD"))

cmd_3dbucket = (
    "3dbucket -prefix "
    + hrf_stats.replace("stats", "stats.HRF")
    + " "
    + hrf_stats
    + f"[{','.join([str(idx) for idx in HRF_idx])}]"
)
run_command(cmd_3dbucket)

if os.path.isfile(hrf_stats.replace("stats", "stats.HRF") + ".gz"):
    run_command(f'gzip -d {hrf_stats.replace("stats", "stats.HRF") + ".gz"}')

IED_all = sorted(list(set(IED_all)))
combinations = list(itertools.product(IED_all, TR_smoothed))
combinations = [(f"{ied}-{tr}s", ied, tr) for ied, tr in combinations]
df_prediction = pd.DataFrame(combinations, columns=["label", "IED", "TR"])

# save tables required by 3dMSS
df_dataTable.to_csv(
    f"{script_dir}/{basis}_{space}_{modulation}.smooth-HRF.table", sep="\t", index=False
)
df_prediction.to_csv(
    f"{script_dir}/{basis}_{space}_{modulation}.HRF.table", sep="\t", index=False
)
