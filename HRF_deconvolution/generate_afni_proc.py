import os
import pandas as pd
import numpy as np
import glob
import re
import json
from decimal import Decimal
import nibabel as nib

root_dir = "/home/djangoc/projects/def-jgotman/djangoc/EP_EEGfMRI_deconvolution"

afni_dir = os.path.join(root_dir, "derivatives/afni")
scratch_path = "/scratch/djangoc/EP_EEGfMRI_deconvolution/derivatives/afni"

SUMA_path = "/scratch/djangoc/EP_EEGfMRI_deconvolution/derivatives/freesurfer_AFNIUAC"

template = f"{afni_dir}/MNI152_2009_template_SSW.nii.gz"

# set cluster HPC options
slurm_account = "def-jgotman"
time = "5:00:00"
cpus = "16"
mem_per_cpu = "4G"
email = "zhengchen.cai@gmail.com"

# set basis function, space, modulation
basis = "TENT"  # TENT
space = "orig"  # tlrc orig
modulation = "dm"
regress_anaticor = False  # AFNI local white matter regression
if regress_anaticor:
    anaticor_radius = str(20)  # regress_anaticor FWHM defined the local region
    suffix = f"{basis}_{space}_{modulation}_anaticor{anaticor_radius}"
else:
    suffix = f"{basis}_{space}_{modulation}"

subj_df = pd.read_csv(os.path.join(root_dir, "participants.tsv"), sep="\t")

subj = subj_df["participant_id"]


for i_sub in subj:
    print(f"Making afni_proc.py file for {i_sub}")
    # create .sh file with name i_sub_SSwarper.sh
    file_name = i_sub + f"_{suffix}.sh"
    if not os.path.isdir(os.path.join(root_dir, "code", "afniproc_script")):
        os.mkdir(os.path.join(root_dir, "code", "afniproc_script"))
    file_path = os.path.join(root_dir, "code", "afniproc_script", file_name)

    # set input data path
    sdir_basic = f"{root_dir}/{i_sub}"
    sdir_ssw = f"{afni_dir}/{i_sub}/anat"
    cavity_dir = f"{afni_dir}/{i_sub}/cavity"
    func_dir = f"{root_dir}/{i_sub}/func"
    sdir_timing = f"{afni_dir}/{i_sub}/events"
    out_dir = f"{scratch_path}/{i_sub}/{suffix}"
    if not os.path.isdir(f"{scratch_path}/{i_sub}/script"):
        os.makedirs(f"{scratch_path}/{i_sub}/script", exist_ok=True)
    script_file = f"{scratch_path}/{i_sub}/script/proc.{i_sub}_{suffix}"
    log_file = f"{scratch_path}/{i_sub}/script/output.{i_sub}_{suffix}"

    # data inputs
    dset_anat_00 = f"{sdir_basic}/anat/{i_sub}_T1w.nii.gz"
    anat_cp = f"{sdir_ssw}/anatSS.{i_sub}.nii"
    anat_aseg_REN = f"{SUMA_path}/{i_sub}/SUMA/aparc.a2009s+aseg_REN_all.nii.gz"
    if regress_anaticor:
        anat_aseg_wm = f"{SUMA_path}/{i_sub}/SUMA/fs_ap_wm.nii.gz"
    cavity = glob.glob(f"{cavity_dir}/*_resection_*.nii")
    dsets_epi = glob.glob(f"{func_dir}/{i_sub}_task-spike_run-*_bold.nii.gz")
    dsets_epi = sorted(dsets_epi, key=lambda x: int(re.findall(r"run-(\d+)", x)[0]))

    # if img < 50 frames, exclude
    goodbold = []
    for r, runfile in enumerate(dsets_epi):
        bold_img = nib.load(runfile)
        if bold_img.shape[3] > 50:
            goodbold.append(r)

    dsets_epi = [dsets_epi[i] for i in goodbold]

    # SSW alignment data
    dsets_NL_warp = [
        f"{sdir_ssw}/anatQQ.{i_sub}.nii",
        f"{sdir_ssw}/anatQQ.{i_sub}.aff12.1D",
        f"{sdir_ssw}/anatQQ.{i_sub}_WARP.nii",
    ]

    # timing files
    timing_files = sorted(glob.glob(os.path.join(sdir_timing, f"*_{modulation}.1D")))

    stim_classes = list(
        map(
            lambda x: "type" + re.search(rf"type(.+?)_{modulation}", x).group(1),
            timing_files,
        )
    )

    dur_all = []

    for i_evt in timing_files:
        dur_tmp = []
        with open(i_evt) as f:
            for line in f:
                parts = line.strip().split()
                for part in parts:
                    subparts = part.split(":")
                    if len(subparts) > 1:
                        dur_tmp.append(subparts[1])
        dur_all.append(round(np.mean(np.array(dur_tmp, dtype=float)), 1))

    # control variables
    nt_rm = 6
    blur_size = 6.5
    cen_motion = 1
    cen_outliers = 0.05
    run_csim = "yes"
    njobs = cpus

    # define basis function based on TR of the EPI
    run_JSON = dsets_epi[0].replace(".nii.gz", ".json")
    with open(run_JSON, "r") as f:
        metadata = json.load(f)

    TR = Decimal(str(metadata["RepetitionTime"]))

    times_offset = str(-1 * nt_rm * TR)

    # use common TENT
    hrf_start = str(0)
    hrf_end = str(17 * Decimal(1.5))  # this allow TR=1.75 and 1.9 to have 25s after 0s
    n_knot = "18"

    # set basis function, stim type
    if modulation == "dm":
        stim_types = " ".join(["AM1"] * len(timing_files))
    elif modulation == "am":
        stim_types = " ".join(["AM2"] * len(timing_files))
    else:
        raise ValueError(f"{i_sub}: event modulation type error")

    if basis in ["CSPLIN", "TENT"]:
        basis_hrf = f"'{basis}({hrf_start},{hrf_end},{n_knot})'"
    elif basis == "dmUBLOCK":
        basis_hrf = []
        for i_dur in dur_all:
            basis_hrf.append(f"'dmUBLOCK(-{i_dur})'")
        basis_hrf = " ".join(basis_hrf)
    else:
        raise ValueError(f"{i_sub}: basis function type error")

    # write afni_proc.py code

    with open(file_path, "w") as f:
        f.write("#!/bin/bash\n")
        # HPC set up
        f.write(f"#SBATCH --account={slurm_account}\n")
        f.write(f"#SBATCH --job-name={i_sub}_afniproc_{suffix}.job\n")
        f.write(f"#SBATCH --output={scratch_path}/{i_sub}_{suffix}.out\n")
        f.write(f"#SBATCH --error={scratch_path}/{i_sub}_{suffix}.err\n")
        f.write(f"#SBATCH --time={time}\n")
        f.write(f"#SBATCH --cpus-per-task={cpus}\n")
        f.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
        f.write(f"#SBATCH --mail-user={email}\n")
        f.write(f"#SBATCH --mail-type=BEGIN\n")
        f.write(f"#SBATCH --mail-type=END\n")
        f.write(f"#SBATCH --mail-type=FAIL\n\n")
        # load modules required by afni and python
        f.write("module load StdEnv/2020  gcc/9.3.0 afni/23.1.08\n")
        # write afni_proc.py line by line
        f.write("afni_proc.py\\\n")
        f.write(
            f"    -subj_id                  {i_sub}                                        \\\n"
        )
        f.write(
            f"    -out_dir                  {out_dir}                                      \\\n"
        )
        f.write(
            f"    -script                   {script_file}                                   \\\n"
        )
        # if no tlrc, remove tlrc block
        if space == "tlrc":
            f.write(
                f"    -blocks                   despike tshift align tlrc volreg mask blur scale regress     \\\n"
            )
        elif space == "orig":
            f.write(
                f"    -blocks                   despike tshift align volreg mask blur scale regress         \\\n"
            )
        else:
            raise ValueError(f"{i_sub}: space definition error")

        f.write(
            f"    -dsets                    {' '.join(dsets_epi)}                          \\\n"
        )
        f.write(
            f"    -tcat_remove_first_trs    {nt_rm}                                        \\\n"
        )
        f.write(
            f"    -radial_correlate_blocks  tcat volreg                                    \\\n"
        )
        f.write(
            f"    -copy_anat                {anat_cp}                                      \\\n"
        )
        f.write(
            f"    -anat_has_skull           no                                             \\\n"
        )
        f.write(
            f"    -anat_follower            anat_w_skull anat {dset_anat_00}               \\\n"
        )
        # two anat_follower_ROI for resampling atals and cavity on BOLD space
        f.write(
            f"    -anat_follower_ROI        FS_REN_epi epi {anat_aseg_REN}                  \\\n"
        )

        for i_cavity in cavity:
            pattern = r"resection_(.*?)_space"
            match = re.search(pattern, os.path.basename(i_cavity))
            cavity_name = "resection_" + match.group(1)
            f.write(
                f"    -anat_follower_ROI        {cavity_name}_epi epi {i_cavity}            \\\n"
            )

        if regress_anaticor:
            f.write(
                f"    -anat_follower_ROI        FS_wm_e epi {anat_aseg_wm}                  \\\n"
            )
            f.write(
                f"    -anat_follower_erode      FS_wm_e                                    \\\n"
            )

        f.write(
            f"    -volreg_align_to          MIN_OUTLIER                                    \\\n"
        )
        f.write(
            f"    -volreg_align_e2a                                                        \\\n"
        )

        # if no tlrc, remove this options
        if space == "tlrc":
            f.write(
                f"    -volreg_tlrc_warp                                                        \\\n"
            )

        f.write(
            f"    -align_opts_aea           -check_flip -cost lpc+ZZ -giant_move -AddEdge  \\\n"
        )
        f.write(
            f"    -align_unifize_epi        local                                          \\\n"
        )

        # if no tlrc, remove these template and tlrc file options
        if space == "tlrc":
            f.write(
                f"    -tlrc_base                {template}                                     \\\n"
            )
            f.write(
                f"    -tlrc_NL_warp                                                            \\\n"
            )
            f.write(
                f"    -tlrc_NL_warped_dsets     {' '.join(dsets_NL_warp)}                      \\\n"
            )

        f.write(
            f"    -mask_epi_anat            yes                                            \\\n"
        )
        # this -clfrac 0.1 controls the clipping threshold see 3dcliplevel and 3dAutomask
        f.write(
            f"    -mask_opts_automask -clfrac 0.1 -dilate 1                                \\\n"
        )
        f.write(
            f"    -blur_in_mask             yes                                            \\\n"
        )
        f.write(
            f"    -blur_size                {blur_size}                                    \\\n"
        )
        f.write(
            f"    -blur_to_fwhm                                                            \\\n"
        )
        f.write(
            f"    -regress_stim_times       {' '.join(timing_files)}                       \\\n"
        )
        f.write(
            f"    -regress_stim_labels      {' '.join(stim_classes)}                       \\\n"
        )
        f.write(
            f"    -regress_stim_times_offset  {times_offset}                               \\\n"
        )

        if basis in ["CSPLIN", "TENT"]:
            f.write(
                f"    -regress_basis            {basis_hrf}                                    \\\n"
            )
        elif basis == "dmUBLOCK":
            f.write(
                f"    -regress_basis_multi      {basis_hrf}                                    \\\n"
            )
        else:
            assert True, f"{i_sub}: basis function type error"

        f.write(
            f"    -regress_stim_types       {stim_types}                                   \\\n"
        )
        f.write(
            f"    -regress_local_times                                                     \\\n"
        )
        f.write(
            f"    -regress_motion_per_run                                                  \\\n"
        )

        # regress_anaticor option
        if regress_anaticor:
            f.write(
                f"    -regress_anaticor_fast                                               \\\n"
            )
            f.write(
                f"    -regress_anaticor_radius  {anaticor_radius}                      \\\n"
            )
            f.write(
                f"    -regress_anaticor_label   FS_wm_e                                 \\\n"
            )

        f.write(
            f"    -regress_censor_motion    {cen_motion}                                   \\\n"
        )
        f.write(
            f"    -regress_censor_outliers  {cen_outliers}                                 \\\n"
        )
        f.write(
            f"    -regress_compute_fitts                                                   \\\n"
        )
        f.write(
            f"    -regress_make_ideal_sum   sum_ideal.1D                                   \\\n"
        )

        if regress_anaticor:
            f.write(
                f"    -regress_make_corr_vols   FS_wm_e                                    \\\n"
            )

        f.write(
            f"    -regress_est_blur_epits                                                  \\\n"
        )
        f.write(
            f"    -regress_est_blur_errts                                                  \\\n"
        )
        f.write(
            f"    -regress_run_clustsim     {run_csim}                                     \\\n"
        )
        f.write(
            f"    -regress_reml_exec                                                       \\\n"
        )
        f.write(
            f"    -regress_opts_reml        -GOFORIT 99                                    \\\n"
        )
        f.write(
            f"    -regress_opts_3dD                                                        \\\n"
        )
        f.write(
            f"        -bout                                                                \\\n"
        )
        f.write(
            f"        -jobs                 {njobs}                                        \\\n"
        )
        f.write(
            f"        -allzero_OK                                                          \\\n"
        )
        f.write(
            f"        -GOFORIT              99                                             \\\n"
        )
        f.write(
            f"    -regress_3dD_stop                                                        \\\n"
        )
        f.write(
            f"    -html_review_style        pythonic                                       \n\n"
        )

        f.write(f"time tcsh -xef {script_file} 2>&1 | tee {log_file}")
