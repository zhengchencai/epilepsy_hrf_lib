import os
import pandas as pd

root_path = "/home/djangoc/projects/def-jgotman/djangoc/EP_EEGfMRI_deconvolution"
scratch_path = "/scratch/djangoc/EP_EEGfMRI_deconvolution"
dir_fs = "/scratch/djangoc/EP_EEGfMRI_deconvolution/derivatives/freesurfer_AFNIUAC"

slurm_account = "def-jgotman"
time = "5:00:00"
cpus = "8"
mem_per_cpu = "4G"
email = "zhengchen.cai@gmail.com"

subj_df = pd.read_csv(os.path.join(root_path, "participants.tsv"), sep="\t")

subj = subj_df["participant_id"]

for i_sub in subj:
    # create .sh file with name i_sub_FSSUMA.sh
    file_name = i_sub + "_FSSUMA.sh"
    file_path = os.path.join(root_path, "code", "FSSUMA_script", file_name)

    dset_anat_00 = f"{root_path}/derivatives/afni/{i_sub}/anat/anatUAC.{i_sub}.nii"

    # write commands to .sh file
    # modified based on https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/tutorials/fs/fs_fsprep.html
    with open(file_path, "w") as f:
        f.write(f"#!/bin/bash\n")

        f.write(f"#SBATCH --account={slurm_account}\n")
        f.write(f"#SBATCH --job-name={i_sub}_FSSUMA.job\n")
        f.write(f"#SBATCH --output={scratch_path}/{i_sub}_FSSUMA.out\n")
        f.write(f"#SBATCH --error={scratch_path}/{i_sub}_FSSUMA.err\n")
        f.write(f"#SBATCH --time={time}\n")
        f.write(f"#SBATCH --cpus-per-task={cpus}\n")
        f.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
        f.write(f"#SBATCH --mail-user={email}\n")
        f.write(f"#SBATCH --mail-type=BEGIN\n")
        f.write(f"#SBATCH --mail-type=END\n")
        f.write(f"#SBATCH --mail-type=FAIL\n\n")

        f.write(f"export FREESURFER_HOME=$HOME/freesurfer\n\n")
        f.write(f"source $FREESURFER_HOME/SetUpFreeSurfer.sh\n\n")

        f.write(f"module load StdEnv/2020  gcc/9.3.0 afni/23.1.08\n\n")

        f.write(f"recon-all \\\n")
        f.write(f"\t-all \\\n")
        f.write(f"\t-3T \\\n")
        f.write(f"\t-sd {dir_fs} \\\n")
        f.write(f"\t-subjid {i_sub} \\\n")
        f.write(f"\t-i {dset_anat_00} \\\n")
        f.write(f"\t-parallel \n\n")

        f.write(f"@SUMA_Make_Spec_FS \\\n")
        f.write(f"\t-fs_setup \\\n")
        f.write(f"\t-NIFTI \\\n")
        f.write(f"\t-sid {i_sub} \\\n")
        f.write(f"\t-fspath {dir_fs}/{i_sub} \n")
