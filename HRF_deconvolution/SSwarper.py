import os
import pandas as pd

root_path = "/home/djangoc/projects/def-jgotman/djangoc/EP_EEGfMRI_deconvolution"
derivatives_dir = os.path.join(root_path, "derivatives")
scratch_path = "/scratch/djangoc/EP_EEGfMRI_deconvolution"
template = f"{derivatives_dir}/MNI152_2009_template_SSW.nii.gz"

slurm_account = "def-jgotman"
time = "3:00:00"
cpus = "8"
mem_per_cpu = "4G"
email = "zhengchen.cai@gmail.com"

subj_df = pd.read_csv(os.path.join(root_path, "participants.tsv"), sep="\t")

subj = subj_df["participant_id"]

for i_sub in subj:
    # create .sh file with name i_sub_SSwarper.sh
    file_name = i_sub + "_SSwarper.sh"
    file_path = os.path.join(root_path, "code", "SSwarper_script", file_name)

    dset_anat_00 = f"{root_path}/{i_sub}/anat/{i_sub}_T1w.nii.gz"

    out_dir = f"{scratch_path}/{i_sub}/anat"

    # write commands to .sh file
    with open(file_path, "w") as f:
        f.write(f"#!/bin/bash\n")

        f.write(f"#SBATCH --account={slurm_account}\n")
        f.write(f"#SBATCH --job-name={i_sub}_SSwarper.job\n")
        f.write(f"#SBATCH --output={scratch_path}/{i_sub}_SSwarper.out\n")
        f.write(f"#SBATCH --error={scratch_path}/{i_sub}_SSwarper.err\n")
        f.write(f"#SBATCH --time={time}\n")
        f.write(f"#SBATCH --cpus-per-task={cpus}\n")
        f.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
        f.write(f"#SBATCH --mail-user={email}\n")
        f.write(f"#SBATCH --mail-type=BEGIN\n")
        f.write(f"#SBATCH --mail-type=END\n")
        f.write(f"#SBATCH --mail-type=FAIL\n\n")

        f.write(f"module load StdEnv/2020  gcc/9.3.0 afni/23.1.08\n\n")

        f.write(f"time @SSwarper \\\n")
        f.write(f"\t-tmp_name_nice \\\n")
        f.write(f"\t-base {template} \\\n")
        f.write(f"\t-subid {i_sub} \\\n")
        f.write(f"\t-input {dset_anat_00} \\\n")
        f.write(f"\t-odir {out_dir}\n")
