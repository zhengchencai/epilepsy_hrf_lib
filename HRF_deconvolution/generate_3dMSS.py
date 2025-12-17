import os
import pandas as pd

root_dir = "/home/djangoc/projects/def-jgotman/djangoc/EP_EEGfMRI_deconvolution"

scratch_path = "/scratch/djangoc/EP_EEGfMRI_deconvolution/derivatives/afni"

# set cluster HPC options
slurm_account = "def-jgotman"
time = "2:00:00"
cpus = "16"
mem_per_cpu = "4G"
email = "zhengchen.cai@gmail.com"

# set basis function, space, modulation
basis = "TENT"  # CSPLIN TENT dmUBLOCK
space = "orig"  # tlrc orig
modulation = "dm"  # amplitude modulated: am
regress_anaticor = False  # AFNI local white matter regression
if regress_anaticor:
    anaticor_radius = str(20)  # regress_anaticor FWHM defined the local region
    suffix = f"{basis}_{space}_{modulation}_anaticor{anaticor_radius}"
else:
    suffix = f"{basis}_{space}_{modulation}"

subj_df = pd.read_csv(os.path.join(root_dir, "participants.tsv"), sep="\t")

subj = subj_df["participant_id"]

for i_sub in subj:
    print(f"Making post_afni.py file for {i_sub}")
    # create .sh file with name i_sub_SSwarper.sh
    file_name = i_sub + f"_postafni_{suffix}.sh"
    if not os.path.isdir(os.path.join(root_dir, "code", "postafni_script")):
        os.mkdir(os.path.join(root_dir, "code", "postafni_script"))
    file_path = os.path.join(root_dir, "code", "postafni_script", file_name)

    sub_dir = os.path.join(scratch_path, i_sub)

    script_dir = os.path.join(sub_dir, "script")
    result_dir = os.path.join(sub_dir, suffix)

    # write afni_proc.py code

    with open(file_path, "w") as f:
        f.write("#!/bin/bash\n")
        # HPC set up
        f.write(f"#SBATCH --account={slurm_account}\n")
        f.write(f"#SBATCH --job-name={i_sub}_postafni_{suffix}.job\n")
        f.write(f"#SBATCH --output={scratch_path}/{i_sub}_postafni_{suffix}.out\n")
        f.write(f"#SBATCH --error={scratch_path}/{i_sub}_postafni_{suffix}.err\n")
        f.write(f"#SBATCH --time={time}\n")
        f.write(f"#SBATCH --cpus-per-task={cpus}\n")
        f.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
        f.write(f"#SBATCH --mail-user={email}\n")
        f.write(f"#SBATCH --mail-type=BEGIN\n")
        f.write(f"#SBATCH --mail-type=END\n")
        f.write(f"#SBATCH --mail-type=FAIL\n\n")
        # load modules required by afni and python
        f.write("module load StdEnv/2020  gcc/9.3.0 afni/23.1.08 scipy-stack\n\n")

        f.write("source /home/djangoc/afni/bin/activate \n\n")

        # write post_afni.py line by line

        if regress_anaticor:
            f.write(
                f"python {root_dir}/code/postAFNI.py -sub {i_sub} -root_dir {root_dir} -afni_dir {scratch_path}  -regress_anaticor \n\n"
            )
        else:
            f.write(
                f"python {root_dir}/code/postAFNI.py -sub {i_sub} -root_dir {root_dir} -afni_dir {scratch_path} \n\n"
            )

        f.write(f"3dMSS -prefix {result_dir}/smoothHRF.{i_sub} -jobs {cpus} \\\n")
        f.write(f"  -lme 's(TR)' \\\n")
        f.write(f"  -qVars 'TR' \\\n")
        f.write(
            f"  -prediction @{script_dir}/{basis}_{space}_{modulation}.HRF.table \\\n"
        )
        f.write(
            f"  -dataTable  @{script_dir}/{basis}_{space}_{modulation}.smooth-HRF.table \n\n"
        )

        f.write(
            f"find {result_dir} -name 'smoothHRF.{i_sub}*.gz' -exec gzip -d {{}} \\;"
        )
