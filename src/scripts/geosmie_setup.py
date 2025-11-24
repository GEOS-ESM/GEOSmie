#!/usr/bin/env python 

"""
  This script is expected to be run from the GEOSAQcGAN/install/bin folder.
  It assumes that the template file geosmie_run.j is available in the folder.

  The script asks the user to provide:
    - an experiment name (exp_name)
    - the group id (group_id), i.e., the NCCS sponsor code to be used in SLURM.
 
  It will then create an experiment directory that has a self-contained and ready
  to use SLURM script geosmie_run.j.
"""

from pathlib import Path
import sys
import os
import subprocess
import shutil
import glob

def print_message():
    mssg = """
    ---------------------------------------------------------------------------------
    This setup script creates a self-contained experiment directory to run
    geosmie.

    The script is interactive and asks the user to provide:

        - an experiment name (exp_name)
        - the run script you want to start from
        - the location of GRASP dust kernels
 
    It will then create an experiment directory that has a self-contained and ready
    to use run script
    ---------------------------------------------------------------------------------
    """
    print(mssg)

def search_reaplace_in_file(loc_filename: str, 
                            target_dir: Path, 
                            dict_words: dict) -> None:
    """
    Take a file template to search and replace collection of words.
    The new file (with the same name) will be created in the target directory.

    Parameters
    ----------
    loc_filename : str
       Local template file name.
    target_dir : Path
       Target directory where the new file will be created.
    dict_words : dict
       Dictionary where the keys are old words and the corresponding values
       are the new words.
    """
    new_filename = target_dir / loc_filename
    shutil.copy(loc_filename, new_filename)

    try:
        with open(new_filename, 'r') as fid:
            file_content = fid.read()

        for key in dict_words:
            file_content = file_content.replace(key, dict_words[key])
            print(f"Successfully replaced '{key}' with '{dict_words[key]}' in '{new_filename}'.")
            print("")

        with open(new_filename, 'w') as file:
            file.write(file_content)

    except FileNotFoundError:
        print(f"Error: File '{new_filename}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


def create_experiment_directory():

    # Get the current directory
    # Will be in the form FULL_PATH/GEOSmie/install/bin
    current_directory = Path.cwd()

    # Determine the source code main directory
    # Will be FULL_PATH/GEOSmie
    source_directory = current_directory.parent.parent

    reference_directory = source_directory.parent

    # Obtain the experiment name
    experiment_name = input("Provide the experiment name (one word):  ")
    experiment_name = experiment_name.strip()

    if not experiment_name:
        print("You need to provide and experiment name")
        sys.exit()

    if len(experiment_name.split()) > 1:
        print(f"The experiment name ({experiment_name}) should be in one word.")
        sys.exit()

    # Create the experiment directory

    experiment_directory = reference_directory / experiment_name
    print(f"The following experiment directory will be created: \n\n\t {experiment_directory}")
    print()

    experiment_directory.mkdir(parents=True, exist_ok=True)

    # Copy the geosparticles configuration files to the experiment directory.
    config_filepath = source_directory / "install/etc/geosparticles"
    shutil.copytree(config_filepath, experiment_directory / config_filepath.name,dirs_exist_ok=True)

    # Copy geosmie scripts to the experiment directory
    config_filepath = current_directory / "geosmie/*"
    for p in glob.glob(str(config_filepath)):
        if os.path.isfile(p):
            shutil.copy(p, experiment_directory)
        elif os.path.isdir(p):
            shutil.copytree(p, experiment_directory / os.path.basename(p),dirs_exist_ok=True)

    # Copy gsf scripts to the experiment directory
    config_filepath = current_directory / "gsf/*"
    for p in glob.glob(str(config_filepath)):
        if os.path.isfile(p):
            shutil.copy(p, experiment_directory)
        elif os.path.isdir(p):
            shutil.copytree(p, experiment_directory / os.path.basename(p),dirs_exist_ok=True)

    # Copy utils scripts to the experiment directory
    config_filepath = current_directory / "utils/*"
    for p in glob.glob(str(config_filepath)):
        if os.path.isfile(p):
            shutil.copy(p, experiment_directory)
        elif os.path.isdir(p):
            shutil.copytree(p, experiment_directory / os.path.basename(p),dirs_exist_ok=True)

    # Get the template script
    nscript = "proc.v2.0.0.csh"
    script_name = input(f"Provide the script name [default: {nscript}]:  ")
    script_name = script_name.strip()
    if not script_name:
        script_name = nscript

    target_dir = experiment_directory
    dict_words = {"@SRCDIR": str(source_directory)}
    search_reaplace_in_file(script_name, target_dir, dict_words)

    # Dust kernels
    dkernel = "/home/pcolarco/geos_aerosols/pcolarco/GEOSmie/kernels"
    if(os.uname().nodename[0:6] == 'bender'):
        dkernel = "/home/colarco/ExtData/chemistry/kernels"
    kernel_dir = input(f"Provide the location of the GRASP dust kernels [default: {dkernel}]:  ")
    kernel_dir = kernel_dir.strip()
    if not kernel_dir:
        kernel_dir = dkernel

    # copy kernels over
    os.symlink(kernel_dir, experiment_directory / os.path.basename(kernel_dir), target_is_directory=False)


    print()
    print("-"*70)
    print(f"The experiment directory was created: \n\n\t {experiment_directory}")
    print()
    print(f"Go to the folder and if necessary edit the file {script_name}.")
    print("-"*70)
    print()

if __name__ == "__main__":
    print_message()
    create_experiment_directory()
