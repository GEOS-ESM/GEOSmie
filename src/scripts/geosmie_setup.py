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
        - the group id (group_id), i.e., the NCCS sponsor code to be used in SLURM.
 
    It will then create an experiment directory that has a self-contained and ready
    to use SLURM script geosmis_run.j.
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

        with open(new_filename, 'w') as file:
            file.write(file_content)

    except FileNotFoundError:
        print(f"Error: File '{new_filename}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def ignore_top_level(dir, files):
    # Ignore the source directory itself
    if dir == source_directory:
        return [source_directory]  # Return a list containing the source directory name
    return [] # Don't ignore anything else


def create_experiment_directory():

    # Get the current directory
    # Will be in the form FULL_PATH/GEOSmie/install/bin
    current_directory = Path.cwd()

    # Determine the source code main directory
    # Will be FULL_PATH/GEOSmie
    source_directory = current_directory.parent.parent.parent

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
    config_filepath = current_directory.parent.parent / "etc/geosparticles"
    shutil.copytree(config_filepath, experiment_directory / config_filepath.name,dirs_exist_ok=True)

    # Copy geosmie scripts to the experiment directory
    config_filepath = current_directory.parent / "geosmie/*"
    for p in glob.glob(str(config_filepath)):
        if os.path.isfile(p):
            shutil.copy(p, experiment_directory)
        elif os.path.isdir(p):
            shutil.copytree(p, experiment_directory / os.path.basename(p),dirs_exist_ok=True)


    # Get the sponsor code id

    result = subprocess.run(["groups"], shell=True, capture_output=True, text=True)

    groups = result.stdout.strip().split()

    print(f"The list of available group ids is: \n\n\t {result.stdout.strip()}")
    print()
    my_group = input(f"Provide the group id do you want to use (default: {groups[0]}): ")
    my_group = my_group.strip()

    if not my_group:
        my_group = groups[0]

    print()
    print(f"Your group is is: {my_group}")
    print()

    if my_group not in groups:
        print(f"You selected and invalid group.")
        print(f"The group {groups[0]} will be used.")
        print(f"You can change the group id in the SLURM script available in the experiment directory")
        print()

    loc_filename = "geosmie_runoptics.j"
    target_dir = experiment_directory
    dict_words = {"@SRCDIR": str(source_directory), "@GROUPID": my_group}
    search_reaplace_in_file(loc_filename, target_dir, dict_words)

    loc_filename2 = "geosmie_runbands.j"
    target_dir = experiment_directory
    dict_words = {"@SRCDIR": str(source_directory), "@GROUPID": my_group}
    search_reaplace_in_file(loc_filename2, target_dir, dict_words)


    print()
    print("-"*70)
    print(f"The experiment directory was created: \n\n\t {experiment_directory}")
    print()
    print(f"Go to the folder and if necessary edit the file {loc_filename} or {loc_filename2}.")
    print("-"*70)
    print()

if __name__ == "__main__":
    print_message()
    create_experiment_directory()
