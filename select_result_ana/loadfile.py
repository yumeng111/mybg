#import ROOT as rt
import ROOT
from ROOT import TFile, TChain
import sys, os, string, re

from array import array
import pwd, subprocess
import optparse, shlex, re
import time
from time import gmtime, strftime
import math
import pickle

dir= '/eos/user/y/yumeng/ana/flaf/anaTuples/localeditMuh072924/Run3_2022/'
"""
mc_samples= []

# Check if the directory exists
if not os.path.exists(dir):
    print(f"The directory {dir} does not exist.")
else:
    # Iterate over all files in the directory
    for file in os.listdir(dir):
        # Construct the full file path
        full_path = os.path.join(dir, file)
        # Check if it is a file (not a directory)
        if os.path.isfile(full_path):
            # Add the file name to mc_samples
            mc_samples.append(file)
        else:
            print(f"Skipping {full_path} as it is not a file.")

    print("Files found:", mc_samples)
"""

# Can modify: List of sample files
mc_samples= [
    'GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-250/nano_0.root',
    'GluGlutoRadiontoHHto2B2Vto2B2JLNu_M-250/nano_0.root',
    'GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-450/nano_0.root',
    'GluGlutoRadiontoHHto2B2Vto2B2JLNu_M-450/nano_0.root',
    'GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-650/nano_0.root',
    'GluGlutoRadiontoHHto2B2Vto2B2JLNu_M-650/nano_0.root',
    'GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-1000/nano_0.root',
    'GluGlutoRadiontoHHto2B2Vto2B2JLNu_M-1000/nano_0.root',
    'GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-3000/nano_0.root',
    'GluGlutoRadiontoHHto2B2Vto2B2JLNu_M-3000/nano_0.root',
    'GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-5000/nano_0.root',
    'GluGlutoRadiontoHHto2B2Vto2B2JLNu_M-5000/nano_0.root',
    'DYto2L_M-50_madgraphMLM/nano_0.root', 
]

# Dictionaries to hold TChains and entry counts
Event_Tree = {}
entry_counts = {}

def categorize_samples(samples):
    categorized_samples = {
        "dl_signal": [],
        "sl_signal": [],
        "dl_bg": [],
        "sl_bg": []
    }

    for sample in samples:
        if 'HHto2B2Vto2B2L2Nu' in sample:
            categorized_samples["dl_signal"].append(sample)
        elif 'HHto2B2Vto2B2JLNu' in sample:
            categorized_samples["sl_signal"].append(sample)
        else:
            if '2L' in sample:
                categorized_samples["dl_bg"].append(sample)
            else:
                categorized_samples["sl_bg"].append(sample)
    
    return categorized_samples

def load_files():
    categorized_samples = categorize_samples(mc_samples)

    for category in ["dl_signal", "sl_signal", "dl_bg", "sl_bg"]:
        for sample in categorized_samples[category]:
            sample_name = sample.split('/')[0]
            file_path = dir + sample

            if sample_name not in Event_Tree:
                Event_Tree[sample_name] = TChain("Events")
                entry_counts[sample_name] = 0

            Event_Tree[sample_name].Add(file_path)
            file = TFile.Open(file_path)
            tree = file.Get("Events")
            if tree:
                entry_counts[sample_name] += tree.GetEntries()
            else:
                print(f"Warning: Events tree not found in {file_path}")

    return categorized_samples

# Load the files and categorize them
categorized_samples = load_files()

# Print categorized samples
print("Categorized Samples:")
for category, samples in categorized_samples.items():
    print(f"{category}:")
    for sample in samples:
        print(f"  {sample}")

# Print Event_Tree contents
print("\nEvent_Tree Contents:")
for sample_name, chain in Event_Tree.items():
    print(f"{sample_name}: {chain.GetEntries()} entries")

# Print entry counts
print("\nEntry Counts:")
for sample_name, count in entry_counts.items():
    print(f"{sample_name}: {count} entries")

#############################################
# Save the categorized samples, chains(name: Event_Tree, which contains the tree Events from each of the ROOT files), and entry counts to a pickle file for use in plotter.py
with open('file_chains.pkl', 'wb') as f:
    pickle.dump((categorized_samples, Event_Tree, entry_counts), f)

print("Files loaded and categorized successfully.")