import ROOT as rt
import sys, os, string, re
#from ROOT import TCanvas, TLegend, TChain, THStack, TH1D, TH1, kGray, kGreen, kRed, kBlue, kYellow, kMagenta, kSpring, kCyan, kTRUE, TLatex, gPad, TPad
from array import array
import math
from math import *
import pickle


with open('file_chains.pkl', 'rb') as f:
    categorized_samples, Event_Tree, entry_counts = pickle.load(f)

def check_samples(samples, sample_type):
    for sample_path in samples:
        sample_name = sample_path.split('/')[0]
        print(f"Processing {sample_type} sample: {sample_name} from path: {sample_path}")
        if sample_name not in Event_Tree:
            print(f"Warning: {sample_name} not found in Event_Tree.")
        else:
            print(f"{sample_name} found in Event_Tree with {entry_counts[sample_name]} entries")

# Check DL signals and backgrounds
dl_signals = categorized_samples["dl_signal"]
dl_backgrounds = categorized_samples["dl_bg"]

print("Checking DL Signals")
check_samples(dl_signals, "DL Signal")

print("\nChecking DL Backgrounds")
check_samples(dl_backgrounds, "DL Background")

# Check SL signals and backgrounds
sl_signals = categorized_samples["sl_signal"]
sl_backgrounds = categorized_samples["sl_bg"]

print("\nChecking SL Signals")
check_samples(sl_signals, "SL Signal")

print("\nChecking SL Backgrounds")
check_samples(sl_backgrounds, "SL Background")
