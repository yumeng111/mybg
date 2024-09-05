import ROOT as rt
import sys, os, string, re
#from ROOT import TCanvas, TLegend, TChain, THStack, TH1D, TH1, kGray, kGreen, kRed, kBlue, kYellow, kMagenta, kSpring, kCyan, kTRUE, TLatex, gPad, TPad
from array import array
import math
from math import *
import pickle

from tdrStyle import *
setTDRStyle()

# Load the categorized samples and chains from the pickle file
with open('file_chains.pkl', 'rb') as f:
    categorized_samples, Event_Tree, entry_counts = pickle.load(f)

#sample_path = f"GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-250/nano_0.root"
def print_selected_events(samples, cut,output_file):
    
    sample_path = "GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-250/nano_0.root"
    if sample_path == "GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-250/nano_0.root":
    #for sample_path in samples:
        sample_name = sample_path.split('/')[0]
        if sample_name not in Event_Tree:
            print(f"Warning: {sample_name} not found in Event_Tree.")
            #continue
        
        tree = Event_Tree[sample_name]
        entries = tree.GetEntries()
        print(f"Processing {sample_name} with {entries} entries...")

        with open(output_file, 'w') as f:
            f.write(f"Processing {sample_name} with {entries} entries...\n")
        # Define the cut as a TFormula
        
            formula = ROOT.TTreeFormula("cut_formula", cut, tree)

            # Loop over the events
            for i in range(entries):
                tree.GetEntry(i)

                # Evaluate the cut formula
                if formula.EvalInstance():
                    # Print all branch values
                    branches = tree.GetListOfBranches()
                    for branch in branches:
                        branch_name = branch.GetName()
                        value = getattr(tree, branch_name)
                        #print(f"{branch_name}: {value}")
                        f.write(f"{branch_name}: {value}\n")

cut = "(lep1_type==1 && lep2_type==-1 && lep2_genLep_kind!=-1)"
dl_signals = categorized_samples["dl_signal"]
dl_backgrounds = categorized_samples["dl_bg"]
output_file = "output.txt"

print_selected_events(dl_signals, cut,output_file)
print_selected_events(dl_backgrounds, cut,output_file)



