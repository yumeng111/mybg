#yumeng
import datetime
import os
import sys
import ROOT
import shutil
import zlib

# #import Common.BaselineSelection as Baseline
# import Common.Utilities as Utilities
# import Common.ReportTools as ReportTools
# import Common.triggerSel as Triggers
# from Common.Setup import Setup
# from Corrections.Corrections import Corrections
# from Corrections.lumi import LumiFilter

# #from Common.Utilities import *

# import AnaProd.HH_bbWW.baseline as AnaBaseline
# import Common.BaselineSelection as CommonBaseline
# from Corrections.Corrections import Corrections

#ROOT.gROOT.ProcessLine(".include "+ os.environ['ANALYSIS_PATH'])
ROOT.gInterpreter.Declare(f'#include "/afs/cern.ch/user/y/yumeng/FLAF/include/BaselineGenSelection.h"')

def check_root_file(file_path):
    if not os.path.exists(file_path):
        print(f"Error: The file '{file_path}' does not exist.")
        return False

    root_file = ROOT.TFile.Open(file_path)
    if not root_file or root_file.IsZombie():
        print(f"Error: The file '{file_path}' could not be opened.")
        return False

    print(f"Success: The file '{file_path}' is valid")
    root_file.Close()
    return True

mass_point = "650"
file_path = f"/eos/cms/store/group/phys_higgs/HLepRare/skim_2024_v1/Run3_2022/GluGlutoRadiontoHHto2B2Vto2B2JLNu_M-{mass_point}/nano_0.root"
check_root_file(file_path)

df = ROOT.RDataFrame("Events", file_path)

df=df.Define("H_to_VV", """GetGenHVVCandidate(event, genLeptons, GenPart_pdgId, GenPart_daughters, GenPart_statusFlags, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenJet_p4, true)""")

df=df.Define("Leg0_Kind","H_to_VV.legs.at(0).leg_kind.at(0)")

df_prompt_muon = df.Filter("Leg0_Kind == Vleg::PromptMuon")
df_tau_muon = df.Filter("Leg0_Kind == Vleg::TauDecayedToMuon")

hist_prompt_muon = df_prompt_muon.Histo1D(("PromptMuon", "Prompt Muon Distribution", 100, 0, 100), "Leg0_Kind")
hist_tau_muon = df_tau_muon.Histo1D(("TauDecayedToMuon", "Tau Decayed to Muon Distribution", 100, 0, 100), "Leg0_Kind")

n_prompt_muon = df_prompt_muon.Count().GetValue()
n_tau_muon = df_tau_muon.Count().GetValue()

hist = ROOT.TH1F("genkind", "GenKind Distribution", 2, 0, 2)
hist.GetXaxis().SetBinLabel(1, "PromptMuon")
hist.GetXaxis().SetBinLabel(2, "TauDecayedToMuon")

hist.SetBinContent(1, n_prompt_muon)
hist.SetBinContent(2, n_tau_muon)

canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
hist.Draw("HIST")

canvas.SaveAs("genkind_distribution.png")
