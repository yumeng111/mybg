#yumeng
import datetime
import os
import sys
import ROOT
import shutil
import zlib

import argparse

import traceback

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

ROOT.gROOT.ProcessLine(".include "+ os.environ['ANALYSIS_PATH'])
header_path_RootExt ="include/RootExt.h"
header_path_GenLepton ="include/GenLepton.h"
header_path_Gen ="include/BaselineGenSelection.h"
header_path_Reco ="include/BaselineRecoSelection.h"
header_path_HHbTag ="include/HHbTagScores.h"
ROOT.gInterpreter.Declare(f'#include "{header_path_RootExt}"')
ROOT.gInterpreter.Declare(f'#include "{header_path_GenLepton}"')
ROOT.gInterpreter.Declare(f'#include "{header_path_Gen}"')
ROOT.gInterpreter.Declare(f'#include "{header_path_Reco}"')

# ROOT.gROOT.ProcessLine(".include "+ os.environ['ANALYSIS_PATH'] + '/include')
# ROOT.gInterpreter.Declare('#include "BaselineGenSelection.h"')

parser = argparse.ArgumentParser()
parser.add_argument('--particleFile', type=str,
                        default=f"{os.environ['ANALYSIS_PATH']}/config/pdg_name_type_charge.txt")

args = parser.parse_args()

ROOT.gInterpreter.ProcessLine(f"ParticleDB::Initialize(\"{args.particleFile}\");")

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

#df = ROOT.RDataFrame("Events", file_path)
df = ROOT.RDataFrame("EventsNotSelected", file_path)

############ nanoGenkind:
# Genjet is not available in eventnotselect:
# for var in ["GenJet", "GenJetAK8", "SubGenJetAK8"]:
#     df = df.Define(f"{var}_idx", f"CreateIndexes({var}_pt.size())")
#     df = df.Define(f"{var}_p4", f"GetP4({var}_pt,{var}_eta,{var}_phi,{var}_mass, {var}_idx)")

df = df.Define("GenJet_p4", "RVecLV()")
df = df.Define("GenPart_daughters", "GetDaughters(GenPart_genPartIdxMother)")
df = df.Define("genLeptons", """reco_tau::gen_truth::GenLepton::fromNanoAOD(GenPart_pt, GenPart_eta,
                                GenPart_phi, GenPart_mass, GenPart_genPartIdxMother, GenPart_pdgId,
                                GenPart_statusFlags, event)""")

df=df.Define("H_to_VV", """GetGenHVVCandidate(event, genLeptons, GenPart_pdgId, GenPart_daughters, GenPart_statusFlags, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenJet_p4, true)""")

df=df.Define("Leg0_Kind","H_to_VV.legs.at(0).leg_kind.at(0)")

df_prompt_muon = df.Filter("Leg0_Kind == Vleg::PromptElectron")
df_prompt_muon = df.Filter("Leg0_Kind == Vleg::PromptMuon")
df_tau_muon = df.Filter("Leg0_Kind == Vleg::TauDecayedToElectron")
df_tau_muon = df.Filter("Leg0_Kind == Vleg::TauDecayedToMuon")
#df_tau_muon = df.Filter("Leg0_Kind == Vleg::Other")

n_prompt_muon = 0
n_tau_muon = 0

###########for debug
try:
    # Count the events for each kind of muon
    n_prompt_muon = df_prompt_muon.Count().GetValue()
except Exception as e:
    print(f"Error counting prompt muons: {e}")

try:
    n_tau_muon = df_tau_muon.Count().GetValue()
except Exception as e:
    print(f"Error counting tau-decayed muons: {e}")

# try:
#     # Count the events for each kind of muon
#     n_prompt_muon = df_prompt_muon.Count().GetValue()
# except Exception as e:
#     print(f"Error counting prompt muons: {e}")
#     traceback.print_exc()

# try:
#     n_tau_muon = df_tau_muon.Count().GetValue()
# except Exception as e:
#     print(f"Error counting tau-decayed muons: {e}")
#     traceback.print_exc()

hist = ROOT.TH1F("genkind", "GenKind Distribution", 2, 0, 2)
hist.GetXaxis().SetBinLabel(1, "PromptMuon")
hist.GetXaxis().SetBinLabel(2, "TauDecayedToMuon")

hist.SetBinContent(1, n_prompt_muon) 
hist.SetBinContent(2, n_tau_muon)
canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
hist.Draw("HIST")

canvas.SaveAs("genkind_distribution.png")
# print("here")
# raise RuntimeError("stop")

##############for future: cutflow
# #For python import:
# #sys.path.append(os.environ['ANALYSIS_PATH'])
# #from Common.BaselineSelection import Initialize, DefineGenObjects, SelectRecoP4, CreateRecoP4
# #Initialize()
# #df = df.Filter("v_ops::pt(Electron_p4) > 10", "e_pt")

# df = df.Filter("Electron_pt > 10", "e_pt")

# report = df.Report()

# df.Count()

# report.Print()

# #with open(f"cutflow_report_mass_point_{mass_point}.txt", "w") as report_file:
# #    report_file.write(report.Print())
