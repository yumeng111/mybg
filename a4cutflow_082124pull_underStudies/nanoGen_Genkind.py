#yumeng
import datetime
import os
import sys
import ROOT
import shutil
import zlib

import argparse
import traceback

ROOT.gROOT.ProcessLine(".include " + os.environ['ANALYSIS_PATH'])
header_path_RootExt = "include/RootExt.h"
header_path_GenLepton = "include/GenLepton.h"
header_path_Gen = "include/BaselineGenSelection.h"
header_path_Reco = "include/BaselineRecoSelection.h"
header_path_HHbTag = "include/HHbTagScores.h"
ROOT.gInterpreter.Declare(f'#include "{header_path_RootExt}"')
ROOT.gInterpreter.Declare(f'#include "{header_path_GenLepton}"')
ROOT.gInterpreter.Declare(f'#include "{header_path_Gen}"')
ROOT.gInterpreter.Declare(f'#include "{header_path_Reco}"')

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

def get_event_counts(df):
    try:
        if df.Count().GetValue() > 0:
            return df.Count().GetValue()
        else:
            return 0
    except Exception as e:
        print(f"Error counting events: {e}")
        #traceback.print_exc()
        #return 0

def process_tree(file_path, tree_name):
    df = ROOT.RDataFrame(tree_name, file_path)

    df = df.Define("GenJet_p4", "RVecLV()")
    df = df.Define("GenPart_daughters", "GetDaughters(GenPart_genPartIdxMother)")
    
    if tree_name == "Events":
        df = df.Define("genLeptons", """reco_tau::gen_truth::GenLepton::fromNanoAOD(GenPart_pt, GenPart_eta,
                                    GenPart_phi, GenPart_mass, GenPart_genPartIdxMother, GenPart_pdgId,
                                    GenPart_statusFlags, event)""")
        df = df.Define("H_to_VV", """GetGenHVVCandidate(event, genLeptons, GenPart_pdgId, GenPart_daughters, 
                                  GenPart_statusFlags, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenJet_p4, true)""")
    else:
        df = df.Define("genLeptons", """reco_tau::gen_truth::GenLepton::fromNanoAOD(GenPart_pt, GenPart_eta,
                                    GenPart_phi, GenPart_mass, GenPart_genPartIdxMother, GenPart_pdgId,
                                    GenPart_statusFlags,-1)""")
        df = df.Define("H_to_VV", """GetGenHVVCandidate(-1,genLeptons, GenPart_pdgId, GenPart_daughters, 
                                  GenPart_statusFlags, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenJet_p4, true)""")
    
    df = df.Define("Leg0_Kind", "H_to_VV.legs.at(0).leg_kind.at(0)")

    counts = {
        "PromptElectron": get_event_counts(df.Filter("Leg0_Kind == Vleg::PromptElectron")),
        "PromptMuon": get_event_counts(df.Filter("Leg0_Kind == Vleg::PromptMuon")),
        "TauDecayedToElectron": get_event_counts(df.Filter("Leg0_Kind == Vleg::TauDecayedToElectron")),
        "TauDecayedToMuon": get_event_counts(df.Filter("Leg0_Kind == Vleg::TauDecayedToMuon")),
        "TauDecayedToHadrons": get_event_counts(df.Filter("Leg0_Kind == Vleg::TauDecayedToHadrons")),
        "NuorJet": get_event_counts(df.Filter("Leg0_Kind == Vleg::Nu||Leg0_Kind == Vleg::Jet"))
    }
    return counts

mass_points = [250, 450, 650, 1000, 3000, 5000]
colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange]

canvas = ROOT.TCanvas("canvas", "Mass Points Comparison", 800, 600)
canvas.SetGrid()

histograms = []

for i, mass_point in enumerate(mass_points):
    file_path = f"/eos/cms/store/group/phys_higgs/HLepRare/skim_2024_v1/Run3_2022/GluGlutoRadiontoHHto2B2Vto2B2JLNu_M-{mass_point}/nano_0.root"
    if check_root_file(file_path):
        df_eventstree = ROOT.RDataFrame("Events", file_path)
        df_notselected = ROOT.RDataFrame("EventsNotSelected", file_path)
        eventscount=get_event_counts(df_eventstree)
        notselectedcount=get_event_counts(df_notselected)
    initial_count = eventscount + notselectedcount

    events_counts = process_tree(file_path, "Events")
    not_selected_counts = process_tree(file_path, "EventsNotSelected")

    total_counts = {key: events_counts.get(key, 0) + not_selected_counts.get(key, 0) for key in events_counts}

    hist = ROOT.TH1F(f"genkind_{mass_point}", f"GenKind Distribution for Events_beforeSkim at Masspoints", 6, 0, 6)
    hist.GetXaxis().SetBinLabel(1, "PromptElectron")
    hist.GetXaxis().SetBinLabel(2, "PromptMuon")
    hist.GetXaxis().SetBinLabel(3, "TauDecayedToElectron")
    hist.GetXaxis().SetBinLabel(4, "TauDecayedToMuon")
    hist.GetXaxis().SetBinLabel(5, "TauDecayedToHadrons")
    hist.GetXaxis().SetBinLabel(6, "NuorJet")

    hist.SetBinContent(1, total_counts["PromptElectron"]/initial_count)
    hist.SetBinContent(2, total_counts["PromptMuon"]/initial_count)
    hist.SetBinContent(3, total_counts["TauDecayedToElectron"]/initial_count)
    hist.SetBinContent(4, total_counts["TauDecayedToMuon"]/initial_count)
    hist.SetBinContent(5, total_counts["TauDecayedToHadrons"]/initial_count)
    hist.SetBinContent(6, total_counts["NuorJet"]/initial_count)

    hist.SetLineColor(colors[i])
    hist.SetLineWidth(2)

    ROOT.gStyle.SetOptStat(0)
    
    histograms.append(hist)

    if i == 0:
        hist.Draw("HIST")
    else:
        hist.Draw("HIST SAME")

legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
for i, mass_point in enumerate(mass_points):
    legend.AddEntry(histograms[i], f"M-{mass_point}", "l")
legend.Draw()

canvas.SaveAs("genkind_beforeskim.png")