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

root_file = ROOT.TFile.Open(file_path)
tree = root_file.Get("EventsNotSelected")
tree.Print()  # Print out all available branches in the tree
root_file.Close()