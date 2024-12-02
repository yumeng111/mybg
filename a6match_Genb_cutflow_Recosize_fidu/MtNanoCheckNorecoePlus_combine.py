#yumeng
import array

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
header_path_AnalysisTools = "include/AnalysisTools.h"
header_path_studyYumengAnalysisTools = "./yumenginclude/studyYumengAnalysisTools.h"
ROOT.gInterpreter.Declare(f'#include "{header_path_RootExt}"')
ROOT.gInterpreter.Declare(f'#include "{header_path_GenLepton}"')
ROOT.gInterpreter.Declare(f'#include "{header_path_Gen}"')
ROOT.gInterpreter.Declare(f'#include "{header_path_Reco}"')
ROOT.gInterpreter.Declare(f'#include "{header_path_AnalysisTools}"') 
ROOT.gInterpreter.Declare(f'#include "{header_path_studyYumengAnalysisTools}"') 

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
        print(f"Error: The file '{file_path}' can not be opened.")
        return False

    print(f"Success: The file '{file_path}' is valid")
    root_file.Close()
    return True

def get_event_counts(df):
    try:
        # if df.Count().GetValue() > 0:
        #     return df.Count().GetValue()
        # else:
        #     return 0
        return df.Count().GetValue()
    except Exception as e:
        print(f"Error counting events: {e}")
        #traceback.print_exc()
        #return 0

def process_tree_Gen(file_path, tree_name):
    df = ROOT.RDataFrame(tree_name, file_path)

    df = df.Define("GenJet_p4", "RVecLV()")
    df = df.Define("GenPart_daughters", "GetDaughters(GenPart_genPartIdxMother)")
    
    if tree_name == "Events":
        df = df.Define("genLeptons", """reco_tau::gen_truth::GenLepton::fromNanoAOD(GenPart_pt, GenPart_eta,
                                    GenPart_phi, GenPart_mass, GenPart_genPartIdxMother, GenPart_pdgId,
                                    GenPart_statusFlags, event)""")
        df = df.Define("H_to_VV", """GetGenHVVCandidate(event, genLeptons, GenPart_pdgId, GenPart_daughters, 
                                  GenPart_statusFlags, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenJet_p4, true)""")
        df = df.Define("H_to_BB", """GetGenHBBCandidate(event, GenPart_pdgId, GenPart_daughters, GenPart_statusFlags, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenJet_p4, true)""")

        df = df.Define(f"{mychfull}_p4", f"yumengGetP4({mychfull}_pt, {mychfull}_eta, {mychfull}_phi, {mychfull}_mass)")

    else:
        df = df.Define("genLeptons", """reco_tau::gen_truth::GenLepton::fromNanoAOD(GenPart_pt, GenPart_eta,
                                    GenPart_phi, GenPart_mass, GenPart_genPartIdxMother, GenPart_pdgId,
                                    GenPart_statusFlags,-1)""")
        df = df.Define("H_to_VV", """GetGenHVVCandidate(-1,genLeptons, GenPart_pdgId, GenPart_daughters, 
                                  GenPart_statusFlags, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenJet_p4, true)""")
        df = df.Define("H_to_BB", """GetGenHBBCandidate(-1, GenPart_pdgId, GenPart_daughters, GenPart_statusFlags, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenJet_p4, true)""")
    
    df = df.Define("Leg0_Kind", "H_to_VV.legs.at(0).leg_kind.at(0)")
    df= df.Define("GenLep0_p4","H_to_VV.legs.at(0).leg_vis_p4.at(0)")
    
    df= df.Define("Genb0_p4","H_to_BB.leg_p4.at(0)")
    df= df.Define("Genb1_p4","H_to_BB.leg_p4.at(1)")

    if mychannel == "sle":
        #df_GenSL=df.Filter("Leg0_Kind==Vleg::PromptElectron|| Leg0_Kind==Vleg::TauDecayedToElectron", "")
        df_GenSL=df.Filter("Leg0_Kind==Vleg::PromptElectron", "")
        genSL_count = get_event_counts(df_GenSL)

        #df_Gen_Accept = df_GenSL.Filter("GenLep0_p4.pt()>10 && abs (GenLep0_p4.eta()) <2.5", "")
        df_Gen_Accept = df_GenSL.Filter("GenLep0_p4.pt()>5 && abs (GenLep0_p4.eta()) <2.5", "")
        gen_accept_count = get_event_counts(df_Gen_Accept)
        
    if mychannel == "slmuon":
        #df_GenSL=df.Filter("Leg0_Kind==Vleg::PromptMuon|| Leg0_Kind==Vleg::TauDecayedToMuon", "")
        df_GenSL=df.Filter("Leg0_Kind==Vleg::PromptMuon", "")
        genSL_count = get_event_counts(df_GenSL)

        #df_Gen_Accept = df_GenSL.Filter("GenLep0_p4.pt()>15 && abs (GenLep0_p4.eta()) <2.4", "")
        df_Gen_Accept = df_GenSL.Filter("GenLep0_p4.pt()>5 && abs (GenLep0_p4.eta()) <2.5", "")
        gen_accept_count = get_event_counts(df_Gen_Accept)
    
    df_Genb_Accept = df_Gen_Accept.Filter("Genb0_p4.pt()>20 && Genb1_p4.pt()>20 && abs (Genb0_p4.eta()) <2.5 && abs (Genb1_p4.eta()) <2.5", "")
    genb_accept_count = get_event_counts(df_Genb_Accept)

    df_GenWq_Accept = df_Genb_Accept.Filter("H_to_VV.legs.at(1).leg_p4.at(0).pt()>20 && H_to_VV.legs.at(1).leg_p4.at(1).pt()>20 && abs (H_to_VV.legs.at(1).leg_p4.at(0).eta()) <5 && abs (H_to_VV.legs.at(1).leg_p4.at(1).eta()) <5", "")
    #df_GenWq_Accept = df_Genb_Accept.Filter("H_to_VV.legs.at(1).leg_p4.at(0).pt()>20 && H_to_VV.legs.at(1).leg_p4.at(1).pt()>20", "")
    #df_GenWq_Accept = df_Genb_Accept.Filter("std::max(H_to_VV.legs.at(1).leg_p4.at(0).pt(), H_to_VV.legs.at(1).leg_p4.at(1).pt()) > 20 && abs (H_to_VV.legs.at(1).leg_p4.at(0).eta()) <5 && abs (H_to_VV.legs.at(1).leg_p4.at(1).eta()) <5", "")
    genWq_accept_count = get_event_counts(df_GenWq_Accept)

    return df_GenSL, df_Gen_Accept, df_Genb_Accept, df_GenWq_Accept, genSL_count, gen_accept_count, genb_accept_count, genWq_accept_count

def RecoHWWCandidateSelection(df):
    def add_filter(df, base_cut, and_cut, label):
        full_cut = f"Sum({base_cut} && {and_cut}) > 0"
        df = df.Filter(full_cut, label)
        return df, f"{base_cut} && {and_cut}"
    
    if mychannel == "sle":
        #find reco idx matched to gen by mindR
        df=df.Filter("nElectron>0", "HasRecoElectron")
        df=df.Define("RecoLepIdx", "FindMatching(GenLep0_p4, Electron_p4, 0.2)")
        df=df.Filter("RecoLepIdx>=0", "Reco_e_MatchGen")
        df = df.Filter("Electron_pt.at(RecoLepIdx)>5", "Reco_e_pt>5")#10
        df = df.Filter("Electron_eta.at(RecoLepIdx)<2.5", "Reco_e_eta<2.5")
        df = df.Filter("Electron_dz.at(RecoLepIdx)<0.1", "Reco_e_dz<0.1")
        df = df.Filter("Electron_dxy.at(RecoLepIdx)<0.05", "Reco_e_dxy<0.05")
        df = df.Filter("Electron_sip3d.at(RecoLepIdx) <= 8", "Reco_e_sip3d <= 8")

        # base_cut = "Electron_pt > 10"
        # df = df.Filter("Sum(Electron_pt > 10) > 0", "Reco_e_pt>10")
        # df, base_cut = add_filter(df, base_cut, "abs(Electron_eta) <2.5", "Reco_e_eta<2.5")
        # df, base_cut = add_filter(df, base_cut, "abs(Electron_dz) <0.1", "Reco_e_dz<0.1")
        # df, base_cut = add_filter(df, base_cut, "abs(Electron_dxy) <0.05", "Reco_e_dxy<0.05")
        # df, base_cut = add_filter(df, base_cut, "Electron_sip3d <= 8", "Reco_e_sip3d <= 8")
        # df, base_cut = add_filter(df, base_cut, "Electron_miniPFRelIso_all < 0.4", "Reco_e_miniPFRelIso_all < 0.4")
        # df, base_cut = add_filter(df, base_cut, "Electron_mvaIso_WP90", "Reco_e_mvaIso_WP90")
        # df, base_cut = add_filter(df, base_cut, "Electron_mvaIso_WP80", "Reco_e_mvaIso_WP80")

    # #double check:
    #     df = df.Filter("Sum(Electron_pt > 10) > 0", "Reco_e_pt>10")
    #     df = df.Filter("Sum(Electron_pt > 10 && abs(Electron_eta) <2.5) > 0", "Reco_e_eta<2.5")
    #     df = df.Filter("Sum(Electron_pt > 10 && abs(Electron_eta) <2.5 && abs(Electron_dz) <0.1) > 0", "Reco_e_dz<0.1")

    if mychannel== "slmuon":
        df=df.Filter("nMuon>0", "HasRecoMuon")
        df=df.Define("RecoLepIdx", "FindMatching(GenLep0_p4, Muon_p4, 0.01)")
        df=df.Filter("RecoLepIdx>=0", "Reco_mu_MatchGen")
        df = df.Filter("Muon_pt.at(RecoLepIdx)>5", "Reco_mu_pt>5")  #15
        df = df.Filter("Muon_eta.at(RecoLepIdx)<2.5", "Reco_mu_eta<2.5") #2.4
        df = df.Filter("Muon_dz.at(RecoLepIdx)<0.1", "Reco_mu_dz<0.1")
        df = df.Filter("Muon_dxy.at(RecoLepIdx)<0.05", "Reco_mu_dxy<0.05")
        df = df.Filter("Muon_sip3d.at(RecoLepIdx) <= 8", "Reco_mu_sip3d <= 8")

        # base_cut = "Muon_pt > 15"
        # df = df.Filter("Sum(Muon_pt > 15) > 0", "Reco_mu_pt>15")
        # df, base_cut = add_filter(df, base_cut, "abs(Muon_eta) <2.4", "Reco_mu_eta<2.4")
        # df, base_cut = add_filter(df, base_cut, "abs(Muon_dz) <0.1", "Reco_mu_dz<0.1")
        # df, base_cut = add_filter(df, base_cut, "abs(Muon_dxy) <0.05", "Reco_mu_dxy<0.05")
        # df, base_cut = add_filter(df, base_cut, "Muon_sip3d <= 8", "Reco_mu_sip3d <= 8")
        # df, base_cut = add_filter(df, base_cut, "Muon_miniPFRelIso_all < 0.4", "Reco_mu_miniPFRelIso_all < 0.4")
        # df, base_cut = add_filter(df, base_cut, "Muon_looseId", "Mu_looseId")
        # df, base_cut = add_filter(df, base_cut, "Muon_tightId", "Muon_tightId")
    return df

#20240910: This definition has not been used.
def RecoHWWJetSelection(df):
    df = df.Define("Jet_Incl", "Jet_pt>20 && abs(Jet_eta) < 2.5 && ( Jet_jetId & 2 )")
    df = df.Define("FatJet_Incl", "FatJet_pt >200 && abs(FatJet_eta) < 2.5 ) && ( FatJet_jetId & 2 ) && (FatJet_msoftdrop > 30) ")
    ######from here not modified
    df = df.Define("Jet_sel", """return RemoveOverlaps(Jet_p4, Jet_Incl,HwwCandidate.getLegP4s(), 0.4);""")
    df = df.Define("FatJet_sel", """return RemoveOverlaps(FatJet_p4, FatJet_Incl,HwwCandidate.getLegP4s(), 0.4);""")
    df = df.Define("Jet_cleaned", " RemoveOverlaps(Jet_p4, Jet_sel,{ {FatJet_p4[FatJet_sel][0], },}, 1, 0.8)")

    df = df.Define("n_eff_Jets", "(FatJet_p4[FatJet_sel].size()*2)+(Jet_p4[Jet_cleaned].size())")
    df = df.Define("n_eff_jets_SL","(is_SL && n_eff_Jets>=3)")
    df = df.Define("n_eff_jets_DL","(!is_SL && n_eff_Jets>=2)")

    return df.Filter(" (n_eff_jets_SL || n_eff_jets_DL)", "Reco bjet candidates")

def Report_cutflow_hist(canvas, mass_point, report, initial_count, genSL_count,  gen_accept_count, genb_accept_count, genWq_accept_count, skim_count, reportName="Cutflow", printOut=False):
    cuts = [c for c in report]
    num_cuts = len(cuts)

    hist_eff = ROOT.TH1D(reportName + f"__M-{mass_point}" + "_Efficiency", reportName+f" {ctype}"+" Efficiency" + "_for Masspoints", num_cuts + 4, 0, num_cuts + 4)
    #hist_eff.GetXaxis().SetBinLabel(1, "Initial_NoSkim")
    hist_eff.GetXaxis().SetBinLabel(1, "Initial_NANOAOD")
    hist_eff.SetBinContent(1, 1.0)
    hist_eff.SetBinError(1, 0.0) 
    hist_eff.SetMinimum(0)

    if ctype == "Relative":
        hist_eff.GetXaxis().SetBinLabel(2, f"GenSL{mychbrief}")
        hist_eff.SetBinContent(2, genSL_count/initial_count)

        hist_eff.GetXaxis().SetBinLabel(3, f"Gen_{mychbrief}_Lep_Accept")
        hist_eff.SetBinContent(3, gen_accept_count/genSL_count)

        hist_eff.GetXaxis().SetBinLabel(4, f"Gen_b_Accept")
        hist_eff.SetBinContent(4, genb_accept_count/gen_accept_count)

        hist_eff.GetXaxis().SetBinLabel(5, f"Gen_Wq_Accept")
        hist_eff.SetBinContent(5, genWq_accept_count/genb_accept_count)

        # hist_eff.GetXaxis().SetBinLabel(6, "Skim")
        # hist_eff.SetBinContent(6,skim_count/genWq_accept_count)

        for c_id, cut in enumerate(cuts):
                p = cut.GetEff() / 100
                n = cut.GetAll()
                binomial_error = ROOT.TMath.Sqrt(p * (1 - p) / n)
                hist_eff.SetBinContent(c_id + 6, cut.GetEff()/100)
                hist_eff.GetXaxis().SetBinLabel(c_id + 6, cut.GetName())
                hist_eff.SetBinError(c_id + 6, binomial_error)
    
    if ctype == "Cumulative":
        hist_eff.GetXaxis().SetBinLabel(2, f"GenSL{mychbrief}")
        hist_eff.SetBinContent(2, genSL_count/initial_count)

        hist_eff.GetXaxis().SetBinLabel(3, f"Gen_{mychbrief}_Lep_Accept")
        hist_eff.SetBinContent(3, gen_accept_count/initial_count)

        hist_eff.GetXaxis().SetBinLabel(4, f"Gen_b_Accept")
        hist_eff.SetBinContent(4, genb_accept_count/initial_count)

        hist_eff.GetXaxis().SetBinLabel(5, f"Gen_Wq_Accept")
        hist_eff.SetBinContent(5, genWq_accept_count/initial_count)

        # hist_eff.GetXaxis().SetBinLabel(6, "Skim")
        # hist_eff.SetBinContent(6,skim_count/initial_count)

        for c_id, cut in enumerate(cuts):
                p = cut.GetPass()/initial_count
                n = initial_count
                binomial_error = ROOT.TMath.Sqrt(p * (1 - p) / n)
                hist_eff.SetBinContent(c_id + 6, cut.GetPass()/initial_count)
                hist_eff.GetXaxis().SetBinLabel(c_id + 6, cut.GetName())
                hist_eff.SetBinError(c_id + 6, binomial_error)
    
    hist_eff.GetXaxis().SetLabelSize(0.03);  
    
    #hist_eff.GetXaxis().LabelsOption("v");
    canvas.SetBottomMargin(0.25);

    return hist_eff

def genDistribution_hist(mass_point,df):
    if var == "pt":
        df = df.Define("GenLep0_pt", "GenLep0_p4.pt()") 
        
        df= df.Define("neutrino0_p4","H_to_VV.legs.at(0).leg_p4.at(1)")
        
        df= df.Define("GenLep0E","GenLep0_p4.E()")
        df= df.Define("neutrino0E","neutrino0_p4.E()") 
        df= df.Define("NuoverNuPlusLep","neutrino0E/(neutrino0E+GenLep0E)")

        df= df.Define("neutrino0_pt","neutrino0_p4.pt()")
        df= df.Define("NuoverNuPlusLepPt","neutrino0_pt/(neutrino0_pt+GenLep0_pt)") 
        
        df = df.Define("VCandidateE", "H_to_VV.legs.at(0).cand_p4.E()")
        df= df.Define("VEMinusLegsE","VCandidateE-(neutrino0E+GenLep0E)")

        df = df.Define("WBoostVec", "ROOT::Math::Boost(H_to_VV.legs.at(0).cand_p4.BoostToCM())")
        df = df.Define ("Nu_p4W", "WBoostVec(neutrino0_p4)")
        df = df.Define ("GenLep_p4W", "WBoostVec(GenLep0_p4)")
        df = df.Define("NuoverNuPlusLepEW","Nu_p4W.E()/(Nu_p4W.E()+GenLep_p4W.E())")

        def linspace(start, stop, num):
            #Function to mimic numpy.
            step = (stop - start) / (num - 1)
            return [start + step * i for i in range(num)]

        bin_edges = (
            # linspace(0, 20, 6) +   
            # linspace(20, 50, 6)[1:] + 
            # linspace(50, 100, 3)[1:] +  
            # linspace(100, 200, 3)[1:] +
            # linspace(200, 1000, 3)[1:] +
            linspace(0, 500, 126) +  # 125 bins between 0 and 500
            linspace(500, 1000, 2)[1:] +  #[1: ]remove duplicate 500
            [2150] #need define oveflow bin upper edge (but not make overflow here), to get all statistics for the full range
        )

        bin_edges_array = array.array('d', bin_edges)

        # # pt_min = df.Min("GenLep0_pt").GetValue()
        # # pt_max = df.Max("GenLep0_pt").GetValue()
        pt_min = 0
        # #pt_max = 2150 for eventstree 5000, >1000 for eventstree 3000
        #variable changes
        # pt_max = 500
        # bins = int((pt_max - pt_min)/10)
        
        bins=20
        pt_max =1

        #model = ROOT.RDF.TH1DModel(f"GenLep0_pT_M-{mass_point}", f"GenLep0 pT for Masspoints ({mychannel})", len(bin_edges_array) - 1, bin_edges_array)
        #model = ROOT.RDF.TH1DModel(f"GenLep0_pT_M-{mass_point}", f"GenLep0 pT for Masspoints ({mychannel})", bins, pt_min, pt_max)
        #model = ROOT.RDF.TH1DModel(f"GenLevelComparison_M-{mass_point}", f"GenLevel Comparison for Masspoints ({mychannel})", bins, pt_min, pt_max)

        #variable changes
        #if want to make overflow: all the bin content out of the range to be drawn in the last bin.
        #Root will have the overflow bin automatically when x filled is bigger than defined range, and hist_pt.GetNbinsX() + 1 is the overflow bin number.
        # hist_pt = df.Histo1D(model, "GenLep0_pt")
        # #df.Foreach(lambda pt: hist_pt.Fill(pt), ["GenLep0_pt"])
        # hist_pt.GetXaxis().SetTitle(f"{mychannel} pT(GeV)")
        # hist_pt.GetXaxis().SetRangeUser(0, 500)

        # hist_pt = df.Histo1D(model, "NuoverNuPlusLep")
        # hist_pt.GetXaxis().SetTitle(f"E(Nu)/E(Nu+Lep) ({mychannel} channel)")
        # hist_pt.GetXaxis().SetRangeUser(0, 1)

        # hist_pt = df.Histo1D(model, "VEMinusLegsE")
        # hist_pt.GetXaxis().SetTitle(f"E(VCandidate)-E(Nu+Lep) ({mychannel} channel)")
        # hist_pt.GetXaxis().SetRangeUser(0, 500)
       
        # hist_pt = df.Histo1D(model, "NuoverNuPlusLepPt")
        # hist_pt.GetXaxis().SetTitle(f"pT(Nu)/(pT(Nu)+pT(Lep)) ({mychannel} channel)")
        # hist_pt.GetXaxis().SetRangeUser(0, 1)

        # hist_pt = df.Histo1D(model, "NuoverNuPlusLepEW")
        # hist_pt.GetXaxis().SetTitle(f"E(Nu)/E(Nu+Lep) (WFrame,{mychannel} channel)")
        # hist_pt.GetXaxis().SetRangeUser(0, 1)

        ###df=df.Filter("HVVCandidatept<700", "")
        hist_pt = df.Histo1D((f"mindR(GenLep, Fatjet){mass_point}", f"mindR(GenLep, Fatjet)({mychannel})", 50, 0, 5), "mindRGenlepFatjet")

        hist_pt.GetXaxis().SetTitle(f"mindR(GenLep, Fatjet), SL{mychbrief} channel)")

        hist_pt.SetBinContent(hist_pt.GetNbinsX(), hist_pt.GetBinContent(hist_pt.GetNbinsX()) + hist_pt.GetBinContent(hist_pt.GetNbinsX() + 1))

        return hist_pt
    
    if var == "eta":
        df = df.Define("GenLep0_eta", "GenLep0_p4.eta()") 

        # eta_min = df.Min("GenLep0_eta").GetValue()
        # eta_max = df.Max("GenLep0_eta").GetValue()
        eta_min = -6 #theta/eta 0/5~180/-5
        eta_max = 6

        bins=10 #normally bins=50, bins = int((eta_max - eta_min)/4)
        model = ROOT.RDF.TH1DModel(f"GenLep0_eta_M-{mass_point}", f"GenLep0 eta for Masspoints ({mychannel})", bins, eta_min, eta_max)
        hist_eta = df.Histo1D(model, "GenLep0_eta")
        #df.Foreach(lambda pt: hist_eta.Fill(eta), ["GenLep0_eta"])
        hist_eta.GetXaxis().SetTitle(f"{mychannel} eta")
        return hist_eta
    
def mindR_GenLepJet_hist(mass_point, df):
    df= df.Define("GenJet0_p4","H_to_VV.legs.at(1).leg_vis_p4.at(0)")
    df= df.Define("GenJet1_p4","H_to_VV.legs.at(1).leg_vis_p4.at(1)")

    #dPhi = ROOT.TVector2.Phi_mpi_pi(phi1 - phi2) #return phi in (-pi,pi)

    df = df.Define("mindR_lepJet", 
               """auto lep0p4V = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(GenLep0_p4.pt(), GenLep0_p4.eta(), GenLep0_p4.phi(), GenLep0_p4.mass());
                  auto jet0p4V = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(GenJet0_p4.pt(), GenJet0_p4.eta(), GenJet0_p4.phi(), GenJet0_p4.mass());
                  auto jet1p4V = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(GenJet1_p4.pt(), GenJet1_p4.eta(), GenJet1_p4.phi(), GenJet1_p4.mass());
                  
                  auto dR_lepjet0 = ROOT::Math::VectorUtil::DeltaR(lep0p4V, jet0p4V);
                  auto dR_lepjet1 = ROOT::Math::VectorUtil::DeltaR(lep0p4V, jet1p4V);
                  
                  return std::min(dR_lepjet0, dR_lepjet1);""")
    hist = df.Histo1D((f"mindR_GenLepJet_{mass_point}", f"mindR_GenLepJet for Masspoints({mychannel})", 50, 0, 6), "mindR_lepJet")
    hist = hist.GetValue()
    hist.GetXaxis().SetTitle("mindR_GenLepJet")
    yrange=0.06
    return hist, yrange

def Eff_RecoGenMatch_hist(mass_point,df):
    df = df.Define("GenLep0pt", "GenLep0_p4.pt()") 
    df = df.Define("VCandidatept", "H_to_VV.legs.at(0).cand_p4.Pt()")    
    df = df.Define("HVVCandidatept", "H_to_VV.cand_p4.Pt()")
    df = df.Define("GenLep0eta", "GenLep0_p4.eta()") 
    df = df.Define("VCandidateeta", "H_to_VV.legs.at(0).cand_p4.eta()")    
    df = df.Define("HVVCandidateeta", "H_to_VV.cand_p4.eta()")

    eta_max=5
    eta_min=-5
    
    particl="GenLep0"
    #particl="VCandidate"
    pt_max=2000
    
    # particl="HVVCandidate"
    # pt_max=3000
    
    pt_min = 0
    # #pt_max = 2150 for eventstree 5000, >1000 for eventstree 3000, 500 final; 3000 for HCand;
    #variable changes
    #bins = int((pt_max - pt_min)/10)  

    if var == "pt":
        var_min = pt_min 
        var_max = pt_max
    elif var == "eta":
        var_min = eta_min
        var_max = eta_max
    
    bins= int((var_max - var_min)/50)

    df=df.Define("RecoLepIdx", f"FindMatching(GenLep0_p4, {mychfull}_p4, 0.02)")
    #df=df.Filter("HVVCandidatept<700", "")
    df=df.Filter("HVVCandidatept>=700", "")
    df_sz=df.Filter("nElectron>0", "")
    df_matcdR=df.Filter("RecoLepIdx>=0", "")

    legend = ROOT.TLegend(0.7, 0.8, 0.9, 0.9)
    #variable changes, getvalue for entry

    h_before = df.Histo1D(("h_before", f"{particl}_{var} Distribution", bins, var_min, var_max), f"{particl}{var}").GetValue()
    h_after = df_sz.Histo1D(("h_after", f"HasRecoE_Eff vs {particl}_{var} (M-{mass_point})", bins, var_min, var_max), f"{particl}{var}").GetValue()
    hist_pt = h_before.Clone("h_ratio")
    hist_pt.SetTitle(f"Has_e_Eff vs {particl}_pt M-{mass_point}")
    hist_pt.GetYaxis().SetTitle(f"Eff (after HasRecoElectron)")
    legend.AddEntry(h_before, "Before Has_e", "l")
    legend.AddEntry(h_after, "After Has_e", "l")
    
    # h_before = df_sz.Histo1D(("h_before", f"{particl}_{var} Distribution", bins, var_min, var_max), f"{particl}{var}").GetValue()
    # h_after = df_matcdR.Histo1D(("h_after", f"Match_Eff vs {particl}_{var} (M-{mass_point})", bins, var_min, var_max), f"{particl}{var}").GetValue()
    # hist_pt = h_before.Clone("h_ratio")
    # hist_pt.SetTitle(f"Match_Eff vs {particl}_pt M-{mass_point}")
    # hist_pt.GetYaxis().SetTitle(f"Eff (after Match_RecoGenLep)")
    # legend.AddEntry(h_before, "Before Match", "l")
    # legend.AddEntry(h_after, "After Match", "l")
    
    h_before.SetBinContent(bins, h_before.GetBinContent(bins) + h_before.GetBinContent(bins + 1))
    h_after.SetBinContent(bins, h_after.GetBinContent(bins) + h_after.GetBinContent(bins + 1))

    # hist_pt.Divide(h_before)
    for bin in range(1, hist_pt.GetNbinsX() + 1):
        num = h_after.GetBinContent(bin)
        denom = h_before.GetBinContent(bin)

        #-1(float('inf')orfloat('nan')) diffirentiate from real 0 and 1 case, and denom == 0 -->nom==0
        if denom == 0:
            hist_pt.SetBinContent(bin, -1)
        else:
            hist_pt.SetBinContent(bin, num / denom)
        #print("bin:", bin, hist_pt.GetBinContent(bin))

    c = ROOT.TCanvas("c", "canvas", 1200, 600)
    c.Divide(2,1)
    c.SetGrid()

    pad1 = c.cd(1)
    ROOT.gStyle.SetOptStat(0)
    # make pads closer
    pad1.SetPad(0.01, 0.01, 0.5, 0.99)
    pad1.SetGrid()
    h_before.SetLineColor(ROOT.kGreen)
    h_after.SetLineColor(ROOT.kRed)
    h_before.Draw()
    h_after.Draw("same")
    # legend = ROOT.TLegend(0.7, 0.8, 0.9, 0.9)
    # legend.AddEntry(h_before, "Before Match", "l")
    # legend.AddEntry(h_after, "After Match", "l")
    legend.Draw()
    h_before.GetXaxis().SetTitle(f"{var}, {mychannel} Channel")
    h_before.GetYaxis().SetTitle(f"Counts")

    pad2 = c.cd(2)
    pad2.SetGrid()
    pad2.SetPad(0.5, 0.01, 0.99, 0.99)
    #variable changes
    #hist_pt.SetTitle(f"Match_Eff vs {particl}_pt M-{mass_point}")
    hist_pt.SetLineColor(ROOT.kBlack)
    hist_pt.Draw()
    
    yrange=1.1
    hist_pt.GetXaxis().SetTitle(f"{var}, {mychannel} Channel")
    #hist_pt.GetYaxis().SetTitle(f"Eff (after Match_RecoGenLep)")
    #hist_pt.GetXaxis().SetRangeUser(0, pt_max)
    hist_pt.GetYaxis().SetRangeUser(0, yrange)

    c.SaveAs(f"plots/{mychbrief}{mass_point}{particl}{var}_Eff_MtNano.png")

    return hist_pt,yrange

#change to checkNoRecoe
def EventsDistribution_hist(mass_point,df):  
    df = df.Define("HVVCandidatept", "H_to_VV.cand_p4.Pt()")
    df = df.Define("FatJet_p4", "yumengGetP4(FatJet_pt, FatJet_eta, FatJet_phi, FatJet_mass)")
    df = df.Define("mindRGenlepFatjet","FindMatchingdR(H_to_VV.legs.at(0).leg_vis_p4.at(0), FatJet_p4, 100)")
    df = df.Define("dRGenlepVh","ROOT::Math::VectorUtil::DeltaR(H_to_VV.legs.at(0).leg_vis_p4.at(0), H_to_VV.legs.at(1).cand_p4)")

    #variable changes:
    pt_min = 0
    pt_max = 1
    bins = int((pt_max - pt_min)/10)

    #variable changes:
    yrange=1
    ###df = df.Filter("mindRGenlepFatjet<6")
    ###df=df.Filter("HVVCandidatept<700", "")

    # hist_pt = df.Histo1D((f"mindR(GenLep, Fatjet){mass_point}", f"mindR(GenLep, Fatjet), SL{mychbrief} channel", 50, 0, 0.8), "mindRGenlepFatjet")
    # hist_pt.GetXaxis().SetTitle(f"mindR(GenLep, Fatjet), SL{mychbrief} channel)")
    # fname="MtNanodRFat" 

    # hist_pt = df.Histo1D((f"dR(GenLep, HardronicW){mass_point}", f"dR(GenLep, HardronicW)(H_pT<700,({mychannel})", 50, 0, 5), "dRGenlepVh")
    # #hist_pt = df.Histo1D((f"dR(GenLep, HardronicW){mass_point}", f"dR(GenLep, HardronicW)(H_pT>=700,({mychannel})", 50, 0, 1.5), "dRGenlepVh")
    # hist_pt.GetXaxis().SetTitle(f"dR(GenLep, HardronicW), SL{mychbrief} channel)")
    # fname="MtNanodRVh"

    # ######checkNoRecoe
    # df=df.Filter("nElectron==0", "")
    # df = df.Define("mindRGenlepFatjet_bin", "mindRGenlepFatjet <= 0.8 ? 0 : 1")
    # yrange=1.2
    # hist_pt = df.Histo1D((f"Non Recolep case{mass_point}", f"Cases of Non-RecoEle for M-3000, 5000", 2, 0, 2),"mindRGenlepFatjet_bin")
    # # Set axis titles and bin labels
    # hist_pt.GetXaxis().SetTitle(f"mindR(GenLep, Fatjet), SL{mychbrief} channel")
    # hist_pt.GetXaxis().SetBinLabel(1, "mindR(Genlep, Fatjets)<=0.8")
    # hist_pt.GetXaxis().SetBinLabel(2, "mindR(Genlep, Fatjets)>0.8")
    # fname="MtNanoCheckNorecoe"

    # df_numerator = df.Filter("mindRGenlepFatjet < 0.8")
    # efficiency = ROOT.TEfficiency(
    #     df_numerator.Histo1D(("numerator", "", 50, 0, 0.8), "mindRGenlepFatjet").GetValue(),
    #     df.Histo1D(("denominator", "", 50, 0, 0.8), "mindRGenlepFatjet").GetValue()
    # )
    # efficiency.SetTitle(f"Efficiency VS mindR(Genlep, Fatjets);mindR(GenLep, Fatjets);Efficiency")
    # efficiency.Draw("AP")   
    # hist_pt=efficiency
    # fname="MtNanoEffdR"

    thresholds = [i * 0.04 for i in range(11)]  # 0.0 to 0.8 in steps of 0.02
    
    total_count = df.Count().GetValue()
    
    efficiencies = []
    for threshold in thresholds:
        filtered_df = df.Filter(f"mindRGenlepFatjet >= {threshold} && nElectron>0")
        passing_count = filtered_df.Count().GetValue()
        efficiency = passing_count / total_count
        efficiencies.append(efficiency)
    
    hist_pt = ROOT.TGraph(len(thresholds))
    for i, (threshold, efficiency) in enumerate(zip(thresholds, efficiencies)):
        hist_pt.SetPoint(i, threshold, efficiency)
    
    # Set graph titles and labels
    hist_pt.SetTitle(f"Has_Reco_E_Efficiency vs mindR(Genlep, Fatjests);mindR threshold;Relative Efficiency")
    hist_pt.GetXaxis().SetLimits(0, 0.4)
    hist_pt.GetYaxis().SetRangeUser(0, 1.1)
    fname="MtNanoEffdRthresh"
   
    #hist_pt = hist_pt.GetValue()

    #hist_pt.SetBinContent(hist_pt.GetNbinsX(), hist_pt.GetBinContent(hist_pt.GetNbinsX()) + hist_pt.GetBinContent(hist_pt.GetNbinsX() + 1))

    return hist_pt, yrange, fname
######################################################################def of main functions starts:
######################################################################common configs:
######################################################################change save paths:
def main(mychannel, var, hist_type):
    #mass_points = [250, 450, 650, 1000, 3000, 5000]
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange]
    mass_points = [250,5000]

    histograms = []
    canvas = ROOT.TCanvas("canvas", "Mass Points Comparison", 800, 600)
    canvas.SetGrid()
    yrange=1
    legend = ROOT.TLegend(0.75, 0.9, 1, 1)
    legend.SetNColumns(2)

    for i, mass_point in enumerate(mass_points):
        file_path = f"/eos/user/y/yumeng/ana/flaf/2022run3nanoaod/GluGlutoRadiontoHHto2B2Vto2B2JLNu_M_hadd/hadd{mass_point}.root"
        if check_root_file(file_path):
            df_eventstree = ROOT.RDataFrame("Events", file_path)
            #df_notselected = ROOT.RDataFrame("EventsNotSelected", file_path)
            eventscount=get_event_counts(df_eventstree)
            #notselectedcount=get_event_counts(df_notselected)
        #initial_count = eventscount + notselectedcount
        initial_count = eventscount

        df_GenSL_Eventstree, df_GenAccept_Eventstree, df_GenbAccept_Eventstree, df_GenWqAccept_Eventstree, genSL_count_eventstree, gen_accept_count_eventstree, genb_accept_count_eventstree, genWq_accept_count_eventstree = process_tree_Gen(file_path, "Events")
        #df_GenSL_Notselected, df_GenAccept_Notselected, df_GenbAccept_Notselected, df_GenWqAccept_Notselected, genSL_count_notselected, gen_accept_count_notselected, genb_accept_count_notselected, genWq_accept_count_notselected = process_tree_Gen(file_path, "EventsNotSelected")

        if hist_type =="cutflow":
            hist=process_cutflow_hist(canvas, mass_point, initial_count, genSL_count_eventstree, gen_accept_count_eventstree, genb_accept_count_eventstree, genWq_accept_count_eventstree, df_GenWqAccept_Eventstree)
        if hist_type =="genDistribution":
            hist, yrange=process_genDistribution_hist(mass_point, df_GenSL_Eventstree, df_GenSL_Notselected)
        if hist_type == "mindR_GenLepJet":
            hist, yrange= mindR_GenLepJet_hist(mass_point, df_GenSL_Eventstree)
        if  hist_type == "Eff_RecoGenMatch":
            hist, yrange= Eff_RecoGenMatch_hist(mass_point, df_GenWqAccept_Eventstree)
        if  hist_type == "EventsDistribution":
            hist, yrange, fname= EventsDistribution_hist(mass_point, df_GenWqAccept_Eventstree)
        
        # # #if normalization and adjust nbinsx+1 in 2 places
        # #integral = hist.Integral()
        # integral = hist.Integral(1, hist.GetNbinsX() + 1)
        # if integral != 0:
        #     for bin_idx in range(1, hist.GetNbinsX() + 1):
        #         bin_content = hist.GetBinContent(bin_idx)
        #         if bin_content > 0:
        #             error = ROOT.TMath.Sqrt(bin_content)
        #         else:
        #             error = 0  # No error for empty bins
        #         hist.SetBinError(bin_idx, error)    
        #     hist.Scale(1.0 / integral)
        # else:
        #     print("Warning: Histogram has zero integral, cannot normalize.")
        # hist.GetYaxis().SetTitle("Normalized Distribution")
        # hist.GetYaxis().SetRangeUser(0, yrange)

        #print(f"debug histogram integral: {hist.Integral()}")
        #bin_dR04 = hist.FindBin(0.4)
        #integral_above_04 = hist.Integral(bin_dR04 + 1, hist.GetNbinsX())
        #print(f"Fraction passing mindR(GenLep, GenJet) > 0.4 for M-{mass_point}: {integral_above_04}")
        
        # overflow_bin_content = hist.GetBinContent(hist.GetNbinsX() + 1)
        # print("Integral of overflow bin:", overflow_bin_content)

        legend.AddEntry(hist, f"M-{mass_point}", "l")
        hist.SetLineColor(colors[i])
        histograms.append(hist)

        # if i == 0:
        #     hist.Draw("HIST E")
        #     ROOT.gStyle.SetOptStat(0)
        # else:
        #     hist.Draw("HIST E SAME")
        #     ROOT.gStyle.SetOptStat(0)

        if i == 0:
            hist.Draw("AL")  # First graph with axes and line
        
        else:
            hist.Draw("L SAME")  # Subsequent graphs overlayed
        
        # latex = ROOT.TLatex()
        # latex.SetTextSize(0.03) 
        # latex.SetTextAlign(22)   # Center align the text

        # for i in range(1, hist.GetNbinsX() + 1):
        #     bin_content = hist.GetBinContent(i)
        #     bin_center = hist.GetXaxis().GetBinCenter(i)
        #     latex.SetTextSize(0.025)
        #     latex.SetTextFont(42)
        #     latex.DrawLatex(bin_center, bin_content + 0.02, f"{bin_content:.2f}")

    legend.Draw()

    if hist_type =="cutflow":
            canvas.SaveAs(f"plots/{mychannel}MatcNano{ctype}_.png")
    if hist_type =="genDistribution": 
            canvas.SaveAs(f"plots/cutplusGen_{var}_Combine_{mychannel}.png")
    if hist_type == "mindR_GenLepJet":
            canvas.SaveAs(f"plots/mindR_GenLepJet_Combine_{mychannel}.png")
    if hist_type =="EventsDistribution": 
            canvas.SaveAs(f"plots/{mychbrief}{fname}_Combine.png")

######################################################################def of processing cutflow hist
#def process_cutflow_hist(canvas, mass_point, initial_count, genSL_count_eventstree, genSL_count_notselected, gen_accept_count_eventstree, gen_accept_count_notselected, genb_accept_count_eventstree, genb_accept_count_notselected, genWq_accept_count_eventstree, genWq_accept_count_notselected, df_GenWqAccept_Eventstree):
def process_cutflow_hist(canvas, mass_point, initial_count, genSL_count_eventstree, gen_accept_count_eventstree, genb_accept_count_eventstree, genWq_accept_count_eventstree, df_GenWqAccept_Eventstree):
    # genSL_count= genSL_count_eventstree + genSL_count_notselected
    # gen_accept_count= gen_accept_count_eventstree + gen_accept_count_notselected
    # genb_accept_count= genb_accept_count_eventstree + genb_accept_count_notselected
    # genWq_accept_count= genWq_accept_count_eventstree + genWq_accept_count_notselected
    genSL_count= genSL_count_eventstree
    gen_accept_count= gen_accept_count_eventstree
    genb_accept_count= genb_accept_count_eventstree
    genWq_accept_count= genWq_accept_count_eventstree
    skim_count = genWq_accept_count_eventstree

    df_cutflow=RecoHWWCandidateSelection(df_GenWqAccept_Eventstree)
    #df_cutflow=RecoHWWJetSelection(df_cutflow)

    report=df_cutflow.Report()

    hist = Report_cutflow_hist(canvas, mass_point, report, initial_count, genSL_count, gen_accept_count, genb_accept_count, genWq_accept_count, skim_count)
    return hist

######################################################################def of processing genDistribution hist
def process_genDistribution_hist(mass_point, df_GenSL_Eventstree, df_GenSL_Notselected):
    #if no need for adding Events and EventsNotSelected, only 2 lines below
    hist1 = genDistribution_hist(mass_point,df_GenSL_Eventstree)
    hist = hist1.GetValue()  # Ensure the histogram is fully evaluated
    
    # hist1_evaluated = hist1.GetValue() 
    # hist2 = genDistribution_hist(mass_point, df_GenSL_Notselected)
    # hist2_evaluated = hist2.GetValue()
    #hist = hist1_evaluated.Clone(f"GenLep0_add_{mass_point}")
    #hist.Add(hist2_evaluated)
    
    if var == "eta":
        yrange=0.5 
    if var == "pt":
        #variable changes
        #yrange=0.36
        #yrange=0.1
        #yrange=1.2
        yrange=1
    return hist, yrange

######################################################################def of processing mindR_GenLepJet hist
#no need def process_mindR_GenLepJet_hist(mass_point, df_GenSL_Eventstree):
    
######################################################################run!
mychannel = "sle"
#mychannel = "slmuon"
if mychannel == "sle":
    mychbrief="e"
    mychfull="Electron"
if mychannel == "slmuon":
    mychbrief="Mu"
    mychfull="Muon"

var = "pt"
#var = "eta"

#hist_type ="cutflow" #!!!NO normalization
#ctype = "Relative"
ctype = "Cumulative"
#hist_type = "genDistribution"
#hist_type ="mindR_GenLepJet"
#hist_type = "Eff_RecoGenMatch"
hist_type = "EventsDistribution"

main(mychannel, var, hist_type)



    

    
    


    






