def PrintDecayChain(df, evtIds, outFile):
    printed_event_count = 0
    max_events_to_print = 10

    if len(evtIds) > 0:
        df = df.Filter(f"static const std::set<ULong64_t> evts = {{ {evtIds} }}; return evts.count(event) > 0;")

    df = df.Range(max_events_to_print)
    
    df = df.Define('printer', f'''PrintDecayChain(event, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags,
                                  GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenPart_status, GenPart_daughters,
                                  "{outFile}")''')
    df.Histo1D("printer").GetValue()

    # for entry in df.AsNumpy(['printer'])['printer']:
    #     if printed_event_count >= max_events_to_print:
    #         break
        
    #     print(f"Processed event number {printed_event_count + 1}")

    #     printed_event_count += 1

if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument('--inFile', type=str)
    parser.add_argument('--outFile', type=str)
    parser.add_argument('--evtIds', type=str, default='')
    parser.add_argument('--particleFile', type=str,
                        default=f"{os.environ['ANALYSIS_PATH']}/config/pdg_name_type_charge.txt")
    args = parser.parse_args()

    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.ProcessLine(".include "+ os.environ['ANALYSIS_PATH'])
    ROOT.gROOT.ProcessLine('#include "include/yumengGenTools.h"')
    ROOT.gInterpreter.ProcessLine(f"ParticleDB::Initialize(\"{args.particleFile}\");")

    if os.path.exists(args.outFile):
        os.remove(args.outFile)
    outDir = os.path.dirname(args.outFile)
    if len(outDir) > 0 and not os.path.exists(outDir):
        os.makedirs(outDir)

    df = ROOT.RDataFrame("Events", args.inFile)
    df = df.Define("GenPart_daughters", "GetDaughters(GenPart_genPartIdxMother)")
    PrintDecayChain(df, args.evtIds, args.outFile)

#python ../../Common/PrintGenChain.py --inFile /eos/cms/store/group/phys_higgs/HLepRare/skim_2024_v1/Run3_2022/GluGlutoRadiontoHHto2B2Vto2B2JLNu_M-650/nano_0.root --outFile ./yumengprintchainM650_677.txt --evtIds 677

#python ../../Common/yumengPrintGenChain.py --inFile /eos/cms/store/group/phys_higgs/HLepRare/skim_2024_v1/Run3_2022/GluGlutoRadiontoHHto2B2Vto2B2JLNu_M-650/nano_0.root --outFile ./yumengConditionprintchainM650.txt