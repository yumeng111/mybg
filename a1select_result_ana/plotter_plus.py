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

# combined plots
def make_combined_plot(signals, backgrounds, channel, var, bin, low, high, varxlabel="", xunit="", prelim=False, setLogX=False, setLogY=True, cut="", plotdir=""):
    colors = [rt.kBlue, rt.kRed, rt.kGreen, rt.kMagenta, rt.kCyan, rt.kYellow]
    c1 = rt.TCanvas("c1", "c1", 800, 800)
    if setLogX: c1.SetLogx()
    if setLogY: c1.SetLogy()
    c1.SetRightMargin(0.03)
    
    legend = rt.TLegend(0.7, 0.7, 0.9, 0.9)
    color_index = 0
    
    # Plot signal chains
    for sample_path in signals:
        sample_name = sample_path.split('/')[0]
        if sample_name not in Event_Tree:
            print(f"Warning: {sample_name} not found in Event_Tree.")
            continue
        hist = rt.TH1D(f'hist_{sample_name}', f'hist_{sample_name}', bin, float(low), float(high))
        Event_Tree[sample_name].Draw(f"{var} >> hist_{sample_name}", cut, "goff")
        entries = hist.GetEntries()
        if entries == 0:
            print(f"Warning: {sample_name} histogram has no entries.")
            continue
        else:
            print(f"{sample_name} found in entry_counts with {entries} entries")
        
        hist.SetLineColor(colors[color_index % len(colors)])
        hist.SetLineWidth(2)
        if color_index == 0:
            hist.Draw("hist")
            hist.GetXaxis().SetTitle(varxlabel)
            if xunit == '':
                hist.GetXaxis().SetTitle(varxlabel)
            else:
                hist.GetXaxis().SetTitle(f"{varxlabel} [{xunit}]")
            hist.GetYaxis().SetTitle("Event Counts")
        else:
            hist.Draw("hist same")
        legend.AddEntry(hist, f"Signal: {sample_name}", "l")
        color_index += 1
    c1.SaveAs("my.pdf")
    
    """"
    # Plot background chains
    for sample_path in backgrounds:
        sample_name = sample_path.split('/')[0]
        if sample_name not in Event_Tree:
            print(f"Warning: {sample_name} not found in Event_Tree.")
            continue
        hist = rt.TH1D(f'hist_{sample_name}', f'hist_{sample_name}', bin, float(low), float(high))
        Event_Tree[sample_name].Draw(f"{var} >> hist_{sample_name}", cut, "goff")
        bgentries = hist.GetEntries()
        if bgentries == 0:
            print(f"Warning: {sample_name} histogram has no entries.")
            continue
        else:
            print(f"{sample_name} found in entry_counts with {bgentries} entries")
        hist.SetLineColor(colors[color_index % len(colors)])
        hist.SetLineWidth(2)
        hist.Draw("hist same")

        legend.AddEntry(hist, sample_name, "l")
        color_index += 1
    """""

    legend.Draw()

    # Draw additional plot elements (CMS labels, lumi, etc.)
    latex = rt.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.6 * c1.GetTopMargin())
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(0.85 * c1.GetTopMargin())
    latex.SetTextFont(62)
    latex.SetTextAlign(11)
    latex.DrawLatex(0.2, 0.84, "CMS")
    latex.SetTextSize(0.7 * c1.GetTopMargin())
    latex.SetTextFont(52)
    latex.SetTextAlign(11)
    if prelim:
        latex.DrawLatex(0.2, 0.78, "Preliminary")

""""
    # Save the plot
    os.makedirs(plotdir, exist_ok=True)
    save = channel + "_" + var + '_' + cut
    c1.SaveAs(os.path.join(plotdir, f"hist_{save}.pdf"))
    save = f"{channel}_{var}_{cut}".replace(" ", "_").replace("(", "").replace(")", "").replace("==", "_eq_").replace("&&", "_and_").replace("||", "_or_")
    c1.SaveAs(os.path.join(plotdir, f"combined_hist_{save}.pdf"))
"""

# Plot directories
plotdir = "plots/"

# Define variables for plotting
var = '(lep1_mass+lep2_mass)'
xunit = 'GeV'
bin = 50
low = 0
high = 0.12

varxlabel = 'm(ll)'
cut = ""

# Create DL plots
dl_signals = categorized_samples["dl_signal"]
dl_backgrounds = categorized_samples["dl_bg"]
make_combined_plot(dl_signals, dl_backgrounds, "dl", var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=True, cut=cut, plotdir=plotdir)

"""
# Create SL plots
sl_signals = categorized_samples["sl_signal"]
sl_backgrounds = categorized_samples["sl_bg"]
varxlabel = 'm(jet, lepton)'
cut = ""
make_combined_plot(sl_signals, sl_backgrounds, "sl", var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=True, cut=cut, plotdir=plotdir)
"""