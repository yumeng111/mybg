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
    #rt.gStyle.SetPalette(rt.kRainBow)  # Set the color palette
    #n_colors = rt.gStyle.GetNumberOfColors() 
    c1 = rt.TCanvas("c1", "c1", 1200, 800)
    if setLogX: c1.SetLogx()
    if setLogY: c1.SetLogy()
    c1.SetRightMargin(0.03)
    
    legend = rt.TLegend(0.55, 0.90, 0.95, 0.99)
    legend.SetNColumns(3);
    legend.SetTextSize(0.025)
    #legend.AddEntry(None, "Signal: GluGlutoRadiontoHHto2B2Vto2B2JLNu", "l")
    legend.SetHeader("Signal: GluGlutoRadiontoHHto2B2Vto2B2L2Nu" if channel == "dl" else "Signal: GluGlutoRadiontoHHto2B2Vto2B2JLNu", "C")
    color_index = 0
    first_hist = True

    histograms = {}
    
    # Plot signal chains
    #sample_path = f"GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-250/nano_0.root"
    for sample_path in signals:
        sample_name = sample_path.split('/')[0]
        mass_point = sample_name.split('_')[-1]
        if sample_name not in Event_Tree:
            print(f"Warning: {sample_name} not found in Event_Tree.")
            continue
        hist = rt.TH1D(f'hist_{sample_name}', f'hist_{sample_name}', bin, float(low), float(high))
        Event_Tree[sample_name].Draw(f"{var} >> hist_{sample_name}", cut, "goff")
        entries = hist.GetEntries()
        
        if sample_name in entry_counts:
            sample_entries = entry_counts[sample_name]
        else:
            print(f"Warning: {sample_name} not found in entry_counts.")
            sample_entries = 1  # Avoid division by zero, set to 1 if not found
        if entries == 0:
            print(f"Warning: {sample_name} histogram has no entries.")
            continue
        else:
            percentage = (entries / sample_entries) * 100
            print(f"{sample_name} {sample_entries} found in entry_counts with {percentage:.2f}%({entries}) entries")

        integral =  hist.Integral()
        if integral != 0:
             hist.Scale(1.0 / integral) 

        #color = rt.gStyle.GetColorPalette(color_index)
        #hist.SetLineColor(color)
        hist.SetLineColor(color_index+1)
        hist.SetLineWidth(2)

        #if var == '(lep1_mass+lep2_mass)': 
        if xunit == '':
            hist.GetXaxis().SetTitle(varxlabel)
        else:
            hist.GetXaxis().SetTitle(f"{varxlabel} [{xunit}]")
        hist.GetYaxis().SetTitle("Normalized Distribution")
        #hist.SetTitle(cut);
        # if (var == '(lep1_genLep_kind)'| var == '(lep2_genLep_kind)'): 
        #     if xunit == '':
        #         hist.GetXaxis().SetTitle(varxlabel)
        #     else:
        #         hist.GetXaxis().SetTitle(f"{varxlabel} [{xunit}]")
        # hist.GetYaxis().SetTitle("Normalized Distribution")
        #     #hist.SetTitle(cut);

        hist.SetMinimum(0)
        hist.SetMaximum(1)
        draw_option = "HIST" if first_hist else "HIST SAME"
        histograms[f'hist_{sample_name}']=hist
        histograms[f'hist_{sample_name}'].Draw(draw_option)

        legend.AddEntry(histograms[f'hist_{sample_name}'], mass_point, "l")
        
        first_hist = False
        color_index += 1
    

    # Plot background chains
    for sample_path in backgrounds:
        sample_name = sample_path.split('/')[0]
        if sample_name not in Event_Tree:
            print(f"Warning: {sample_name} not found in Event_Tree.")
            continue
        hist = rt.TH1D(f'hist_{sample_name}', f'hist_{sample_name}', bin, float(low), float(high))
        Event_Tree[sample_name].Draw(f"{var} >> hist_{sample_name}", cut, "goff")
        bgentries = hist.GetEntries()

        if sample_name in entry_counts:
            sample_entries = entry_counts[sample_name]
        else:
            print(f"Warning: {sample_name} not found in entry_counts.")
            sample_entries = 1  # Avoid division by zero, set to 1 if not found
        if bgentries == 0:
            print(f"Warning: {sample_name} histogram has no entries.")
            continue
        else:
            percentage = (bgentries / sample_entries) * 100
            print(f"{sample_name} {sample_entries} found in entry_counts with {percentage:.2f}%({bgentries}) entries")

        integral =  hist.Integral()
        if integral != 0:
            hist.Scale(1.0 / integral) 

        hist.SetLineColor(color_index+1)
        hist.SetLineWidth(2)

        hist.SetMarkerStyle(color_index+17)
        hist.SetMarkerColor(color_index + 1)
        hist.SetMarkerSize(2)

        histograms[f'hist_{sample_name}']=hist
        histograms[f'hist_{sample_name}'].Draw("P HIST SAME")

        legend.AddEntry(histograms[f'hist_{sample_name}'], sample_name, "lp")
        color_index += 1

    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.SetFillStyle(0)
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
    
    c1.Update()
    
    os.makedirs(plotdir, exist_ok=True)
    save = channel + "_" + var + '_' + cut
    c1.SaveAs(os.path.join(plotdir, f"histcomb_{save}.pdf"))

# Plot directories
plotdir = "plots/"

#####var
"""
var = '(lep1_mass+lep2_mass)'
bin = 50
low = 0
xunit = 'GeV'

# Create DL plots
dl_signals = categorized_samples["dl_signal"]
dl_backgrounds = categorized_samples["dl_bg"]

high = 0.25
varxlabel = 'm(ll)'
cut = ""
make_combined_plot(dl_signals, dl_backgrounds, "dl", var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

varxlabel = 'm(mumu)'
cut = "(lep1_type==1 && lep2_type==1)"
make_combined_plot(dl_signals, dl_backgrounds, "dl", var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

high = 0.12
varxlabel = 'm(me)'
cut = "(lep1_type==1 && lep2_type==0)"
make_combined_plot(dl_signals, dl_backgrounds, "dl", var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

high = 0.0015
varxlabel = 'm(ee)'
cut = "(lep1_type==0 && lep2_type==0)"
make_combined_plot(dl_signals, dl_backgrounds, "dl", var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

# Create SL plots
sl_signals = categorized_samples["sl_signal"]
sl_backgrounds = categorized_samples["sl_bg"]

high = 0.25
varxlabel = 'm(ll)'
cut = ""
make_combined_plot(sl_signals, dl_backgrounds, "sl", var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

varxlabel = 'm(mumu)'
cut = "(lep1_type==1 && lep2_type==1)"
make_combined_plot(sl_signals, dl_backgrounds, "sl", var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

high = 0.12
varxlabel = 'm(me)'
cut = "(lep1_type==1 && lep2_type==0)"
make_combined_plot(sl_signals, dl_backgrounds, "sl", var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

high = 0.0015
varxlabel = 'm(ee)'
cut = "(lep1_type==0 && lep2_type==0)"
make_combined_plot(sl_signals, dl_backgrounds, "sl", var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)
"""

"""
#######var
var = '(lep1_genLep_kind)'
xunit=""
bin = 7
low = -1
high = 6

# Create DL plots
dl_signals = categorized_samples["dl_signal"]
dl_backgrounds = categorized_samples["dl_bg"]
channel="dl"

cut = "(lep1_type==1 && lep2_type==1)"
varxlabel = 'mumu: lep1_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(dl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

cut = "(lep1_type==1 && lep2_type==0)"
varxlabel = 'mu,e: lep1_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(dl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

cut = "(lep1_type==0 && lep2_type==0)"
varxlabel = 'ee: lep1_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(dl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

# Create SL plots
sl_signals = categorized_samples["sl_signal"]
sl_backgrounds = categorized_samples["sl_bg"]
channel="sl"

cut = "(lep1_type==1 && lep2_type==-1)"
varxlabel = 'mu: lep1_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(sl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

cut = "(lep1_type==0 && lep2_type==-1)"
varxlabel = 'e: lep1_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(sl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

cut = "(lep1_type==1 && lep2_type==1)"
varxlabel = 'mumu: lep1_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(sl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

cut = "(lep1_type==1 && lep2_type==0)"
varxlabel = 'mu,e: lep1_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(sl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

cut = "(lep1_type==0 && lep2_type==0)"
varxlabel = 'ee: lep1_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(sl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

var = '(lep2_genLep_kind)'
xunit=""
bin = 7
low = -1
high = 6

# Create DL plots
dl_signals = categorized_samples["dl_signal"]
dl_backgrounds = categorized_samples["dl_bg"]
channel="dl"

cut = "(lep1_type==1 && lep2_type==1)"
varxlabel = 'mumu: lep2_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(dl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

cut = "(lep1_type==1 && lep2_type==0)"
varxlabel = 'mu,e: lep2_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(dl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

cut = "(lep1_type==0 && lep2_type==0)"
varxlabel = 'ee: lep2_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(dl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

# Create SL plots
sl_signals = categorized_samples["sl_signal"]
sl_backgrounds = categorized_samples["sl_bg"]
channel="sl"

cut = "(lep1_type==1 && lep2_type==-1)"
varxlabel = 'mu: lep2_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(sl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

cut = "(lep1_type==0 && lep2_type==-1)"
varxlabel = 'e: lep2_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(sl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

cut = "(lep1_type==1 && lep2_type==1)"
varxlabel = 'mumu: lep2_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(sl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

cut = "(lep1_type==1 && lep2_type==0)"
varxlabel = 'mu,e: lep2_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(sl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)

cut = "(lep1_type==0 && lep2_type==0)"
varxlabel = 'ee: lep2_genLep_kind'
print("\n", channel, varxlabel)
make_combined_plot(sl_signals, dl_backgrounds, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=False, cut=cut, plotdir=plotdir)
"""


