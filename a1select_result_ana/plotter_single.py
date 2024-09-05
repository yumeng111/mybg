import ROOT as rt
import sys, os, string, re
#from ROOT import TCanvas, TLegend, TChain, THStack, TH1D, TH1, kGray, kGreen, kRed, kBlue, kYellow, kMagenta, kSpring, kCyan, kTRUE, TLatex, gPad, TPad
from array import array
import math
from math import *

from tdrStyle import *
setTDRStyle()

file_name = '/eos/user/y/yumeng/ana/flaf/anaTuples/v7/Run3_2022/GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-250/nano_0.root'
#or put below inside make_plot
file = rt.TFile.Open(file_name)
if not file or file.IsZombie():
    print(f"Error: Cannot open file {file_name}")
tree_name = 'Events'
tree = file.Get(tree_name)
if not tree:
        print(f"Error: Tree {tree_name} not found in file {file_name}")
#else:
        #browser = rt.TBrowser()
        #browser.Add(tree)
        #tree.StartViewer()

#def make_plot(channel, var, bin, low, high = none, xlabel, xunits, prelim, setLogX, setLogY, passedSelection="", plotdir = ""):
def make_plot(file_name, channel, var, bin, low=None, high=None, varxlabel="", xunit="", prelim=False, setLogX=False, setLogY=True, cut="", plotdir=""):
    
    if low is None or high is None:
        hist_temp = rt.TH1D('hist_temp', 'hist_temp', 100, -1e6, 1e6)  # Initial large range
        tree.Draw(f"{var} >> hist_temp", cut, "goff")
        low = hist_temp.GetXaxis().GetXmin()
        high = hist_temp.GetXaxis().GetXmax()
    
    hist = rt.TH1D('hist', 'hist', bin, float(low), float(high))
    #hist.Sumw2()

    # Draw the tree variable into the hist with cut
    tree.Draw(f"{var} >> hist", cut, "goff")

    c1 = rt.TCanvas("c1", "c1", 800, 800)
    
    if setLogX: c1.SetLogx()
    if setLogY: c1.SetLogy()
    c1.SetRightMargin(0.03)

    hist.SetLineColor(rt.kBlue)
    hist.SetLineWidth(2)
    hist.Draw("hist")

    # Set histogram axis titles
    if xunit == '':
        hist.GetXaxis().SetTitle(varxlabel)
    else:
        hist.GetXaxis().SetTitle(f"{varxlabel} [{xunit}]")
    #hist.GetYaxis().SetTitle("Events / bin")
    hist.GetYaxis().SetTitle("Event Counts")

    # Draw additional plot elements (CMS labels, lumi, etc.)
    latex = rt.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.6 * c1.GetTopMargin())
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    #latex.DrawLatex(0.92, 0.94, "35.9 fb^{-1} (13 TeV)")
    latex.SetTextSize(0.85 * c1.GetTopMargin())
    latex.SetTextFont(62)
    latex.SetTextAlign(11)
    latex.DrawLatex(0.2, 0.84, "CMS")
    latex.SetTextSize(0.7 * c1.GetTopMargin())
    latex.SetTextFont(52)
    latex.SetTextAlign(11)
    if prelim:
        latex.DrawLatex(0.2, 0.78, "Preliminary")

    # Save the plot
    os.makedirs(plotdir, exist_ok=True)
    save = channel + "_" + var + '_' + cut
    #save = f"{channel}_{var}_{cut}".replace(" ", "_").replace("(", "").replace(")", "").replace("==", "_eq_").replace("&&", "_and_").replace("||", "_or_")
    c1.SaveAs(os.path.join(plotdir, f"hist_{save}.pdf"))
    #c1.SaveAs(os.path.join(plotdir, f"hist_{save}.png"))

plotdir = "plots/"

channel = "dl"
var = '(lep1_mass+lep2_mass)'
xunit = 'GeV'
bin = 50
low =0
high =0.2

varxlabel = 'm(ll)'
cut = ""

make_plot(file_name, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=True, cut=cut, plotdir=plotdir)

varxlabel = 'm(me)'
cut = "(lep1_type==1 && lep2_type==0)"
make_plot(file_name, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=True, cut=cut, plotdir=plotdir)

varxlabel = 'm(mumu)'
cut = "(lep1_type==1 && lep2_type==1)"
make_plot(file_name, channel, var, bin, low =0, high = 0.3, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=True, cut=cut, plotdir=plotdir)

varxlabel = 'm(ee)'
cut = "(lep1_type==0 && lep2_type==0)"
make_plot(file_name, channel, var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=True, cut=cut, plotdir=plotdir)
