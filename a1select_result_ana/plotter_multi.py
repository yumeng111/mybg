import ROOT as rt
import os
import pickle
from tdrStyle import setTDRStyle

setTDRStyle()

# Load the categorized samples and chains from the pickle file
with open('file_chains.pkl', 'rb') as f:
    categorized_samples, Event_Tree, entry_counts = pickle.load(f)

# Define a function to create combined plots
def make_combined_plot(signals, channel, var, bin, low, high, varxlabel="", xunit="", prelim=False, setLogX=False, setLogY=True, cut="", plotdir=""):
    colors = [rt.kBlue, rt.kRed, rt.kGreen, rt.kMagenta, rt.kCyan, rt.kYellow]
    c1 = rt.TCanvas("c1", "c1", 800, 800)
    if setLogX: c1.SetLogx()
    if setLogY: c1.SetLogy()
    c1.SetRightMargin(0.03)
    
    legend = rt.TLegend(0.7, 0.7, 0.9, 0.9)
    color_index = 0
    first_hist = True
    
    for M in [250, 450, 650, 1000, 3000, 5000]:
        sample_path = f"GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-{M}/nano_0.root"
        sample_name = sample_path.split('/')[0]
        print(f"Processing signal sample: {sample_name} from path: {sample_path}")
        if sample_name not in Event_Tree:
            print(f"Warning: {sample_name} not found in Event_Tree.")
        hist = rt.TH1D(f'hist_{sample_name}', f'hist_{sample_name}', bin, float(low), float(high))
        Event_Tree[sample_name].Draw(f"{var} >> hist_{sample_name}", cut, "goff")
        entries = hist.GetEntries()
        if entries == 0:
            print(f"Warning: {sample_name} histogram has no entries.")
        else:
            print(f"{sample_name} found in entry_counts with {entries} entries")
            
        integral = hist.Integral()
        if integral != 0:
            hist.Scale(1.0 / integral)   

        hist.SetLineColor(colors[color_index])
        hist.SetLineWidth(2)
        draw_option = "hist" if first_hist else "hist same"
        hist.Draw(draw_option)
        first_hist = False
            
        legend.AddEntry(hist, f"Signal: {sample_name}", "l")
        color_index += 1   

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

    # Update the canvas and save the plot
        c1.Update()
        c1.SaveAs(f"dl_signal_M-{M}.pdf")

# Plot directories
plotdir = "plots/"

# Define variables for plotting
var = '(lep1_mass+lep2_mass)'
xunit = 'GeV'
bin = 50
low = 0
high = 0.2

varxlabel = 'm(ll)'
cut = ""

# Create DL plots
dl_signals = categorized_samples["dl_signal"]
print("Creating DL plots")
print(f"DL Signals: {dl_signals}")
make_combined_plot(dl_signals, "dl", var, bin, low, high, varxlabel=varxlabel, xunit=xunit, prelim=False, setLogX=False, setLogY=True, cut=cut, plotdir=plotdir)
