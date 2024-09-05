import ROOT
import numpy as np

# Create sample data using numpy arrays
x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=np.int32)
y = np.array([11, 12, 13, 14, 15, 16, 17, 18, 19, 20], dtype=np.int32)

# Create a ROOT RDataFrame with the number of entries equal to the length of the arrays
n = len(x)
df = ROOT.RDataFrame(n)

# Define x and y columns directly as arrays
df = df.Define("x_vec", f"ROOT::VecOps::RVec<int>({{{','.join(map(str, x))}}})")
df = df.Define("y_vec", f"ROOT::VecOps::RVec<int>({{{','.join(map(str, y))}}})")

# Define x and y using the rdfentry_ index to access elements in x_vec and y_vec
df = df.Define("x", "x_vec[rdfentry_]").Define("y", "y_vec[rdfentry_]")

# Apply filters to the individual elements
df = df.Filter("x > 6", "x > 3")
df = df.Filter("y < 18", "y < 18")
df = df.Filter("x % 2 == 0", "x is even")

# Generate the report
report = df.Report()

# Print the report
cuts = [c for c in report]
for c in cuts:
    #print(f"Cut name: {c.GetName()}, Events passing: {c.GetAll()}")
    print(f"Cut name: {c.GetName()}, Events passing: {c.GetPass()}")

# Function to save the report as a histogram
def SaveReport(report, reportName="Report", printOut=False):
    cuts = [c for c in report]
    if printOut:
        for c in cuts:
            print(f"{c.GetName()}: {c.GetAll()}")
            #print(f"{c.GetName()}: {c.GetPass()}")
    hist = ROOT.TH1D(reportName, reportName, len(cuts) + 1, 0, len(cuts) + 1)
    for i, c in enumerate(cuts):
        hist.GetXaxis().SetBinLabel(i + 1, c.GetName())
        #hist.SetBinContent(i + 1, c.GetAll())
        hist.SetBinContent(i + 1, c.GetPass())
    return hist

# Create a ROOT file and save the histogram
outputRootFile = ROOT.TFile("report_output.root", "RECREATE")
rep = SaveReport(report, reportName="Report", printOut=True)

# Debug: Check if the histogram is created successfully
if rep:
    print("Histogram 'rep' created successfully")
    outputRootFile.WriteTObject(rep, "Report", "Overwrite")
else:
    print("Error: Histogram 'rep' is None")
outputRootFile.Close()

# Open the ROOT file and retrieve the histogram
inputRootFile = ROOT.TFile("report_output.root", "READ")
rep_retrieved = inputRootFile.Get("Report")

# Debug: Check if the histogram is retrieved successfully
if rep_retrieved:
    print("Histogram 'rep' retrieved successfully")
    canvas = ROOT.TCanvas("canvas", "Cut Flow Report", 800, 600)
    rep_retrieved.Draw()
    canvas.SaveAs("cutflow_report.png")
else:
    print("Error: Retrieved histogram 'rep' is None")

inputRootFile.Close()
