import ROOT

# Create a new ROOT file
#root_file = ROOT.TFile("test.root", "RECREATE")

# Create a TTree
tree = ROOT.TTree("tree", "example tree")

# Define the branches
a = ROOT.std.vector('int')()
b = ROOT.std.vector('int')()
c = ROOT.std.vector('int')()

tree.Branch("a", a)
tree.Branch("b", b)
tree.Branch("c", c)

# Fill the tree with the specified data
rows = [(1, 2, 3), (7, 8, 9)]

for row in rows:
    a.clear()
    b.clear()
    c.clear()
    a.push_back(row[0])
    b.push_back(row[1])
    c.push_back(row[2])
    tree.Fill()

# Write the tree to the file and close the file
#root_file.Write()
#root_file.Close()

# Reopen the file and get the tree
#file = ROOT.TFile("test.root", "READ")
#tree = file.Get("tree")

# Print the content of the tree entries
print("\nTree contents:")
#entries = tree.GetEntries()
#for i in range(entries):
#    tree.GetEntry(i)
tree.Scan()

# Draw the histogram
hist = ROOT.TH1F("hist", "Histogram of b", 10, 0, 10)
tree.Draw("c >> hist", "b>0", "goff")
canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
hist.Draw()
canvas.SaveAs("histogram.png")

# Print the contents of the histogram
print("Histogram contents:")
for bin in range(1, hist.GetNbinsX() + 1):
    print(f"Bin {bin}: {hist.GetBinContent(bin)}")

# Clean up
#file.Close()