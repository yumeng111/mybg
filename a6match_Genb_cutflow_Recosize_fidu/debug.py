import ROOT

def main():
    # Create a histogram with 3 bins
    hist = ROOT.TH1F("hist", "Example Histogram", 3, 0.0, 3.0)

    # Set the contents of the first three bins
    hist.SetBinContent(1, 3)  # First bin
    hist.SetBinContent(2, 2)  # Second bin
    hist.SetBinContent(3, 1)  # Third bin
    hist.SetBinContent(4, 4)

    # Access the underflow bin (bin 0)
    print("BinContent of underflow bin (0):", hist.GetBinContent(0))

    # Access the first bin (bin 1)
    print("BinContent of first bin (1):", hist.GetBinContent(1))
    print("hist.GetNbinsX():", hist.GetNbinsX())
    print("GetBinContent(hist.GetNbinsX() + 1):",hist.GetBinContent(hist.GetNbinsX() + 1))

if __name__ == "__main__":
    main()