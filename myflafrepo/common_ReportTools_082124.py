"""
import ROOT
def SaveReport(report, reoprtName="Report",printOut=False):
    cuts = [c for c in report]
    hist = ROOT.TH1D(reoprtName,reoprtName, len(cuts)+1, 0, len(cuts)+1)
    hist.GetXaxis().SetBinLabel(1, "Initial")
    hist.SetBinContent(1, cuts[0].GetAll())
    for c_id, cut in enumerate(cuts):
        hist.SetBinContent(c_id+2, cut.GetPass())
        hist.GetXaxis().SetBinLabel(c_id+2, cut.GetName())
        if(printOut):
            print(f"for the cut {cut.GetName()} there are {cut.GetPass()} events passed over {cut.GetAll()}, resulting in an efficiency of {cut.GetEff()}")
    return hist

#yumeng
"""
import ROOT

def SaveReport(report, reportName="Report", printOut=False):
    cuts = [c for c in report]
    num_cuts = len(cuts)

    # Histogram for the event numeber
    hist_pass = ROOT.TH1D(reportName + "_Number", reportName + " Number", num_cuts + 1, 0, num_cuts + 1)
    hist_pass.GetXaxis().SetBinLabel(1, "Initial")
    hist_pass.SetBinContent(1, cuts[0].GetAll())

    # Histogram for the efficiency
    hist_eff = ROOT.TH1D(reportName + "_Efficiency", reportName + " Efficiency", num_cuts + 1, 0, num_cuts + 1)
    hist_eff.GetXaxis().SetBinLabel(1, "Initial")
    hist_eff.SetBinContent(1, 1.0)
    hist_eff.SetBinError(1, 0.0) 

    for c_id, cut in enumerate(cuts):
        hist_pass.SetMinimum(0)
        hist_eff.SetMinimum(0)

        # Set the number of events passing the cut
        hist_pass.SetBinContent(c_id + 2, cut.GetPass())
        hist_pass.GetXaxis().SetBinLabel(c_id + 2, cut.GetName())

        p = cut.GetEff() / 100
        n = cut.GetAll()
        binomial_error = ROOT.TMath.Sqrt(p * (1 - p) / n)

        # Set the efficiency
        hist_eff.SetBinContent(c_id + 2, cut.GetEff()/100)
        hist_eff.GetXaxis().SetBinLabel(c_id + 2, cut.GetName())
        hist_eff.SetBinError(c_id + 2, binomial_error)


        if printOut:
            print(f"For the cut {cut.GetName()}, {cut.GetPass()} events passed out of {cut.GetAll()}, "
                  f"resulting in an efficiency of {cut.GetEff()/100}")

    return hist_pass, hist_eff

