from ROOT import TFile,TTree,TCanvas,TH1F,TH1D,TH2F,TLatex,TLine,TMath,TGraphAsymmErrors,TEfficiency,TGraph,TGraphErrors,TLegend,gROOT,gPad,gStyle, kBlack, kWhite, TF1, TPaveStats
import langaus
from stripBox import getStripBox
import optparse
import myStyle
import myFunctions as mf
import os
import numpy as np
from array import array 

gROOT.SetBatch( True )
gStyle.SetOptFit(1011)
colors = myStyle.GetColors(True)

## Defining Style
myStyle.ForceStyle()


class HistoInfo:
    def __init__(self, inHistoName, f, outHistoName, yMax=150,
                 xlabel="", ylabel="Amplitude [mV]",
                 sensor="", center_position = 0.0):
        self.inHistoName = inHistoName
        self.f = f
        self.outHistoName = outHistoName
        self.yMax = yMax
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.sensor = sensor
        self.center_position = center_position
        if("KOJI" in sensor):
            self.th2 = self.getTH2(f, inHistoName, sensor).RebinX(2)
            print(" (!!) Using RebinX(2) to handle low stat bins!")
        else:
            self.th2 = self.getTH2(f, inHistoName, sensor)
        self.th1 = self.getTH1(outHistoName)

    def getTH2(self, f, name, sensor, axis='zx'):
        th3 = f.Get(name)
        th2 = th3.Project3D(axis)

        return th2

    def getTH1(self, hname):
        htitle = ";%s;%s"%(self.xlabel, self.ylabel)
        nxbin = self.th2.GetXaxis().GetNbins()
        xmin, xmax = mf.get_shifted_limits(self.th2, self.center_position)

        # Create and define th1 default style
        th1 = TH1D(hname, htitle, nxbin, xmin, xmax)
        # th1.SetStats(0)
        th1.SetMinimum(0.1)
        th1.SetMaximum(self.yMax)
        # th1.SetLineWidth(3)
        # th1.SetLineColor(kBlack)
        # # th1.GetXaxis().SetRangeUser(-xlength,xlength)

        return th1


# Construct the argument parser
parser = optparse.OptionParser("usage: %prog [options]\n")
parser.add_option('-f', dest='ftbf_dataset', default = "", help="FTBF Dataset")
parser.add_option('-l', dest='laser_dataset', default = "", help="Laser Dataset")

options, args = parser.parse_args()

ftbf_loc = options.ftbf_dataset
laser_loc = options.laser_dataset

ftbf_dataset = [i for i in ftbf_loc.split('/') if "HPK" in i][0]
sensor_Geometry = myStyle.GetGeometry(ftbf_dataset)

outdir = "../Paper_plots/"
if not os.path.exists(outdir):
    myStyle.CreateFolder("../", "Paper_plots/")
outdir_fits = "../Paper_plots/tr_vs_attn_fits/"
if not os.path.exists(outdir_fits):
    myStyle.CreateFolder("../", "Paper_plots/tr_vs_attn_fits/")

laser_bias_voltage = "205"
attenuations = ["91P5", "92", "92P5", "93", "93P5", "94"]
attenuations_val = [91.5, 92, 92.5, 93, 93.5, 94]
tr_vs_attn = []
trErr_vs_attn = []

canvas2 = TCanvas("canvas2", "Histogram Fit", 1200, 900)
for attenuation in attenuations:
    rootfile = TFile(f'{laser_loc}/HPK_W5_17_2_{laser_bias_voltage}V_{attenuation}attn/HPK_W5_17_2_{laser_bias_voltage}V_{attenuation}attn_Analyze.root')
    hist = rootfile.Get("weighted2_timeDiff_tracker")
    myMean = hist.GetMean()
    myRMS = hist.GetRMS()
    gaussian = TF1("gaussian", "gaus")
    gaussian.SetRange(myMean-1.5*myRMS,myMean+1.5*myRMS)
    hist.Fit(gaussian, "R")
    hist.Draw("hist")
    gaussian.Draw("same")
    canvas2.SaveAs(f'{outdir_fits}fit_{laser_bias_voltage}V_{attenuation}attn.png')
    canvas2.Clear()
    tr = gaussian.GetParameter(2)*1000
    tr_err = gaussian.GetParError(2)*1000
    tr_vs_attn.append(tr)
    trErr_vs_attn.append(tr_err)
    rootfile.Close()

ftbf_rootfile = TFile(f'{ftbf_loc}/{ftbf_dataset}_Analyze.root')
ftbf_hist = ftbf_rootfile.Get("weighted2_timeDiff_tracker")
myMean = ftbf_hist.GetMean()
myRMS = ftbf_hist.GetRMS()
ftbf_gaussian = TF1("gaussian", "gaus")
ftbf_gaussian.SetRange(myMean-1.5*myRMS,myMean+1.5*myRMS)
ftbf_hist.Fit(ftbf_gaussian, "R")
ftbf_hist.Draw("hist")
ftbf_gaussian.Draw("same")
canvas2.SaveAs(f'{outdir_fits}fit_ftbf_{ftbf_dataset}.png')
canvas2.Close()
ftbf_rootfile.Close()

ftbf_tr = ftbf_gaussian.GetParameter(2)*1000
ftbf_tr_PS = np.sqrt(ftbf_tr**2 - 100)
ftbf_tr_err = ftbf_gaussian.GetParError(2)*1000
ftbf_tr_err_PS = ftbf_tr_err*(ftbf_tr/ftbf_tr_PS)


canvas = TCanvas("cv","cv",1000,800)
canvas.SetGrid(0,1)
gStyle.SetOptStat(0)
# Define hist for axes style
htemp = TH1F()
htemp.SetStats(0)
htemp.SetLineWidth(3)
htemp.Draw("AXIS")

graph = TGraphErrors(len(attenuations_val), array('d', attenuations_val), array('d', tr_vs_attn), array('d', [0]*len(attenuations_val)), array('d', trErr_vs_attn))
graph.SetLineColor(800+10)
graph.SetMarkerStyle(20)
graph.SetMarkerColor(800+10)
graph.SetMarkerSize(1)
graph.SetTitle("")
graph.GetXaxis().SetTitle("Attenuation [%]")
graph.GetYaxis().SetTitle("Time Resolution [ps]")
graph.GetYaxis().SetRangeUser(20.0, 37.0)
graph.GetXaxis().SetRangeUser(90,95)
graph.Draw("AP same")

band = TGraphAsymmErrors(2, array('d', [91.5, 94]), array('d', [ftbf_tr_PS, ftbf_tr_PS]), array('d', [0, 0]), array('d', [0, 0]), array('d', [ftbf_tr_err_PS]*2), array('d', [ftbf_tr_err_PS]*2))
band.SetFillColor(4)
band.SetFillStyle(3001)  # 3001 is a style for a semi-transparent fill
band.Draw("3 same")

line = TGraph(2, array('d', [91.5, 94]), array('d', [ftbf_tr_PS, ftbf_tr_PS]))
line.SetLineColor(4)
line.SetLineStyle(4)
line.SetLineWidth(2)
line.Draw("L same")

legend = TLegend(0.55, 0.2, 0.9, 0.3)
legend.AddEntry(line, "120 GeV protons", "l")
legend.AddEntry(graph, "Laser", "ep")
legend.SetBorderSize(1)
legend.SetLineColor(kBlack)
legend.SetTextFont(myStyle.GetFont())
legend.SetTextSize(myStyle.GetSize()-4)
legend.Draw("same")

htemp.Draw("AXIS same")

# myStyle.BeamInfo()
myStyle.SensorInfoSmart(ftbf_dataset,isPaperPlot=True)

canvas.SaveAs(f'{outdir}tr_vs_attn.pdf')
canvas.Clear()

