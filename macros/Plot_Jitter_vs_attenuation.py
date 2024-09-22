from ROOT import TFile,TTree,TCanvas,TH1F,TH1D,TH2F,TLatex,TLine,TMath,TGraphAsymmErrors,TEfficiency,TGraph,TGraphErrors,TLegend,gROOT,gPad,gStyle, kBlack, kWhite, kOrange,TF1, TPaveStats
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
outdir_fits = "../Paper_plots/jitter_vs_attn_fits/"
if not os.path.exists(outdir_fits):
    myStyle.CreateFolder("../", "Paper_plots/jitter_vs_attn_fits/")

position_center = 0
laser_bias_voltage = "205"
attenuations = ["91P5", "92", "92P5", "93", "93P5", "94"]
attenuations_val = [91.5, 92, 92.5, 93, 93.5, 94]
jitter_vs_attn = []
jitterErr_vs_attn = []

print("Setting up Langaus")
fit = langaus.LanGausFit()
print("Setup Langaus")

canvas2 = TCanvas("canvas2", "Histogram Fit", 1200, 900)
for attenuation in attenuations:
    rootfile = TFile(f'{laser_loc}/HPK_W5_17_2_{laser_bias_voltage}V_{attenuation}attn/HPK_W5_17_2_{laser_bias_voltage}V_{attenuation}attn_Analyze.root')
    # hist = rootfile.Get("weighted_jitter_hist")
    info_entry = HistoInfo("weighted_jitter_vs_xy", rootfile, "jitter_vs_x", ylabel="Weighted Jitter [mV]", sensor=f'HPK_W5_17_2_{laser_bias_voltage}V_{attenuation}attn', center_position=position_center)
    midgapbin = info_entry.th1.GetXaxis().FindBin(-0.25)
    tmpHist = info_entry.th2.ProjectionY("py",midgapbin,midgapbin)
    myRMS = tmpHist.GetRMS()
    myMean = tmpHist.GetMean()
    myLanGausFunction = fit.fit(tmpHist, fitrange=(myMean-1.5*myRMS, myMean+3*myRMS))
    tmpHist.Draw("hist")
    myLanGausFunction.Draw("same")
    canvas2.SaveAs(f'{outdir_fits}fit_{laser_bias_voltage}V_{attenuation}attn.png')
    canvas2.Clear()
    jitter = myLanGausFunction.GetParameter(1)
    jitter_err = myLanGausFunction.GetParError(1)
    jitter_vs_attn.append(jitter)
    jitterErr_vs_attn.append(jitter_err)
    rootfile.Close()

ftbf_rootfile = TFile(f'{ftbf_loc}/{ftbf_dataset}_Analyze.root')
# ftbf_hist = ftbf_rootfile.Get("weighted2_jitter_Overall")

info_entry2 = HistoInfo("weighted2_jitter_vs_xy", ftbf_rootfile, "jitter_vs_x", ylabel="Weighted Jitter [mV]", sensor=ftbf_dataset, center_position=position_center)
midgapbin = info_entry2.th1.GetXaxis().FindBin(-0.25)
tmpHist2 = info_entry2.th2.ProjectionY("py",midgapbin,midgapbin)
myRMS = tmpHist2.GetRMS()
myMean = tmpHist2.GetMean()
ftbf_myLanGausFunction = fit.fit(tmpHist2, fitrange=(myMean-1.5*myRMS, myMean+3*myRMS))
tmpHist2.Draw("hist")
ftbf_myLanGausFunction.Draw("same")
canvas2.SaveAs(f'{outdir_fits}fit_ftbf_{ftbf_dataset}.png')
canvas2.Close()
ftbf_rootfile.Close()

ftbf_jitter = ftbf_myLanGausFunction.GetParameter(1)
# ftbf_jitter_PS = np.sqrt(ftbf_jitter**2 - 100)
ftbf_jitter_err = ftbf_myLanGausFunction.GetParError(1)
# ftbf_jitter_err_PS = ftbf_jitter_err*(ftbf_jitter/ftbf_jitter_PS)


canvas = TCanvas("cv","cv",1000,800)
canvas.SetGrid(0,1)
gStyle.SetOptStat(0)
# Define hist for axes style
htemp = TH1F()
htemp.SetStats(0)
htemp.SetLineWidth(3)
htemp.Draw("AXIS")

graph = TGraphErrors(len(attenuations_val), array('d', attenuations_val), array('d', jitter_vs_attn), array('d', [0]*len(attenuations_val)), array('d', jitterErr_vs_attn))
graph.SetLineColor(kOrange+10)
graph.SetMarkerStyle(20)
graph.SetMarkerColor(kOrange+10)
graph.SetMarkerSize(1)
graph.SetTitle("")
graph.GetXaxis().SetTitle("Attenuation [%]")
graph.GetYaxis().SetTitle("Weighted jitter [ps]")
graph.GetYaxis().SetRangeUser(15.0, 35.0)
graph.GetXaxis().SetRangeUser(90,95)
graph.Draw("AP same")

band = TGraphAsymmErrors(2, array('d', [91.5, 94]), array('d', [ftbf_jitter, ftbf_jitter]), array('d', [0, 0]), array('d', [0, 0]), array('d', [ftbf_jitter_err]*2), array('d', [ftbf_jitter_err]*2))
band.SetFillColor(4)
band.SetFillStyle(3001)  # 3001 is a style for a semi-transparent fill
band.Draw("3 same")

line = TGraph(2, array('d', [91.5, 94]), array('d', [ftbf_jitter, ftbf_jitter]))
line.SetLineColor(4)
line.SetLineStyle(4)
line.SetLineWidth(2)
line.Draw("L same")

ftbf_jitter = 1.33 * ftbf_jitter
scaledband = TGraphAsymmErrors(2, array('d', [91.5, 94]), array('d', [ftbf_jitter, ftbf_jitter]), array('d', [0, 0]), array('d', [0, 0]), array('d', [ftbf_jitter_err]*2), array('d', [ftbf_jitter_err]*2))
scaledband.SetFillColor(1)
scaledband.SetFillStyle(3001)  # 3001 is a style for a semi-transparent fill
scaledband.Draw("3 same")

scaledline = TGraph(2, array('d', [91.5, 94]), array('d', [ftbf_jitter, ftbf_jitter]))
scaledline.SetLineColor(1)
scaledline.SetLineStyle(4)
scaledline.SetLineWidth(2)
scaledline.Draw("L same")

legend = TLegend(0.15, 0.7, 0.65, 0.9)
legend.AddEntry(line, "120 GeV protons", "l")
legend.AddEntry(scaledline, "120 GeV protons (Scaled)", "l")
legend.AddEntry(graph, "Laser", "ep")
legend.SetBorderSize(1)
legend.SetLineColor(kBlack)
legend.SetTextFont(myStyle.GetFont())
legend.SetTextSize(myStyle.GetSize()-4)
legend.Draw("same")

htemp.Draw("AXIS same")

# myStyle.BeamInfo()
myStyle.SensorInfoSmart(ftbf_dataset,isPaperPlot=True)

canvas.SaveAs(f'{outdir}jitter_vs_attn.pdf')
canvas.Clear()

