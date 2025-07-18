from ROOT import TFile,TTree,TCanvas,TH1D,TH1F,TH2F,TLatex,TMath,TEfficiency,TGraphAsymmErrors,TLegend,TPaveText,gROOT,gStyle, gPad, kWhite, kBlack
import ROOT
import math
import os
from stripBox import getStripBox
import optparse
import myStyle
import langaus
import myFunctions as mf

gROOT.SetBatch( True )
gStyle.SetOptFit(1011)
colors = myStyle.GetColors(True)

## Defining Style
myStyle.ForceStyle()

organized_mode=True

def shift_histogram(hist, name):
    xmin = -0.5
    xmax = 0.05
    hist_cloned = hist.Clone(name)
    for i in range(1, hist_cloned.GetXaxis().GetNbins()+1):
        hist_cloned.SetBinContent(i, 0)
    offset = hist_cloned.GetXaxis().FindBin(-1*position_center) - hist_cloned.GetXaxis().FindBin(xmin)
    for i in range(hist.GetXaxis().FindBin(xmin), hist.GetXaxis().FindBin(xmax)):
        val = hist.GetBinContent(i)
        hist_cloned.SetBinContent(i+offset, val)
    return hist_cloned


# Construct the argument parser
parser = optparse.OptionParser("usage: %prog [options]\n")
parser.add_option('-D', dest='Dataset', default = "", help="Dataset, which determines filepath")
parser.add_option('-x', dest='xlength', default = 0.49, help="X-Length")
parser.add_option('-y', dest='ylength', default = 100, help="Y-Length")
options, args = parser.parse_args()
position_center = 0.25

dataset = options.Dataset
map = {'HPK_W2_3_2_50T_1P0_500P_50M_E240_180V': 'LeCroy_W2_3_2_198V_99P9attn', 'HPK_W4_17_2_50T_1P0_500P_50M_C240_204V': 'LeCroy_W4_17_2_222V_99P9attn', 'HPK_W5_17_2_50T_1P0_500P_50M_E600_190V':'HPK_W5_17_2_205V_94attn', 'HPK_W9_15_2_20T_1P0_500P_50M_E600_114V':'LeCroy_W9_15_2_121V_92P3attn'}
laser_dataset = map[dataset]
xlength = float(options.xlength)
ylength = float(options.ylength)


outdir=""
if organized_mode: 
    outdir = myStyle.getOutputDir(dataset)
    inputfile = TFile("%s/Jitter/JitterVsX.root"%(outdir))
    inputfile2 = TFile("%s/Resolution_Time/TimeDiffVsX.root"%(outdir))
else: 
    inputfile = TFile("../test/myoutputfile.root") 
    inputfile2 = TFile("../test/myoutputfile.root") 

outdirSave = "../Paper_plots/"
if not os.path.exists(outdirSave):
    myStyle.CreateFolder("../", "Paper_plots/")

inputfile_Analyze = TFile("%s%s_Analyze.root"%(outdir,dataset))
colors = myStyle.GetColors(True)

sensor_Geometry = myStyle.GetGeometry(dataset)

sensor = sensor_Geometry['sensor']
pitch  = sensor_Geometry['pitch']

canvas = TCanvas("canvas", "All TRs vs X", 1000, 800)

jitter = inputfile.Get("jitter_vs_x")
tr = inputfile2.Get("Time_DiffW2Tracker")
landau = jitter.Clone("landau_vs_x")

for number in range(1, jitter.GetXaxis().GetNbins()+1):
    j = jitter.GetBinContent(number)
    j_error = jitter.GetBinError(number)
    x_position = jitter.GetXaxis().GetBinCenter(number)
    bin_number_tr = tr.FindBin(x_position)
    t = tr.GetBinContent(bin_number_tr)
    t_error = tr.GetBinError(bin_number_tr)
    # t = tr.GetBinContent(number)
    if (t>j):
        l = math.sqrt(t*t - j*j)
        l_error = math.sqrt((t**2 * t_error**2 / l**2) + (j**2 * j_error**2 / l**2))
    else:
        l = 0
        l_error = 0
    landau.SetBinContent(number, l)
    landau.SetBinError(number,l_error)

allHistos = [tr, jitter]
names = ["tr_vs_x","jitter_vs_x"]
colors_new = [colors[2], colors[0], colors[4]]
labels = ["Total time resolution (t)", "Weighted jitter (j)"]

# Define hist for axes style
htemp = TH1F("htemp", "", 1, -xlength, xlength)
htemp.SetStats(0)
htemp.GetXaxis().SetTitle("Track x position [mm]")
htemp.GetYaxis().SetRangeUser(0.0, ylength)
htemp.GetYaxis().SetTitle("Time resolution [ps]")
htemp.SetLineColor(colors[2])

for i,info_entry in enumerate(allHistos):
    hist = info_entry
    ymin = hist.GetMinimum()
    ymax = hist.GetMaximum()
    haxis = htemp.Clone()
    haxis.SetMinimum(ymin)
    haxis.SetMaximum(ymax)
    # Define and draw gray bars in the background (Position of metallic sections)
    boxes = getStripBox(inputfile_Analyze, ymin=ymin, ymax=ylength, strips=True,
                        shift=position_center, pitch=pitch/1000.)

htemp.Draw("AXIS")
for box in boxes:
    box.Draw()
gPad.RedrawAxis("g")
# Draw all amplitude per channel vs X
# Define legend
pcenter = myStyle.GetPadCenter()
pmargin = myStyle.GetMargin()
legend = TLegend(pcenter-0.30, 1-pmargin-0.23, pcenter+0.30, 1-pmargin-0.05)
legend.SetLineColor(kBlack)
# legend.SetFillColor(kWhite)
legend.SetTextFont(myStyle.GetFont())
legend.SetTextSize(myStyle.GetSize()-4)
legend.SetNColumns(1)

legendBot = TLegend(pcenter-0.30, 1-pmargin-0.3, pcenter+0.30, 1-pmargin-0.23)
legendBot.SetNColumns(2)
legendBot.SetLineColor(kBlack)
legendBot.SetTextFont(myStyle.GetFont())
legendBot.SetTextSize(myStyle.GetSize()-4)

all_histos_shifted = []
for i,j in zip(allHistos,names):
    all_histos_shifted.append(shift_histogram(i,j))
for i, (hist, label, tmpcolor) in enumerate(zip(all_histos_shifted, labels, colors_new)):
    hist.GetXaxis().SetRangeUser(-0.25, 0.25)
    hist.SetLineColor(tmpcolor)
    hist.SetLineWidth(2)
    hist.SetLineStyle(7)
    hist.SetStats(0)
    hist.Draw("hist E same")
    # legend.AddEntry(hist, label, "lp")

inputfileLaser = TFile("%s%sAllTRs.root"%("/uscms/home/dshekar/nobackup/laser_analysis/TestbeamReco/output/",laser_dataset+'/'))
tr_laser = inputfileLaser.Get("weighted2_time_diffTracker")
scaled_jitter_laser = inputfileLaser.Get("scaled_jitter_vs_x")
landau_laser = inputfileLaser.Get("landau_vs_x")
laser_allHistos = [tr_laser, scaled_jitter_laser]
laser_all_histos_shifted = []
for i,j in zip(laser_allHistos,names):
    laser_all_histos_shifted.append(shift_histogram(i,j))
for i, (hist, label, tmpcolor) in enumerate(zip(laser_all_histos_shifted, labels, colors_new)):
    hist.GetXaxis().SetRangeUser(-0.2, 0.2)
    hist.SetLineColor(tmpcolor)
    hist.SetLineWidth(2)
    # hist.SetStats(0)
    hist.Draw("hist E same")
    legend.AddEntry(hist, label, "lp")

legend.Draw()
ftbf_tmphist = allHistos[-1].Clone()
ftbf_tmphist.SetLineColor(kBlack)
ftbf_tmphist.SetLineWidth(2)
ftbf_tmphist.SetLineStyle(7)
laser_tmphist = laser_allHistos[-1].Clone()
laser_tmphist.SetLineColor(kBlack)
laser_tmphist.SetLineWidth(2)
legendBot.AddEntry(ftbf_tmphist, "120 GeV protons")
legendBot.AddEntry(laser_tmphist, "Laser")
legendBot.Draw()
legendBox = TPaveText(pcenter-0.30, 1-pmargin-0.3, pcenter+0.30, 1-pmargin-0.05, "NDC")
legendBox.SetBorderSize(1)
legendBox.SetLineColor(kBlack)
legendBox.SetFillColor(0)
legendBox.SetFillColorAlpha(0, 0.0)
legendBox.Draw("same")

# myStyle.BeamInfo()
myStyle.SensorInfoSmart(dataset,isPaperPlot=True)

htemp.Draw("AXIS same")
canvas.Update()
canvas.SaveAs(outdirSave+"All_TRs_noLandau_"+dataset+".pdf")