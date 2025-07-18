from ROOT import TFile,TTree,TCanvas,TH1D,TH1F,TH2F,TLatex,TMath,TEfficiency,TGraphAsymmErrors,TLegend,TPaveText,gROOT,gStyle, gPad, kWhite, kBlack
import ROOT
import math
import os
from stripBox import getStripBox
import optparse
import myStyle
import langaus
import myFunctions as mf
import numpy as np

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
if "W9" in dataset:
    ylength=150

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
jitter_new_ftbf = jitter.Clone("jitterNew_vs_x")
tr = inputfile2.Get("Time_DiffW2Tracker")
landau_50umSi = 22
if "W9" in dataset:
    landau_50umSi = 35

for number in range(1, jitter.GetXaxis().GetNbins()+1):
    x_position = jitter.GetXaxis().GetBinCenter(number)
    bin_number_tr = tr.FindBin(x_position)
    t = tr.GetBinContent(bin_number_tr)
    t_error = tr.GetBinError(bin_number_tr)
    # t = tr.GetBinContent(number)
    if (t>landau_50umSi):
        l = math.sqrt(t*t - landau_50umSi*landau_50umSi)
    else:
        l = 0
    jitter_new_ftbf.SetBinContent(number, l)

allHistos = [tr, jitter_new_ftbf]
names = ["tr_vs_x","jitter_vs_x"]
colors_new = [colors[2], colors[0], colors[4]]
labels = ["120 GeV total time resolution (t)", "120 GeV jitter (#sqrt{t^{2} - 22^{2}})"]

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
    legend.AddEntry(hist, label, "lp")

inputfileLaser = TFile("%s%sAllTRs.root"%("/uscms/home/dshekar/nobackup/laser_analysis/TestbeamReco/output/",laser_dataset+'/'))
labels=["Laser jitter = total time resolution"]
tr_laser = inputfileLaser.Get("weighted2_time_diffTracker")
jitter_laser = inputfileLaser.Get("weighted2_time_diffTracker")
laser_allHistos = [tr_laser]
laser_all_histos_shifted = []
for i,j in zip(laser_allHistos,names):
    laser_all_histos_shifted.append(shift_histogram(i,j))
for i, (hist, label, tmpcolor) in enumerate(zip(laser_all_histos_shifted, labels, colors_new)):
    hist.GetXaxis().SetRangeUser(-0.2, 0.2)
    hist.SetLineColor(tmpcolor)
    hist.SetLineWidth(2)
    # hist.SetStats(0)
    hist.Draw("hist E same")
    tr_hist_laser = hist.Clone("ratio_hist")
    legend.AddEntry(hist, label, "lp")

legend.Draw()
ftbf_tmphist = allHistos[-1].Clone()
ftbf_tmphist.SetLineColor(kBlack)
ftbf_tmphist.SetLineWidth(2)
ftbf_tmphist.SetLineStyle(7)
laser_tmphist = laser_allHistos[-1].Clone()
laser_tmphist.SetLineColor(kBlack)
laser_tmphist.SetLineWidth(2)
legendBox = TPaveText(pcenter-0.30, 1-pmargin-0.23, pcenter+0.30, 1-pmargin-0.05)
legendBox.SetBorderSize(1)
legendBox.SetLineColor(kBlack)
legendBox.SetFillColor(0)
legendBox.SetFillColorAlpha(0, 0.0)
legendBox.Draw("same")
# myStyle.BeamInfo()
myStyle.SensorInfoSmart(dataset,isPaperPlot=True)
htemp.Draw("AXIS same")
canvas.Update()
canvas.SaveAs(outdirSave+"PaperCommentsAll_TRs_"+str(landau_50umSi)+"landau_"+dataset+".pdf")


# CALCULATE RATIO PLOTS OF NOISE AND THEN JITTERS
noiseInputfileLaser = TFile(f'/uscms/home/dshekar/nobackup/laser_analysis/TestbeamReco/output/{laser_dataset}/NoiseStudy/PlotNoiseOverallVsX.root')
Noise_laser = noiseInputfileLaser.Get("baselineRMS_vs_x")
laser_noise_shifted = shift_histogram(Noise_laser,"laser_noise_vs_x")

ftbf_noise_file = TFile("%s/Noise/NoiseVsX.root"%(outdir))
noise = ftbf_noise_file.Get("Noise")
ftbf_noise_shifted = shift_histogram(noise,"ftbf_noise_vs_x")


jitter_new_ftbf = all_histos_shifted[1]
x_min = -0.25
x_max = 0.25
step_size = 0.05
x_values = np.round(np.arange(x_min, x_max + step_size, step_size), 2)  # Generate x values with step size
print("x_bin centers = ", x_values)
bin_edges = np.append(x_values - step_size / 2, x_values[-1] + step_size / 2)
ratio_jitter = ROOT.TH1F("ratio_jitter", "Ratio of jitter (laser/ftbf)", len(bin_edges) - 1, bin_edges)
ratio_noise = ROOT.TH1F("ratio_noise", "Ratio of noise (laser/ftbf)", len(bin_edges) - 1, bin_edges)


# Loop through the bins of the first histogram
for bin_center in x_values:
    bin_index_laser = tr_hist_laser.GetXaxis().FindBin(bin_center)
    jitter_content_laser = tr_hist_laser.GetBinContent(bin_index_laser)
    noise_content_laser = laser_noise_shifted.GetBinContent(bin_index_laser)

    bin_index_ftbf = jitter_new_ftbf.GetXaxis().FindBin(bin_center)
    jitter_content_ftbf = jitter_new_ftbf.GetBinContent(bin_index_ftbf)
    noise_content_ftbf = ftbf_noise_shifted.GetBinContent(bin_index_ftbf)

    # Avoid division by zero
    if jitter_content_ftbf != 0:
        ratio = jitter_content_laser / jitter_content_ftbf
        ratio_jitter.Fill(bin_center, ratio)  # Fill with bin center and ratio
    else:
        print(f"Bin at x = {bin_center}: Division by zero (jitter_new_ftbf content is 0)")
    if noise_content_ftbf != 0:
        ratio2 = noise_content_laser / noise_content_ftbf
        ratio_noise.Fill(bin_center, ratio2)  # Fill with bin center and ratio
    else:
        print(f"Bin at x = {bin_center}: Division by zero (noise_ftbf content is 0)")
    
# Style the ratio histograms
ratio_jitter.SetLineColor(1)  # Black for jitter ratio
ratio_jitter.SetLineWidth(2)  # Set line width
ratio_jitter.SetMarkerStyle(20)  # Add markers for jitter ratio
ratio_jitter.SetStats(0)  # Disable statistics box
ratio_jitter.GetXaxis().SetTitle("Track X position [mm]")  # Set x-axis title
ratio_jitter.GetYaxis().SetTitle("Ratio")  # Set y-axis title

ratio_noise.SetLineColor(2)  # Red for noise ratio
ratio_noise.SetLineWidth(2)  # Set line width
ratio_noise.SetMarkerStyle(21)  # Add markers for noise ratio
ratio_noise.SetStats(0)  # Disable statistics box

# Create a new canvas for the ratio histograms
ratio_canvas_manual = ROOT.TCanvas("ratio_canvas_manual", "Ratio Histogram (Manual)", 800, 600)

# Adjust canvas margins to prevent axis titles and legend from being cut off
ratio_canvas_manual.SetLeftMargin(0.15)  # Increase left margin for y-axis title
ratio_canvas_manual.SetBottomMargin(0.15)  # Increase bottom margin for x-axis title
ratio_canvas_manual.SetRightMargin(0.1)  # Add some space on the right for the legend

# Set the x-axis range for the histograms
ymax = 1.5
if "W9" in dataset:
    ymax = 2.5
ratio_jitter.GetXaxis().SetRangeUser(-0.3, 0.3)  # Example: Set x-axis range
ratio_jitter.GetYaxis().SetRangeUser(0, ymax)  # Example: Set x-axis range
ratio_noise.GetXaxis().SetRangeUser(-0.3, 0.3)
ratio_noise.GetYaxis().SetRangeUser(0, ymax)

# Draw the histograms
ratio_jitter.Draw("hist")  # Draw jitter ratio first
ratio_noise.Draw("hist SAME")  # Draw noise ratio on the same canvas
# Add a legend
legend2 = ROOT.TLegend(0.3, 0.2, 0.7, 0.4)  # Define legend2 position (bottom-right corner)
legend2.SetLineColor(ROOT.kBlack)  # Set the border color of the legend2 box
legend2.SetTextFont(myStyle.GetFont())  # Set the font for the legend2 text
legend2.SetTextColor(ROOT.kBlack)  # Ensure the text color is visible (black)
legend2.SetTextFont(62)  # or 63 for pixel-based
legend2.SetTextSize(0.03)
# legend2.SetHeader("Ratios (laser/beam)", "C")  # Add a title to the legend2
legend2.AddEntry(ratio_jitter, "Jitter", "l")  # Add jitter ratio to legend2
legend2.AddEntry(ratio_noise, "Noise", "l")  # Add noise ratio to legend2
legend2.SetBorderSize(0)  # Remove border around legend2
legend2.SetTextSize(0.03)  # Set text size
legend2.Draw()  # Draw the legend2

# Force the canvas to update and display everything
ratio_canvas_manual.Update()

# Save the ratio histogram to a file
ratio_canvas_manual.SaveAs(outdirSave + f'PaperComments_manual_ratio_histogram_{landau_50umSi}landau_{dataset}.png')

# Close the file
ftbf_noise_file.Close()
inputfile.Close()
inputfile_Analyze.Close()
noiseInputfileLaser.Close()
inputfile2.Close()
inputfileLaser.Close()
ftbf_noise_file.Close()