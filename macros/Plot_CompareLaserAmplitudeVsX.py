from ROOT import TFile,TTree,TCanvas,TH1D,TH1F,TH2F,TLatex,TMath,TEfficiency,TGraphAsymmErrors,TLegend,TPaveText,gROOT,gStyle, gPad, kWhite, kBlack
import ROOT
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
        self.th2 = self.getTH2(f, inHistoName, sensor)
        self.th1 = self.getTH1(outHistoName)

    def getTH2(self, f, name, sensor, axis='zx'):
        th3 = f.Get(name)
        th2 = th3.Project3D(axis)

        # # Rebin low statistics sensors
        # if sensor=="BNL2020":
        #     th2.RebinX(5)
        # elif sensor=="BNL2021":
        #     th2.RebinX(10)
        # # if "1cm_500up_300uw" in sensor:
        # #     th2.RebinX(2)

        return th2
    
    def getName(self, hname):
        return self.outHistoName

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
parser.add_option('-x','--xlength', dest='xlength', default = 2.5, help="Limit x-axis in final plot")
parser.add_option('-y','--ylength', dest='ylength', default = 160, help="Max Amp value in final plot")
parser.add_option('-D', dest='Dataset', default = "", help="Dataset, which determines filepath")
parser.add_option('-d', dest='debugMode', action='store_true', default = False, help="Run debug mode")
parser.add_option('-t', dest='useTight', action='store_true', default = False, help="Use tight cut for pass")
options, args = parser.parse_args()

dataset = options.Dataset
outdir = myStyle.getOutputDir(dataset)
inputfile = TFile("%s%s_Analyze.root"%(outdir,dataset))

sensor_Geometry = myStyle.GetGeometry(dataset)

sensor = sensor_Geometry['sensor']
pitch = sensor_Geometry['pitch'] # For pads
# pitch = 0 # For strips

xlength = float(options.xlength)
ylength = float(options.ylength)
debugMode = options.debugMode

is_tight = options.useTight

# Get position of the central channel in the "x" direction
# position_center = mf.get_central_channel_position(inputfile, "x")
position_center = 0.25

outdirOld = myStyle.GetPlotsDir(outdir, "AmplitudeCompare/")

outdirSave = "../Paper_plots/"
if not os.path.exists(outdirSave):
    myStyle.CreateFolder("../", "Paper_plots/")

# Save list with histograms to draw
list_overall_htitles = [
    # [hist_input_name, short_output_name, y_axis_title]
    ["amplitude_vs_xy", "Amplitude", "MPV signal amplitude [mV]"],
    ["amplitudeNoSum_vs_xy", "AmplitudeNoSum", "MPV signal amplitude [mV]"]
]

# Amplitude per channel
list_channel_htitles = []
# indices = mf.get_existing_indices(inputfile, "amplitude_vs_xy_channel")
indices = ["03","02"]
for index in indices:
    channel_element = ["amplitude_vs_xy_channel%s"%index, "Amplitude_ch%s"%index, "MPV signal amplitude [mV]"]
    list_channel_htitles.append(channel_element)
    ncol = int(index[1])+1 if index[0] != "0" else 2

# Use tight cut histograms
if (is_tight):
    print(" >> Using tight cuts!")
    for titles in list_overall_htitles:
        titles[0]+= "_tight"

# List with histograms using HistoInfo class
histoInfo_overall, histoInfo_channel = [], []
for titles in (list_overall_htitles + list_channel_htitles):
    hname, outname, ytitle = titles
    if not (inputfile.Get(hname)):
        print(" >> Histogram %s does not exist! Skipping."%hname)
        continue
    info_obj = HistoInfo(hname, inputfile, outname, yMax=ylength, ylabel=ytitle,
                         sensor=dataset, center_position=position_center)
    if titles in list_overall_htitles:
        histoInfo_overall.append(info_obj)
    elif titles in list_channel_htitles:
        histoInfo_channel.append(info_obj)

canvas = TCanvas("cv","cv",1000,800)
canvas.SetGrid(0,1)
gStyle.SetOptStat(0)

if debugMode:
    outdir_q = myStyle.CreateFolder(outdir, "Amp_vs_X_fits0/")

# Get total number of bins in x-axis to loop over (all hists have the same number, in principle)
all_histoInfos = histoInfo_overall + histoInfo_channel
nbins = all_histoInfos[0].th2.GetXaxis().GetNbins()

midgap_bins = [all_histoInfos[0].th1.GetXaxis().FindBin(-0.25)]

plot_xlimit = abs(inputfile.Get("stripBoxInfo00").GetMean(1) - position_center)
if ("pad" not in dataset) and ("500x500" not in dataset):
    plot_xlimit-= pitch/(2 * 1000.)

print("Setting up Langaus")
fit = langaus.LanGausFit()
print("Setup Langaus")

# Loop over X bins
for i in range(all_histoInfos[0].th1.GetXaxis().FindBin(-0.25), all_histoInfos[0].th1.GetXaxis().FindBin(0.30)):
    for info_entry in all_histoInfos:
        totalEvents = info_entry.th2.GetEntries()
        tmpHist = info_entry.th2.ProjectionY("py",i,i)
        myRMS = tmpHist.GetRMS()
        myMean = tmpHist.GetMean()
        nEvents = tmpHist.GetEntries()
        value = myMean

        # Define minimum of bin's entries to be fitted
        minEvtsCut = 0.5*totalEvents/nbins
        if("20T" in dataset):
            minEvtsCut = 0.3*totalEvents/nbins
        if(is_tight):
            minEvtsCut = 0
        if (i == 1):
            msg_nentries = "%s: nEvents > %.2f "%(info_entry.inHistoName, minEvtsCut)
            msg_nentries+= "(Total events: %i)"%(totalEvents)
            print(msg_nentries)

        #Do fit
        if(nEvents > minEvtsCut):
            tmpHist.Rebin(2)
            # if (myMean > 50):
            #     tmpHist.Rebin(5)
            # else:
            #     tmpHist.Rebin(10)

            myLanGausFunction = fit.fit(tmpHist, fitrange=(myMean-1.5*myRMS, myMean+3*myRMS))
            myMPV = myLanGausFunction.GetParameter(1)
            value = myMPV

            # For Debugging
            if (debugMode):
                tmpHist.Draw("hist")
                myLanGausFunction.Draw("same")
                canvas.SaveAs("%sq_%s%i.gif"%(outdir_q, info_entry.outHistoName, i))
                bin_center = info_entry.th1.GetXaxis().GetBinCenter(i)
                msg_amp = "Bin: %i (x center = %.3f)"%(i, bin_center)
                # msg_amp+= " -> Amplitude: %.3f mV"%(value)
                msg_amp+= " -> Entries: %.3f"%(tmpHist.GetEntries())
                print(msg_amp)
        else:
            value = 0.0

        value = value if (value>0.0) else 0.0

        # Fill only when inside limits
        if not mf.is_inside_limits(i, info_entry.th1, xmax=plot_xlimit):
            continue

        info_entry.th1.SetBinContent(i, value)

# Define output file
output_path = "%sCompare_AmplitudeVsX"%(outdirSave)
if (is_tight):
    output_path+= "_tight"
output_path+= ".root"

outputfile = TFile(output_path,"RECREATE")

# Define hist for axes style
# htemp = TH1F("htemp", "", 1, -xlength, xlength)
htemp = TH1F("htemp", "", 1, -0.49, 0.49)
htemp.SetStats(0)
htemp.GetXaxis().SetTitle("Track x position [mm]")
htemp.GetYaxis().SetRangeUser(0.0, ylength)
htemp.GetYaxis().SetTitle("MPV signal amplitude [mV]")
htemp.SetLineColor(colors[2])

# Draw overall amplitude vs X
for i,info_entry in enumerate(histoInfo_overall):
    hist = info_entry.th1
    hist.SetLineColor(colors[0])
    hist.SetLineWidth(2)

    ymin = hist.GetMinimum()
    ymax = hist.GetMaximum()

    haxis = htemp.Clone()
    haxis.SetMinimum(ymin)
    haxis.SetMaximum(ymax)

    # Define and draw gray bars in the background (Position of metallic sections)
    boxes = getStripBox(inputfile, ymin=ymin, ymax=ymax, strips=True,
                        shift=position_center, pitch=pitch/1000.)

    hist.Write()


# Draw all amplitude per channel vs X
# Define legend
pcenter = myStyle.GetPadCenter()
pmargin = myStyle.GetMargin()
legend = TLegend(pcenter-0.30, 1-pmargin-0.15, pcenter+0.30, 1-pmargin-0.05)
legend.SetLineColor(kBlack)
# legend.SetFillColor(kWhite)
legend.SetTextFont(myStyle.GetFont())
legend.SetTextSize(myStyle.GetSize()-4)
legend.SetNColumns(ncol)

legendBot = TLegend(pcenter-0.30, 1-pmargin-0.25, pcenter+0.30, 1-pmargin-0.15)
legendBot.SetNColumns(2)
legendBot.SetLineColor(kBlack)
legendBot.SetTextFont(myStyle.GetFont())
legendBot.SetTextSize(myStyle.GetSize()-4)

htemp.Draw("AXIS")
for box in boxes:
    box.Draw()
gPad.RedrawAxis("g")
map = {'02': [2],
        '03': [1]}
for i,info_entry in enumerate(histoInfo_channel):
    hist = info_entry.th1
    hist.SetLineColor(colors[i])
    hist.SetLineWidth(2)
    hist.Draw("hist same")

    idx = indices[i]
    ltitle = "Pad %s"%(idx) if "10" in indices else "Strip %i"%(map[str(idx)][0])
    legend.AddEntry(hist, ltitle, "lp")

    hist.Write()


inputfileLaser = TFile("%s%sPlotAmplitudeVsX.root"%("/uscms/home/snanda/nobackup/LaserDanushUpdated/TestbeamReco/output/","LeCroy_W2_3_2_198V_99P9attn/"))
histRightCh = inputfileLaser.Get("amplitude_vs_x_channel02")
histRightCh.SetLineColor(colors[0])
histRightCh.SetLineWidth(2)
histRightCh.SetLineStyle(7)
histLeftCh = inputfileLaser.Get("amplitude_vs_x_channel01")
histLeftCh.SetLineColor(colors[1])
histLeftCh.SetLineWidth(2)
histLeftCh.SetLineStyle(7)
histRightCh.Draw("SAME")
histLeftCh.Draw("SAME")

legend.Draw()
ftbf_tmphist = histoInfo_overall[-1].th1.Clone()
ftbf_tmphist.SetLineColor(kBlack)
laser_tmphist = histRightCh.Clone()
laser_tmphist.SetLineColor(kBlack)
laser_tmphist.SetLineStyle(7)
legendBot.AddEntry(ftbf_tmphist, "120 GeV protons")
legendBot.AddEntry(laser_tmphist, "Laser")
legendBot.Draw()
legendBox = TPaveText(pcenter-0.30, 1-pmargin-0.25, pcenter+0.30, 1-pmargin-0.05, "NDC")
legendBox.SetBorderSize(1)
legendBox.SetLineColor(kBlack)
legendBox.SetFillColor(0)
legendBox.SetFillColorAlpha(0, 0.0)
legendBox.Draw("same")

# myStyle.BeamInfo()
myStyle.SensorInfoSmart(dataset,isPaperPlot=True)

htemp.Draw("AXIS same")
save_path = "%sAmplitudeAllChannels_vs_x"%(outdirSave)
canvas.SaveAs("%s.pdf"%save_path)

canvas.Clear()
outputfile.Close()
