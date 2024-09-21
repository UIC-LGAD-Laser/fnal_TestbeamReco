from ROOT import TFile,TTree,TCanvas,TH1F,TH1D,TH2F,TLatex,TMath,TEfficiency,TGraphAsymmErrors,TLegend,gROOT,gPad,gStyle, kBlack, kWhite, TF1, TPaveStats
import langaus
from stripBox import getStripBox
import optparse
import myStyle
import myFunctions as mf
import os

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
parser.add_option('-x', dest='xlength', default = 150, help="X-axis max val")

options, args = parser.parse_args()

ftbf_loc = options.ftbf_dataset
laser_loc = options.laser_dataset
xlength = float(options.xlength)

ftbf_dataset = [i for i in ftbf_loc.split('/') if "HPK" in i]
sensor_Geometry = myStyle.GetGeometry(ftbf_dataset[0])

sensor = sensor_Geometry['sensor']
pitch = sensor_Geometry['pitch']

canvas = TCanvas("cv","cv",1000,800)
canvas.SetGrid(0,1)
gStyle.SetOptStat(0)

outdir = "../Paper_plots/"
if not os.path.exists(outdir):
    myStyle.CreateFolder("../", "Paper_plots/")

print("Setting up Langaus")
fit = langaus.LanGausFit()
print("Setup Langaus")

rootfiles = [TFile(f'{ftbf_loc}/Amplitude/MidGapAmp_distribution.root'), TFile(f'{laser_loc}/MidGapAmp_distribution.root')]
fit_type = ["landaugausfunction", "gaussian"]
labels = ["120 GeV protons", "Laser"]
colors=[4,800+8]

# Define hist for axes style
htemp = TH1F()
htemp.SetStats(0)
htemp.SetLineWidth(3)
htemp.Draw("AXIS")

legend = TLegend(0.55, 0.5, 0.9, 0.6)

count = 0
for file_iter,fit,color,label in zip(rootfiles,fit_type,colors,labels):
    hist = file_iter.Get("py")
    hist.SetStats(0)
    hist.SetLineColor(color)
    if count > 0: suffix = "same" 
    else: suffix = ""
    hist.Draw(f'hist {suffix}')
    fit_func = file_iter.Get(fit)
    fit_func.SetLineColor(color)
    fit_func.Draw("same")
    legend.AddEntry(hist, label, "l")
    count+=1

legend.SetBorderSize(1)
legend.SetLineColor(kBlack)
legend.SetTextFont(myStyle.GetFont())
legend.SetTextSize(myStyle.GetSize()-4)
legend.Draw("same")

htemp.Draw("AXIS same")
htemp.SetTitle("")
# htemp.GetXaxis().SetRangeUser(0.0, xlength)
htemp.GetXaxis().SetTitle("Amplitude [mV]")
htemp.GetYaxis().SetTitle("Frequency")

# myStyle.BeamInfo()
myStyle.SensorInfoSmart(ftbf_dataset[0],isPaperPlot=True)

canvas.SaveAs(f'{outdir}AmplitudeDistribution.pdf')

canvas.Clear()

for outputfile in rootfiles:
    outputfile.Close()
