from ROOT import TFile,TTree,TCanvas,TH1D,TH1F,TH2D,TH2F,TLatex,TMath,TColor,TLegend,TEfficiency,TGraphAsymmErrors,gROOT,gPad,TF1,gStyle,kBlack,kWhite,TH1
import ROOT
import os
from stripBox import getStripBox
import optparse
import myStyle

gROOT.SetBatch( True )
gStyle.SetOptFit(1011)

## Defining Style
myStyle.ForceStyle()

class HistoInfo:
    def __init__(self, plane, f, outHistoName, doFits=True, yMax=30.0, title="", xlabel="", ylabel="Position resolution [#mum]", sensor=""):
        self.plane = plane
        self.f = f
        self.outHistoName = outHistoName
        self.doFits = doFits
        self.yMax = yMax
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.th2 = self.getTH2(f, plane, sensor)
        self.th1 = self.getTH1(self.th2, outHistoName, self.shift(), self.fine_tuning(plane, sensor), plane)
        # self.sensor = sensor

    def getTH2(self, f, plane, sensor):
        th2 = f.Get("deltaX_vs_Xtrack_vs_Ytrack").Project3D(plane)
        # th2_temp = TH2D(outHist,"",42,-0.210,0.210,th2.GetYaxis().GetNbins(),th2.GetYaxis().GetXmin(),th2.GetYaxis().GetXmax())
        # for i in range(th2.GetXaxis().FindBin(-0.210+centerShift),th2.GetXaxis().FindBin(0.210+centerShift)+1,1):
        #     th2_temp.Fill()
        if sensor=="BNL2020" and plane=="zx": th2.RebinX(7)
        elif sensor=="BNL2020" and plane=="zy": th2.RebinX(8)
        elif sensor=="BNL2021" and plane=="zx": th2.RebinX(9)
        elif sensor=="BNL2021" and plane=="zy": th2.RebinX(8)
        return th2

    def getTH1(self, th2, name, centerShift, fine_value, plane):
        if plane=="zx" : th1_temp = TH1D(name,"",th2.GetXaxis().GetNbins(),th2.GetXaxis().GetXmin()-centerShift-fine_value,th2.GetXaxis().GetXmax()-centerShift-fine_value)
        elif plane=="zy": th1_temp = TH1D(name,"",th2.GetXaxis().GetNbins(),th2.GetXaxis().GetXmin()-fine_value,th2.GetXaxis().GetXmax()-fine_value)
        return th1_temp

    def shift(self):
        return (self.f.Get("stripBoxInfo02").GetMean(1)+self.f.Get("stripBoxInfo03").GetMean(1))/2.

    def fine_tuning(self, plane, sensor):
        value = 0.0
        if sensor=="BNL2020" and plane=="zx": value = 0.0075
        return value

# Construct the argument parser
parser = optparse.OptionParser("usage: %prog [options]\n")
parser.add_option('-f','--file', dest='file', default = "myoutputfile.root", help="File name (or path from ../test/)")
parser.add_option('-s','--sensor', dest='sensor', default = "BNL2020", help="Type of sensor (BNL, HPK, ...)")
parser.add_option('-b','--biasvolt', dest='biasvolt', default = 220, help="Bias Voltage value in [V]")
options, args = parser.parse_args()

file = options.file
sensor = options.sensor
bias = options.biasvolt

inputfile = TFile("../test/"+file,"READ")

all_histoInfos = [
    HistoInfo("zx",   inputfile, "x_track", True,  35.0, "", "Track x position [mm]","Mean #Deltax [#mum]",sensor),
    HistoInfo("zy",    inputfile, "y_track",  True,  35.0, "", "Track y position [mm]","Mean #Deltax [#mum]",sensor),
]

canvas = TCanvas("cv","cv",1000,800)
canvas.SetGrid(0,1)
TH1.SetDefaultSumw2()
gStyle.SetOptStat(0)

print("Finished setting up langaus fit class")

#loop over X bins
for i in range(0, all_histoInfos[0].th2.GetXaxis().GetNbins()+1):
    ##For Debugging
    #if not (i==46 and j==5):
    #    continue

    for info in all_histoInfos:
        tmpHist = info.th2.ProjectionY("py",i,i)
        myMean = tmpHist.GetMean()
        myRMS = tmpHist.GetRMS()
        myRMSError = tmpHist.GetRMSError()
        nEvents = tmpHist.GetEntries()
        fitlow = myMean - 1.5*myRMS
        fithigh = myMean + 1.5*myRMS
        value = myRMS
        error = myRMSError

        #Do fit 
        if(nEvents > 50):
            if(info.doFits):
                tmpHist.Rebin(2)
                
                fit = TF1('fit','gaus',fitlow,fithigh)
                tmpHist.Fit(fit,"Q", "", fitlow, fithigh)
                myMPV = fit.GetParameter(1)
                myMPVError = fit.GetParError(1)
                mySigma = fit.GetParameter(2)
                mySigmaError = fit.GetParError(2)
                value = 1000.0* myMPV
                error = 1000.0* myMPVError
                # value = 1000.0*mySigma
                # error = 1000.0*mySigmaError
            
                ##For Debugging
                # tmpHist.Draw("hist")
                # fit.Draw("same")
                # canvas.SetLogy()
                # canvas.SaveAs("q_"+str(i)+".gif")
                
                #print ("Bin : " + str(i) + " -> " + str(value) + " +/- " + str(error))
            else:
                value *= 1000.0
                error *= 1000.0
        else:
            value = 0.0
            error = 0.0

        # Removing telescope contribution
        # if value>6.0:
        #     error = error*value/TMath.Sqrt(value*value - 6*6)
        #     value = TMath.Sqrt(value*value - 6*6)
        # else:
        #     value = 0.0 # 20.0 to check if there are strange resolution values
        #     error = 0.0
        # if i<=info.th1.FindBin(-0.2) and sensor=="BNL2020":
        #     value = 0.0
        #     error = 0.0

        info.th1.SetBinContent(i,value)
        info.th1.SetBinError(i,error)
                        
# Plot 2D histograms
outputfile = TFile("PlotXRecoDiffVsX.root","RECREATE")
for info in all_histoInfos:
    if info.plane=="zx":
        if sensor=="BNL2020": xlimit_high = 0.32
        else: xlimit_high = 0.43
        xlimit_low = -xlimit_high
    elif info.plane=="zy":
        if sensor=="BNL2020":
            xlimit_high = 11.8
            xlimit_low = 9.6
        else:
            xlimit_high = 12.0
            xlimit_low = 9.0
    htemp = TH1F("htemp","",1,xlimit_low,xlimit_high)
    htemp.SetStats(0)
    htemp.SetMinimum(-15.0)
    htemp.SetMaximum(10.0)
    htemp.GetXaxis().SetTitle(info.xlabel)
    htemp.GetYaxis().SetTitle(info.ylabel)
    info.th1.SetLineColor(kBlack)
    htemp.Draw("AXIS")

    ymin = info.th1.GetMinimum()
    ymax = info.yMax

    boxes = getStripBox(inputfile,ymin,ymax,False,18,True,info.shift())
    for box in boxes:
        box.Draw()

    gPad.RedrawAxis("g")

    htemp.Draw("AXIS same")
    # info.th1.Draw("AXIS same")
    info.th1.Draw("hist e same")

    myStyle.BeamInfo()
    myStyle.SensorInfo(sensor, bias)

    canvas.SaveAs("PositionRes_vs_"+info.outHistoName+"_Mean.gif")
    canvas.SaveAs("PositionRes_vs_"+info.outHistoName+"_Mean.pdf")
    info.th1.Write()
    htemp.Delete()

outputfile.Close()

