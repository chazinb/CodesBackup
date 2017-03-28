import sys
import subprocess
import ROOT 



file_r = ROOT.TFile("StatisticVariation_r_Signal.root", "RECREATE")
h_r=ROOT.TH1F('h_r', 'Statistic Variation of r', 200, 0, 2)
theR = []
#combinecommand = "combine -options -blah blah blah"
combinecommand = "combine -M Asymptotic"
for i in range(0, 200):
   command = combinecommand + " T2tt_350_225_MT2ll_metPfType1_" + str(i) + "_makeVariationsOfSignal_Weigthed.txt" 
   #command = combinecommand + " Iteration_" + str(i) + ".root"
   output_aux = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
   output = output_aux.communicate()[0]
   ListOfFiles = output.split("\n")
   for line in output.split("\n"):
      #if line.find("Expected 50.0%: r < ") >= 0:
      if not line.find("Expected 50.0%: r < "):
         #print line
         chops = line.split()
         theR.append(float(chops[4]))
         h_r.Fill(float(chops[4]))
         #print (chops[4])
  # print (theR[i])
h_r.Write()
#canvas = ROOT.TCanvas("", "", 800, 600)
#cd()
#canvas.Draw()
#h_r.Draw()
#canvas.SaveAs("StatisticVariation_r_Signal.png")

file_r.Close()

 
