
# HPK 500um Pitch sensors
# -----------------------

# HPK_500P=("HPK_W9_15_2_20T_1P0_500P_50M_E600_114V" "HPK_W2_3_2_50T_1P0_500P_50M_E240_180V" "HPK_W4_17_2_50T_1P0_500P_50M_C240_204V")
# source setup.sh
# ./configure
# make clean
# make -j4

# # HPK_500P=("HPK_W5_17_2_50T_1P0_500P_50M_E600_190V")


# for sensor in "${HPK_500P[@]}"; do
#     printf "\nRunning over ${sensor} sensor\n"
#     # cd ../test
#     # ./MyAnalysis -A InitialAnalyzer -D ${sensor}
#     # cd ../macros
#     # python FindStripCenters.py -D ${sensor}
#     # python FindDelayCorrections.py -D ${sensor}
#     # python FindInputHistos4YReco.py -D ${sensor} -I
#     # cd ../test
#     # ./MyAnalysis -A Analyze -D ${sensor}
#     cd ../macros

#     # Paper plots
#     # python Plot_RisetimeVsX.py          -D ${sensor} -x 1.9 -d
#     # python Plot_NoiseVsX.py          -D ${sensor} -x 1.9 -d
#     # python Plot_NoiseOverallVsX.py          -D ${sensor} -x 1.9 -d
#     # python Plot_JitterVsX.py            -D ${sensor} -x 1.9 -d
#     # python Plot_AmplitudeVsX.py         -D ${sensor} -x 1.9 -y 200.0 -d
#     # python Plot_AmplitudeVsXY.py        -D ${sensor} -z 0.0 -Z 200.0 -d
#     # python Plot_Resolution1D.py         -D ${sensor} -c
#     # python Plot_Efficiency.py           -D ${sensor} -x 1.9 -d
#     python Plot_ResolutionXRecoVsX.py   -D ${sensor} -x 1.9 -d
#     # python Plot_ResolutionTimeVsX.py    -D ${sensor} -x 1.9 -y 150 -d
#     # python Plot_ResolutionTimeVsXForSimulation.py    -D ${sensor} -x 1.9 -y 150 -d

#     # python Plot_RisetimeVsX.py          -D ${sensor} -t -x 1.9 -d
#     # python Plot_JitterVsX.py            -D ${sensor} -t -x 1.9 -d
#     # python Plot_AmplitudeVsX.py         -D ${sensor} -t -x 1.9 -y 200.0 -d
#     # python Plot_AmplitudeVsXY.py        -D ${sensor} -t -z 0.0 -Z 200.0 -d
#     # python Plot_Resolution1D.py         -D ${sensor} -t -d
#     # python Plot_Efficiency.py           -D ${sensor} -t -x 1.9 -d
#     # python Plot_ResolutionXRecoVsX.py   -D ${sensor} -t -x 1.9 -d
#     # python Plot_ResolutionCombinedPosMethod1.py   -D ${sensor} -t -x 1.9 -d
#     # python Plot_ResolutionTimeVsX.py    -D ${sensor} -t -x 1.9 -y 100 -d
# done



# cd ../../TestbeamReco/test/
# python3 ../macros/PlotAmplitudeVsX.py -D LeCroy_W2_3_2_198V_99P9attn
# python ../macros/Paper_CompareXRes.py -D LeCroy_W2_3_2_198V_99P9attn -x 0.8
# python ../macros/Paper_CompareXRes.py -D LeCroy_W4_17_2_222V_99P9attn -x 0.8
# python ../macros/Paper_CompareXRes.py -D HPK_W5_17_2_205V_94attn -x 0.8
# python ../macros/Paper_CompareXRes.py -D LeCroy_W9_15_2_121V_92P3attn -x 0.8
# cd ../../fnal_TestbeamReco/test/
# python3 ../macros/Plot_CompareLaserAmplitudeVsX.py -D HPK_W2_3_2_50T_1P0_500P_50M_E240_180V
# python3 ../macros/PlotAmplitudeDistribution.py -f ../output/HPK_W2_3_2_50T_1P0_500P_50M_E240_180V/ -l ../../TestbeamReco/output/LeCroy_W2_3_2_198V_99P9attn/ -x 150
# # Use for plotting ftbf results + noise-scaled results from laser
# python3 ../macros/Plot_CompareAllTRs.py -D HPK_W2_3_2_50T_1P0_500P_50M_E240_180V -y 80
# python3 ../macros/Plot_CompareAllTRs.py -D HPK_W4_17_2_50T_1P0_500P_50M_C240_204V -y 80
# python3 ../macros/Plot_CompareAllTRs.py -D HPK_W9_15_2_20T_1P0_500P_50M_E600_114V -y 150
# python3 ../macros/Plot_CompareAllTRs_noLandau.py -D HPK_W2_3_2_50T_1P0_500P_50M_E240_180V -y 80
# python3 ../macros/Plot_CompareAllTRs_noLandau.py -D HPK_W4_17_2_50T_1P0_500P_50M_C240_204V -y 80
# python3 ../macros/Plot_CompareAllTRs_noLandau.py -D HPK_W9_15_2_20T_1P0_500P_50M_E600_114V -y 150
# # Use for plotting ftbf results + noise-scaled results from laser
# python3 ../macros/Plot_CompareAllTRsSimSFscaled.py -D HPK_W4_17_2_50T_1P0_500P_50M_C240_204V -y 80
# python3 ../macros/Plot_CompareAllTRsSimSFscaled.py -D HPK_W2_3_2_50T_1P0_500P_50M_E240_180V -y 80
# python3 ../macros/Plot_CompareAllTRsSimSFscaled.py -D HPK_W9_15_2_20T_1P0_500P_50M_E600_114V -y 150
python3 ../macros/Plot_CompareAllTRsSimSFscaled_noLandau.py -D HPK_W4_17_2_50T_1P0_500P_50M_C240_204V -y 80
python3 ../macros/Plot_CompareAllTRsSimSFscaled_noLandau.py -D HPK_W2_3_2_50T_1P0_500P_50M_E240_180V -y 80
python3 ../macros/Plot_CompareAllTRsSimSFscaled_noLandau.py -D HPK_W9_15_2_20T_1P0_500P_50M_E600_114V -y 150

# python3 ../macros/Plot_CompareLaserPosResVsX.py -D HPK_W2_3_2_50T_1P0_500P_50M_E240_180V -y 60
# python3 ../macros/Plot_CompareLaserPosResVsX.py -D HPK_W4_17_2_50T_1P0_500P_50M_C240_204V -y 60
# python3 ../macros/Plot_CompareLaserPosResVsX.py -D HPK_W9_15_2_20T_1P0_500P_50M_E600_114V -y 60
# python3 ../macros/Plot_CompareLaserRisetimeVsX.py -D HPK_W2_3_2_50T_1P0_500P_50M_E240_180V -y 900
# python3 ../macros/Plot_CompareLaserRisetimeVsX.py -D HPK_W4_17_2_50T_1P0_500P_50M_C240_204V -y 900
# python3 ../macros/Plot_CompareLaserRisetimeVsX.py -D HPK_W9_15_2_20T_1P0_500P_50M_E600_114V -y 900
# python3 ../macros/Plot_CompareLaserNoiseVsX.py -D HPK_W2_3_2_50T_1P0_500P_50M_E240_180V -y 5
# python3 ../macros/Plot_CompareLaserNoiseVsX.py -D HPK_W4_17_2_50T_1P0_500P_50M_C240_204V -y 5
# python3 ../macros/Plot_CompareLaserNoiseVsX.py -D HPK_W9_15_2_20T_1P0_500P_50M_E600_114V -y 5
# python3 ../macros/Plot_TR_vs_attenuation.py -f ../output/HPK_W5_17_2_50T_1P0_500P_50M_E600_190V/ -l ../../TestbeamReco/output/
# python3 ../macros/Plot_Jitter_vs_attenuation.py -f ../output/HPK_W5_17_2_50T_1P0_500P_50M_E600_190V/ -l ../../TestbeamReco/output/
# python3 ../macros/Plot_CompareAllTRs.py -D HPK_W5_17_2_50T_1P0_500P_50M_E600_190V -y 80
# python3 ../macros/Plot_CompareLaserRisetimeVsX.py -D HPK_W5_17_2_50T_1P0_500P_50M_E600_190V -y 1000
# python3 ../macros/Plot_CompareLaserPosResVsX.py -D HPK_W5_17_2_50T_1P0_500P_50M_E600_190V -y 60
# python3 ../macros/Plot_CompareLaserNoiseVsX.py -D HPK_W5_17_2_50T_1P0_500P_50M_E600_190V -y 5
# 

# Codes written to study some results for addressing comments by Ryan.
# python3 ../macros/PaperComments_CompareAmpVsX.py -D HPK_W2_3_2_50T_1P0_500P_50M_E240_180V