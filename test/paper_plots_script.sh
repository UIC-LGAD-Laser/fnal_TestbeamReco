cd ../../test

# HPK 500um Pitch sensors
# -----------------------

# HPK_500P=("HPK_W9_15_4_20T_0P5_500P_50M_E600_110V")
HPK_500P=("HPK_W4_17_2_50T_1P0_500P_50M_C240_204V" "HPK_W5_17_2_50T_1P0_500P_50M_E600_190V")

for sensor in "${HPK_500P[@]}"; do
    printf "\nRunning over ${sensor} sensor\n"
    cd ../test
    ./MyAnalysis -A InitialAnalyzer -D ${sensor}
    cd ../macros
    python FindDelayCorrections.py -D ${sensor}
    python FindInputHistos4YReco.py -D ${sensor} -I
    cd ../test
    ./MyAnalysis -A Analyze -D ${sensor}
    cd ../macros

done


# cd ../../TestbeamReco/test/
# python3 ../macros/PlotAmplitudeVsX.py -D LeCroy_W2_3_2_198V_99P9attn
# python ../macros/Paper_CompareXRes.py -D LeCroy_W2_3_2_198V_99P9attn -x 0.8
# python ../macros/Paper_CompareXRes.py -D LeCroy_W4_17_2_222V_99P9attn -x 0.8
# python ../macros/Paper_CompareXRes.py -D HPK_W5_17_2_205V_94attn -x 0.8
# python ../macros/Paper_CompareXRes.py -D LeCroy_W9_15_2_121V_92P3attn -x 0.8
# cd ../../fnal_TestbeamReco/test/
# python3 ../macros/Plot_CompareLaserAmplitudeVsX.py -D HPK_W2_3_2_50T_1P0_500P_50M_E240_180V
# python3 ../macros/PlotAmplitudeDistribution.py -f ../output/HPK_W2_3_2_50T_1P0_500P_50M_E240_180V/ -l ../../TestbeamReco/output/LeCroy_W2_3_2_198V_99P9attn/ -x 150
# python3 ../macros/Plot_TR_vs_attenuation.py -f ../output/HPK_W5_17_2_50T_1P0_500P_50M_E600_190V/ -l ../../TestbeamReco/output/
# python3 ../macros/Plot_Jitter_vs_attenuation.py -f ../output/HPK_W5_17_2_50T_1P0_500P_50M_E600_190V/ -l ../../TestbeamReco/output/
# python3 ../macros/Plot_CompareAllTRs.py -D HPK_W2_3_2_50T_1P0_500P_50M_E240_180V -y 80
# python3 ../macros/Plot_CompareAllTRs.py -D HPK_W4_17_2_50T_1P0_500P_50M_C240_204V -y 80
# python3 ../macros/Plot_CompareAllTRs.py -D HPK_W5_17_2_50T_1P0_500P_50M_E600_190V -y 80
# python3 ../macros/Plot_CompareAllTRs.py -D HPK_W9_15_2_20T_1P0_500P_50M_E600_114V -y 150
# python3 ../macros/Plot_CompareLaserPosResVsX.py -D HPK_W2_3_2_50T_1P0_500P_50M_E240_180V -y 60
# python3 ../macros/Plot_CompareLaserPosResVsX.py -D HPK_W4_17_2_50T_1P0_500P_50M_C240_204V -y 60
# python3 ../macros/Plot_CompareLaserPosResVsX.py -D HPK_W5_17_2_50T_1P0_500P_50M_E600_190V -y 60
# python3 ../macros/Plot_CompareLaserPosResVsX.py -D HPK_W9_15_2_20T_1P0_500P_50M_E600_114V -y 60
# python3 ../macros/Plot_CompareLaserRisetimeVsX.py -D HPK_W2_3_2_50T_1P0_500P_50M_E240_180V -y 900
# python3 ../macros/Plot_CompareLaserRisetimeVsX.py -D HPK_W4_17_2_50T_1P0_500P_50M_C240_204V -y 900
# python3 ../macros/Plot_CompareLaserRisetimeVsX.py -D HPK_W5_17_2_50T_1P0_500P_50M_E600_190V -y 1000
# python3 ../macros/Plot_CompareLaserRisetimeVsX.py -D HPK_W9_15_2_20T_1P0_500P_50M_E600_114V -y 900
# python3 ../macros/Plot_CompareLaserNoiseVsX.py -D HPK_W2_3_2_50T_1P0_500P_50M_E240_180V -y 5
# python3 ../macros/Plot_CompareLaserNoiseVsX.py -D HPK_W4_17_2_50T_1P0_500P_50M_C240_204V -y 5
# python3 ../macros/Plot_CompareLaserNoiseVsX.py -D HPK_W5_17_2_50T_1P0_500P_50M_E600_190V -y 5
# python3 ../macros/Plot_CompareLaserNoiseVsX.py -D HPK_W9_15_2_20T_1P0_500P_50M_E600_114V -y 5
