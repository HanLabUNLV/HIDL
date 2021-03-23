#grep -A 2 '$alpha'  A_model_fit_* | grep '\[1' > ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  C_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  D_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  E_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  F_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  G_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  H_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  I_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  K_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  L_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  M_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  N_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  P_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  Q_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  R_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  S_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  T_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  V_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  W_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#grep -A 2 '$alpha'  Y_model_fit_* | grep '\[1' >> ../AAC_uv_results_linker.txt
#
#grep -A 1 -e "aic.c\|AICc" A_* | grep "\[1\|AICc" > ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" C_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" D_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" E_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" F_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" G_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" H_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" I_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" K_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" L_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" M_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" N_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" P_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" Q_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" R_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" S_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" T_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" V_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" W_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc
#grep -A 1 -e "aic.c\|AICc" Y_* | grep "\[1\|AICc" >> ../AAC_uv_results_linker.aicc

grep -A 1 -e "global" A_* | grep -v global > ../AAC_uv_results_linker.optima
grep -A 1 -e "global" C_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" D_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" E_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" F_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" G_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" H_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" I_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" K_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" L_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" M_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" N_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" P_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" Q_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" R_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" S_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" T_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" V_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" W_* | grep -v global >> ../AAC_uv_results_linker.optima
grep -A 1 -e "global" Y_* | grep -v global >> ../AAC_uv_results_linker.optima

