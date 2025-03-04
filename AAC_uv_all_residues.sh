 #!/bin/bash
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % A 2>>ouch_A_residues.log"
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % C 2>>ouch_C_residues.log"
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % D 2>>ouch_D_residues.log"
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % E 2>>ouch_E_residues.log"
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % H 2>>ouch_H_residues.log"
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % K 2>>ouch_K_residues.log"
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % L 2>>ouch_L_residues.log"
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % M 2>>ouch_M_residues.log"
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % P 2>>ouch_P_residues.log"
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % T 2>>ouch_T_residues.log"
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % V 2>>ouch_V_residues.log"
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % R 2>>ouch_R_residues.log"
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % F 2>>ouch_F_residues.log"
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % G 2>>ouch_G_residues.log"
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % N 2>>ouch_N_residues.log"
#find DomainFastasFeatures -name '*.AAC.tsv' | parallel --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % I 2>>ouch_I_residues.log"
find DomainFastasFeatures -name '*.AAC.tsv' | parallel -j 20 --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % Q 2>>/local/han_lab/logs/ouch_Q_residues.log"
find DomainFastasFeatures -name '*.AAC.tsv' | parallel -j 20 --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % S 2>>/local/han_lab/logs/ouch_S_residues.log"
find DomainFastasFeatures -name '*.AAC.tsv' | parallel -j 20 --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % W 2>>/local/han_lab/logs/ouch_W_residues.log"
find DomainFastasFeatures -name '*.AAC.tsv' | parallel -j 20 --max-args 1 -I% "Rscript ouch_parallel_mv_AAC_domain.R % Y 2>>/local/han_lab/logs/ouch_Y_residues.log"

