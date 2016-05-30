
%% BEWARE : must be in the folder of .mex files to run the following :

mex stats_Psigma_PK__Pkalpha.c fill_zk_partial__order1.c fill_zk_partial__order2.c get_zk_partial.c -lm
mex stats_Psigma_PK_meansigmaK__alphabetagamma.c get_zk_partial.c fill_zk_partial__order1.c fill_zk_partial__order2.c  -lm
mex Zk_MaxEnt_PsigmaK.c -lm
mex Zk_MaxEnt_PsigmaK_partialorder1.c  get_zk_partial.c fill_zk_partial__order1.c -lm
mex Zk_MaxEnt_PsigmaK_partialorder2.c fill_zk_partial__order2.c -lm  