% Script to generate plots from stored output data
% Note: must first run tests to generate output data; this only calls
% plotting functions for existing data


% HSD plots
u = eps('single');
% plotdatabc('opdata_lsir_hsd_3','opdata_glsir_hsd_3','opdata_glsir_bd_hsd_3','hsd3_3',1e3,u,1,15)
% plotdatabc('opdata_lsir_hsd_4','opdata_glsir_hsd_4','opdata_glsir_bd_hsd_4','hsd3_4',1e4,u,0,15)
% plotdatabc('opdata_lsir_hsd_5','opdata_glsir_hsd_5','opdata_glsir_bd_hsd_5','hsd3_5',1e5,u,0,15)
% plotdatabc('opdata_lsir_hsd_6','opdata_glsir_hsd_6','opdata_glsir_bd_hsd_6','hsd3_6',1e6,u,0,15)

plotdatab('opdata_lsir_hsd_3','opdata_glsir_hsd_3','hsd3_3',1e3,u,1,15)
plotdatab('opdata_lsir_hsd_4','opdata_glsir_hsd_4','hsd3_4',1e4,u,0,15)
plotdatab('opdata_lsir_hsd_5','opdata_glsir_hsd_5','hsd3_5',1e5,u,0,15)
plotdatab('opdata_lsir_hsd_6','opdata_glsir_hsd_6','hsd3_6',1e6,u,0,15)
keyboard
% SDQ plots
u = eps('double');
plotdatab('opdata_lsir_sdq_4','opdata_glsir_sdq_4','sdq_4',1e4,u,1,15)
plotdatab('opdata_lsir_sdq_6','opdata_glsir_sdq_6','sdq_6',1e6,u,0,15)
plotdatab('opdata_lsir_sdq_7','opdata_glsir_sdq_7','sdq_7',1e7,u,0,15)
plotdatab('opdata_lsir_sdq_13','opdata_glsir_sdq_13','sdq_13',1e13,u,0,15)

% HDQ plots
plotdatab('opdata_lsir_hdq_3','opdata_glsir_hdq_3','hdq_3',1e3,u,1,15)
plotdatab('opdata_lsir_hdq_4','opdata_glsir_hdq_4','hdq_4',1e4,u,0,15)
plotdatab('opdata_lsir_hdq_8','opdata_glsir_hdq_8','hdq_8',1e8,u,0,15)
plotdatab('opdata_lsir_hdq_11','opdata_glsir_hdq_11','hdq_11',1e11,u,0,15)



