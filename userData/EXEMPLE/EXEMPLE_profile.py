# Note: all the variables in list format have to contain the same number of elements. 
#       For exemple if you choose 2 experiements, sentence_to_use has to be 2 elements long
#       even if you leave it blank (sentence_to_use = ['','']). Same happens for other list 
#       variables.
#       Some path and file names can be left blank if they don't exist. Check the
#       userData/PROFILE_EXEMPLE/PROFILE_EXEMPLE_profile.py.

# Figures / PDF / GIF settings 
#---------------------------------------------------#
fig_format = 'png' # Figures format ('png' / 'pdf')

create_pdf = 'y'          # ('y' = yes / 'n' = no)
sentence_to_use = ['',''] # Notes to write in pdf, can be left blank

create_gif = 'n' # ('y' = yes / 'n' = no)
timeInterv = 250 # Time in ms between each frame of the gif file, must be integer
#---------------------------------------------------#


# files you want to be loaded and their path location
#---------------------------------------------------#
dataDirOutput = ["/ccc/store/cont003/gen2212/p519don/IGCM_OUT/IPSLCM5A2/PROD/paleo/C30MaTotV1-3X/",
                 "/ccc/store/cont003/gen2212/p25sepul/IGCM_OUT/IPSLCM5A2/PROD/piControl/NORIVER-00/"] # Where model output files are located
filePrefix    = ["C30MaTotV1-3X_SE_4805_4854_1M",
                 "NORIVER-00_SE_2000_2009_1M"] # [filePrefix]_grid_[X].nc with [X] = T, U, V or W
                                               # [filePrefix]_diad_T.nc
                                               # ...

dataDirMask = ["/ccc/work/cont003/gen2212/p519don/",
               "/ccc/work/cont003/gen2212/p25sepul/OHC_BSF/"] # Where mesh mask file is located
maskFile    = ["C30MaTMP_mesh_mask.nc",
               "NORIVER-00_mesh_mask.nc"]

dataDirBathy = ["/ccc/work/cont003/gen2212/p519don/BC_PALEOIPSL/NEMO/30MaTot/",
                ""] # Where bathy and subabsin files are located
bathyFile    = ["bathyORCA2.RupelianTotalV1.nc",
                ""]
subBasinFile = ["subbasins_rupelianTot.nc",
                ""]
#---------------------------------------------------#

# Manual or automatic colormap limits
# 'y' for manual, 'n' for automatic
#---------------------------------------------------#
manual_lim_bathy   = 'n'; min_lim_bathy   = 0; max_lim_bathy   = 100; step_bathy   = 1
manual_lim_sss     = 'n'; min_lim_sss     = 0; max_lim_sss     = 100; step_sss     = 1
manual_lim_zosalin = 'n'; min_lim_zosalin = 0; max_lim_zosalin = 100; step_zosalin = 1
manual_lim_sst     = 'n'; min_lim_sst     = 0; max_lim_sst     = 100; step_sst     = 1
manual_lim_zotemp  = 'n'; min_lim_zotemp  = 0; max_lim_zotemp  = 100; step_zotemp  = 1
manual_lim_zostrf  = 'n'; min_lim_zostrf  = 0; max_lim_zostrf  = 100; step_zostrf  = 1
manual_lim_bstrf   = 'n'; min_lim_bstrf   = 0; max_lim_bstrf   = 100; step_bstrf   = 1
manual_lim_omlnh   = 'n'; min_lim_omlnh   = 0; max_lim_omlnh   = 100; step_omlnh   = 1
manual_lim_omlsh   = 'n'; min_lim_omlsh   = 0; max_lim_omlsh   = 100; step_omlsh   = 1
manual_lim_intpp   = 'n'; min_lim_intpp   = 0; max_lim_intpp   = 100; step_intpp   = 1
manual_lim_epc100  = 'n'; min_lim_epc100  = 0; max_lim_epc100  = 100; step_epc100  = 1
manual_lim_po4     = 'n'; min_lim_po4     = 0; max_lim_po4     = 100; step_po4     = 1
manual_lim_zopo4   = 'n'; min_lim_zopo4   = 0; max_lim_zopo4   = 100; step_zopo4   = 1
manual_lim_no3     = 'n'; min_lim_no3     = 0; max_lim_no3     = 100; step_no3     = 1
manual_lim_zono3   = 'n'; min_lim_zono3   = 0; max_lim_zono3   = 100; step_zono3   = 1
manual_lim_o2      = 'n'; min_lim_o2      = 0; max_lim_o2      = 100; step_o2      = 1
manual_lim_zoo2    = 'n'; min_lim_zoo2    = 0; max_lim_zoo2    = 100; step_zoo2    = 1
#---------------------------------------------------#
