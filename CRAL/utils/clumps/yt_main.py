import yt_funcs

'''
simdir_1="/u/sraghu/sim/SHUB/b128/star"
simdir_2="/u/sraghu/sim/SHUB/b128/smbh2"
simdir_3="/u/sraghu/sim/SHUB/b128/smbhstar"
simdir_4="/u/sraghu/sim/SHUB/b128/smbhfeed"
outnum=15
'''

'''
yt_funcs.create_halo_catalog(simdir_1, outnum, 0, 0)
yt_funcs.halomass_vs_rvir(simdir_1, "star", outnum)
'''

'''
#compare hmf
simdir_arr = [simdir_1, simdir_2, simdir_3, simdir_4]
simtyp_arr  = [ "star", "sink_nostar", "sink_star", "sink_feed"]
fitfunc_arr = [1,2,3]
for simdir in simdir_arr: 
  yt_funcs.create_halo_catalog(simdir, outnum, 0, 0)
yt_funcs.compare_hmf(simdir_arr, simtyp_arr, outnum, fitfunc_arr)
'''

'''
simdir_arr = [simdir_1, simdir_3, simdir_4]
simtyp_arr  = [ "star",  "sink_star", "sink_feed"]
#compare sfr
for simdir in simdir_arr:
  yt_funcs.write_starfile(simdir, outnum,0)
yt_funcs.compare_sfr(simdir_arr, simtyp_arr ,outnum )
'''

'''
simdir_arr = [simdir_1, simdir_3, simdir_4]
simtyp_arr  = [ "star",  "sink_star", "sink_feed"]
log_M_min = 3
log_M_max = 13
delta_log_M = 0.1
yt_funcs.compare_star_mass_function(simdir_arr, simtyp_arr, outnum, log_M_min, log_M_max, delta_log_M)
'''

'''
#simdir_arr = [simdir_3, simdir_4]
#simtyp_arr  = [ "sink_star", "sink_feed"]
yt_funcs.create_halo_catalog(simdir_4, outnum, 0, 0)
yt_funcs.sink_vs_star_masstot_halo(simdir_4, "sinkfeed", outnum)
'''

simdir="/u/sraghu/sim/SHUB/b128/smbh"
simtyp="smbhstar"
outnum=15
#yt_funcs.plot_proj_gasdens(simdir, outnum, simtyp)
#yt_funcs.plot_proj_gastemp(simdir, outnum, simtyp)
#yt_funcs.plot_proj_DMdens(simdir, outnum, simtyp)
#yt_funcs.create_halo_catalog(simdir, outnum, 0, 0)
#yt_funcs.plot_proj_halo_gasdens(simdir, outnum, simtyp,1)

yt_funcs.write_halos(simdir, outnum)