import yt_funcs
import cfuncs

'''
#1 drag 
Dir_sim_gasdrag ='/ptmp/mpa/sraghu/sink_gasdrag' 

simdir_0        = "/ptmp/mpa/sraghu/star"
simdir_1        = Dir_sim_gasdrag + "/star_sink"
simdir_2        = Dir_sim_gasdrag + "/star_sn_sink"
simdir_3        = Dir_sim_gasdrag + "/star_sink_agn"
simdir_4        = Dir_sim_gasdrag + "/star_sn_sink_agn"
simtyp_arr     = ['star', 'star_sink', 'star_sn_sink', 'star_sink_agn', 'star_sn_sink_agn']
simdir_arr     = [simdir_0, simdir_1, simdir_2, simdir_3, simdir_4 ]
cbar_arr       = ['m', 'r', 'g', 'b', 'y' ]
linestyle_arr  = ['-', '-', '-', '-', '-']
#mbar_arr       = ['*'     , '^'     , 's'     ,  'h'     , 's']
sfr_outnum_arr = [11, 10,16,10, 16]
halonum_arr    = [-1, -1,-1,-1, -1]
outnum_1       = [ ]
#outnum_2       = [ 12, 10, 8 ]
#outnum_3       = [ 12, 10, 8 ]
#outnum_4       = [ 12, 10, 8 ]
#outnum_arr     = [ outnum_1, outnum_2, outnum_3, outnum_4]
'''

'''
#1 drag vs no drag 
Dir_sim_nodrag ='/ptmp/mpa/sraghu/sink_nodrag' 
simdir_2       = Dir_sim_nodrag + "/star_sn_sink"
simdir_3       = Dir_sim_nodrag + "/star_sink_agn"
simdir_4       = Dir_sim_nodrag + "/star_sn_sink_agn"
Dir_sim_gasdrag ='/ptmp/mpa/sraghu/sinktest_gasdrag' 
#Dir_sim_gasdrag = Dir_sim_nodrag
simdir_5        = Dir_sim_gasdrag + "/star_sn_sink"
simdir_6        = Dir_sim_gasdrag + "/star_sink_agn"
simdir_7        = Dir_sim_gasdrag + "/star_sn_sink_agn"
simtyp_arr     = ['star_sn_sink(nodrag)',  'star_sn_sink(gasdrag)',   'star_sink_agn(nodrag)',  'star_sink_agn(gasdrag)',  'star_sn_sink_agn(nodrag)',  'star_sn_sink_agn(gasdrag)']
simdir_arr     = [simdir_2              ,  simdir_5               , simdir_3                 , simdir_6                 ,   simdir_4                 ,  simdir_7                   ]
cbar_arr       = ['r'                   ,  'r'                    , 'g'                      , 'g'                      ,   'b'                      ,  'b'                        ]
linestyle_arr  = ['-'                   ,  '--'                   , '-'                      , '--'                     ,  '-'                       , '--']
#mbar_arr       = ['*'     , '^'     , 's'     ,  'h'     , 's'     ]
#sfr_outnum_arr = [11, 10,16,12, 16]
#halonum_arr    = [-1, -1,-1,-1, -1]
outnum_1       = [ ]
#outnum_2       = [ 12, 10, 8 ]
#outnum_3       = [ 12, 10, 8 ]
#outnum_4       = [ 12, 10, 8 ]
#outnum_arr     = [ outnum_1, outnum_2, outnum_3, outnum_4]
'''


'''
Dir_sim='/ptmp/mpa/sraghu/sink_nodrag/' 
print("Dir_sim :",Dir_sim)
simdir_1       = "/u/sraghu/sim/SHUB/b256/sn_eq"
simdir_2       = "/u/sraghu/sim/SHUB/b256/sn_neq"
simdir_3       = "/ptmp/mpa/sraghu/tests_nosink/sn_neq_18"
simtyp_arr     = ['star+sn_eq_lowres', 'star+sn_neq_lowres', 'star+sn_neq_highres' ]
simdir_arr     = [simdir_1, simdir_2, simdir_3]
cbar_arr       = ['r'     , 'g'     , 'b']
linestyle_arr  = [ '-'    , '-'     , '-']
mbar_arr       = ['*'     , '^'     , '^']
sfr_outnum_arr = [13,13,13]
halonum_arr    = [-1, -1, -1]
outnum_1       = []
outnum_2       = []
outnum_arr     = []
'''

'''
simdir_1       = "/ptmp/mpa/sraghu/sink_nodrag/star_sn_sink"
simdir_2       = "/ptmp/mpa/sraghu/sink_gasdrag/star_sn_sink"
simdir_3       = "/u/sraghu/sim/SHUB/b256/gasdrag/star_sn_sink"
simdir_4       = "/u/sraghu/sim/SHUB/b256/nodrag/star_sn_sink"

simtyp_arr     = ['star+sn+sink(bursty)',  'star+sn+sink(smooth)']
simdir_arr     = [simdir_2, simdir_3]
cbar_arr       = ['r'     , 'b'     ]
linestyle_arr  = [ '-'    , '-'     ]
mbar_arr       = ['*'     , '^'     ]
#sfr_outnum_arr = [17,16]
sfr_outnum_arr = [14,13]
outnum_max   = [-1,-1]
halonum_arr    = [-1, -1]
outnum_1       = [14]
outnum_2       = [13]
outnum_arr     = [outnum_1, outnum_2]
'''


simdir_1       = "/u/sraghu/sim/SHUB/b256/edd_test/star_sink_agn"
simdir_2       = "/ptmp/mpa/sraghu/smooth/msink_1e7"

simtyp_arr     = ['star+sn+sink+drag']

simdir_arr     = [simdir_1]
cbar_arr       = ['g'      ]
linestyle_arr  = [ '-'    , '-'     ]
mbar_arr       = ['*'     , '*'     ]
sfr_outnum_arr = [14,14]
halonum_arr    = [-1, -1]

outnum_max   = [14,14]

outnum_1       = [14]
outnum_2       = [14]
outnum_arr     = [outnum_1]


'''
simdir_1       = "/u/sraghu/sim/SHUB/b256/star_sn_sink"
simtyp_arr     = ['star_sn_sink']
simdir_arr     = [simdir_1]
cbar_arr       = ['b'  ]
linestyle_arr  = [ '-' ]
mbar_arr       = ['*' ]
outnum_1       = [14]
outnum_arr     = [outnum_1]
halonum_arr    = [-1]
sfr_outnum_arr = [14]
outnum_max     = [-1]
'''

#prop_arr      = ['sfr', 'hmf' , 'dmdens_proj' , 'dens_proj', 'temp_proj' , 'phaseplot']#,  'globalsink_vs_z',  'haloprops']
prop_arr = ['globalsink_vs_z']
Nout     = len(outnum_1)


for prop in prop_arr :
  print('######################## prop :', prop, "##########################")
  if(prop=='sfr') :
    #compare SFR
    plotanal=False
    plotdata=True
    overwrite=True
    physical = True
    fullbox = True
    #write star files
    iout=0
    for simdir in simdir_arr:
      outnum = sfr_outnum_arr[iout]
      if(halonum_arr[iout] > 0) :
        halo_attrb = yt_funcs.give_haloprops(simdir, outnum, halonum_arr[iout], physical, False)
        center = halo_attrb[0]
        radius = halo_attrb[1]
        mass   = halo_attrb[2] 
        fullbox = False     
      else :
        center = [-1,-1,-1]
        radius = -1  
        mass = -1e9
      print("halo center, radius, mass : ", center, radius, mass/1e9)  
      yt_funcs.write_starfile(simdir, outnum, overwrite, radius, center)
      iout = iout+1
    # plot sfr  
    yt_funcs.compare_sfr(simdir_arr, simtyp_arr ,sfr_outnum_arr, cbar_arr, linestyle_arr ,plotanal, plotdata, fullbox )
  if(prop=='hmf') :
    # write halo catalog 
    overwrite=False
    for iout in range(0,Nout) :
      isim=0
      hmf_outnum_arr = []
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]  
        hmf_outnum_arr.append(outnum)
        yt_funcs.create_halo_catalog(simdir, outnum, 0, overwrite)
        isim = isim + 1
      fitfunc_arr = []
      print('hmf_outnum_arr : ', hmf_outnum_arr)
      yt_funcs.compare_hmf(simdir_arr, simtyp_arr, cbar_arr, mbar_arr, hmf_outnum_arr, fitfunc_arr)
  if(prop=='dmdens_proj') :
    boxlen_comov=6.77
    overwrite=False
    physical = False
    fullbox = True
    projdir = 'y'
    
    for iout in range(0,Nout) :
      isim=0
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]
        simtyp = simtyp_arr[isim]
        if(fullbox or halonum_arr[iout] < 0) :
          radius = -1
          center = [0.5, 0.5, 0.5]
        else :
          halo_attrb = yt_funcs.give_haloprops(simdir, outnum, halonum_arr[iout], physical, False)
          center = halo_attrb[0]
          radius = halo_attrb[1]
          aexp   = halo_attrb[3]
          #if(physical) :
          #  radius = radius / aexp  
        yt_funcs.plot_proj_dmdens(simdir, outnum, simtyp, radius, center, projdir, boxlen_comov)
        isim = isim + 1
  if(prop=='dens_proj') :
    boxlen_comov=6.77
    overwrite=False
    physical = False
    fullbox = False
    projdir = 'y'
    
    for iout in range(0,Nout) :
      isim=0
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]
        simtyp = simtyp_arr[isim]
        if(fullbox or halonum_arr[iout] < 0) :
          radius = -1
          center = [0.5, 0.5, 0.5]
        else :
          halo_attrb = yt_funcs.give_haloprops(simdir, outnum, halonum_arr[iout], physical, False)
          center = halo_attrb[0]
          radius = halo_attrb[1]
          aexp   = halo_attrb[3]
          #if(physical) :
          #  radius = radius / aexp  
        yt_funcs.plot_proj_gasdens(simdir, outnum, simtyp, radius, center, projdir, boxlen_comov)
        isim = isim + 1
  if(prop=='temp_proj') :
    boxlen_comov=6.77
    overwrite=False
    physical = False
    fullbox = True
    projdir = 'z'
    
    for iout in range(0,Nout) :
      isim=0
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]
        simtyp = simtyp_arr[isim]
        if(fullbox or halonum_arr[iout] < 0) :
          radius = -1
          center = [0.5, 0.5, 0.5]
        else :
          halo_attrb = yt_funcs.give_haloprops(simdir, outnum, halonum_arr[iout], physical, False)
          center = halo_attrb[0]
          radius = halo_attrb[1]
          aexp   = halo_attrb[3]
          #if(physical) :
          #  radius = radius / aexp  
        yt_funcs.plot_proj_gastemp(simdir, outnum, simtyp, radius, center, projdir, boxlen_comov)
        isim = isim + 1
  if(prop=='metalicity_proj') :
    boxlen_comov=6.77
    overwrite=False
    for iout in range(0,Nout) :
      isim=0
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]
        simtyp = simtyp_arr[isim]
        yt_funcs.plot_proj_gastemp(simdir, outnum, simtyp, boxlen_comov, overwrite)
        isim = isim + 1
  if(prop=='phaseplot') :
    boxlen_comov=6.77
    overwrite=False
    for iout in range(0,Nout) :
      isim=0
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]
        simtyp = simtyp_arr[isim]
        yt_funcs.tn_phaseplot(simdir, outnum, simtyp, boxlen_comov, overwrite)
        isim = isim + 1
  if(prop=='globalsink_vs_z')  :
    boxlen_comov = 6.77
    #cfuncs.compare_sinkmasstot(simdir_arr, simtyp_arr, cbar_arr, linestyle_arr, boxlen_comov)
    #cfuncs.plot_sink_pos_vs_z(simdir_arr, cbar_arr, simtyp_arr, 4, boxlen_comov, True, 100)
    #cfuncs.plot_sink_pos_vs_z(simdir_arr, cbar_arr, simtyp_arr, boxlen_comov, 100)
    #cfuncs.plot_sink_relpos_vs_z(simdir_arr, cbar_arr, simtyp_arr, 9, boxlen_comov, 100, outnum_max)
    #cfuncs.hist_z_sinkspawn(simdir_arr, simtyp_arr, cbar_arr)
    cfuncs.sink_m_dmdt(simdir_arr, 1, simtyp_arr, cbar_arr, boxlen_comov)
  if(prop=='haloprops') :
    overwrite_sink = True
    overwrite_star = False
    plotsinks      = True
    boxlen_comov = 6.77
   
    for iout in range(0,Nout) :
      isim = 0
      haloprop_outnum_arr = []
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout] 
        haloprop_outnum_arr.append(outnum)
        isim = isim + 1
      print('haloprop_outnum_arr :', haloprop_outnum_arr) 
      cfuncs.haloprops(simdir_arr, simtyp_arr, haloprop_outnum_arr, cbar_arr, linestyle_arr, 20, boxlen_comov, overwrite_star, overwrite_sink, plotsinks)



   








