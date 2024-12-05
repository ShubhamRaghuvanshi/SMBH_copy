import cfuncs
import yt_funcs


simdir              = "/ptmp/mpa/sraghu/eq/star_sn_sink"
simtyp              = "eq"
outnum              = 9
boxlen_comov        = 6.77

#yt_funcs.location_phasepoints(simdir, outnum, True, 1e5 , 1e10, 0.1, 1000.0)


halonum             = 1
plotsinks           = False
plotstars           = True
plotcircle          = True #fix this
overwrite_densbin   = True
overwrite_sink      = False
overwrite_star      = False
fullbox             = False

dir                 = 'z'
prop                = 'temp'

if(fullbox) :
  cfuncs.plot_halodens(simdir, simtyp,outnum, halonum, plotsinks, plotcircle, plotstars, overwrite_densbin, overwrite_sink, overwrite_star, boxlen_comov, dir, prop, fullbox) 
else :   
  halonum_arr = [1]
  for halonum in halonum_arr :
    cfuncs.plot_halodens(simdir, simtyp,outnum, halonum, plotsinks, plotcircle, plotstars, overwrite_densbin, overwrite_sink, overwrite_star, boxlen_comov, dir, prop, fullbox)
