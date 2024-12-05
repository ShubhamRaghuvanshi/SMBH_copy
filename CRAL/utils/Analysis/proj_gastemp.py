import yt
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
from yt.units import Mpc
import os
from pylab import *

simdir_1="/u/sraghu/sim/SHUB/b256/star_sink"
simdir_2="/u/sraghu/sim/SHUB/b256/star_sn_sink"
simdir_3="/u/sraghu/sim/SHUB/b256/star_sink_agn"
simdir_4="/u/sraghu/sim/SHUB/b256/star_sn_sink_agn"

simdir_arr = [simdir_1, simdir_2, simdir_3, simdir_4]
simtyp_arr = ["star+sink", "star+SN+sink", "star+sink+AGN", "star+SN+sink+AGN"]

isim = 0 
for simdir in simdir_arr : 
  if(simdir == simdir_1 ) :
    outnum_arr = ["08", "07", "06"]
  elif(simdir == simdir_2) :
    outnum_arr = ["09", "08", "06"]
  elif(simdir == simdir_3) :
    outnum_arr = ["08", "07", "06"]
  elif(simdir == simdir_4) :
    outnum_arr = ["09", "08", "06"]
  else :
    outnum_arr = []  

  simtyp  = simtyp_arr[isim]
  for outnum in outnum_arr :
    outdir  = simdir + "/output_000" + outnum
    ##
    ds  = yt.load(outdir)
    Z   = ds.current_redshift
    Z_int = round(Z)
    p   = yt.ProjectionPlot(ds, "z", ("gas", "temperature"), weight_field= ("gas", "density") ,width=(6.77,"Mpccm/h")) 
    p.set_cmap(("gas", "temperature"), "twilight_shifted")
    p.set_xlabel("cMpc/h")
    p.set_ylabel("cMpc/h")
    p.set_zlim( ("gas", "temperature"), 1e1, 1e5 )
    p.annotate_title(simtyp)

    #v, c = ds.find_max(("gas", "density"))
    #c = [0.009694367650967489, 0.001431078281960231, 0.13738351506818214]
    #p.set_center((c[0], c[1]))
    #p.zoom(4)

    z_str   = f"{Z:.2f}"
    imgname = "gastemp_proj_" + simtyp+ "_z"  + str(Z_int) +".png"
    p.save(imgname)
    ##

    image       = Image.open(imgname)
    draw        = ImageDraw.Draw(image)
    myfont      = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 42)
    text_z      = f"z:{Z:.2f}"
    draw.text((730, 36), text_z, font=myfont, fill=(255, 255, 255))
    imgname_z =   "gastemp_proj_" + simtyp+ "_z"  + str(Z_int) +".png"
    print("Saving plot: ", imgname_z)
    image.save(imgname_z)
  isim = isim + 1
