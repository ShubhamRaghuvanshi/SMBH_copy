import yt 
from yt.extensions.astro_analysis.halo_analysis import HaloCatalog
import yt.visualization.eps_writer as eps
import os
import numpy as np
import matplotlib.pyplot as plt
from halo_mass_function import *
from sfr_spectrum import * 
from PIL import Image, ImageDraw, ImageFont

# HOP(HaloFinderOverDensity
finder_kwargs_hop = {
    "threshold":80,  
    "ptype":"DM"
}

# FOF
finder_kwargs_fof = {
    "padding": 0.02,       
    "link": 0.2, 
    "ptype": "DM"
}

# ROCKSTAR
finder_kwargs_rockstar = {
    "num_particles_min": 100,  
    "halo_finder_path": "/path/to/rockstar_executable"
}

FinderMethods = ['hop', 'fof', 'rockstar' ]
FinderKwargs  = [finder_kwargs_hop, finder_kwargs_fof, finder_kwargs_rockstar ] 

def create_halo_catalog(simdir, outnum, MethId, overwrite) :
  outnum_char=str(outnum).zfill(5)
  outdir = simdir+'/output_' + outnum_char
  halocatalog_file = outdir + '/info_'+ outnum_char + '/info_' + outnum_char +'.0.h5'

  if(overwrite) :
    command = f'rm {sfr_file}'
    print("eceuting", command)
    os.system(command)   

  if( not os.path.exists(halocatalog_file) ) :
    ds = yt.load(outdir)
    halo_catalog = HaloCatalog(data_ds=ds, finder_method=FinderMethods[MethId], finder_kwargs=FinderKwargs[MethId], output_dir = outdir)
    halo_catalog.create()
  else :
    print("already exists :", halocatalog_file)

log_M_min = 9
log_M_max = 13
delta_log_M = 0.1
label_arr = ["PS", "JN", "ST", "WR", "TK" ]
cbar = ['r', 'g', 'b', 'k']
mbar = ["o", "^", "s", "h", "D" ]

# simply plot one hmf
def plot_hmf(simdir, outnum, fit_func): 

  outnum_char      = str(outnum).zfill(5)
  outdir           = simdir+'/output_' + outnum_char
  halocatalog_file = outdir + '/info_'+ outnum_char + '/info_' + outnum_char +'.0.h5'

  ds               = yt.load(halocatalog_file)
  Z                = ds.current_redshift
  ad               = ds.all_data()
    
  NumPart          = np.array( ad['halos', 'particle_number'] )  
  HMass            = np.array( ad["halos", "particle_mass"])
  halo_masses      = []
  
  for ihalo in range (0,len(NumPart)):
    if( NumPart[ihalo] > 60  ):
      halo_masses.append(HMass[ihalo])   

  Bins        = 10**np.arange(log_M_min, log_M_max + delta_log_M, delta_log_M)
  hist, edges = np.histogram(halo_masses, bins=Bins)
  bin_centers = (edges[1:] + edges[:-1]) / 2
  bin_widths  = edges[1:] - edges[:-1]
  dn_dlogm    = hist
  err         = np.sqrt(hist)

  hmf    = HaloMassFcn(halos_ds=ds, log_mass_min=math.log10(np.min(halo_masses)), log_mass_max=math.log10(np.max(halo_masses)), 
                       fitting_function=fit_func, num_sigma_bins= len(bin_widths) +1)
  n_anal = hmf.dndM_dM_analytic[:-1]/0.703
  n_anal = n_anal*1000.0*np.log10(bin_widths)
  m_anal = hmf.masses_analytic[:-1]

  anal_label = label_arr[fit_func -1]
  if(fit_func == 1):
      plt.title("DM HMF compared with Press-Schechter[PS]")
      funcname = 'ps'
  elif (fit_func == 2):
      plt.title("DM HMF compared with Jenkins[JN]")
      funcname = 'jn'
  elif (fit_func == 3):
      plt.title("DM HMF compared with Sheth-Tormen[ST]")
      funcname = 'st'        
  elif (fit_func == 4):
      plt.title("DM HMF compared with Warren[WR]")
      funcname = 'wr'
  elif (fit_func == 5):
      plt.title("DM HMF compared with Tinker[TK]")
      funcname = 'tk'
  else:
      ValueError("Invalid value for fitting function")
  

  plt.xscale('log')
  plt.yscale('log')
  plt.errorbar(bin_centers, dn_dlogm, yerr=err, fmt="o", capsize=3, label=f'[sim]Z={Z:.2f}')
  plt.plot(m_anal, n_anal, color='red',label=anal_label)
  plt.legend(loc='upper right')
  plt.savefig("hmf.png")
  
# compare hmf from with same outnum at same redshift
def compare_hmf(simdir_arr, simname_arr ,outnum ,fitfunc_arr): 

  outnum_char      = str(outnum).zfill(5)

  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel("Halo Mass (Msun)")
  plt.ylabel("Number of Halos")
  iout=0
  for simdir in simdir_arr :

    outdir           = simdir+'/output_' + outnum_char
    halocatalog_file = outdir + '/info_'+ outnum_char + '/info_' + outnum_char +'.0.h5'

    ds               = yt.load(halocatalog_file)
    Z                = ds.current_redshift
    ad               = ds.all_data()  
    NumPart          = np.array( ad['halos', 'particle_number'] )  
    HMass            = np.array( ad["halos", "particle_mass"])
    halo_masses      = []

    if (iout == 0 ):
      title = "DM HMF" + f' [Z={Z:.2f}]'
      plt.title(title)

    for ihalo in range (0,len(NumPart)):
      if( NumPart[ihalo] > 60  ):
        halo_masses.append(HMass[ihalo])   

    Bins        = 10**np.arange(log_M_min, log_M_max + delta_log_M, delta_log_M)
    hist, edges = np.histogram(halo_masses, bins=Bins)
    bin_centers = (edges[1:] + edges[:-1]) / 2
    bin_widths  = edges[1:] - edges[:-1]
    dn_dlogm    = hist
    err         = np.sqrt(hist)  
    labl        = simname_arr[iout]
    plt.errorbar(bin_centers, dn_dlogm, yerr=err, fmt=mbar[iout], color=cbar[iout],capsize=3, label=labl)
    iout = iout + 1

  # now plot the analytical hmf
  ifunc=0
  for fit_func in fitfunc_arr : 
    hmf    = HaloMassFcn(halos_ds=ds, log_mass_min=math.log10(np.min(halo_masses)), log_mass_max=math.log10(np.max(halo_masses)), 
                       fitting_function=fit_func, num_sigma_bins= len(bin_widths) +1)
    n_anal = hmf.dndM_dM_analytic[:-1]/0.703
    n_anal = n_anal*1000.0*np.log10(bin_widths)
    m_anal = hmf.masses_analytic[:-1]
    anal_label = label_arr[fit_func -1]
    plt.plot(m_anal, n_anal, color=cbar[ifunc],label=anal_label)
    plt.legend()
    ifunc = ifunc + 1
  plt.tight_layout()  
  imgname="comp_hmf.png"
  print("Saving : ", imgname)
  plt.savefig(imgname)


#star functions 
def write_starfile(simdir, outnum, overwrite):
  outnum_char      = str(outnum).zfill(5)
  outdir           = simdir+'/output_' + outnum_char
  sfr_file         = outdir + "/SFR_" +  outnum_char + ".out" 
  print("outdir : ", outdir)

  if(overwrite) :
    command = f'rm {sfr_file}'
    print("eceuting", command)
    os.system(command)   

  if( not os.path.exists(sfr_file) ) :     
    ds         = yt.load(outdir)
    ad         = ds.all_data()
    mass       = ad[("star", "particle_mass")].in_units('Msun')
    age        = ad[("star", "age")].in_units('Myr')
    ct         = ad[("star", "particle_birth_time")]

    threshold  = ds.quan(0.1, "Myr")
    mass_old   = mass[age > threshold]
    ct_old     = ct[age > threshold]
    sp         = ds.sphere("c", (1000.0, "Mpc"))

    sfr  = StarFormationRate(ds, star_mass=mass_old, star_creation_time=ct_old, volume=sp.volume())
    #sfr = StarFormationRate(ds, star_mass=mass_old)
    #sfr = StarFormationRate(ds, data_source=ad, star_creation_time=ct_old)
    print("Writing : ", sfr_file)
    #print( "Number of stars particles :", len(star_mass))
    sfr.write_out(name=sfr_file)
  else :
    print(sfr_file, "already exists")  

def read_star_prop(simdir, outnum, propnum):
  outnum_char = str(outnum).zfill(5)
  outdir = simdir+'/output_' + outnum_char
  outfile = outdir+ "/SFR_"+ outnum_char + ".out"
  cmass_arr = []
  with open(outfile, 'r') as clumpfile:
    numline=0
    for line in clumpfile:
      numline = numline + 1
      if(numline > 1):
        cmass = float(line.split()[propnum])      
        cmass_arr.append(cmass) 
  return cmass_arr

# compare hmf from with same outnum at same redshift
def compare_sfr(simdir_arr, simname_arr ,outnum ): 

  outnum_char      = str(outnum).zfill(5)

  #plt.xscale('log')
  plt.yscale('log')
  plt.xlabel("Redshift(Z)")
  plt.ylabel(r"SFR $M_\odot$/yr/Mpcc ")
  
  mass_thres = 0.0 
  tage_thres = 0.0
  iout=0
  middlename=''
  for simdir in simdir_arr :
    outdir      = simdir+'/output_' + outnum_char
    sfr_file    = outdir + '/SFR_'+ outnum_char +'.out'
    Redshift    = []
    SFR         = []
    
    with open(sfr_file, 'r') as file:
      numline     = 0
      for line in file:
        numline = numline + 1
        if(numline > 1):
          tage = float(line.split()[0]) 
          mass = float(line.split()[5]) 
          if(mass > mass_thres and tage > tage_thres):
            sfrcc = float(line.split()[4])
            z     = float(line.split()[2])
            Redshift.append(z)
            SFR.append(sfrcc)
    plt.plot(Redshift, SFR, linestyle="-"  ,color=cbar[iout], label=simname_arr[iout])  
    
    middlename = middlename + "_" + simname_arr[iout]
    iout = iout + 1    

  Nanal = len(Redshift)
  Z_anal   =   []
  SFR_anal =   []

  for i in range(0, Nanal) :
    z   = Redshift[i]
    sfr = 0.015 * (1 + z)**2.7 / (1 + ((1+z)/2.9)**5.6)
    SFR_anal.append(sfr)
  
  plt.plot(Redshift, SFR_anal, linestyle="-"  ,color='k', label="Madau, Dickinson")

  #imgname = "sfr" + middlename + ".png"  
  imgname = "sfr_cmp.png"
  plt.legend()
  print("Saving : ", imgname)
  plt.tight_layout()
  plt.savefig(imgname)


def prop_diffsq(simdir_1, simdir_2, outnum):
  outnum_char      = str(outnum).zfill(5)
  outdir_1      = simdir_1+'/output_' + outnum_char
  outdir_2      = simdir_2+'/output_' + outnum_char

  ds_1 = yt.load(outdir_1)
  ds_2 = yt.load(outdir_2)
  ad_1 = ds_1.all_data()
  ad_2 = ds_2.all_data()

  Z_1 = ds_1.current_redshift
  Z_2 = ds_2.current_redshift

  prop_1 = np.array(ad_1['deposit', 'DM_density'])
  prop_2 = np.array(ad_2['deposit', 'DM_density'])

  print("Data lengths : ", len(prop_1), len(prop_2))
  if( len(prop_1) == len(prop_2) and Z_1 == Z_2) :
    sum_sq = 0  
    for i in range(0, len(prop_1)):
      diff = prop_1[i] - prop_2[i]
      sum_sq = sum_sq + diff*diff
    print("Sum of squares of difference : ", sum_sq)  
  else : 
    print( len(prop_1) ,len(prop_2), Z_1, Z_2 )
    raise ValueError("Apples and oranges")

def sink_vs_star_masstot_halo(simdir, simname, outnum):
  outnum_char      = str(outnum).zfill(5)
  outdir           = simdir+'/output_' + outnum_char
  halocatalog_file = outdir + '/info_'+ outnum_char + '/info_' + outnum_char +'.0.h5'

  ds       = yt.load(outdir)
  halo_ds  = yt.load(halocatalog_file)

  ad       = ds.all_data()
  halo_ad  = halo_ds.all_data()

  Z        = ds.current_redshift

  NumPart   = np.array( halo_ad['halos', 'particle_number'] )  
  Rvir      = np.array( halo_ad['halos', 'virial_radius'].in_units("kpc"))
  halo_posx = np.array( halo_ad['halos', 'particle_position_x'] )
  halo_posy = np.array( halo_ad['halos', 'particle_position_y'] )
  halo_posz = np.array( halo_ad['halos', 'particle_position_z'] )

  sinkmass_arr_halo = []
  starmass_arr_halo = []

  print("Nhalo tot: ",len(NumPart))
  for ihalo in range (0,len(NumPart)):
    if( NumPart[ihalo] > 60  ):
      center    = [halo_posx[ihalo], halo_posy[ihalo], halo_posz[ihalo]]
      sp        = ds.sphere(center, (40.0*Rvir[ihalo], "kpc"))
      star_mass = np.array(sp[("star", "particle_mass")].in_units('Msun'))
      sink_mass = np.array(sp[("sink", "particle_mass")].in_units('Msun'))
      sinkmasstot = 0.0 
      for ipart in range(0, len(sink_mass) ) :
        sinkmasstot = sinkmasstot + sink_mass[ipart] 
      starmasstot = 0.0 
      for ipart in range(0, len(star_mass) ) :
        starmasstot = starmasstot + star_mass[ipart] 
      if(starmasstot > 0 and sinkmasstot > 0) : 
        sinkmass_arr_halo.append(sinkmasstot)
        starmass_arr_halo.append(starmasstot)
      #print(sinkmasstot, starmasstot)

  print("plotting for Nhalo : ",len(sinkmass_arr_halo), len(sinkmass_arr_halo))
  plt.scatter( starmass_arr_halo, sinkmass_arr_halo, marker='o', label = simname)
  
  plt.title("BH and star mass in Haloes for Z={:.2f}".format(Z))
  plt.xlabel("Mstar[Msun]")
  plt.ylabel("Mbh[Msun]")
  plt.xscale("log")
  plt.yscale("log")
  plt.legend()
  plt.tight_layout()
  imgname = "sink_vs_star_mass_perhalo_" + simname +  ".png"
  print("Saving : ", imgname)
  plt.savefig(imgname)


def halomass_vs_rvir(simdir, simname, outnum):
  outnum_char      = str(outnum).zfill(5)
  outdir           = simdir+'/output_' + outnum_char
  halocatalog_file = outdir + '/info_'+ outnum_char + '/info_' + outnum_char +'.0.h5'

  ds       = yt.load(outdir)
  halo_ds  = yt.load(halocatalog_file)

  ad       = ds.all_data()
  halo_ad  = halo_ds.all_data()

  Z        = ds.current_redshift

  NumPart   = np.array( halo_ad['halos', 'particle_number'] )
  Rvir      = np.array( halo_ad['halos', 'virial_radius'].in_units("kpc"))
  halo_mass = np.array( halo_ad['halos', 'particle_mass'].in_units("Msun") )

  halomass_arr = []
  Rvir_arr     = []   

  print("Nhalo tot: ",len(NumPart))
  for ihalo in range (0,len(NumPart)):
    if( NumPart[ihalo] > 60  ):
      halomass_arr.append(halo_mass[ihalo])
      Rvir_arr.append( Rvir[ihalo] )  
  
  plt.scatter(Rvir_arr, halomass_arr, color='red',marker='o')
  plt.title("Z={:.2f}".format(Z))
  plt.xlabel("Rvir[Kpc]")
  plt.ylabel("Mhalo[Msun]")
  #plt.xscale("log")
  plt.yscale("log")
  #plt.legend()
  plt.tight_layout()
  imgname = "halomass_vs_rvir.png" 
  plt.savefig(imgname)


def compare_star_mass_function(simdir_arr, simname_arr, outnum, Mmin_log, Mmax_log, dM_log): 

  outnum_char  = str(outnum).zfill(5)
  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel("Mstar[Msun]")
  plt.ylabel("Number")

  isim =0 
  for simdir in simdir_arr:
    outdir       = simdir+'/output_' + outnum_char
    ds           = yt.load(outdir)
    Z            = ds.current_redshift
    if(isim == 0 ) :
      plt.title("Star Mass Function for  Z={:.2f}".format(Z))

    ad           = ds.all_data()
    star_masses  = np.array( ad["star", "particle_mass"].in_units('Msun'))
    print("Found Nstars: ", len(star_masses))  

    Bins         = 10**np.arange(Mmin_log, Mmax_log + dM_log, dM_log)
    hist, edges  = np.histogram(star_masses, bins=Bins)
    bin_centers  = (edges[1:] + edges[:-1]) / 2
    bin_widths   = edges[1:] - edges[:-1]
    dn_dlogm     = hist
    err          = np.sqrt(hist)
    labl         = simname_arr[isim] 
    plt.errorbar(bin_centers, dn_dlogm, yerr=err, fmt=mbar[isim], color=cbar[isim],capsize=3, label=labl)
    isim = isim + 1

  plt.legend() 
  plt.tight_layout()
  figname = "star_mass_func.png"
  print("saving : ", figname)
  plt.savefig(figname)
       

def compare_sink_mass_function(simdir_arr, simname_arr, outnum, Mmin_log, Mmax_log, dM_log): 

  outnum_char  = str(outnum).zfill(5)
  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel("Msink[Msun]")
  plt.ylabel("Number")

  isim =0 
  for simdir in simdir_arr:
    outdir       = simdir+'/output_' + outnum_char
    ds           = yt.load(outdir)
    Z            = ds.current_redshift
    if(isim == 0 ) :
      plt.title("Sink Mass Function for  Z={:.2f}".format(Z))

    ad           = ds.all_data()
    star_masses  = np.array( ad["sink", "particle_mass"].in_units('Msun'))
    print("Found Nstars: ", len(star_masses))  

    Bins         = 10**np.arange(Mmin_log, Mmax_log + dM_log, dM_log)
    hist, edges  = np.histogram(star_masses, bins=Bins)
    bin_centers  = (edges[1:] + edges[:-1]) / 2
    bin_widths   = edges[1:] - edges[:-1]
    dn_dlogm     = hist
    err          = np.sqrt(hist)
    labl         = simname_arr[isim] 
    plt.errorbar(bin_centers, dn_dlogm, yerr=err, fmt=mbar[isim], color=cbar[isim],capsize=3, label=labl)
    isim = isim + 1

  plt.legend() 
  plt.tight_layout()
  figname = "sink_mass_func.png"
  print("saving : ", figname)
  plt.savefig(figname)

def plot_proj_gasdens(simdir, outnum, simtyp): 
  outnum_char  = str(outnum).zfill(5)
  outdir       = simdir+'/output_' + outnum_char
  imgname      =  simdir  + "/temp_gasdens_proj_"+ simtyp + "_" + outnum_char  +".png"
  imgname_z    =  simdir +  "/gasdens_proj_" + simtyp + "_"+outnum_char +  ".png"

  ds = yt.load(outdir)
  Z  = ds.current_redshift
  p  = yt.ProjectionPlot(ds, "z", ("gas", "density"), weight_field=("gas", "density"),width=(10.0,"Mpccm/h"))
  p.save(imgname)

  image = Image.open(imgname)
  draw = ImageDraw.Draw(image)
  myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 42)
  text = f"Z:{Z:.2f}"
  text_color = (255, 0, 0)
  draw.text((730, 36), text, font=myfont, fill=(0, 0, 0))
  print("Saving plot: ", imgname_z)
  image.save(imgname_z)


def plot_proj_gastemp(simdir, outnum, simtyp): 
  outnum_char  = str(outnum).zfill(5)
  outdir       = simdir+'/output_' + outnum_char
  imgname      =  simdir  + "/temp_gastemp_proj_"+ simtyp + "_" + outnum_char  +".png"
  imgname_z    =  simdir +  "/gastemp_proj_" + simtyp + "_"+outnum_char +  ".png"

  ds = yt.load(outdir)
  Z  = ds.current_redshift

  p = yt.ProjectionPlot(ds, "z", ("gas", "temperature"), weight_field=("gas", "density") ,width=(10.0,"Mpccm/h")) 
  p.set_cmap(("gas", "temperature"), "afmhot")
  p.set_zlim( ("gas", "temperature"), 1e2, 1e6 )
  p.save(imgname)

  image = Image.open(imgname)
  draw = ImageDraw.Draw(image)
  myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 42)
  text = f"Z:{Z:.2f}"
  text_color = (255, 0, 0)
  draw.text((730, 36), text, font=myfont, fill=(0, 0, 0))
  print("Saving plot: ", imgname_z)
  image.save(imgname_z)

def plot_proj_DMdens(simdir, outnum, simtyp): 
  outnum_char  = str(outnum).zfill(5)
  outdir       =  simdir + "/output_" + outnum_char
  imgname      =  simdir + "/temp_dmdens_proj_"+ simtyp + "_" + outnum_char  +".png"
  imgname_z    =  simdir + "/dmdens_proj_" + simtyp + "_"+outnum_char +  ".png"

  ds = yt.load(outdir)
  Z  = ds.current_redshift
  p = yt.ProjectionPlot(ds, "z", ('deposit', 'DM_density'), weight_field=("deposit", "DM_density") ,width=(10.0,"Mpccm/h")) 
  p.set_cmap(('deposit', 'DM_density'), "cividis")
  p.save(imgname)

  image = Image.open(imgname)
  draw = ImageDraw.Draw(image)
  myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 42)
  text = f"Z:{Z:.2f}"
  text_color = (255, 0, 0)
  draw.text((730, 36), text, font=myfont, fill=(0, 0, 0))
  print("Saving plot: ", imgname_z)
  image.save(imgname_z)

def plot_proj_halo_gasdens(simdir, outnum, simtyp, numplot): 
  outnum_char      = str(outnum).zfill(5)
  outdir           =  simdir + "/output_" + outnum_char
  halocatalog_file = outdir + '/info_'+ outnum_char + '/info_' + outnum_char +'.0.h5'

  ds        = yt.load(outdir)
  halo_ds   = yt.load(halocatalog_file)

  ad        = ds.all_data()
  halo_ad   = halo_ds.all_data()

  Z         = ds.current_redshift

  #NumPart   = np.array( halo_ad['halos', 'particle_number'] )
  halo_mass = np.array( halo_ad['halos', 'particle_mass'].in_units("Msun") )
  Rvir      = np.array( halo_ad['halos', 'virial_radius'].in_units("Mpccm/h"))
  halo_posx = np.array( halo_ad['halos', 'particle_position_x'])
  halo_posy = np.array( halo_ad['halos', 'particle_position_y'] )
  halo_posz = np.array( halo_ad['halos', 'particle_position_z'])

  sink_posx = np.array( ad['sink', 'particle_position_x'].in_units("Mpccm/h"))
  sink_posy = np.array( ad['sink', 'particle_position_y'].in_units("Mpccm/h") )
  sink_posz = np.array( ad['sink', 'particle_position_z'].in_units("Mpccm/h"))


  for isink in range(0, len(sink_posx)) :
    sep = 0 
    for idim in range(0,3):
      sep = sep + (halo_posx[0] - sink_posx[isink]) * (halo_posx[0] - sink_posx[isink])
    sep = math.sqrt(sep) 
    if(sep < 5.0*Rvir[0]): 
      print("********", isink+1, sep*1000.0, Rvir[0]*1000.0)


  for iplot in range(0, numplot) :

    halocenter     = [halo_posx[iplot], halo_posy[iplot], halo_posz[iplot]]
    #Rvir[iplot] = 0.001
    #left_edge  = [halo_posx[iplot] - Rvir[iplot], halo_posy[iplot] - Rvir[iplot], halo_posz[iplot] - Rvir[iplot] ]
    #right_edge = [halo_posx[iplot] + Rvir[iplot], halo_posy[iplot] + Rvir[iplot], halo_posz[iplot] + Rvir[iplot] ]
    #left_edge = [halo_posx[iplot] - 0.1, halo_posy[iplot] - 0.1, halo_posz[iplot] -0.1]
    #right_edge = [halo_posx[iplot] + 0.1, halo_posy[iplot] + 0.1, halo_posz[iplot] +0.1]
    #center = ds.domain_center.copy()
    #cds        = ds.region(halocenter, left_edge, right_edge)
    #print(center)

    print(halocenter)
    print(Rvir[iplot])
    print(len(sink_posx))
    #print(Rvir[iplot])
    #print(halo_ad['halos', 'particle_position_x'][iplot])
    #print(halo_ad['halos', 'virial_radius'][iplot])
    #print(left_edge)
    #print(right_edge)
    #cds = ds.region(ds.domain_center, ds.domain_left_edge, ds.domain_right_edge)

    L = [1, 1, 0]  # vector normal to cutting plane
    north_vector = [-1, 1, 0]

    print("******* plotting halo of mass", halo_mass[iplot])

    imgname      =  simdir + "/temp_halo_gasdens_proj_"+ simtyp + "_" + outnum_char+ "_"+ str(iplot) +".png"
    imgname_z    =  simdir + "/halo_gasdens_proj_" + simtyp + "_"+outnum_char+ "_" + str(iplot) +".png"

    #p  = yt.ProjectionPlot(ds, L, ('gas', 'density'), weight_field=('gas', 'density'), width=(5.0*Rvir[iplot],"kpccm/h"),north_vector=north_vector)
    p  = yt.ProjectionPlot(ds, "z", ('gas', 'density'), weight_field=('gas', 'density'), center=halocenter, width=(20.0*Rvir[iplot], "Mpccm/h"))
    #p  = yt.ParticlePlot(ds, ("sink", "particle_position_x"), ("sink", "particle_position_y"), ("sink", "particle_mass"), center=halocenter, width=(5.0*Rvir[iplot], "Mpccm/h"))
    p.annotate_particles(width=(20.0*Rvir[iplot], "Mpccm/h"),ptype='sink', p_size=5.0)
    #p.set_unit(("sink", "particle_mass"), "Msun")
    #eps_fig = eps.single_plot(p)
    #eps_fig.circle(radius=0.2, loc=(0.5, 0.5))
    #eps_fig.save_fig("zoom-circle", format="eps")

    #p.set_width((10, "Mpccm/h"))
    #p.set_xlim(halo_posz[iplot] - Rvir[iplot], halo_posz[iplot] + Rvir[iplot])
    #p.set_ylim(halo_posz[iplot] - Rvir[iplot], halo_posz[iplot] + Rvir[iplot])
    #p.set_zlim(ad, halo_posz[iplot] - Rvir[iplot], halo_posz[iplot] + Rvir[iplot])

    #p.set_center((halo_posx[iplot], halo_posy[iplot]))
    p.save(imgname)

    image = Image.open(imgname)
    draw = ImageDraw.Draw(image)
    myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 42)
    text = f"Z:{Z:.2f}"
    text_color = (255, 0, 0)
    draw.text((730, 36), text, font=myfont, fill=(0, 0, 0))
    #width, height = image.size
    #radius = (width /Rvir[iplot] ) * Rvir[iplot]
    #imgcenter = (int(halocenter[2] * width), int(halocenter[0] * height))
    #draw.ellipse([(imgcenter[0] - radius, imgcenter[1] - radius),
    #          (imgcenter[0] + radius, imgcenter[1] + radius)],
    #        outline="black", width=2)
    print("Saving plot: ", imgname_z)
    image.save(imgname_z)

  #center    = [halo_posx[ihalo], halo_posy[ihalo], halo_posz[ihalo]]
  #sp        = ds.sphere(center, (40.0*Rvir[ihalo], "kpc"))
  #star_mass = np.array(sp[("star", "particle_mass")].in_units('Msun'))
  #sink_mass = np.array(sp[("sink", "particle_mass")].in_units('Msun'))

def write_halos(simdir, outnum):
  outnum_char       = str(outnum).zfill(5)
  outdir           = simdir+'/output_' + outnum_char
  halocatalog_file  = outdir + '/info_'+ outnum_char + '/info_' + outnum_char +'.0.h5'
  halo_ds           = yt.load(halocatalog_file)
  halo_ad           = halo_ds.all_data()

  #NumPart   = np.array( halo_ad['halos', 'particle_number'] )
  halo_mass = np.array( halo_ad['halos', 'particle_mass'].in_units("Msun") )
  Rvir      = np.array( halo_ad['halos', 'virial_radius'].in_units("kpccm/h"))
  halo_posx = np.array( halo_ad['halos', 'particle_position_x'])
  halo_posy = np.array( halo_ad['halos', 'particle_position_y'] )
  halo_posz = np.array( halo_ad['halos', 'particle_position_z'])
  NumPart   = np.array( halo_ad['halos', 'particle_number'] )
  nhalo     = len(halo_mass)

  outfile   = outdir + '/info_'+ outnum_char + "/halodata.csv"

  print("writing : ", nhalo, "halos in : ", outfile)
  with open(outfile, 'w') as file:
    for ihalo in range(0, nhalo):
      line = f"{ihalo+1}\t{halo_posx[ihalo]:.6f}\t{halo_posy[ihalo]:.6f}\t{halo_posz[ihalo]:.6f}\t{Rvir[ihalo]:.6f}\t{halo_mass[ihalo]:.6f}\t{NumPart[ihalo]:.6f}\n"
      file.write(line)
      if(ihalo < 10):
        print(ihalo+1, halo_posx[ihalo], halo_posy[ihalo], halo_posz[ihalo], Rvir[ihalo], halo_mass[ihalo],NumPart[ihalo])



