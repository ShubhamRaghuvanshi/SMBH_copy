import yt 
from yt.extensions.astro_analysis.halo_analysis import HaloCatalog
import yt.visualization.eps_writer as eps
import os
import numpy as np
import matplotlib.pyplot as plt
from halo_mass_function import *
from sfr_spectrum import * 
from PIL import Image, ImageDraw, ImageFont
from matplotlib.ticker import FuncFormatter

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

  if( not os.path.exists(halocatalog_file) or overwrite) :
    print("creating halo catalog: ", halocatalog_file)
    print("finder method : ", FinderMethods[MethId])
    print("finder kwargs : ", FinderKwargs[MethId])
    ds = yt.load(outdir)
    halo_catalog = HaloCatalog(data_ds=ds, finder_method=FinderMethods[MethId], finder_kwargs=FinderKwargs[MethId], output_dir = outdir)
    halo_catalog.create()
  else :
    print("halo catalog already exists :", halocatalog_file)

log_M_min = 7
log_M_max = 14
delta_log_M = 0.1
label_arr = ["PS", "JN", "ST", "WR", "TK" ]
cbar = ['m', 'r', 'g', 'b', 'y', 'm', 'c', 'k']

# compare hmf from with same outnum at same redshift
def compare_hmf(simdir_arr, simname_arr, cbar_arr, mbar_arr, outnum_arr, fitfunc_arr): 

  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel("Halo Mass [Msun]")
  plt.ylabel("Number of Halos")
  iout=0
  for simdir in simdir_arr :
    outnum           = outnum_arr[iout]
    outnum_char      = str(outnum).zfill(5)
    outdir           = simdir+'/output_' + outnum_char
    halocatalog_file = outdir + '/info_'+ outnum_char + '/info_' + outnum_char +'.0.h5'
    print("Analyzing : ", halocatalog_file)

    ds               = yt.load(halocatalog_file)
    Z                = ds.current_redshift
    ad               = ds.all_data()  
    NumPart          = np.array( ad['halos', 'particle_number'] )  
    HMass            = np.array( ad["halos", "particle_mass"])
    halo_masses      = []

    
    print("z :", f"{Z:.2f}")
    if (iout == 0 ):
      title = "DM Halo Mass function at" + f' z={Z:.2f}'
      plt.title(title)
      Zint  = round(Z)
      Zint2 = round(Z)

    if(Zint != round(Z)) :
      print("output reshift dont match :", outdir)
      raise ValueError("Compare Halo Mass Functions onlt at the same redshift")

    for ihalo in range (0,len(NumPart)):
      if( NumPart[ihalo] > 10  ):
        halo_masses.append(HMass[ihalo])   

    Bins        = 10**np.arange(log_M_min, log_M_max + delta_log_M, delta_log_M)
    hist, edges = np.histogram(halo_masses, bins=Bins)
    bin_centers = (edges[1:] + edges[:-1]) / 2
    bin_widths  = edges[1:] - edges[:-1]
    dn_dlogm    = hist
    err         = np.sqrt(hist)  
    labl        = simname_arr[iout]
    print("label :", labl)
    plt.errorbar(bin_centers, dn_dlogm, yerr=err, fmt=mbar_arr[iout], color=cbar_arr[iout],capsize=3, label=labl)
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
    ifunc = ifunc + 1
 
  plt.legend()
  plt.tight_layout()  
  imgname='hmf_cmp_z' + str(Zint2) + '.png' 
  print('Saving : ', imgname)
  plt.savefig(imgname)
  plt.clf()


#star functions 
def write_starfile(simdir, outnum, overwrite, radius, center):
  outnum_char      = str(outnum).zfill(5)
  outdir           = simdir+'/output_' + outnum_char
  sfr_file         = outdir + "/SFR_" +  outnum_char + ".out" 
  print("sfr_outdir : ", outdir)

  #if(overwrite) :
  #  command = f'rm {sfr_file}'
  #  print("executing", command)
  #  os.system(command)   

  if( not os.path.exists(sfr_file) or overwrite ) :     
    ds         = yt.load(outdir)
    ad         = ds.all_data()
    mass       = ad[("star", "particle_mass")].in_units('Msun')
    age        = ad[("star", "age")].in_units('Myr')
    ct         = ad[("star", "particle_birth_time")]
    Z        = ds.current_redshift
   
    threshold  = ds.quan(10.0, "Myr")
    mass_old   = mass[age > threshold]
    ct_old     = ct[age > threshold]
    if(radius < 0) :
      sp        = ds.all_data()
    else :
      sp        =   ds.sphere(center, (radius, "Mpccm/h"))

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


# compare sfr from with same outnum at same redshift
def compare_sfr(simdir_arr, simname_arr ,outnum_arr, cbar_arr, linestyle_arr ,plotanal, plotdata, fullbox ): 

  #plt.xscale('log')
  
  fig, ax1 = plt.subplots()
  ax1.set_yscale('log')
  #ax1.set_xscale('log')
  ax1.set_xlabel("Redshift(z)")
  ax1.set_ylabel(r"SFR  [$M_\odot$/yr/Mpcc] ")

  if(fullbox ) :
    ax1.set_title("global SFRD vs redshift")
  else :
    ax1.set_title("SFRD vs redshift for one halo")

  mass_thres = 0.0 
  tage_thres = 0.0
  iout=0
  nout = len(simdir_arr)
  for simdir in simdir_arr :
    outnum      = outnum_arr[iout]  
    outnum_char = str(outnum).zfill(5)
    outdir      = simdir+'/output_' + outnum_char
    sfr_file    = outdir + '/SFR_'+ outnum_char +'.out'
    
    Redshift    = []
    SFR         = []
    tuniv       = []  
    mstar       = []

    with open(sfr_file, 'r') as file:
      numline     = 0
      for line in file:
        numline = numline + 1
        if(numline > 1):
          tage = float(line.split()[0]) 
          z    = float(line.split()[2])
          mass = float(line.split()[5])
          if(tage > tage_thres and z > 8.0 and mass>0): 
            sfrcc = float(line.split()[4])
            aexp  =  1.0/(1.0+z)
            Redshift.append(z)
            SFR.append(sfrcc)
            tuniv.append(aexp*13.6)
            mstar.append(mass)
            #print(13.6/(1.0+z))
    ax1.plot(Redshift, SFR, linestyle=linestyle_arr[iout], color=cbar_arr[iout], label=simname_arr[iout])  
    #ax2 = ax1.secondary_xaxis('top', functions=(lambda x: 13.6*x/(1.0+x), lambda x: x/(13.6-x ) ))
    iout = iout + 1    

  if(plotanal and fullbox) :
    Z_anal   =   [8,10,12,14,16]
    SFR_anal =   []
    SFR_aniket_bursty=[-1.96, -2.33, -2.77,  ]
    for i in range(0, len(Z_anal)) :
      z   = Z_anal[i]
      sfr = 0.015 * (1 + z)**2.7 / (1 + ((1+z)/2.9)**5.6)
      SFR_anal.append(sfr)
    ax1.plot(Z_anal, SFR_anal, linestyle="-"  ,color='k', label="Madau, Dickinson")


  if(plotdata and fullbox) :
    #z_arr     = [8   ,  9   , 10.5 , 12   , 13.25, 16]
    #logsfrd   = [-2.3, -2.61, -2.75, -3.23, -3.53, -3.59]
    #err_p = [0, 0.18, 0, 0.29, 0, 0.33]
    #err_m = [0, 0.16, 0, 0.27, 0, 2.83]

    z_arr             = [5.91 , 7.00 , 7.80 , 7.90 , 8.01 , 9.00 , 10.00, 11.00, 12.01, 12.01, 13.75, 17.00 ]
    logsfrd           = [-1.65, -1.92, -1.07, -2.07, -2.19, -2.59, -2.85, -2.70, -3.37, -3.30, -3.77, -3.62 ]
    err_p             = [0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0]
    err_m             = [0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0] 
    logsfrd_ab_bursty = [-1.49, -1.78, -1.82, -1.83, -1.83 , -2.12, -2.29, -2.56, -2.79, -2.79, -3.17, -4.23 ]
    logsfrd_ab_smooth = [-0.86, -1.10, -1.50, -1.53, -1.49 , -1.87, -2.02, -2.24, -2.38, -2.38, -2.83, -3.90 ]

    sfrd_p    = []
    sfrd_m    = []
    sfrd      = []  
    sfrd_ab_bursty = []  
    sfrd_ab_smooth = []
    Z_data    = []

    for id in range(0,len(logsfrd)) :
      if(z_arr[id]> 7.95):
        x  = 10** (logsfrd[id])
        xp = 10** ( logsfrd[id] + err_p[id] )
        xm = 10** ( logsfrd[id] - err_m[id] )
        Z_data.append(z_arr[id])
        sfrd.append(x)
        sfrd_ab_bursty.append(10**(logsfrd_ab_bursty[id]))
        sfrd_ab_smooth.append(10**(logsfrd_ab_smooth[id]))
        sfrd_p.append(xp-x)
        sfrd_m.append(x-xm)
    ax1.errorbar(Z_data, sfrd, yerr=[sfrd_m, sfrd_p] ,fmt='o' ,color='grey', label="Donnan(2022)+Harikane(2023)")
    #ax1.plot(Z_data, sfrd_ab_bursty, marker='*', linestyle="--"  ,color='r', label="SPICE, bursty")
    ax1.plot(Z_data, sfrd_ab_smooth, marker='*', linestyle="--"  ,color='grey', label="SPICE, smooth")

  imgname = "sfr_cmp.png"
  plt.legend()
  print("Saving : ", imgname)
  plt.tight_layout()
  plt.savefig(imgname)
  plt.clf()

def plot_proj_dmdens(simdir, outnum, simtyp, radius, center, projdir, boxlen_comov): 
  outnum_char  = str(outnum).zfill(5)
  outdir       = simdir+'/output_' + outnum_char
  print("outdir :",outdir)
  imgname      = "temp_dmdens_proj_"+ simtyp + "_" + outnum_char  +".png"
  
  ds = yt.load(outdir)
  Z  = ds.current_redshift
  Zint2 = round(Z)
  imgname_z    =  "dmdens_proj_" + simtyp + "_z"+ str(Zint2) +  ".png"

  if(radius <= 0) :
    radius = boxlen_comov/2.0

  p  = yt.ProjectionPlot(ds, projdir, ("deposit", "DM_density"), weight_field=("deposit", "DM_density"),width=(2.0*radius,"Mpccm/h"))

  if(radius > 0) :
    if(projdir == 'z')   :
      p.set_center(    (center[0], center[1])  )
    elif(projdir == 'x') :
      p.set_center(    (center[1], center[2])  )
    elif(projdir == 'y') :
      p.set_center(    (center[2], center[0])  )
    else :
      raise ValueError("Invalid direction for projection")              
  p.set_cmap(("deposit", "DM_density"), "viridis")
  p.set_xlabel("cMpc/h")
  p.set_ylabel("cMpc/h")
  p.set_zlim( ("deposit", "DM_density"), 1e-26, 1e-21 )
  p.annotate_title(simtyp)
  p.save(imgname)
  
  image = Image.open(imgname)
  draw = ImageDraw.Draw(image)
  myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 42)
  text = f"Z:{Z:.2f}"
  text_color = (255, 0, 0)
  draw.text((730, 36), text, font=myfont, fill=(0, 0, 0))
  print("Saving plot: ", imgname_z)
  image.save(imgname_z)

def plot_proj_gasdens(simdir, outnum, simtyp, radius, center, projdir, boxlen_comov): 
  outnum_char  = str(outnum).zfill(5)
  outdir       = simdir+'/output_' + outnum_char
  print("outdir :",outdir)
  imgname      = "temp_gasdens_proj_"+ simtyp + "_" + outnum_char  +".png"
  
  ds = yt.load(outdir)
  Z  = ds.current_redshift
  Zint2 = round(Z)
  imgname_z    =  "gasdens_proj_" + simtyp + "_z"+ str(Zint2) +  ".png"

  if(radius <= 0) :
    radius = boxlen_comov/2.0

  p  = yt.ProjectionPlot(ds, projdir, ("gas", "density"), weight_field=("gas", "density"),width=(2.0*radius,"Mpccm/h"))

  if(radius > 0) :
    if(projdir == 'z')   :
      p.set_center(    (center[0], center[1])  )
    elif(projdir == 'x') :
      p.set_center(    (center[1], center[2])  )
    elif(projdir == 'y') :
      p.set_center(    (center[2], center[0])  )
    else :
      raise ValueError("Invalid direction for projection")              
  p.set_cmap(("gas", "density"), "viridis")
  p.set_xlabel("cMpc/h")
  p.set_ylabel("cMpc/h")
  p.set_zlim( ("gas", "density"), 1e-27, 1e-23 )
  p.annotate_title(simtyp)
  p.save(imgname)
  
  image = Image.open(imgname)
  draw = ImageDraw.Draw(image)
  myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 42)
  text = f"Z:{Z:.2f}"
  text_color = (255, 0, 0)
  draw.text((730, 36), text, font=myfont, fill=(0, 0, 0))
  print("Saving plot: ", imgname_z)
  image.save(imgname_z)

def plot_proj_gastemp(simdir, outnum, simtyp, radius, center, projdir, boxlen_comov): 
  outnum_char  = str(outnum).zfill(5)
  outdir       = simdir+'/output_' + outnum_char
  print("outdir :",outdir)
  imgname      = "temp_gastemp_proj_"+ simtyp + "_" + outnum_char  +".png"
  
  ds = yt.load(outdir)
  Z  = ds.current_redshift
  Zint2 = round(Z)
  imgname_z    =  "gastemp_proj_" + simtyp + "_z"+ str(Zint2) +  ".png"

  if(radius <= 0) :
    radius = boxlen_comov/2.0

  p  = yt.ProjectionPlot(ds, projdir, ("gas", "temperature"), weight_field=("gas", "density"),width=(2.0*radius,"Mpccm/h"))

  if(radius > 0) :
    if(projdir == 'z')   :
      p.set_center(    (center[0], center[1])  )
    elif(projdir == 'x') :
      p.set_center(    (center[1], center[2])  )
    elif(projdir == 'y') :
      p.set_center(    (center[2], center[0])  )
    else :
      raise ValueError("Invalid direction for projection")              
  p.set_cmap(("gas", "temperature"), "twilight_shifted")
  p.set_xlabel("cMpc/h")
  p.set_ylabel("cMpc/h")
  p.set_zlim( ("gas", "temperature"), 1e1, 1e6 )
  p.annotate_title(simtyp)
  p.save(imgname)
  
  image = Image.open(imgname)
  draw = ImageDraw.Draw(image)
  myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 42)
  text = f"Z:{Z:.2f}"
  text_color = (255, 0, 0)
  draw.text((730, 36), text, font=myfont, fill=(0, 0, 0))

  print("Saving plot: ", imgname_z)
  image.save(imgname_z)

def tn_phaseplot(simdir, outnum, simtyp, boxlen_comov, overwrite): 
  outnum_char  = str(outnum).zfill(5)
  outdir       = simdir+'/output_' + outnum_char
  print("outdir :",outdir)
  imgname      = "temp_tn_"+ simtyp + "_" + outnum_char  +".png"
  
  ds = yt.load(outdir)
  ad = ds.all_data()
  Z  = ds.current_redshift
  Zint2 = round(Z,2)
  imgname_z    =  "tn_" + simtyp + "_z"+ str(Zint2) +  ".png"

  #family = np.array(ad['all', 'particle_family'])
  #alldens = np.array(ad['ramses', 'Density'])
  #alltemp = np.array(ad['ramses', 'Pressure']) 

  #print("len : ", len(family),  len(alldens), len(alltemp))
  #for ip in range(0, len(family)) :
  #  if(family[ip] == 2.0) :
  #    print("family :", family[ip], alldens[ip], alltemp[ip])

  if( not os.path.exists(imgname_z) or overwrite) :
    plot = yt.PhasePlot(ad, (('gas', 'number_density')), ("gas", "temperature"), ("gas", "mass"))
    plot.set_unit(("gas", "mass"), "Msun")
    plot.set_unit(("gas", "number_density"), "1/cm**3")
    plot.set_ylim( 1, 1e9 )
    plot.set_xlim( 1e-7,1e4 )
    plot.set_zlim( ("gas", "mass"), 1e0, 1e10 )
    plot.annotate_title(simtyp)
    plot.save(imgname)
  
    image = Image.open(imgname)
    draw = ImageDraw.Draw(image)
    myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 42)
    text = f"Z:{Z:.2f}"
    text_color = (255, 0, 0)
    draw.text((730, 36), text, font=myfont, fill=(0, 0, 0))
    print("Saving plot: ", imgname_z)
    image.save(imgname_z)
  else :
    print(imgname_z, "already exists")  

def location_phasepoints(simdir, outnum, overwrite, T1, T2, n1, n2):

  file_name='tn_xyz.txt'
  if(overwrite or not os.path.exists(file_name)) :
    outnum_char  = str(outnum).zfill(5)
    outdir       = simdir+'/output_' + outnum_char  
    ds           = yt.load(outdir)

    my_sphere    =   ds.sphere([0.17, 0.13, 0.44], (0.4, "Mpccm/h"))
    #my_sphere    = ds.all_data()
    Z            = ds.current_redshift

    temp = my_sphere["gas", "temperature"]
    nden = my_sphere["gas", "number_density"].in_units("1/cm**3")
    posx = np.array(my_sphere["gas", "x"].in_units("Mpccm/h"))
    posy = np.array(my_sphere["gas", "y"].in_units("Mpccm/h"))
    posz = np.array(my_sphere["gas", "z"].in_units("Mpccm/h"))
    mass = np.array(my_sphere['gas', 'cell_mass'].in_units("Msun")  )
    print("proceesing gas cells....")

    iline=0
    gasmass=  0.0
    with open(file_name, 'w') as f:
      for itemp in range(0, len(temp)) :
        if(temp[itemp] > T1 and temp[itemp] < T2 and nden[itemp] > n1 and nden[itemp] < n2 ) :
          print("temp, dens :", temp[itemp], nden[itemp], math.log10(temp[itemp]), math.log10(nden[itemp]) )
          x = posx[itemp]
          y = posy[itemp]
          z = posz[itemp]
          gasmass = gasmass + mass[itemp]
          f.write(f"{x}\t{y}\t{z}\n")
          if(iline < 5) :
            print(x,y,z)
            iline = iline + 1 
    print("done... total gas mass : ", gasmass)   
  else :
    print(file_name, "exists, do you want to overwrite ? ")   



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
    p  = yt.ProjectionPlot(ds, "z", ('gas', 'density'), weight_field=('gas', 'density'), center=halocenter, width=(10.0*Rvir[iplot], "Mpccm/h"))
    #p  = yt.ParticlePlot(ds, ("sink", "particle_position_x"), ("sink", "particle_position_y"), ("sink", "particle_mass"), center=halocenter, width=(5.0*Rvir[iplot], "Mpccm/h"))
    p.annotate_particles(width=(10.0*Rvir[iplot], "Mpccm/h"),ptype='sink', p_size=5.0)
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
  outnum_char      = str(outnum).zfill(5)
  outdir           = simdir+'/output_' + outnum_char
  halocatalog_file = outdir + '/info_'+ outnum_char + '/info_' + outnum_char +'.0.h5'
  halo_ds          = yt.load(halocatalog_file)
  halo_ad          = halo_ds.all_data()
  Zcurr            = halo_ds.current_redshift

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
      line = f"{ihalo+1}\t{halo_posx[ihalo]:.6f}\t{halo_posy[ihalo]:.6f}\t{halo_posz[ihalo]:.6f}\t{Rvir[ihalo]:.6f}\t{halo_mass[ihalo]:.6f}\t{NumPart[ihalo]:.6f}\t{Zcurr:.4f}\n"
      file.write(line)
      #if(ihalo < 10):
      #  print(ihalo+1, halo_posx[ihalo], halo_posy[ihalo], halo_posz[ihalo], Rvir[ihalo], halo_mass[ihalo],NumPart[ihalo])

def volumefrac(simdir, outnum, lbox, t_frac):
  Mpc  = 3.0856e24
  h    = 0.67739997
  Mpch = Mpc/h  

  outnum_char = str(outnum).zfill(5)
  outdir      = simdir + "/output_" +  outnum_char
  if(os.path.exists(outdir)) :
    print("outdir : ", outnum_char)
    ds    = yt.load(outdir)
    ad    = ds.all_data()
    Z     = ds.current_redshift
    aexp  = 1.0/(1.0+Z)
    temp  = t_frac*np.array( ad[('gas', 'temperature')]  )
    cvol  = np.array( ad[('gas', 'cell_volume')] )
    ncell = len(temp)

    volfrac = 0.0 
    for icell in range(0, ncell) :
      celldx  = pow(cvol[icell], 1.0/3.0)/(Mpch*lbox*aexp)
      if( temp[icell] > 1e6 ) :
        volfrac = celldx*celldx*celldx + volfrac
      #print(- math.log2(celldx), volfrac)
    print(outnum_char, Z, volfrac)
    return[Z, volfrac]
  else :
    raise ValueError("My error: File" , outdir ,"does not exist")








