import os 
import struct
import math 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import yt_funcs
import yt
from scipy import stats
from PIL import Image, ImageDraw, ImageFont
plt.rcParams['text.usetex'] = True
plt.rc('text', usetex=True)

#constants 
pi = 3.14159265
Kilo = 1e3
million = 1e6
Mpc = 3.0856e24 #cm
pc = 3.0856e18
Kpc = 3.0856e21
Msun = 1.988e33 #g
yr = 31536000.0 #sec
Myr = yr*million
G_cgs = 6.674e-8
mH = 1.6735575e-24
kB = 1.3806490e-16
mP = 1.6726219e-24
X = 0.76
G_cgs = 6.674e-8

h = 0.6773999
Omega_m = 0.30889999
Omega_b = 0.04488
Omega_r = 0.0
Omega_l = 0.0
Omega_k = 0.0
H0 = 100.0 # h(km/s)/Mpc

H0_cgs = H0 * h * 1e5 / Mpc
H0_IGyr = H0_cgs * 1e9* yr

#cgs  units
unit_M            = Msun
unit_L            = Mpc
unit_T            = 1.0
unit_vel          = unit_L/unit_T
unit_E            = unit_M*unit_vel*unit_vel
unit_rho          = unit_M/(unit_L*unit_L*unit_L)

#
G                 = G_cgs*unit_M*unit_T*unit_T/( unit_L*unit_L*unit_L )
H0_sys            = (H0*Kilo/Mpc)*unit_T*100.0
rho_crit0         = (3.0*H0_sys*H0_sys)/(8.0*pi*G)

#params

mass_halo_AGN     = 1e9  
mass_clump_AGN    = 1e9  
n_star            = 0.1

status=0

#code units 
Z         = -1.0
scale_M   = -1.0
scale_L   = -1.0
scale_T   = -1.0
scale_d   = -1.0
scale_nH  = -1.0
scale_vel = -1.0
scale_T2  = -1.0
simtime   = -1.0 
def give_units(aexp, boxlength):
  global Z
  global scale_M
  global scale_L
  global scale_T
  global scale_d
  global scale_nH
  global scale_vel
  global scale_T2

  Lx          = boxlength # h^-1 mpc
  Ly          = Lx
  Lz          = Lx
  V           = Lx*Ly*Lz
  M_V         = rho_crit0*Omega_m*V

  scale_L     = (aexp*boxlength/h)*unit_L
  scale_M     = (M_V/h)*unit_M
  scale_d     = (rho_crit0*Omega_m*h*h/( aexp*aexp*aexp )) *unit_rho
  scale_T     = aexp*aexp/(h*H0*1e5/Mpc)  
                # aexp**2 / (h0*1d5/Mpc2cm)

  scale_nH    = X*scale_d/mH
  #scale_vel   = (scale_L/scale_T)*unit_vel
  scale_vel   = (scale_L/scale_T)
  scale_T2    = mH/kB * scale_vel*scale_vel


def read_units(simdir, outnum):
  global Z
  global scale_M
  global scale_L
  global scale_T
  global scale_d
  global scale_nH
  global scale_vel
  global scale_T2
  global simtime

  outnumn_char = str(outnum).zfill(5)
  outdir = simdir+'/output_' + outnumn_char
  infofile = outdir+'/info_'+outnumn_char+'.txt' 
  if(os.path.exists(infofile)) :
    aexp    = None
    scale_L = None 
    scale_d = None
    scale_T = None 
    simtime = None
    with open(infofile , 'r') as file:
      for line in file:
        if 'time' in line:
          simtime = float(line.split()[2])
        if 'aexp' in line:
          aexp = float(line.split()[2])
        if 'unit_l' in line:
          scale_L = float(line.split()[2])
        if 'unit_d' in line:
          scale_d = float(line.split()[2])
        if 'unit_t' in line:
          scale_T = float(line.split()[2])
          break
    if(aexp == None or scale_L == None or scale_d == None or scale_T == None or simtime == None) :
      raise ValueError(' units not read ')
    Z           = 1.0/aexp - 1.0
    scale_nH    = X/mH * scale_d
    scale_M     = scale_d*scale_L*scale_L*scale_L
    scale_vel   = scale_L/scale_T
    scale_T2    = mH/kB * scale_vel*scale_vel
  else : 
    print("read_unit says: ", outdir, "doesn't exist")

def merge_clump_props(simdir, outnum, NTask ) :
  outnumn_char = str(outnum).zfill(5)
  outdir=simdir+"/output_"+outnumn_char
  merged_file=outdir+"/clump_"+outnumn_char+".txt" 

  if( not os.path.exists(merged_file) ) :
    with open(merged_file, 'w') as output:  
      for thistask in range(1, NTask+1):
        thistask_char = str(thistask).zfill(5)
        filename = outdir + "/clump_"+ outnumn_char +  ".txt"  + thistask_char
        linenum=0
        with open(filename , 'r') as datafile:
          for line in datafile:
            if(thistask==1 and linenum ==0):
              output.write(line )  
            linenum = linenum+1
            if(linenum > 1):
              output.write(line )
  else :
    print("merged clump file already exists")

def merge_halo_props(simdir, outnum, NTask ) :
  outnumn_char = str(outnum).zfill(5)
  outdir=simdir+"/output_"+outnumn_char
  merged_file=outdir+"/halo_"+outnumn_char+".txt" 

  if( not os.path.exists(merged_file) ) :
    with open(merged_file, 'w') as output:  
      for thistask in range(1, NTask+1):
        thistask_char = str(thistask).zfill(5)
        filename = outdir + "/halo_"+ outnumn_char +  ".txt"  + thistask_char
        linenum=0
        with open(filename , 'r') as datafile:
          for line in datafile:
            if(thistask==1 and linenum ==0):
              output.write(line )  
            linenum = linenum+1
            if(linenum > 1):
              output.write(line )
  else :
    print("merged halo file already exists")


def read_clump_prop(simdir, outnum, propnum):
  outnum_char = str(outnum).zfill(5)
  outdir = simdir+'/output_' + outnum_char
  outfile = outdir+ "/clump_"+ outnum_char + ".txt"
  cmass_arr = []
  with open(outfile, 'r') as clumpfile:
    numline=0
    for line in clumpfile:
      numline = numline + 1
      if(numline > 1):
        cmass = float(line.split()[propnum])      
        cmass_arr.append(cmass) 
  return cmass_arr

def read_halo_prop(simdir, outnum, propnum):
  outnum_char = str(outnum).zfill(5)
  outdir = simdir+'/output_' + outnum_char
  outfile = outdir+ "/halo_"+ outnum_char + ".txt"
  hmass_arr = []
  with open(outfile, 'r') as clumpfile:
    numline=0
    for line in clumpfile:
      numline = numline + 1
      if(numline > 1):
        hmass = float(line.split()[propnum])      
        hmass_arr.append(hmass) 
  return hmass_arr

def read_sink_prop(simdir, outnum, propnum):
  outnum_char = str(outnum).zfill(5)
  outdir = simdir+'/output_' + outnum_char
  outfile = outdir+ "/sink_"+ outnum_char + ".csv"
  hmass_arr = []
  with open(outfile, 'r') as clumpfile:
    numline=0
    for line in clumpfile:
      numline = numline + 1
      if(numline > 2):
        hmass = float(line.split(',')[propnum])      
        hmass_arr.append(hmass) 
  return hmass_arr

def sink_prop_return(simdir, sinknum, propnum):
  global status
  outdir=simdir+"/sinkdata" 
  sinknum_char = str(sinknum).zfill(5)
  outfile=outdir+"/sink_"+sinknum_char+".csv"
  prop_arr = []
  if( os.path.exists(outfile) ) :
    status=1
    #print(outfile) 
    with open(outfile, 'r') as sinkfile:
      for line in sinkfile:
        propval = float(line.split(',')[propnum])      
        prop_arr.append(propval) 
  else :
    print("warning : ",outfile," doesn't exist, using default value of property")    
    status=0  
  return prop_arr


def give_massive_sinkid(simdir):

  sinkdir   = simdir+'/sinkdata/'
  csv_files = [file for file in os.listdir(sinkdir) if file.endswith('.csv')]
  Nfiles = len(csv_files)

  Nsink = Nfiles
  massive_id = 1
  msinkmax = 0
  status =0
  for isink in range (1, Nsink+1) :
    sinkmass   = sink_prop_return(simdir, isink, 1)
    Ntime      = len(sinkmass)
    if(sinkmass[Ntime -1 ] > msinkmax ):
      msinkmax = sinkmass[Ntime -1 ]
      massive_id = isink  

  print("massive sink : ", massive_id)
  return massive_id

cbar = ['r', 'g', 'b', 'k']
mbar = ["o", "^", "s", "h", "D" ]

def compare_sinkmasses(simdir_arr, label_arr, sinkid_arr) :
  cbar = ['r', 'g', 'b', 'y', 'm', 'c', 'k']
  plt.xlabel('Redshift')
  plt.ylabel('Msink[Msun]')
  plt.title('Total Sink Mass')
  plt.yscale('log')

  isink =0
  imgname = "sinkmass" 
  for simdir in simdir_arr :
    sinkid   = sinkid_arr[isink]
    plabel   = label_arr[isink]
    sinkmass = sink_prop_return(simdir, sinkid, 1)
    aexp_arr = sink_prop_return(simdir, sinkid, 10)

    Ntime        = len(aexp_arr)
    Z_arr        = []

    for itime in range(0, Ntime):
      aexp = aexp_arr[itime]
      z    = 1.0/aexp -1.0
      give_units(aexp)
      sinkmass[itime]   = (scale_M/Msun) * sinkmass[itime]
      Z_arr.append(z)
    plt.plot(Z_arr, sinkmass,  label=plabel, color=cbar[isink] ,marker='o', linestyle='-')
    imgname = imgname + "_" + plabel
    isink   = isink + 1
  plt.legend(loc='best')
  plt.gca().invert_xaxis()
  plt.tight_layout()
  imgname = imgname + ".png"
  print("Saving :" , imgname)
  plt.savefig(imgname)


def compare_sinkmasstot(simdir_arr, label_arr, cbar_arr, linestyle_arr, boxlen_comov) :

  fig1, ax1 = plt.subplots()
  fig2, ax2 = plt.subplots()
  fig3, ax3 = plt.subplots()
  fig4, ax4 = plt.subplots()
  fig5, ax5 = plt.subplots()
  fig6, ax6 = plt.subplots()
  fig7, ax7 = plt.subplots()
  fig8, ax8 = plt.subplots()
  fig9, ax9 = plt.subplots()

  isim =0
  for simdir in simdir_arr :
    print("simdir :", simdir)
    sinkid     = 1
    sinkid_mostmassvie = 2
    msinkfirst = sink_prop_return(simdir, sinkid, 1)
    msink_mostmassive = sink_prop_return(simdir, sinkid_mostmassvie, 1)
    msinktot   = sink_prop_return(simdir, sinkid, 11)
    mstartot   = sink_prop_return(simdir, sinkid, 12)
    Nsinktot   = sink_prop_return(simdir, sinkid, 0)
    aexp_arr   = sink_prop_return(simdir, sinkid, 10)
    sink_dmBHdt = sink_prop_return(simdir, sinkid, 9)

    if(status) :
      Ntime               = len(aexp_arr)
      Z_arr               = []
      msinkfirst_clipped  = []
      msink_mostmassive_clipped = []
      msinktot_clipped    = []
      mstartot_clipped    = []
      msinktomstar        = []
      Nsinktot_clipped    = []
      nsinktot_clipped    = []
      rhosinktot_clipped  = []
      sink_dmBHdt_clipped = []
      for itime in range(0, Ntime):
        aexp = aexp_arr[itime]
        z    = 1.0/aexp -1.0
        if(z > 8.0) :
          vol = pow ( boxlen_comov*aexp, 3 )
          give_units(aexp,boxlen_comov)
          #print(itime, scale_M/Msun)
          msinkfirst_clipped.append (  (scale_M/Msun) * msinkfirst[itime]    )
          msink_mostmassive_clipped.append (  (scale_M/Msun) * msink_mostmassive[itime]    )
          msinktot_clipped.append   (  (scale_M/Msun) * msinktot[itime]      )
          mstartot_clipped.append   (  (scale_M/Msun) * mstartot[itime]      )
          msinktomstar.append       (  msinktot[itime]/mstartot[itime]       )
          Nsinktot_clipped.append   (          Nsinktot[itime]               )
          nsinktot_clipped.append   (        Nsinktot[itime]/vol             )
          rhosinktot_clipped.append ( (scale_M/Msun) * (msinktot[itime]/vol) )
          sink_dmBHdt_clipped.append(  (scale_M/Msun) *( Myr/scale_T/1e6 ) * sink_dmBHdt[itime]   )
          Z_arr.append(z)
      ax1.plot(Z_arr, msinktot_clipped         , label=label_arr[isim], color=cbar_arr[isim], linestyle=linestyle_arr[isim])
      ax2.plot(Z_arr, mstartot_clipped         , label=label_arr[isim], color=cbar_arr[isim], linestyle=linestyle_arr[isim])
      ax3.plot(Z_arr, msinktomstar             , label=label_arr[isim], color=cbar_arr[isim], linestyle=linestyle_arr[isim])
      ax4.plot(Z_arr, Nsinktot_clipped         , label=label_arr[isim], color=cbar_arr[isim], linestyle=linestyle_arr[isim])
      ax5.plot(Z_arr, nsinktot_clipped         , label=label_arr[isim], color=cbar_arr[isim], linestyle=linestyle_arr[isim])
      ax6.plot(Z_arr, rhosinktot_clipped       , label=label_arr[isim], color=cbar_arr[isim], linestyle=linestyle_arr[isim])
      ax7.plot(Z_arr, msinkfirst_clipped       , label=label_arr[isim], color=cbar_arr[isim], linestyle=linestyle_arr[isim])
      ax9.plot(Z_arr, msink_mostmassive_clipped, label=label_arr[isim], color=cbar_arr[isim], linestyle=linestyle_arr[isim])
      ax8.plot(Z_arr, sink_dmBHdt_clipped      , label=label_arr[isim], color=cbar_arr[isim], linestyle=linestyle_arr[isim])
      isim   = isim + 1
    else :
      print("simdir :", simdir)
      raise ValueError('no sink in simulation box')

  imgname_1 = "msinktot_cmp.png"
  imgname_2 = "mstartot_cmp.png"
  imgname_3 = "msinktomstar_cmp.png"
  imgname_4 = "Nsinktot_cmp.png"
  imgname_5 = "nsinktot_cmp.png"
  imgname_6 = "rhosinktot_cmp.png"
  imgname_7 = "msinkfirst_cmp.png"
  imgname_8 = "dmdtsinkfirst_cmp.png"
  imgname_9 = "msink_mostmassive_cmp.png"

  ax1.set_xlabel("Redshift")  
  ax1.set_ylabel("Msinktot[Msun]")
  ax1.set_title("Total sink mass in the simulation box")
  ax1.set_yscale('log')
  ax1.legend(loc='best')
  fig1.gca().invert_xaxis()
  fig1.tight_layout()

  ax2.set_xlabel("Redshift")  
  ax2.set_ylabel("Mstartot[Msun]")
  ax2.set_title("Total star mass in the simulation box")
  ax2.set_yscale('log')
  ax2.legend(loc='best')
  fig2.gca().invert_xaxis()
  fig2.tight_layout()

  ax3.set_xlabel("Redshift")  
  ax3.set_ylabel("Msinktot/Mstartot")
  ax3.set_title("Ratio of total sink mass to total star mass in simulation box")
  ax3.set_yscale('log')
  ax3.legend(loc='best')
  fig3.gca().invert_xaxis()
  fig3.tight_layout()

  ax4.set_xlabel("Redshift")  
  ax4.set_ylabel("Nsinktot")
  ax4.set_title("Number of sinks in simulation box")
  ax4.set_yscale('log')
  ax4.legend(loc='best')
  fig4.gca().invert_xaxis()
  fig4.tight_layout()

  ax5.set_xlabel("Redshift")  
  ax5.set_ylabel("nsink[per cubic Mpc/h]")
  ax5.set_title("Number density of sinks in simulation box")
  ax5.set_yscale('log')
  ax5.legend(loc='best')
  fig5.gca().invert_xaxis()
  fig5.tight_layout()

  ax6.set_xlabel("Redshift")  
  ax6.set_ylabel("rhosink[per cubic Mpc/h]")
  ax6.set_title("Mass density of sinks in simulation box")
  ax6.set_yscale('log')
  ax6.legend(loc='best')
  fig6.gca().invert_xaxis()
  fig6.tight_layout()

  ax7.set_xlabel("Redshift")  
  ax7.set_ylabel("Msink_first[Msun]")
  ax7.set_title("Mass of first sink")
  ax7.set_yscale('log')
  ax7.legend(loc='best')
  fig7.gca().invert_xaxis()
  fig7.tight_layout()

  ax8.set_xlabel("Redshift")  
  ax8.set_ylabel("dMsink/dt[Msun/yr]")
  ax8.set_title("dmdt of first sink")
  ax8.set_yscale('log')
  ax8.legend(loc='best')
  fig8.gca().invert_xaxis()
  fig8.tight_layout()

  ax9.set_xlabel("Redshift")  
  ax9.set_ylabel("Msink_first[Msun]")
  ax9.set_title("Mass of most massive sink")
  ax9.set_yscale('log')
  ax9.legend(loc='best')
  fig9.gca().invert_xaxis()
  fig9.tight_layout()


  print("Saving : ", imgname_1)
  fig1.savefig(imgname_1)
  fig1.clf()

  print("Saving : ", imgname_2)
  fig2.savefig(imgname_2)
  fig2.clf()
 
  print("Saving : ", imgname_3)
  fig3.savefig(imgname_3)
  fig3.clf()

  print("Saving : ", imgname_4)
  fig4.savefig(imgname_4)
  fig4.clf()
 
  print("Saving : ", imgname_5)
  fig5.savefig(imgname_5)
  fig5.clf()
 
  print("Saving : ", imgname_6)
  fig6.savefig(imgname_6)
  fig6.clf()

  print("Saving : ", imgname_7)
  fig7.savefig(imgname_7)
  fig7.clf()

  print("Saving : ", imgname_8)
  fig8.savefig(imgname_8)
  fig8.clf()

  print("Saving : ", imgname_9)
  fig9.savefig(imgname_9)
  fig9.clf()


def give_outnum_from_aexp(simdir, aexp_out) :

  Z_out = 1.0/aexp_out - 1.0
  outnum_dzmin = 1
  dzmin = 9999.0
  z_prev = -1
  for outnum in range(1, 20) :
    outnum_char = str(outnum).zfill(5)
    outdir      = simdir + "/output_" + str(outnum_char)
    if(os.path.exists(outdir)) :
      read_units(simdir, outnum)
      #print("z_out :", Z_out, Z, abs(Z - Z_out ), outnum)
      if( abs(Z - Z_out) < 0.03 ) :
        print("found exact outnum for the given redshift : ", outnum)
        outnum_dzmin = outnum
        break
      if(dzmin > abs(Z - Z_out)) :
        dzmin = abs(Z - Z_out)
        outnum_dzmin = outnum
        z_prev = Z
    else :
      print(outdir, "is the last snapshot")
      break    
  return [outnum_dzmin, z_prev]


def plot_sink_pos_vs_z(simdir_arr, cbar_arr, label_arr, sinkid, boxlen_comov, relhalo, halonum_max) :

  fig1, ax1 = plt.subplots()
  fig2, ax2 = plt.subplots()
  fig3, ax3 = plt.subplots()

  isim=0
  for simdir in simdir_arr :
    xsink       = sink_prop_return(simdir, sinkid, 2)
    ysink       = sink_prop_return(simdir, sinkid, 3)
    zsink       = sink_prop_return(simdir, sinkid, 4)
    msink       = sink_prop_return(simdir, sinkid, 1)
    sink_dmBHdt = sink_prop_return(simdir, sinkid, 9)
    aexp_arr    = sink_prop_return(simdir, sinkid, 10)    

    if(relhalo) :
      OZ        = give_outnum_from_aexp(simdir, aexp_arr[0])
      outnum    = OZ[0]
      print("outnum with closest redshift to spawn of sink :", outnum)
      outnum_char    = str(outnum).zfill(5)
      outdir        = simdir+'/output_' + outnum_char
      halofile       = outdir + '/info_'+ outnum_char + "/halodata.csv"

      if  (not os.path.exists(halofile) ) :
        print("Halo catalog not found, creating one.......................") 
        yt_funcs.create_halo_catalog(simdir, outnum, 0, 0)
        yt_funcs.write_halos(simdir, outnum)
      else :
        print(halofile, "already exists")  
      ##
      ihalo = 0
      drmin = 99999
      halonum_drmin = 0
      with open(halofile, 'r') as file:
        for line in file:
          if(ihalo < halonum_max-1) :
            dx = xsink[0] - float(line.split()[1])
            dy = ysink[0] - float(line.split()[2])
            dz = zsink[0] - float(line.split()[3])
            dr2  = dx*dx + dy*dy + dz*dz
            if(drmin > dr2) :
              drmin         = dr2
              halonum_drmin = ihalo
              halo_posx     = float(line.split()[1])
              halo_posy     = float(line.split()[2])
              halo_posz     = float(line.split()[3])  
              halo_rvir     = float(line.split()[4])
              halo_mass     = float(line.split()[5])
              rvir          = (halo_rvir/1000.0) /boxlen_comov
          else :
            break    
          ihalo = ihalo + 1
      print("nearest halo center :", math.sqrt(dr2), math.sqrt(dr2)/rvir,  halo_mass/1e9, ihalo)    
    else : 
      halo_posx = 0.0
      halo_posy = 0.0
      halo_posz = 0.0  
      halo_rvir = -1
      rvir      = 1
      halo_mass = -1

    dr_arr          = []
    z_arr           = []
    msink_arr       = []
    sink_dmBHdt_arr = []
    itime    = 0
    for aexp in aexp_arr :
      give_units(aexp, boxlen_comov)     
      z = 1.0/aexp - 1.0
      if(z > 7.0) :
        dx    = xsink[itime] - halo_posx
        dy    = ysink[itime] - halo_posy
        dz    = zsink[itime] - halo_posz
        dr = math.sqrt( dx*dx + dy*dy + dz*dz )
        dr_arr.append(dr/rvir )
        msink_arr.append(msink[itime] * (scale_M/Msun) * 0.677)
        sink_dmBHdt_arr.append( (scale_M/Msun)*(Myr/scale_T)*sink_dmBHdt[itime] )
        z_arr.append(z) 
      itime = itime + 1  

    ax1.plot(z_arr, dr_arr         , label=label_arr[isim], color=cbar_arr[isim], linestyle='-')
    ax2.plot(z_arr, msink_arr      , label=label_arr[isim], color=cbar_arr[isim], linestyle='-')
    ax3.plot(z_arr, sink_dmBHdt_arr, label=label_arr[isim], color=cbar_arr[isim], linestyle='-')  
    isim = isim + 1 

  ax1.set_xlabel("redshift")
  if(relhalo) :
    ax1.set_ylabel(r"$ \left| (\vec{r}_{sink} - \vec{r}_{halo}) \right| / Rvir_{halo} $") 
    ax1.set_title(f"sink position relative to its nearest halo center, sinkid : {sinkid}") 
  else :
    ax1.set_ylabel("sink position in code units")
    ax1.set_title(f"sink position in code units, sinkid : {sinkid}")
  ax1.legend()
  fig1.tight_layout()  
  fig1.gca().invert_xaxis()
  figname = "sink_dr_" + str(sinkid) + ".png"
  print("saving :", figname)
  fig1.savefig(figname)

  ax2.set_xlabel("redshift")
  ax2.set_ylabel(" Msink [Msun] ") 
  ax2.set_title(f"sink mass, sinkid : {sinkid}") 
  ax2.legend()
  fig2.tight_layout()  
  fig2.gca().invert_xaxis()
  figname = "sink_mass_" + str(sinkid) + ".png"
  print("saving :", figname)
  fig2.savefig(figname)

  ax3.set_xlabel("Redshift")  
  ax3.set_ylabel("dMsink/dt[Msun/yr]")
  ax3.set_title(f"sink mass accretion rate, sinkid : {sinkid}")
  ax3.set_yscale('log')
  ax3.legend(loc='best')
  fig3.gca().invert_xaxis()
  fig3.tight_layout()
  figname = "sink_dmdt_" + str(sinkid) + ".png"
  print("saving :", figname)
  fig3.savefig(figname)

def give_nearest_ihalo(simdir, outnum, posx, posy, posz, boxlen_comov, halonum_max):
  outnum_char    = str(outnum).zfill(5)
  outdir         = simdir+'/output_' + outnum_char
  halofile       = outdir + '/info_'+ outnum_char + "/halodata.csv"

  if  (not os.path.exists(halofile) ) :
    print("Halo catalog not found, creating one.......................") 
    yt_funcs.create_halo_catalog(simdir, outnum, 0, 0)
    yt_funcs.write_halos(simdir, outnum)
  #else :
    #print(halofile, "already exists")  
  ##
  ihalo = 0
  drmin = 99999
  ihalo_drmin = 0
  with open(halofile, 'r') as file:
    for line in file:
      if(ihalo < halonum_max-1) :
        dx = posx - float(line.split()[1])
        dy = posy - float(line.split()[2])
        dz = posz - float(line.split()[3])
        dr2  = dx*dx + dy*dy + dz*dz
        if(drmin > dr2) :
          drmin         = dr2
          halo_posx   = float(line.split()[1])
          halo_posy   = float(line.split()[2])
          halo_posz   = float(line.split()[3])  
          halo_rvir   = float(line.split()[4])/1000.0/boxlen_comov
          halo_mass   = float(line.split()[5])
          ihalo_drmin = ihalo
      else :
        break    
      ihalo = ihalo + 1
  print("nearest halo center :", ihalo_drmin, ihalo, math.sqrt(drmin)/halo_rvir)
  return [halo_posx, halo_posy, halo_posz, halo_rvir, halo_mass, ihalo_drmin]

def give_haloprops(simdir, outnum, halonum, boxlen_comov) :
  outnum_char = str(outnum).zfill(5)
  outdir      = simdir+'/output_' + outnum_char
  halofile    = outdir + '/info_'+ outnum_char + "/halodata.csv"
 
  if(os.path.exists(halofile)) :
    ihalo = 0
    with open(halofile, 'r') as file:
      for line in file:
        if(ihalo < halonum-1) :
          ihalo = ihalo + 1
        else :
          break
    halo_posx   = float(line.split()[1])
    halo_posy   = float(line.split()[2])
    halo_posz   = float(line.split()[3])  
    halo_rvir   = float(line.split()[4])
    halo_mass   = float(line.split()[5])
    redshift    = float(line.split()[7])
    rvir        = halo_rvir/1000.0/boxlen_comov #Mpccm/h
    return [halo_posx, halo_posy, halo_posz, rvir, halo_mass]
  else :
    print("give_haloprops says :", halofile, "doesn't exist")
    return [-1,-1,-1,-1,-1]
  
def plot_sink_relpos_vs_z(simdir_arr, cbar_arr, label_arr, sinkid, boxlen_comov, halonum_max, outnum_max) :

  fig1, ax1 = plt.subplots()
  fig2, ax2 = plt.subplots()
  fig3, ax3 = plt.subplots()

  isim=0
  for simdir in simdir_arr :
    xsink         = sink_prop_return(simdir, sinkid, 2)
    ysink         = sink_prop_return(simdir, sinkid, 3)
    zsink         = sink_prop_return(simdir, sinkid, 4)
    msink         = sink_prop_return(simdir, sinkid, 1)
    sink_dmBHdt   = sink_prop_return(simdir, sinkid, 9)
    aexp_arr      = sink_prop_return(simdir, sinkid, 10)    

    ##
    OZ        = give_outnum_from_aexp(simdir, aexp_arr[0])
    outnum    = OZ[0]
    z_prev    = OZ[1]
    read_units(simdir, outnum+1)
    z_next = Z
    print("outnum with closest redshift to spawn of sink :", 1.0/aexp_arr[0] - 1.0, z_prev, z_next, outnum)
    halo_attrb = give_nearest_ihalo(simdir, outnum, xsink[0], ysink[0], zsink[0], boxlen_comov, halonum_max)
    #halo_attrb    = give_haloprops(simdir, outnum, nearest_ihalo, boxlen_comov)
    halo_posx     = halo_attrb[0]
    halo_posy     = halo_attrb[1]
    halo_posz     = halo_attrb[2]
    rvir          = halo_attrb[3]
    ##

    dr_arr          = []
    z_arr           = []
    msink_arr       = []
    sink_dmBHdt_arr = []
    z_snap_arr      = []
    dr_snap_arr     = []
    itime    = 0

    if(outnum_max[isim] > 0) :
      z_snap_arr.append( 1.0/aexp_arr[0] - 1.0)
      dx_snap  = xsink[0] - halo_posx
      dy_snap  = ysink[0] - halo_posy
      dz_snap  = zsink[0] - halo_posz  
      dr_snap  = dx_snap*dx_snap + dy_snap*dy_snap + dz_snap*dz_snap
      dr_snap_arr.append(math.sqrt(dr_snap)/rvir)
    
    for aexp in aexp_arr :
      give_units(aexp, boxlen_comov)     
      z = 1.0/aexp - 1.0
      if(z > 7.0) :
        if(z- z_next  < z_prev -z and outnum < outnum_max[isim]) :
          z_snap_arr.append(z)
          z_prev = z_next 
          outnum = outnum + 1
          read_units(simdir, outnum+1)
          z_next = Z
          print("new halo center : ", z_prev, z, z_next, z_next-z, z-z_prev, outnum)
          halo_attrb = give_nearest_ihalo(simdir, outnum, xsink[itime], ysink[itime], zsink[itime], boxlen_comov, halonum_max)
          #halo_attrb    = give_haloprops(simdir, outnum, nearest_ihalo, boxlen_comov)
          halo_posx     = halo_attrb[0]
          halo_posy     = halo_attrb[1]
          halo_posz     = halo_attrb[2]
          rvir          = halo_attrb[3]
          dx_snap  = xsink[itime] - halo_posx
          dy_snap  = ysink[itime] - halo_posy
          dz_snap  = zsink[itime] - halo_posz  
          dr_snap  = dx_snap*dx_snap + dy_snap*dy_snap + dz_snap*dz_snap
          dr_snap_arr.append(math.sqrt(dr_snap)/rvir)
          
        dx    = xsink[itime] - halo_posx
        dy    = ysink[itime] - halo_posy
        dz    = zsink[itime] - halo_posz
        dr = math.sqrt( dx*dx + dy*dy + dz*dz )
        dr_arr.append(dr/rvir )
        give_units(aexp, boxlen_comov)
        #print("dr :",scale_M, Msun)
        msink_arr.append(msink[itime] * scale_M/Msun)
        sink_dmBHdt_arr.append( sink_dmBHdt[itime] *(scale_M/Msun)*(Myr/scale_T) )
        z_arr.append(z) 
      itime = itime + 1  

    ax1.plot(z_arr, dr_arr          , label=label_arr[isim], color=cbar_arr[isim], linestyle='-')
    if(outnum_max[isim] > 0) :
      ax1.scatter(z_snap_arr, dr_snap_arr, color=cbar_arr[isim], marker='o')
    ax2.plot(z_arr, msink_arr       , label=label_arr[isim], color=cbar_arr[isim], linestyle='-')
    ax3.plot(z_arr, sink_dmBHdt_arr , label=label_arr[isim], color=cbar_arr[isim], linestyle='-')  
    isim = isim + 1 

  ax1.set_xlabel("redshift")
  if(outnum_max[0] > 0) :
    ax1.set_ylabel(r"$ \left| (\vec{R}_{sink}(z) - \vec{R}_{halo}(z)) \right| / R_{vir}(z) $") 
    ax1.set_title(f" (R_sink(z) - R_halo(z))/R_vir(z) ..... [for sinkid : {sinkid}]")
  else :
    ax1.set_ylabel(r"$ \left| (\vec{R}_{sink}(z) - \vec{R}_{halo, initial}) \right| / R_{vir, initial} $") 
    ax1.set_title(f"(R_sink(z) - R_halo_initial)/R_vir_initial ..... [for sinkid : {sinkid}]")

  ax1.legend()
  fig1.tight_layout()  
  fig1.gca().invert_xaxis()
  figname = "sink_dr_" + str(sinkid) + ".png"
  print("saving :", figname)
  fig1.savefig(figname)

  ax2.set_xlabel("redshift")
  ax2.set_ylabel(" Msink [Msun] ") 
  ax2.set_title(f"sink mass, sinkid : {sinkid}") 
  ax2.legend()
  fig2.tight_layout()  
  fig2.gca().invert_xaxis()
  figname = "sink_mass_" + str(sinkid) + ".png"
  print("saving :", figname)
  fig2.savefig(figname)

  ax3.set_xlabel("Redshift")  
  ax3.set_ylabel("dMsink/dt[Msun/yr]")
  ax3.set_title(f"sink mass accretion rate, sinkid : {sinkid}")
  ax3.set_yscale('log')
  ax3.legend(loc='best')
  fig3.gca().invert_xaxis()
  fig3.tight_layout()
  figname = "sink_dmdt_" + str(sinkid) + ".png"
  print("saving :", figname)
  fig3.savefig(figname)

def hist_z_sinkspawn(simdir_arr, label_arr, cbar_arr) :
  fig1, ax1 = plt.subplots() 
  isim =0
  for simdir in simdir_arr :
    sinkdir   = simdir+'/sinkdata/'
    csv_files = [file for file in os.listdir(sinkdir) if file.endswith('.csv')]
    Nfiles = len(csv_files)
    nsink = Nfiles
    spawn_z_arr = []
    for isink in range(1, nsink+1):     
      outnum_char = str(isink).zfill(5)
      sinkfile = simdir+'/sinkdata/sink_' + outnum_char + ".csv"
      with open(sinkfile, 'r') as file : 
        first_line = file.readline()
        aexp = float(first_line.split(',')[10]) 
        zspawn = 1.0/aexp - 1.0
        spawn_z_arr.append(zspawn)
        print("isink, zspawn :", isink, zspawn)
    ax1.hist(spawn_z_arr, bins=40, color=cbar_arr[isim], label=label_arr[isim]  ,edgecolor=cbar_arr[isim], histtype='step')
    isim = isim + 1

  imgname_1 = "z_sink_spawn.png"
  ax1.set_xlabel("z_sink_spawn")  
  ax1.set_ylabel("dNsink")
  ax1.set_title("redshift of sink spawn")
  #ax1.set_xscale('log')
  #ax1.set_yscale('log')
  ax1.legend(loc='best')
  fig1.gca().invert_xaxis()
  fig1.tight_layout()
  print("Saving : ", imgname_1)
  fig1.savefig(imgname_1)

def sink_m_dmdt(simdir_arr, sinkid, label_arr, cbar_arr, boxlen_comov) :
  fig1, ax1 = plt.subplots()
  fig2, ax2 = plt.subplots() 
  fig3, ax3 = plt.subplots() 

  isim =0
  print("Age of universe : ", 1.0/H0_IGyr)
  for simdir in simdir_arr :
    sinkdir   = simdir+'/sinkdata/'
    msink         = sink_prop_return(simdir, sinkid, 1)
    sink_dmBHdt   = sink_prop_return(simdir, sinkid, 9)
    aexp_arr      = sink_prop_return(simdir, sinkid, 13)
    
    z_arr           = []
    msink_arr       = []
    dmdt_arr        = []
    dmdtdirect_arr = []

    dmbhoverdt_arr  = []
    dmEddoverdt_arr = [] 
    itime = 0
    Ntime = len(aexp_arr)    
    for aexp in aexp_arr :
      give_units(aexp, boxlen_comov)     
      z = 1.0/aexp - 1.0
      if(z > 7.0 and itime < Ntime - 1) :
        H1 = Omega_m/(aexp*aexp*aexp) + Omega_r/(aexp*aexp*aexp*aexp) + Omega_l - Omega_k/(aexp*aexp)  
        H1 = H0_IGyr * math.sqrt(H1)
        dt = (aexp_arr[itime+1] - aexp_arr[itime]) / (H1 * aexp_arr[itime]) #Gyr
        dt = dt * 1e9
        dmdt = ( msink[itime+1] - msink[itime] )/dt
        z_arr.append(z)
        msink_arr.append(msink[itime] * scale_M/Msun)
        dmdt_arr.append(dmdt * scale_M/Msun)   
        dmdtdirect_arr.append( sink_dmBHdt[itime] * (scale_M/Msun) * (yr/scale_T) )
        #print("aexp : ", aexp , 1.0/aexp -1.0, dmdt/sink_dmBHdt[itime])
      itime = itime + 1
    ax1.plot(z_arr, msink_arr     , label=label_arr[isim], color=cbar_arr[isim], linestyle='-')
    ax2.plot(z_arr, dmdt_arr      , label=label_arr[isim], color=cbar_arr[isim], linestyle='-')
    ax3.plot(z_arr, dmdtdirect_arr, label=label_arr[isim], color=cbar_arr[isim], linestyle='-')
    isim = isim + 1

  ax1.set_xlabel("redshift")
  ax1.set_ylabel(" Msink [Msun] ") 
  ax1.set_title(f"sink mass, sinkid : {sinkid}") 
  ax1.legend()
  fig1.tight_layout()  
  fig1.gca().invert_xaxis()
  figname = "sink_mass_" + str(sinkid) + ".png"
  print("saving :", figname)
  fig1.savefig(figname)

  ax2.set_xlabel("Redshift")  
  ax2.set_ylabel("dMsink/dt[Msun/yr]")
#  ax2.set_ylim(1e-20, 100)
  ax2.set_title(f"sink mass accretion rate, sinkid : {sinkid}")
#  ax2.set_yscale('log')
  ax2.legend(loc='best')
  fig2.gca().invert_xaxis()
  fig2.tight_layout()
  figname = "sink_dmdt_" + str(sinkid) + ".png"
  print("saving :", figname)
  fig2.savefig(figname)

  ax3.set_xlabel("Redshift")  
  ax3.set_ylabel("dMsink/dt[Msun/yr]")
#  ax3.set_ylim(1e-20, 100)
  ax3.set_title(f"sink mass accretion rate, sinkid : {sinkid}")
#  ax3.set_yscale('log')
  ax3.legend(loc='best')
  fig3.gca().invert_xaxis()
  fig3.tight_layout()
  figname = "sink_dmdtdirect_" + str(sinkid) + ".png"
  print("saving :", figname)
  fig3.savefig(figname)


#give file with most number of lines
def give_longest_file(simdir) :
  max_lines    = 0
  max_lines_id = 0

  sinkdir   = simdir+'/sinkdata/'
  csv_files = [file for file in os.listdir(sinkdir) if file.endswith('.csv')]
  Nfiles = len(csv_files)

  for iout in range(1, Nfiles+1):     
    outnum_char = str(iout).zfill(5)
    sinkfile = simdir+'/sinkdata/sink_' + outnum_char + ".csv"
    numlines = 0 
    with open(sinkfile, 'r') as file : 
      for line in file:
        numlines = numlines + 1
    #print(iout, numlines)    
    if(max_lines < numlines) :
      max_lines = numlines
      max_lines_id = iout 
  print("longest :", max_lines_id, max_lines)
  return max_lines_id 

Nx = -1
Ny = -1
Nz = -1
def read_amr2map_binary(file_path):
    global Nx
    global Ny
    global Nz
    with open(file_path, 'rb') as file:
        record_length1 = np.fromfile(file, dtype=np.int32, count=1)[0]
        header_bytes = file.read(record_length1)
        cnt = record_length1/8
        #print('record length', record_length1, cnt)
        narr  = np.frombuffer(header_bytes, dtype=np.float64, count=int(cnt))
        print('[t, xxmax-xxmin, yymax-yymin, zzmax-zzmin] : ', narr)
        record_length2 = np.fromfile(file, dtype=np.int32, count=1)[0]
        if( record_length1 != record_length2 ) :
            raise ValueError("Invalid dimensions read from the file.")


        record_length1 = np.fromfile(file, dtype=np.int32, count=1)[0]
        header_bytes = file.read(record_length1)
        cnt = record_length1/4
        #print('record length', record_length1, cnt)
        arrsize  = np.frombuffer(header_bytes, dtype=np.int32, count=int(cnt))
        print('[imax-imin+1,jmax-jmin+1] : ',arrsize)
        datalength = arrsize[0]*arrsize[1]*4
        Nx = arrsize[0]
        Ny = arrsize[1]
        Nz = Nx
        record_length2 = np.fromfile(file, dtype=np.int32, count=1)[0]
        if( record_length1 != record_length2 ) :
            raise ValueError("Invalid dimensions read from the file.")

        record_length1 = np.fromfile(file, dtype=np.int32, count=1)[0]
        if( record_length1 != datalength ) :
            print("Data length and record length don't match: ",record_length1, datalength)
            raise ValueError("Invalid dimensions read from the file.")
        #print("data record length : ", record_length1, datalength)
        header_bytes = file.read(datalength)
        cnt = datalength/4
        #print('record length', record_length1, cnt)
        toto  = np.frombuffer(header_bytes, dtype=np.float32, count=int(cnt))
        record_length2 = np.fromfile(file, dtype=np.int32, count=1)[0]
        if( record_length1 != record_length2 ) :
            print("Invavild record lengths: ",record_length1, record_length2)
            raise ValueError("Invalid dimensions read from the file.")

        #toto  = np.frombuffer(header_bytes, dtype=np.float32, count=int(cnt))
    return toto       

def plot_dens(simdir, outnum, binname, plotsinks, boxlen_comov): 
  read_units(simdir, outnum)
  print("scale_M : ", scale_L)
  print("scale_L : ", scale_M)
  print("scale_T : ", scale_T)
  print("Z       : ", Z)

  #Nx = 128
  #Ny = 128 
  #Nz = 128
 
  scale_V   = scale_L*scale_L*scale_L
  binpath = simdir + '/' + binname
  print("reading : ", binpath) 
  dens_arr  = read_amr2map_binary(binpath)
 
  if(Nx < 0 or Ny < 0 or Nz < 0):
    print("Nx, Ny, Nz", Nx, Ny, Nz)
    raise ValueError("Invalid grid size.")

  #dens_arr  = dens_arr*scale_M*Nx*Ny*Nz/scale_V
  dens_arr = dens_arr*scale_d
  #Sx = 32768
  #Sy = Sx
  dens_reshape = dens_arr.reshape(Nx, Ny)
  #dens_trans   = np.transpose(dens_reshape)
  dens_trans   = dens_reshape
  extent       = [0, boxlen_comov, 0, boxlen_comov ]

  #plt.xlabel(r'x(Mpch^{-1}(1+z)^{-1})')
  #plt.ylabel(r'y(Mpch^{-1}(1+z)^{-1})')
  #norm = colors.SymLogNorm(linthresh=5, vmin=min(dens_slice), vmax=max(dens_slice))
  plt.xlabel('x')
  plt.ylabel('y')
  plt.figure(figsize=(10, 8))
  

  if(plotsinks == 1) :
    nsink=0
    with open('sinkdata.csv', 'r') as file:
      for line in file:
        nsink=nsink+1
    print("nsink: ", nsink)
    xsink = np.zeros(nsink)
    ysink = np.zeros(nsink)
    zsink = np.zeros(nsink)   
    nsink=0
    with open('sinkdata.csv', 'r') as file:
      for line in file:
        xsink[nsink] = float(line.split(',')[2]) 
        ysink[nsink] = float(line.split(',')[3]) 
        zsink[nsink] = float(line.split(',')[4]) 
        print(nsink, "xsink, ysink", xsink[nsink], ysink[nsink])
        nsink=nsink+1
    plt.scatter(xsink, ysink, marker='o', color='black', s=1)

  plt.imshow(dens_trans, cmap='viridis', norm=LogNorm(),extent = extent, origin='lower' ,label=f'Z = {Z:.2f}' )
  plt.text(0.95, 0.95, f'Z = {Z:.2f}', ha='right', va='top', color='black',fontweight='bold', fontsize=14, transform=plt.gca().transAxes)
  cbar = plt.colorbar()
  cbar.set_label(r'density$ [gm/cm^3 ] $ ')
  figname =  'dens.png'
  
  plt.savefig(figname, dpi=300)

def plot_halodens(simdir, simtyp, outnum, halonum, plotsinks, plotcircle, plotsfr, overwrite_densbin, overwrite_sink, overwrite_star, boxlen_comov, dir, var, fullbox):
  read_units(simdir, outnum)
#  print("scale_M : ", scale_M)
#  print("scale_L : ", scale_L)
#  print("scale_T : ", scale_T)
#  print("scale_T2 : ", scale_T2)
  print("Z          : ", Z)
  Zint = round(Z,2)
  scale_vel = scale_L/scale_T
  scale_E = scale_M*scale_vel*scale_vel
#  print("scale_E : ", scale_E)

  aexp = 1.0/(1.0 + Z)
  Mpch = 3.0856e24 / 0.67739997
  boxlen_comov =  scale_L/(aexp*Mpch)
  print("Lbox_mpc/h : ", boxlen_comov)
  print("")
  comoving = False
  if(comoving) :
    aexp = 1.0

  if(var == 'dens') :
    scale_var = scale_d
    plotcolor = 'viridis'
    typ       =  1
    cbarlabel = r'density$ [ \textrm{gm}/\textrm{cm}^3 ] $ '
    val_low   = 1e-27
    val_high  = 1e-23  
  elif(var == 'temp') :
    gamma     = 1.4
    scale_var = (gamma - 1.0)*scale_T2
    print("scale_T2 : ", scale_var)
    plotcolor = 'twilight_shifted'
    typ       =  5
    cbarlabel = 'Temperature[K]'
    val_low   = 10.0
    val_high  = 1e6  
  elif(var == 'sfr') :
    scale_var =  (scale_M/scale_T) * (yr/Msun)
    plotcolor = 'viridis'
    typ       =  6
    cbarlabel = 'sfr'
    val_low   = 1e-6
    val_high  = 1e4 

  outnum_char    = str(outnum).zfill(5)
  outdir        = simdir+'/output_' + outnum_char
  halofile       = outdir + '/info_'+ outnum_char + "/halodata.csv"


  if(fullbox == True) :
    zoomin     = 1.0
    halo_posx  = 0.5
    halo_posy  = 0.5
    halo_posz  = 0.5  
    rvir       = 0.5/zoomin
    xmin       = halo_posx - rvir
    xmax       = halo_posx + rvir
    ymin       = halo_posy - rvir
    ymax       = halo_posy + rvir
    zmin       = halo_posz - rvir
    zmax       = halo_posz + rvir
    plotcircle = False
  else :
    ## this doens't work 
    if  (not os.path.exists(halofile) ) :
      print("Halo catalog not found, creating one.......................") 
      yt_funcs.create_halo_catalog(simdir, outnum, 0, 0)
      yt_funcs.write_halos(simdir, outnum)
    else :
      print(halofile, "already exists")  
    ##
    ihalo = 0
    with open(halofile, 'r') as file:
      for line in file:
        if(ihalo < halonum-1) :
          ihalo = ihalo + 1
        else :
          break
    halo_posx = float(line.split()[1])
    halo_posy = float(line.split()[2])
    halo_posz = float(line.split()[3])  
    halo_rvir = float(line.split()[4])
    halo_mass = float(line.split()[5])
    rvir    =   aexp*(halo_rvir/1000.0) /boxlen_comov
    print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA Halo Rvir[kpc] : ", rvir*1000.0*boxlen_comov, halo_mass/1e10)

    zoomout  = 4.0
    xmin = max(halo_posx - zoomout*rvir,0)
    xmax = min(halo_posx + zoomout*rvir,1)
    ymin = max(halo_posy - zoomout*rvir,0)
    ymax = min(halo_posy + zoomout*rvir,1)
    zmin = max(halo_posz - zoomout*rvir,0)
    zmax = min(halo_posz + zoomout*rvir,1)

    xbox = xmax-xmin
    ybox = ymax-ymin
    zbox = zmax-zmin

    rbox = -1
    if(xbox < ybox and xbox < zbox):
      rbox = min( halo_posx - xmin, xmax - halo_posx )
    if(ybox < zbox and ybox < zbox):
      rbox = min( halo_posy - ymin, ymax - halo_posy )
    if(zbox < ybox and zbox < xbox):
      rbox = min( halo_posz - zmin, zmax - halo_posz )

    if(rbox > 0):  
      xmin = halo_posx - rbox
      xmax = halo_posx + rbox
      ymin = halo_posy - rbox
      ymax = halo_posy + rbox
      zmin = halo_posz - rbox
      zmax = halo_posz + rbox
      xbox = xmax-xmin
      ybox = ymax-ymin
      zbox = zmax-zmin
    else :
      rbox = min( halo_posx - xmin, xmax - halo_posx )   

  extent  = [(xmin-halo_posx)*boxlen_comov, (xmax-halo_posx)*boxlen_comov, (ymin-halo_posy)*boxlen_comov, (ymax-halo_posy)*boxlen_comov ]
  xbox = xmax-xmin
  ybox = ymax-ymin
  zbox = zmax-zmin

  lmax = 13
  nx = int ( (2**lmax)*xbox - 1)
  ny = nx

  print("nx, ny : ", nx, ny)

  densbin = outdir + '/info_'+ outnum_char
  if ( not os.path.exists(densbin) ) :
    command = "mkdir " + densbin
    os.system(command)

  densbin = outdir + '/info_'+ outnum_char + "/halodens_" + str(halonum) + dir
  if ( (not os.path.exists(densbin) ) or (os.path.exists(densbin) and overwrite_densbin == True) ) :
    command = ("./amr2map -inp " + outdir + " -out " + densbin + " -dir " + dir + " -xmi " + str(xmin) + " -xma " + str(xmax) + " -ymi " + 
               str(ymin) + " -yma " + str(ymax) + " -zmi " + str(zmin) + " -zma " + str(zmax) + " -typ " + str(typ) + " -lma " + str(lmax) + 
               " -nx " + str(nx) + " -ny " + str(ny))
    #print("command : ", command)
    print("writing : ", densbin)
    os.system(command)
  else :
    print(densbin, "already exists")  

  dens_arr  = read_amr2map_binary(densbin)

  if(Nx < 0 or Ny < 0 or Nz < 0):
    print("Nx, Ny, Nz", Nx, Ny, Nz)
    raise ValueError("Invalid grid size.")

  dens_arr = dens_arr*scale_var
  print(var, " min, max : ", np.min(dens_arr), np.max(dens_arr) )
  dens_reshape = dens_arr.reshape(Nx, Ny)
  dens_trans   = dens_reshape

  #dens_trans_clipped = np.clip(dens_trans, val_low, val_high)
  dens_trans_clipped = dens_trans
  plt.imshow(dens_trans_clipped, cmap=plotcolor, norm=LogNorm(vmin=val_low, vmax=val_high),extent = extent, origin='lower' ,label=f'Z = {Z:.2f}' )
  plt.text(0.95, 0.95, f'Z = {Z:.2f}', ha='right', va='top', color='white',fontweight='bold', fontsize=14, transform=plt.gca().transAxes)
  plt.title(simtyp)
  plt.subplots_adjust(left=0.0, right=1.0, top=0.99, bottom=0.1) 
  #plt.gca().set_aspect('equal', adjustable='box')
  cbar = plt.colorbar(pad=0.001)
  cbar.set_label(cbarlabel)
  #cbar.ax.yaxis.labelpad = -30  # Adjust this value as needed

  plot_somepoints = False
  if(plot_somepoints):
    print("plotting some points ", halo_posx, halo_posy, halo_posz)
    somex = []
    somey = []
    somez = []
    with open('tn_xyz.txt', 'r') as file:
      for line in file :
        posx = float(line.split('\t')[0])/boxlen_comov
        posy = float(line.split('\t')[1])/boxlen_comov
        posz = float(line.split('\t')[2])/boxlen_comov
        if( posx>xmin and posx<xmax and posy>ymin and posy<ymax and posz>zmin and posz<zmax):
          px  = (posx- halo_posx)*boxlen_comov
          py  = (posy- halo_posy)*boxlen_comov
          pz  = (posz- halo_posz)*boxlen_comov
          somex.append(px)
          somey.append(py)
          somez.append(pz)
    plt.scatter(somex, somey, marker='.', color="red" ,s=1)  

  plotstars = False
  if(plotsfr == True) :
    xstar = []
    ystar = []
    zstar = []
    tstar = []
    mstar = []
    nstar=0
    starfile = outdir + "/starpos.csv"
    if ( (not os.path.exists(starfile) ) or (os.path.exists(starfile) and overwrite_star == True) ) :  
      command = "./read_star -inp " + outdir + " -out " + starfile
      print("writing : ", starfile)  
      os.system(command)
    else :
      print(starfile, "already exists")

    with open(starfile, 'r') as file:
      for line in file:
        mass = float(line.split(',')[1])
        posx = float(line.split(',')[2]) 
        posy = float(line.split(',')[3]) 
        posz = float(line.split(',')[4]) 
        time  = float(line.split(',')[8])
        if( posx>xmin and posx<xmax and posy>ymin and posy<ymax and posz>zmin and posz<zmax and time > 0.0) :
          star_posx = (posx - halo_posx)*boxlen_comov 
          star_posy = (posy - halo_posy)*boxlen_comov
          star_posz = (posz - halo_posz)*boxlen_comov
          xstar.append( star_posx )
          ystar.append( star_posy ) 
          zstar.append( star_posz )
          tstar.append(    time   )
          mstar.append(    mass   )
          nstar=nstar+1
    if(nstar > 0 and plotstars == True) :   
      print('ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ')
      if(dir == 'z') :  
        plt.scatter(xstar, ystar, marker='.', color="yellow" ,s=1)
      elif (dir == 'x') :
        plt.scatter(ystar, zstar, marker='.', color="yellow" ,s=1)
      elif (dir == 'y') :
        plt.scatter(xstar, zstar, marker='.', color="yellow" ,s=1)

  if(plotsinks == True) :
    xsink = []
    ysink = []
    zsink = []
    nsink=0
    sinkfile = outdir + "/sinkpos.csv"
    if ( (not os.path.exists(sinkfile) ) or (os.path.exists(sinkfile) and overwrite_sink == True) ) :  
      command = "./read_sink -inp " + outdir + " -out " + sinkfile  
      print("writing : ", sinkfile)
      os.system(command)
    else :
      print(sinkfile, "already exists")

    with open(sinkfile, 'r') as file:
      for line in file:
        posx = float(line.split(',')[2]) 
        posy = float(line.split(',')[3]) 
        posz = float(line.split(',')[4]) 
        if( posx>xmin and posx<xmax and posy>ymin and posy<ymax and posz>zmin and posz<zmax ) :
          sink_posx = (posx - halo_posx)*boxlen_comov 
          sink_posy = (posy - halo_posy)*boxlen_comov 
          sink_posz = (posz - halo_posz)*boxlen_comov
          xsink.append( sink_posx )
          ysink.append( sink_posy ) 
          zsink.append( sink_posz )
          print("SINK id in this projection :", float(line.split(',')[0]))
          nsink=nsink+1
    if(nsink > 0) :    
      if(dir == 'z') :  
        plt.scatter(xsink, ysink, marker='o', color='black', s=2)
      elif (dir == 'x') :
        plt.scatter(ysink, zsink, marker='o', color='black', s=2)  
      elif (dir == 'y') :
        plt.scatter(xsink, zsink, marker='o', color='black', s=2)  
    else :
      print("No sink in the given region")
  
  if(plotcircle == True ):
    # rvir is in box units
    rgal =   rvir*1000.0*boxlen_comov
    rcirc =  rvir*boxlen_comov #rcirc =1 => xbox=rvir
    xtext = 0.5 + 1.1*rvir/xbox 
    ytext = 0.5 - 1.1*rvir/xbox 
    if(comoving) :
      text_rgal=f'{rgal:.2f} kpccm'
    else :  
      text_rgal=f'{rgal:.2f} kpc'
    plt.text(xtext, ytext, text_rgal, ha='right', va='top', color='black',fontweight='bold', fontsize=14, transform=plt.gca().transAxes)
    cen_x = 0
    cen_y = 0
    print("circle : ",  rcirc, rvir, xbox)
    circle = plt.Circle((cen_x, cen_y), rcirc , edgecolor='red', facecolor='none', linewidth=2)
    plt.gca().add_patch(circle)

    # rvir is in box units
    rgal =   rbox*1000.0*boxlen_comov
    rcirc =  rbox*boxlen_comov #rcirc =1 => xbox=rvir
    xtext = 0.5 + 0.9*rbox/xbox 
    ytext = 0.5 - 0.9*rbox/xbox 
    if(comoving) :
      text_rgal=f'{rgal:.2f} kpccm'
    else :  
      text_rgal=f'{rgal:.2f} kpc'
    plt.text(xtext, ytext, text_rgal, ha='right', va='top', color='black',fontweight='bold', fontsize=14, transform=plt.gca().transAxes)
    cen_x = 0
    cen_y = 0
    circle = plt.Circle((cen_x, cen_y), rcirc , edgecolor='red', facecolor='none', linewidth=2)
    print("circle : ",  rcirc, rbox, xbox)
    plt.gca().add_patch(circle)

  if(dir == 'z') :
    plt.xlabel(r' $\textrm{X} \left[  \textrm{Mpccm/h} \right]$ ')
    plt.ylabel(r' $\textrm{Y} \left[  \textrm{Mpccm/h} \right] $')
  elif (dir == 'x') :
    plt.xlabel(r' $\textrm{Y} \left[  \textrm{Mpccm/h} \right]$ ')
    plt.ylabel(r' $\textrm{Z} \left[  \textrm{Mpccm/h} \right] $')
  elif (dir == 'y') :
    plt.xlabel(r' $\textrm{X} \left[  \textrm{Mpccm/h} \right]$ ')
    plt.ylabel(r' $\textrm{Z} \left[  \textrm{Mpccm/h} \right] $')

  if(fullbox) :   
    figname =  var+"_"+simtyp + "_z"+ str(Zint) + '_fullbox.png'
  else :
    figname =  var+"_"+simtyp + "_z"+ str(Zint) + '_' + str(halonum) + '.png'  
  print("saving :", figname, plotsfr)  
  plt.savefig(figname, dpi=200,bbox_inches='tight')
  plt.clf()

  if(plotsfr) :
    nbins = 64
    sfr_arr =  np.zeros(nbins)
    time_arr = np.zeros(nbins)
    tbig = max(tstar)*(1.0 + 0.00001)  #Myr
    dtbin = tbig/float(nbins)
    print("dtbin, tage_max :", dtbin, min(tstar), max(tstar))
    # starage in Myr
    ipart = 0
    for ipart in range(0, len(tstar)) :
      starage = tstar[ipart]
      if(starage >= 0) : 
        ibin       = int(starage/tbig*float(nbins))
        sfr_arr[ibin]  = sfr_arr[ibin] + mstar[ipart]*scale_M/Msun
    for ibin in range (0, nbins) :
      time_arr[ibin] =  ( float(ibin) + 0.5 ) * dtbin
      sfr_arr[ibin]  = sfr_arr[ibin]/dtbin/1e6 

    plt.xlabel("Time[Myr]") 
    plt.ylabel("SFR[Msun/yr]")
    plt.plot(time_arr, sfr_arr,  label='' ,marker='o', linestyle='-')
    if(fullbox) :   
      figname = "sfr_" + simtyp + "_z"+ str(Zint) + '_fullbox.png'
    else :
      figname = "sfr_"+simtyp + "_z"+ str(Zint) + '_' + str(halonum) + '.png'  
    print("Saving : ", figname) 
    plt.savefig(figname, dpi=200,bbox_inches='tight')
    plt.clf()

  #plt.savefig('output_figure.png', bbox_inches='tight')

def sink_star_halo_mass(simdir, outnum, overwrite_sink, overwrite_star, overwrite_halo, plotsinkvsstar ,boxlen_comov):
  read_units(simdir, outnum)
  print("scale_M : ", scale_M/Msun)
  print("scale_L : ", scale_L)
  print("scale_T : ", scale_T)
  print("Z       : ", Z)

  outnum_char    = str(outnum).zfill(5)
  outdir        = simdir+'/output_' + outnum_char
  halofile       = outdir + '/info_'+ outnum_char + "/halodata.csv"

  if  (not os.path.exists(halofile) or (os.path.exists(halofile) and overwrite_halo == True ) ): 
    yt_funcs.write_halos(simdir, outnum)
  else :
    print(halofile, "already exists")  

  # read stars
  xstar = []
  ystar = []
  zstar = []
  mstar = []
  nstar=0
  starfile = outdir + "/starpos.csv"
  if ( (not os.path.exists(starfile) ) or (os.path.exists(starfile) and overwrite_star == True) ) :  
    print("writing :", starfile)
    command = "./read_star -inp " + outdir + " -out " + starfile  
    os.system(command)
  else :
    print(starfile, "already exists")

  with open(starfile, 'r') as file:
    for line in file:
      mass  = float(line.split(',')[1])
      posx  = float(line.split(',')[2]) 
      posy  = float(line.split(',')[3]) 
      posz  = float(line.split(',')[4]) 
      #mstar.append ( mass * scale_M/ Msun/52827.663491655)
      mstar.append ( mass * scale_M/ Msun)
      xstar.append( posx )
      ystar.append( posy ) 
      zstar.append( posz )
      nstar=nstar+1
  
  plt.hist(mstar, bins=40, color='blue', edgecolor='black', histtype='step')
  plt.title('Star mass function')
  plt.xlabel('Mass[Msun]')
  plt.ylabel('Number')
  plt.xscale("log")
  plt.yscale("log")
  figname = "star_imf.png"
  print("saving : ", figname)
  plt.savefig(figname)
  plt.clf()

  # read sinks
  xsink = []
  ysink = []
  zsink = []
  msink = []
  nsink = 0
  sinkfile = outdir + "/sinkpos.csv"
  if ( (not os.path.exists(sinkfile) ) or (os.path.exists(sinkfile) and overwrite_sink == True) ) :  
    print("writing :", sinkfile)
    command = "./read_sink -inp " + outdir + " -out " + sinkfile  
    os.system(command)
  else :
    print(sinkfile, "already exists")

  with open(sinkfile, 'r') as file:
    for line in file:
      mass  = float(line.split(',')[1])
      posx  = float(line.split(',')[2]) 
      posy  = float(line.split(',')[3]) 
      posz  = float(line.split(',')[4]) 
      msink.append ( mass * scale_M/ Msun)
      xsink.append( posx )
      ysink.append( posy ) 
      zsink.append( posz )
      nsink=nsink+1

  print("nstar, nsink : ", nstar, nsink)

  ihalo = 0
  mstar_halo = []
  msink_halo = []
  mhost_halo = []
 
  with open(halofile, 'r') as file:
    for line in file:
      halo_posx = float(line.split()[1])
      halo_posy = float(line.split()[2])
      halo_posz = float(line.split()[3])  
      halo_rvir = float(line.split()[4])
      halo_mass = float(line.split()[5])
      numpart   = float(line.split()[6])

      if(numpart > 60) :
        rvir  = (halo_rvir/1000.0) /boxlen_comov
        rvir2 = rvir*rvir 
        mstartot = 0 
        msinktot = 0 
        nsinktot = 0

        for isink in range(0,nsink): 
          sink_dx = xsink[isink] - halo_posx
          sink_dy = ysink[isink] - halo_posy
          sink_dz = zsink[isink] - halo_posz
          sink_dr = sink_dx*sink_dx + sink_dy*sink_dy + sink_dz*sink_dz
          if(sink_dr <= 0.33*rvir2)  :
            msinktot = msinktot + msink[isink] 
            nsinktot = nsinktot + 1

        if(nsinktot > 0 and msinktot> 0 ) :
          print(nsinktot, "halo with sink : ", format(halo_mass, "8.6E"), format(halo_mass/numpart, "8.6E") )  
          mhost_halo.append(halo_mass)
          msink_halo.append(msinktot)

        if(nsinktot > 0 and msinktot> 0 and plotsinkvsstar) :
          print(nsinktot, "halo with sink : ", format(halo_mass, "8.6E") )  
          for istar in range(0,nstar): 
            star_dx = xstar[istar] - halo_posx
            star_dy = ystar[istar] - halo_posy
            star_dz = zstar[istar] - halo_posz
            star_dr = star_dx*star_dx + star_dy*star_dy + star_dz*star_dz
            if(star_dr <= rvir2)  :
              mstartot = mstartot + mstar[istar] 
          mstar_halo.append(mstartot)
              
        ihalo = ihalo + 1

  if(plotsinkvsstar == True) :
    plt.scatter(mstar_halo, msink_halo, marker='o', color='black', s=2)
    plt.title("BH and star mass in Haloes for Z={:.2f}".format(Z))
    plt.xlabel("Mstar[Msun]")
    plt.ylabel("Mbh[Msun]")
    plt.xscale("log")
    plt.yscale("log")
    #plt.legend()
    plt.tight_layout()
    imgname = "sink_vs_star_mass_perhalo.png"
    print("Saving : ", imgname)
    plt.savefig(imgname)
    plt.clf() 

  plt.scatter(mhost_halo, msink_halo, marker='o', color='black', s=2)
  plt.title("BHand halomass Z={:.2f}".format(Z))
  plt.xlabel("Mhalo[Msun]")
  plt.ylabel("Msink[Msun]")
  plt.xscale("log")
  plt.yscale("log")
  #plt.legend()
  plt.tight_layout()
  imgname = "sink_vs_halo_mass_perhalo.png"
  print("Saving : ", imgname)
  plt.savefig(imgname)

def haloprops(simdir_arr, simtyp_arr, outnum_arr, cbar_arr, linestyle_arr, halonum, boxlen_comov, overwrite_star, overwrite_sink, plotsinks):

  fb = Omega_b/Omega_m
  #fb = 1.0

  fig1, ax1 = plt.subplots()
  fig2, ax2 = plt.subplots()
  fig3, ax3 = plt.subplots()
  fig4, ax4 = plt.subplots()
  if(plotsinks) :
    fig5, ax5 = plt.subplots()
    fig6, ax6 = plt.subplots()
    fig7, ax7 = plt.subplots()
    fig8, ax8 = plt.subplots()
    fig9, ax9 = plt.subplots()
    fig10, ax10 = plt.subplots()
    fig11, ax11 = plt.subplots()
    fig12, ax12 = plt.subplots()

  isim=0
  for simdir in simdir_arr:
    outnum = outnum_arr[isim]
    read_units(simdir, outnum)
    print("simdir :", simdir)
    print("Z       : ", Z)
    if(isim==0):
      Zint = round(Z)
      Zint2 = round(Z,2)
    else :
      if(round(Z) != Zint) :
        raise ValueError("Can't compare different redshifts")  
    isim = isim + 1    

  isim=0
  for simdir in simdir_arr:
    simtyp = simtyp_arr[isim]
    outnum = outnum_arr[isim]
    read_units(simdir, outnum)      
    plt.title("Z={:.2f}".format(Z))
    Zint = round(Z)
    aexp = 1.0/(1.0 + Z)

    outnum_char = str(outnum).zfill(5)
    outdir      = simdir+'/output_' + outnum_char
    halofile    = outdir + '/info_'+ outnum_char + "/halodata.csv"

    #read halos
    if  (not os.path.exists(halofile) ) : 
      print('outdir :', outdir)
      print('creating halodata.csv, this needs halo catalog')
      yt_funcs.create_halo_catalog(simdir, outnum, 0, 0)
      yt_funcs.write_halos(simdir, outnum)
    else :
      print(halofile, "already exists")  

    Nhalo = 0 
    with open(halofile, 'r') as file:
      for line in file:
        Nhalo  = Nhalo + 1
    Nhalo  = min(Nhalo, halonum)    

    if  (not os.path.exists(halofile) ) :
      raise ValueError( halofile, "not found")

    #read stars 
    starfile = outdir + "/starpos.csv"
    if ( (not os.path.exists(starfile) ) or (overwrite_star == True) ) :  
      print("writing :", starfile)
      if(os.path.exists(starfile)) :
        command = "rm " + starfile
        os.system(command)
      command = "./read_star -inp " + outdir + " -out " + starfile  
      os.system(command)
    else :
      print(starfile, "already exists")

    #raise ValueError( starfile, "is written")  
    xstar  = []
    ystar  = []
    zstar  = []
    vxstar = []
    vystar = []
    vzstar = []
    mstar  = []
    tstar  = []
    nstar=0
    with open(starfile, 'r') as file:
      for line in file:
        mass  = float(line.split(',')[1])
        posx  = float(line.split(',')[2]) 
        posy  = float(line.split(',')[3]) 
        posz  = float(line.split(',')[4]) 
        velx  = float(line.split(',')[5]) 
        vely  = float(line.split(',')[6]) 
        velz  = float(line.split(',')[7]) 
        time  = float(line.split(',')[8])

        mstar.append( mass * scale_M/ Msun)
        xstar.append( posx )
        ystar.append( posy ) 
        zstar.append( posz )
        vxstar.append( velx )
        vystar.append( vely ) 
        vzstar.append( velz )
        tstar.append( time )
        nstar=nstar+1
    #print( "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ : ", np.min(tstar), np.max(tstar) )
    #print( "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ : ", np.min(mstar), np.max(mstar) )

    #read sinks
    nsink=0
    if(plotsinks == True) :
      sinkfile = outdir + "/sinkpos.csv"
      if ( (not os.path.exists(sinkfile) ) or ( overwrite_sink == True) ) :  
        print("writing :", sinkfile)
        if(os.path.exists(sinkfile)) :
          command = "rm " + sinkfile
          os.system(command)
        command = "./read_sink -inp " + outdir + " -out " + sinkfile  
        os.system(command)
      else :
        print(sinkfile, "already exists")
      xsink = []
      ysink = []
      zsink = []
      dmdbh = []
      dmedd = []
      tsink = []
      msink = []
      sinktostarmass = []
      nsink = 0
      with open(sinkfile, 'r') as file:
        for line in file:
          mass  = float(line.split(',')[1])
          posx  = float(line.split(',')[2]) 
          posy  = float(line.split(',')[3]) 
          posz  = float(line.split(',')[4]) 
          velx  = float(line.split(',')[5]) 
          vely  = float(line.split(',')[6]) 
          velz  = float(line.split(',')[7]) 
          time  = float(line.split(',')[8]) 
          dmdt  = float(line.split(',')[9])
          dedd  = float(line.split(',')[10])
        
#  write(346,'(I10,12(A1,ES20.10))')id(isink),',',sinkmass(isink),&
#          ',',sinkpos(isink,1),',',sinkpos(isink,2),',',sinkpos(isink,3),&
#          ',',sinkvel(isink,1),',',sinkvel(isink,2),',',sinkvel(isink,3),&
#          ',',t-birth(isink),',',dMsmbh(isink), ',' ,dMEd_coarse(isink), ',', bhspin(isink,1), ',', bhspin(isink,2) 

          msink.append ( mass * scale_M/ Msun)
          xsink.append( posx )
          ysink.append( posy ) 
          zsink.append( posz )
          dmdbh.append(dmdt * (scale_M/ Msun) * (yr/scale_T) )
          dmedd.append(dedd)
          nsink=nsink+1
      ax9.hist(np.log10(msink) , bins=16, density=True, color=cbar_arr[isim], label=simtyp_arr[isim], edgecolor=cbar_arr[isim], histtype='step')    
      #ax9.hist(binned_x, binned_xy[1], color=cbar_arr[isim], linestyle=linestyle_arr[isim], label=simtyp_arr[isim])
    print("nstar, nsink, Nhalo : ", nstar, nsink, Nhalo)
  
    ihalo = 0
    mhost_halo = []
    #star related
    mstar_halo = []
    mstar_to_mhalo = []
    sfr_halo = []
    rhost_halo = []

    #sink related 
    if(plotsinks):
      msink_halo = []
      Nsink_halo = []
      drsink_halo = []
      dmdbh_halo = []
      dmedd_halo = []
      dbhtoedd_halo = []
      vsigma_halo = []
      mstar_halo_withsink = []
      mhost_halo_withsink = []

    with open(halofile, 'r') as file:
      for line in file:
        if(ihalo < halonum and ihalo < Nhalo)  :
          #print("halo : ", ihalo, Nhalo,end=' ')
          print(f"\rhalo: {ihalo}, Nhalo: {Nhalo}", end='')
          halo_posx = float(line.split()[1])
          halo_posy = float(line.split()[2])
          halo_posz = float(line.split()[3])  
          halo_rvir = float(line.split()[4])
          halo_mass = float(line.split()[5])
          numpart   = float(line.split()[6])
          rvir = (halo_rvir/1000.0) /boxlen_comov
          rvir2 = rvir*rvir 

          if(numpart > 10) :
            #calculate star props
            nstartot = 0.0
            mstartot = 0.0
            vsigma    = 0.0
            vmean     = 0.0
            mass_bin = 0.0
            for istar in range(0,nstar): 
              star_dx = xstar[istar] - halo_posx
              star_dy = ystar[istar] - halo_posy
              star_dz = zstar[istar] - halo_posz
              star_dvx = vxstar[istar]
              star_dvy = vystar[istar]
              star_dvz = vzstar[istar]
              star_dr = star_dx*star_dx + star_dy*star_dy + star_dz*star_dz
              star_dv = star_dvx*star_dvx + star_dvy*star_dvy + star_dvz*star_dvz
              if(star_dr <= rvir2) :
                nstartot = nstartot + 1.0
                mstartot = mstartot + mstar[istar]
                vsigma = vsigma + star_dv*mstar[istar]
                ##sfr, check this
                tage       = abs( (tstar[istar])    )  #Myr
                if(tage < 10.0) :
                  mass_bin = mass_bin + mstar[istar]
                ##  
            if(nstartot > 0 ) :
              mstar_halo.append(mstartot)
              mstar_to_mhalo.append(mstartot/(halo_mass*fb))
              mhost_halo.append(halo_mass)
              rhost_halo.append(halo_rvir)
              sfr_halo.append( mass_bin/(10.0*1e6) ) #Msun/yr
              #print("sfr :", mass_bin/(10.0*1e6))
      
            #calculate sink props
            if(plotsinks):      
              msinktot = 0 
              nsinktot = 0
              drsinktot = 0
              dmdbhtot = 0
              dmdedtot = 0
              for isink in range(0,nsink): 
                sink_dx = xsink[isink] - halo_posx
                sink_dy = ysink[isink] - halo_posy
                sink_dz = zsink[isink] - halo_posz
                sink_dr2 = sink_dx*sink_dx + sink_dy*sink_dy + sink_dz*sink_dz
                #print("isink :", isink, nsink, sink_dr, rvir2)
                rgalaxy =  (50.0 / 1000.0) * (aexp / boxlen_comov)
                rgalaxy2 = rgalaxy*rgalaxy
                if(sink_dr2 <= (0.5)*(0.5)*rvir2) :
                  #print("     sink_dr :", math.sqrt(sink_dr2), rgalaxy, math.sqrt(sink_dr2)/rgalaxy)
                  #print("isink tot :", isink, nsink)
                  nsinktot = nsinktot + 1.0
                  drsinktot = drsinktot + math.sqrt(sink_dr2)
                  msinktot = msinktot + msink[isink]
                  dmdbhtot = dmdbhtot + dmdbh[isink]
                  dmdedtot = dmdedtot + dmedd[isink]

              if(nstartot > 0 and nsinktot > 0) :
                mstar_halo_withsink.append(mstartot)
                mhost_halo_withsink.append(halo_mass)
                msink_halo.append(msinktot)
                Nsink_halo.append(nsinktot)
                drsink_halo.append( aexp*(boxlen_comov*1000.0/0.607)*drsinktot/float(nsinktot))
                sinktostarmass.append(msinktot/mstartot)
                dmdbhtot = dmdbhtot/nsinktot
                dmdedtot = dmdedtot/nsinktot
                dbhtoedd_halo.append(dmdbhtot)
                vsigma = math.sqrt(vsigma)/mstartot
                vsigma_halo.append(vsigma)
        ihalo = ihalo + 1
      print('')

    ##
    #ax1.scatter(mhost_halo, mstar_halo    , marker='o', color=cbar[isim], s=20)
    binned_xy  = bindata(mhost_halo, mstar_halo, 8, True, 1)
    binned_x   = binned_xy[0]
    ax1.plot(binned_x, binned_xy[1], color=cbar_arr[isim], linestyle=linestyle_arr[isim], label=simtyp_arr[isim])
    ax1.fill_between(binned_x, binned_xy[3], binned_xy[2], color=cbar_arr[isim], alpha=0.25)
    ##

    ##
    #ax2.scatter(mhost_halo, mstar_to_mhalo, marker='o', color=cbar[isim], s=20) 
    binned_xy  = bindata(mhost_halo, mstar_to_mhalo, 8, True, 1)
    binned_x   = binned_xy[0]
    ax2.plot(binned_x, binned_xy[1], color=cbar_arr[isim] , linestyle=linestyle_arr[isim], label=simtyp_arr[isim])
    ax2.fill_between(binned_x, binned_xy[3], binned_xy[2], color=cbar_arr[isim], alpha=0.25)
    ##
   
    ## 
    #ax3.scatter(mstar_halo, sfr_halo, marker='o', color=cbar[isim], s=20)
    binned_xy  = bindata(mstar_halo, sfr_halo, 8, True, 1)
    binned_x   = binned_xy[0]
    ax3.plot(binned_x, binned_xy[1], color=cbar_arr[isim] , linestyle=linestyle_arr[isim], label=simtyp_arr[isim])
    ax3.fill_between(binned_x, binned_xy[3], binned_xy[2], color=cbar_arr[isim], alpha=0.25)
    #print("$$$$$$$$$: ", min(binned_xy[2]), min(mstar_halo), max(mstar_halo))
    #print(binned_x)
    #print(binned_xy[1])
    #print(binned_xy[2])
    #print(binned_xy[3])

    ##

    ##
    #ax4.scatter(mhost_halo, sfr_halo, marker='o', color=cbar[isim], s=20)
    binned_xy  = bindata(mhost_halo, sfr_halo, 8, True, 1)
    binned_x   = binned_xy[0]
    ax4.plot(binned_x, binned_xy[1], color=cbar_arr[isim] , linestyle=linestyle_arr[isim], label=simtyp_arr[isim])
    ax4.fill_between(binned_x, binned_xy[3], binned_xy[2], color=cbar_arr[isim], alpha=0.25)
    print("ax4 :", np.min(mhost_halo)/1e9, np.max(mhost_halo)/1e9)
    print("binned :", np.min(binned_xy[1]), np.max(binned_xy[1]))
    print("binned :", np.min(binned_x)/1e9, np.max(binned_x)/1e9)
    ##

    if(plotsinks):
      ax5.scatter(mhost_halo_withsink, msink_halo    , marker='o', color=cbar_arr[isim], s=20,label=simtyp_arr[isim])
      ax12.scatter(mstar_halo_withsink, msink_halo    , marker='o', color=cbar_arr[isim], s=20,label=simtyp_arr[isim])
      #binned_xy  = bindata(mhost_halo_withsink, msink_halo, 8, True, 2)
      #binned_x   = binned_xy[0]
      #ax5.plot(binned_x, binned_xy[1], color=cbar[isim], linestyle='-', label=simtyp)
      #ax5.fill_between(binned_x, binned_xy[3], binned_xy[2], color=cbar[isim], alpha=0.25)

      #mstar_halo_arr = [pow(10,8.5), pow(10,12)]
      #mBH_halo_arr = [pow(10, 5.8), pow(10, 9.81)] 
      mstar_halo_arr = [pow(10,7.4), pow(10,11.2)]  
      mBH_halo_arr = [pow(10, 4), pow(10, 8)]
      ax12.plot( mstar_halo_arr ,mBH_halo_arr, linestyle='--', color='gray', label='Reines, Volonteri 2015')


      ax6.scatter(mhost_halo_withsink, sinktostarmass, marker='o', color=cbar_arr[isim], s=20,label=simtyp_arr[isim])
      #binned_xy  = bindata(mhost_halo_withsink, sinktostarmass, 8, True, 2)
      #binned_x   = binned_xy[0]
      #ax6.plot(binned_x, binned_xy[1], color=cbar[isim], linestyle='-', label=simtyp)
      #ax6.fill_between(binned_x, binned_xy[3], binned_xy[2], color=cbar[isim], alpha=0.25)

      ax7.scatter(mhost_halo_withsink, dbhtoedd_halo , marker='o', color=cbar_arr[isim], s=20,label=simtyp_arr[isim])
      #binned_xy  = bindata(mhost_halo_withsink, dbhtoedd_halo, 8, True, 2)
      #binned_x   = binned_xy[0]
      #ax7.plot(binned_x, binned_xy[1], color=cbar[isim], linestyle='-', label=simtyp)
      #ax7.fill_between(binned_x, binned_xy[3], binned_xy[2], color=cbar[isim], alpha=0.25)

      ax8.scatter(msink_halo         , vsigma_halo   , marker='o', color=cbar_arr[isim], s=20, label=simtyp_arr[isim])

      ax10.scatter(mhost_halo_withsink, Nsink_halo    , marker='o', color=cbar_arr[isim], s=20,label=simtyp_arr[isim])

      #ax11.scatter(mhost_halo_withsink, drsink_halo    , marker='o', color=cbar_arr[isim], s=20,label=simtyp_arr[isim])
      binned_xy  = bindata(mhost_halo_withsink, drsink_halo, 8, True, 1)
      binned_x   = binned_xy[0]
      ax11.plot(binned_x, binned_xy[1], color=cbar_arr[isim], linestyle=linestyle_arr[isim], label=simtyp_arr[isim])
      ax11.fill_between(binned_x, binned_xy[3], binned_xy[2], color=cbar_arr[isim], alpha=0.25)
     
     
    isim = isim + 1
     
  imgname_1 = "mstar_vs_mhalo_z" + str(Zint) + ".png"
  imgname_2 = "mstar_to_mhalo_z" + str(Zint) + ".png"
  imgname_3 = "sfr_vs_mstar_z" + str(Zint) + ".png"
  imgname_4 = "sfr_vs_mhalo_z" + str(Zint) + ".png"
  
  ax1.set_xlabel("Mhalo[Msun]")  
  ax1.set_ylabel("Mstar[Msun]")
  ax1.set_title("Total star mass per halo at z={:.2f}".format(Z))
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  ax1.legend(loc='best')
  ax1.set_xlim(1e8, 2e10)
  fig1.tight_layout()

  npoints = 4
  ipoint =0  
  aexp = 1.0/(1.0 + Z)
  tuniv = aexp * 13.6
  tuniv=0
  print("ZZZZZZZZZ ==== ", Z)
  logMstar_min = 4.0
  logMstar_max = 9.0
  logMstar_arr = np.linspace(logMstar_min, logMstar_max, npoints)
  logsfr_arr_p   = np.zeros(npoints)
  logsfr_arr_m   = np.zeros(npoints)
  logsfr_arr     = np.zeros(npoints)
  for logMstar in logMstar_arr : 
    logsfr_arr[ipoint]   = (0.80        - 0.017              )*logMstar - (6.487        - 0.039             ) #Iyer 2018
    logsfr_arr_p[ipoint] = (0.80 +0.029 - 0.017 +0.010*tuniv )*logMstar - (6.487 -0.282 - 0.039 -0.008*tuniv) #Iyer 2018
    logsfr_arr_m[ipoint] = (0.80 -0.029 - 0.017 -0.010*tuniv )*logMstar - (6.487 +0.282 - 0.039 +0.008*tuniv) #Iyer 2018
    ipoint = ipoint + 1    
  ipoint = 0   
  for ipoint in range(0, npoints) : 
    logsfr_arr[ipoint]   = 10**logsfr_arr[ipoint]
    logsfr_arr_p[ipoint] = 10**logsfr_arr_p[ipoint]
    logsfr_arr_m[ipoint] = 10**logsfr_arr_m[ipoint]  
    logMstar_arr[ipoint] = 10**logMstar_arr[ipoint]
  ax3.plot( logMstar_arr ,logsfr_arr, linestyle='--', color='gray', label='Iyer 2018')

#  ax3.plot( logMstar_arr ,logsfr_arr_p, linestyle='--', color='gray', label='')
#  ax3.plot( logMstar_arr ,logsfr_arr_m, linestyle='--', color='gray', label='')
#  ax3.fill_between(logMstar_arr, logsfr_arr_p, logsfr_arr_m, color='gray', alpha=0.1, label='Iyer 2018')

  sfms_obs = np.loadtxt('SFMS_obs.txt')
  redz_obs = sfms_obs[:,0]          #Redshift
  mstar_obs = 10**sfms_obs[:,1]     #Stellar mass in Msun
  mstar_obs_up = 10**(sfms_obs[:,1]+sfms_obs[:,2]) - 10**sfms_obs[:,1] #Upper limit stellar mass in Msun
  mstar_obs_lo = 10**sfms_obs[:,1] - 10**(sfms_obs[:,1] - sfms_obs[:,3])  #Lower limit stellar mass in Msun
  sfr_obs = sfms_obs[:,4]     #SFR in Msun/yr
  sfr_obs_up = sfms_obs[:,5]  #SFR upper lim in Msun/yr
  sfr_obs_lo = sfms_obs[:,6]  #SFR lower lim in Msun/yr

  for ii in range(0, len(redz_obs)) :
    if( abs (redz_obs[ii] -Z) < 0.8) :
      ax3.errorbar(mstar_obs[ii],sfr_obs[ii],yerr=[[sfr_obs_lo[ii]],[sfr_obs_up[ii]]],\
                xerr=[[mstar_obs_lo[ii]],[mstar_obs_up[ii]]],\
                mfc='grey',color='grey',fmt='oC2')


  mstomh_obs = np.loadtxt('mstartomhalo_obs.txt')
  z_obs   = mstomh_obs[:,0]          #Redshift
  mh_obs     = mstomh_obs[:,1]     #halo mass in Msun
  msmh_obs     = mstomh_obs[:,2]     #halo mass in Msun

  for ii in range(0, len(mh_obs)) :
    if( abs (redz_obs[ii] -Z) < 0.8) :
      ax2.scatter(mh_obs, msmh_obs, color='grey', marker='o')

  ax2.set_xlabel("Mhalo[Msun]")  
  ax2.set_ylabel("Mstar/(Mhalo*fb)")
  ax2.set_title("ratio of star mass to halo mass at z={:.2f}".format(Z))
  ax2.set_xscale('log')
  ax2.set_yscale('log')
  ax2.legend(loc='best')
  ax2.set_xlim(1e7, 5e11)
  fig2.tight_layout()


  #ax3.set_ylim(5e-3, 100)
  ax3.set_xlabel("Mstar[Msun]")  
  ax3.set_ylabel("SFR[Msun/yr]")
  ax3.set_title("SFR vs star mass per halo at z={:.2f}".format(Z))
  ax3.set_xscale('log')
  ax3.set_yscale('log')
  ax3.set_xlim(1e4, 1e11)
  ax3.legend(loc='best')
  fig3.tight_layout()

  ax4.set_xlabel("Mhalo[Msun]")  
  ax4.set_ylabel("SFR[Msun/yr]")
  ax4.set_title("SFR vs Halo mass at z={:.2f}".format(Z))
  ax4.set_xscale('log')
  ax4.set_yscale('log')
  ax4.set_xlim(1e8, 2e10)
  ax4.legend(loc='best')
  fig4.tight_layout()


  #ax4.set_xlabel("Mhalo[Msun]")  
  #ax4.set_ylabel("Rhalo[Kpc]")
  #ax4.set_title("Halo mass vs virial radius at Z={:.2f}".format(Z))
  #ax4.set_xscale('log')
  #ax4.set_yscale('log')
  #ax4.legend()
  #fig4.tight_layout()

  print("Saving : ", imgname_1)
  fig1.savefig(imgname_1)
  fig1.clf()
  print("Saving : ", imgname_2)
  fig2.savefig(imgname_2)
  fig2.clf()
  print("Saving : ", imgname_3)
  fig3.savefig(imgname_3)
  fig3.clf()
  print("Saving : ", imgname_4)
  fig4.savefig(imgname_4)
  fig4.clf()

  if(plotsinks):
    imgname_5 = "msink_vs_mhalo_z" + str(Zint) + ".png"
    imgname_12 = "msink_vs_mstarhalo_z" + str(Zint) + ".png"
    imgname_10 = "Nsink_vs_mhalo_z" + str(Zint) + ".png"
    imgname_11 = "drsink_vs_mhalo_z" + str(Zint) + ".png"

    imgname_6 = "msinktomstar_vs_mhalo_z" + str(Zint) + ".png"
    imgname_7 = "dmbhtodedd_vs_mhalo_z" + str(Zint) + ".png"
    imgname_8 = "m_sigma" + str(Zint) + ".png"
    imgname_9 = "sink_imf" + str(Zint) + ".png"


    ax5.set_xlabel("Mhalo[Msun]")  
    ax5.set_ylabel("Msink[Msun]")
    ax5.set_title("Total sink mass per Halo at z={:.2f}".format(Z))
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ax5.legend(loc='best')
    ax5.set_xlim(1e8, 2e10)
    fig5.tight_layout()

    ax12.set_xlabel("Mstellar[Msun]")  
    ax12.set_ylabel("Msink[Msun]")
    ax12.set_title("local Msink vs Mstar per Halo at z={:.2f}".format(Z))
    ax12.set_xscale('log')
    ax12.set_yscale('log')
    ax12.legend(loc='best')
    ax12.set_xlim(1e4, 1e12)
    fig12.tight_layout()


    ax10.set_xlabel("Mhalo[Msun]")  
    ax10.set_ylabel("Nsink")
    ax10.set_title("Number of sinks per Halo at z={:.2f}".format(Z))
    ax10.set_xscale('log')
    ax10.set_yscale('log')
    ax10.legend(loc='best')
    ax10.set_xlim(1e8, 2e10)
    fig10.tight_layout()

    ax11.set_xlabel("Mhalo[Msun]")  
    ax11.set_ylabel("drsink[Kpc]")
    ax11.set_title("Average BH distance from center of Halo at z={:.2f}".format(Z))
    ax11.set_xscale('log')
    ax11.set_yscale('log')
    ax11.legend(loc='best')
    ax11.set_xlim(1e8, 2e10)
    fig11.tight_layout()


    ax6.set_xlabel("Mhalo[Msun]")  
    ax6.set_ylabel("Msink/Mstar")
    ax6.set_title("BH to Star mass ratio per Halo at z={:.2f}".format(Z))
    ax6.set_xscale('log')
    ax6.set_yscale('log')
    ax6.set_xlim(1e8, 2e10)
    ax6.legend(loc='best')
    fig6.tight_layout()

    ax7.set_xlabel("Mhalo[Msun]")  
    ax7.set_ylabel("dMBH/dMEdd")
    ax7.set_title("Average BH accretion rate to Eddington rate ratio per Halo at z={:.2f}".format(Z))
    ax7.set_xscale('log')
    ax7.set_yscale('log')
    ax7.set_xlim(1e8, 2e10)
    ax7.legend(loc='best')
    fig7.tight_layout()

    ax8.set_xlabel("Msink[Msun]")  
    ax8.set_ylabel("star velocity sigma")
    ax8.set_title("M-Sigma relation for BH at z={:.2f}".format(Z))
    ax8.set_xscale('log')
    ax8.set_yscale('log')
    ax8.legend(loc='best')
    fig8.tight_layout()

    ax9.set_xlabel("log10(Msink[Msun])")  
    ax9.set_ylabel("")
    ax9.set_title("Sink mass function at z={:.2f}".format(Z))
    ax9.set_xscale('log')
    ax9.set_yscale('log')
    ax9.legend(loc='best')
    fig9.tight_layout()

    print("Saving : ", imgname_5)  
    fig5.savefig(imgname_5)
    fig5.clf()
    print("Saving : ", imgname_6)
    fig6.savefig(imgname_6)
    fig6.clf()
    print("Saving : ", imgname_7)
    fig7.savefig(imgname_7)
    fig7.clf()
    print("Saving : ", imgname_8)
    fig8.savefig(imgname_8)
    fig8.clf()
    print("Saving : ", imgname_9)
    fig9.savefig(imgname_9)
    fig9.clf()
    print("Saving : ", imgname_10)
    fig10.savefig(imgname_10)
    fig10.clf()
    print("Saving : ", imgname_11)
    fig11.savefig(imgname_11)
    fig11.clf()
    print("Saving : ", imgname_12)
    fig12.savefig(imgname_12)
    fig12.clf()


#outstat = 1 for mean and std deviation
#outstat = 2 for median and range   
def bindata(xarr, yarr, nbin, logdata, outstat):

  if( len(xarr) != len(yarr) ) :
    raise ValueError("Data sizes must be same to bin ", len(xarr), len(yarr))
  
  ndata    = len(xarr) 
  perc_p   = 84.0/100.0
  perc_m   = 1.0 - perc_p

  x_binned = []
  y_binned = []
  y_p      = []
  y_m      = []
  xmin     = np.min(xarr)  * (1.0 -  1e-6)
  xmax     = np.max(xarr)  * (1.0 +  1e-6)  
  if(logdata) :
    dx = pow ( (xmax/xmin), 1.0/nbin  )
  else : 
    dx = (xmax - xmin)/nbin  

  npoint_tot = 0
  for ibin in range(0, nbin) :
    ipf = float(ibin)
    if(logdata) :
      x_left  = xmin*pow(dx, ipf)
      x_right = xmin*pow(dx, ipf + 1.0)
      x_mid   = math.sqrt(x_left*x_right)
    else :   
      x_left  = xmin + ipf*dx
      x_right = xmin + (ipf+1.0)*dx
      x_mid   = (x_left + x_right)/2.0
    y_sum   = 0.0
    y_ordered = []
    npoint   = 0
    for i in range(0, ndata) :
      if( xarr[i] > x_left and xarr[i] <= x_right ) :
        npoint = npoint + 1
        y_sum = y_sum + yarr[i]
        y_ordered.append( yarr[i] )
    npoint_tot = npoint_tot + npoint    
    if(npoint > 0) :
      x_binned.append (x_mid)
      if(outstat == 1) :
        y_mean  = y_sum/float(npoint)
        y_binned.append ( y_mean )
        std_dev = 0.0
        for ip in range(0, npoint) :
          std_dev = std_dev + (y_ordered[ip] - y_mean)*(y_ordered[ip] - y_mean)
        std_dev = math.sqrt(std_dev)/float(npoint)
        y_p.append(y_mean + 1.0*std_dev)
        y_m.append(y_mean - 1.0*std_dev)  
      else :  
        for ip in range(0, npoint-1) :
          for jp in range(ip+1, npoint) :
            if( y_ordered[ip] > y_ordered[jp] ) :
              ytemp         = y_ordered[ip]
              y_ordered[ip] = y_ordered[jp]
              y_ordered[jp] = ytemp
        if( npoint%2 == 1 ) :
          midpoint = int ( (npoint + 1)/2 ) -1
          y_binned.append( y_ordered[ midpoint ] )
        else :
          midpoint = int(npoint/2) -1
          y_binned.append( (y_ordered[midpoint] + y_ordered[midpoint+1])/2.0 )
        #print(ibin, x_left, x_right, y_median, npoint, midpoint)
        #print(y_ordered)

        nleft  = math.floor( (npoint-1)*perc_p )
        nright = math.ceil(  (npoint-1)*perc_p )
        p_frac  =  round ( (npoint-1)*perc_p%1, 3 )
        y_p.append(  y_ordered[nleft]  + p_frac * ( y_ordered[nright] - y_ordered[nleft] )  )
        y_plus = y_ordered[nleft]  + p_frac * ( y_ordered[nright] - y_ordered[nleft] )
        nleft  = math.floor( (npoint-1)*perc_m )
        nright = math.ceil(  (npoint-1)*perc_m )
        p_frac  =  round ( (npoint-1)*perc_m%1, 3 )
        y_m.append(  y_ordered[nleft]  + p_frac * ( y_ordered[nright] - y_ordered[nleft] )  )
        y_minus =    y_ordered[nleft]  + p_frac * ( y_ordered[nright] - y_ordered[nleft] )

  #print("Binned data length :", npoint_tot, ndata, len(x_binned))
  if(npoint_tot != ndata) :
    print("Binning might not have happened properly")

  nbin = len(x_binned)  
  for ibin in range(0, nbin-1) :
    for jbin in range(ibin+1, nbin) : 
      if(x_binned[ibin] > x_binned[jbin]) :
        xtemp          = x_binned[ibin]
        x_binned[ibin] = x_binned[jbin]
        x_binned[jbin] = xtemp

        ytemp          = y_binned[ibin]
        y_binned[ibin]    = y_binned[jbin]
        y_binned[jbin]    = ytemp

        ytemp          = y_p[ibin]
        y_p[ibin]      = y_p[jbin]
        y_p[jbin]      = ytemp

        ytemp          = y_m[ibin]
        y_m[ibin]      = y_m[jbin]
        y_m[jbin]      = ytemp

  #print("x_binned : ", x_binned)
  #print("y_binned : ", y_binned)
  #print("y_p : ", y_p)
  #print("y_m :", y_m)
  #print("##################################################################################")

  x_out  = []
  y_out  = []
  ym_out = []
  yp_out = []

  for ibin in range(0, nbin-1) :
    if(y_binned[ibin] > 1e-7) :
      x_out.append(x_binned[ibin])
      y_out.append(y_binned[ibin])
      ym_out.append(y_m[ibin])
      yp_out.append(y_p[ibin])

  return [x_out, y_out, ym_out, yp_out]   

def read_mfblog(simdir, overwrite) :
  print("simdir :", simdir)
  mfblog_all = simdir  + "/mfblog_all.txt" 
  if  (not  os.path.exists(mfblog_all) or overwrite ) :
    if(os.path.exists(mfblog_all) and overwrite) :
      command = "rm " + mfblog_all
      print('command :', command)
      os.system(command)

    command = "touch " + mfblog_all
    print('command :', command)
    os.system(command)

    nlogfile = 1 
    for ifile in range(1,5) :
      simlog = simdir + "/simulation_" + str(ifile) + ".log"
      mfblog = simdir + "/mfblog_"     + str(ifile) + ".txt"
      if(os.path.exists(simlog)) :
        command = "grep -rH 'MFB z= ' " +  simlog + " > " + mfblog
        print('command :', command)
        print('writing :',  mfblog)
        os.system(command)
        command = "cat " + mfblog + " >> " + mfblog_all
        print('command :', command)
        os.system(command)
        nlogfile = nlogfile + 1
      else :
        simlog = simdir + "/simulation" + ".log"
        mfblog = simdir + "/mfblog"     + ".txt"
        command = "grep -rH 'MFB z= ' " +  simlog + " > " + mfblog
        print('command :', command)
        print('writing :',  mfblog)
        os.system(command)
        command = "cat " + mfblog + " >> " + mfblog_all
        print('command :', command)
        os.system(command)
        nlogfile = nlogfile + 1
        break
  else : 
    print(mfblog_all, 'exists')

  z_arr  = []
  N_arr  = []
  nH_arr = []
  T_arr  = []
  Z_arr  = []
  dx_arr = []
  E_arr  = []

  with open(mfblog_all, 'r') as file:  
    lines  = file.readlines()
    Nlines = len(lines)
    for line in lines:
      file_path, data = line.split(':')
      parts           = data.split()
      tag             = parts[0] + parts[1]
      if(float(parts[2]) > 8.4) :
        z_arr.append  (  float(parts[2])      )
        N_arr.append  (  float(parts[4])      )
        nH_arr.append (  10**float(parts[6])      )
        T_arr.append  (  10**float(parts[7])      )
        Z_arr.append  (  10**float(parts[8])  )
        dx_arr.append (  10**float(parts[10]) )
        E_arr.append  (  10**float(parts[10]) )
  return [nH_arr, T_arr]

def plot_mfblog(simdir_arr, simtyp_arr, cbar_arr, overwrite) :

  fig1, ax1 = plt.subplots()
  fig2, ax2 = plt.subplots()
  fig3, ax3 = plt.subplots()
  fig4, ax4 = plt.subplots()
  fig5, ax5 = plt.subplots()
  fig6, ax6 = plt.subplots()

  isim  = 0 
  for simdir in simdir_arr :
    print("simdir :", simdir)
    mfblog_all = simdir  + "/mfblog_all.txt" 
    if  (not  os.path.exists(mfblog_all) or overwrite ) :
      if(os.path.exists(mfblog_all) and overwrite) :
        command = "rm " + mfblog_all
        print('command :', command)
        os.system(command)

      command = "touch " + mfblog_all
      print('command :', command)
      os.system(command)

      nlogfile = 1 
      for ifile in range(1,5) :
        simlog = simdir + "/simulation_" + str(ifile) + ".log"
        mfblog = simdir + "/mfblog_"     + str(ifile) + ".txt"
        if(os.path.exists(simlog)) :
          command = "grep -rH 'MFB z= ' " +  simlog + " > " + mfblog
          print('command :', command)
          print('writing :',  mfblog)
          os.system(command)
          command = "cat " + mfblog + " >> " + mfblog_all
          print('command :', command)
          os.system(command)
          nlogfile = nlogfile + 1
        else :
          simlog = simdir + "/simulation" + ".log"
          mfblog = simdir + "/mfblog"     + ".txt"
          command = "grep -rH 'MFB z= ' " +  simlog + " > " + mfblog
          print('command :', command)
          print('writing :',  mfblog)
          os.system(command)
          command = "cat " + mfblog + " >> " + mfblog_all
          print('command :', command)
          os.system(command)
          nlogfile = nlogfile + 1
          break
    else : 
      print(mfblog_all, 'exists')

    with open(mfblog_all, 'r') as file:  
      lines  = file.readlines()
      Nlines = len(lines)
      z_arr  = []
      N_arr  = []
      nH_arr = []
      T_arr  = []
      Z_arr  = []
      dx_arr = []
      E_arr  = []
      for line in lines:
        file_path, data = line.split(':')
        parts           = data.split()
        tag             = parts[0] + parts[1]
        if(float(parts[2]) > 8.4) :
          z_arr.append  (  float(parts[2])      )
          N_arr.append  (  float(parts[4])      )
          nH_arr.append (  float(parts[6])      )
          T_arr.append  (  float(parts[7])      )
          Z_arr.append  (  10**float(parts[8])  )
          dx_arr.append (  10**float(parts[10]) )
          E_arr.append  (  10**float(parts[10]) )
      ax1.hist(z_arr , bins=256, density=True, color=cbar_arr[isim], label=simtyp_arr[isim], edgecolor=cbar_arr[isim], histtype='step')
      ax2.hist(N_arr , bins=256, density=True, color=cbar_arr[isim], label=simtyp_arr[isim], edgecolor=cbar_arr[isim], histtype='step')
      ax3.hist(nH_arr, bins=256, density=True, color=cbar_arr[isim], label=simtyp_arr[isim], edgecolor=cbar_arr[isim], histtype='step')
      ax4.hist(T_arr , bins=256, density=True, color=cbar_arr[isim], label=simtyp_arr[isim], edgecolor=cbar_arr[isim], histtype='step')
      ax5.hist(Z_arr , bins=256, density=True, color=cbar_arr[isim], label=simtyp_arr[isim], edgecolor=cbar_arr[isim], histtype='step')
      ax6.hist(E_arr , bins=256, density=True, color=cbar_arr[isim], label=simtyp_arr[isim], edgecolor=cbar_arr[isim], histtype='step')
      isim = isim + 1
        #print("iline :", parts[0], parts[1] , parts[2], parts[4], parts[5], parts[6], parts[7], parts[8], parts[9], parts[10]  )

  imgname_1 = "hist_z.png"
  imgname_2 = "hist_NSN.png"
  imgname_3 = "hist_nH.png"
  imgname_4 = "hist_T.png"
  imgname_5 = "hist_Z.png"
  imgname_6 = "hist_E.png"
  
  ax1.set_xlabel("Redshift(z)")  
  ax1.set_ylabel("")
  ax1.set_title("hist_z")
  ax1.set_yscale('log')
  ax1.legend(loc='best')
  fig1.gca().invert_xaxis()
  fig1.tight_layout()

  ax2.set_xlabel("Number of SN events per cell")  
  ax2.set_ylabel("")
  ax2.set_title("hist_NSN")
  ax2.set_yscale('log')
  ax2.legend(loc='best')
  fig2.tight_layout()

  ax3.set_xlabel("log10(nH[cc])")  
  ax3.set_ylabel("")
  ax3.set_title("hist_nH")
  #ax3.set_yscale('log')
  ax3.legend(loc='best')
  fig3.tight_layout()

  ax4.set_xlabel("log10(T[K])")  
  ax4.set_ylabel("")
  ax4.set_title("hist_T")
  ax4.set_yscale('log')
  ax4.legend(loc='best')
  fig4.tight_layout()

  ax5.set_xlabel("Metallicity [Z]")  
  ax5.set_ylabel("")
  ax5.set_title("hist_Z")
  ax5.set_yscale('log')
  ax5.legend(loc='best')
  fig5.tight_layout()

  ax6.set_xlabel("Energy released per SN event [erg]")  
  ax6.set_ylabel("")
  ax6.set_title("hist_E")
  ax6.set_yscale('log')
  ax6.legend(loc='best')
  fig6.tight_layout()

  print("Saving : ", imgname_1)
  fig1.savefig(imgname_1)
  fig1.clf()

  print("Saving : ", imgname_2)
  fig2.savefig(imgname_2)
  fig2.clf()
 
  print("Saving : ", imgname_3)
  fig3.savefig(imgname_3)
  fig3.clf()

  print("Saving : ", imgname_4)
  fig4.savefig(imgname_4)
  fig4.clf()
 
  print("Saving : ", imgname_5)
  fig5.savefig(imgname_5)
  fig5.clf()
 
  print("Saving : ", imgname_6)
  fig6.savefig(imgname_6)
  fig6.clf()

def plot_halohist(simdir_arr, simtyp_arr, outnum_arr, cbar_arr, halonum) :

  fig1, ax1 = plt.subplots()
  fig2, ax2 = plt.subplots()

  isim  = 0 
  for simdir in simdir_arr :
    outnum_char = str(outnum_arr[isim]).zfill(5)
    outdir = simdir+'/output_' + outnum_char
    if(halonum == 0 ) :
      outfile = outdir + '/info_'+ outnum_char + '/TnH.csv'
      imgname_1 = 'hist_halo_nH.png'
      imgname_2 = 'hist_halo_T.png'
    elif(halonum > 0) : 
      outfile = outdir + '/info_'+ outnum_char + '/TnH_'+ str(halonum) + '.csv' 
      imgname_1 = 'hist_halo_nH_' + str(halonum) +'.png'
      imgname_2 = 'hist_halo_T_'  + str(halonum) +'.png'
    else :
      raise ValueError('incorrect value for halonum')

    with open(outfile, 'r') as file:  
      lines  = file.readlines()
      Nlines = len(lines)
      log_nH = np.zeros(Nlines)
      log_T  = np.zeros(Nlines)   
      iline  = 0
      for line in lines:
        log_nH[iline] = float(line.split()[0])
        log_T[iline]  = float(line.split()[1])
        if(iline < 50) :
          print(simtyp_arr[isim], iline, log_nH[iline], log_T[iline])
        if(iline >= Nlines)  :
          #print('T_arr :',  T_arr)
          break
        iline         = iline  + 1
      ax1.hist(log_nH, bins=256, color=cbar_arr[isim], label=simtyp_arr[isim], edgecolor=cbar_arr[isim], histtype='step')
      ax2.hist(log_T , bins=256, color=cbar_arr[isim], label=simtyp_arr[isim], edgecolor=cbar_arr[isim], histtype='step')
      isim = isim + 1

  ax1.set_xlabel("log nH")  
  ax1.set_ylabel("")
  ax1.set_title("hist_nH for one halo")
  ax1.set_yscale('log')
  ax1.legend(loc='best')
  fig1.tight_layout()

  ax2.set_xlabel("log T")  
  ax2.set_ylabel("")
  ax2.set_title("hist_T for one halo")
  ax2.set_yscale('log')
  ax2.legend(loc='best')
  fig2.tight_layout()

  print("Saving : ", imgname_1)
  fig1.savefig(imgname_1)
  fig1.clf()

  print("Saving : ", imgname_2)
  fig2.savefig(imgname_2)
  fig2.clf()
 
def read_sfrlog(simdir,  overwrite) :
  print("simdir :", simdir)
  mfblog_all = simdir  + "/sfrlog_all.txt" 
  if  (not  os.path.exists(mfblog_all) or overwrite ) :
    if(os.path.exists(mfblog_all) and overwrite) :
      command = "rm " + mfblog_all
      print('command :', command)
      os.system(command)

    command = "touch " + mfblog_all
    print('command :', command)
    os.system(command)

    nlogfile = 1 
    for ifile in range(1,5) :
      simlog = simdir + "/simulation_" + str(ifile) + ".log"
      mfblog = simdir + "/sfrlog_"     + str(ifile) + ".txt"
      if(os.path.exists(simlog)) :
        command = "grep -rH 'SF z=' " +  simlog + " > " + mfblog
        print('command :', command)
        print('writing :',  mfblog)
        os.system(command)
        command = "cat " + mfblog + " >> " + mfblog_all
        print('command :', command)
        os.system(command)
        nlogfile = nlogfile + 1
      else :
        simlog = simdir + "/simulation" + ".log"
        mfblog = simdir + "/sfrlog"     + ".txt"
        command = "grep -rH 'SF z=' " +  simlog + " > " + mfblog
        print('command :', command)
        print('writing :',  mfblog)
        os.system(command)
        command = "cat " + mfblog + " >> " + mfblog_all
        print('command :', command)
        os.system(command)
        nlogfile = nlogfile + 1
        break
  else : 
    print(mfblog_all, 'exists')
      
  z_arr  = []
  N_arr  = []
  nH_arr = []
  T_arr  = []
  Z_arr  = []
  dx_arr = []
  E_arr  = []

  with open(mfblog_all, 'r') as file: 
    lines  = file.readlines()     
    for line in lines:
      #SF z=   7.61010  nH=   1.650 T=  -0.904 N=     1 eps_ff=   0.13134 alpha=   0.00000
      parts = line.split()
      if(float(parts[2]) > 8.4) :
        z_arr.append  (  float(parts[2])      )
        nH_arr.append (  10**float(parts[4])      )
        T_arr.append  (  10**float(parts[6])      )
  return[nH_arr, T_arr]

def read_starfiles(simdir, outnum ,overwrite) :
  nproc = 200
  outnum_char = str(outnum).zfill(5)
  outdir = simdir+'/output_' + outnum_char
  if( os.path.exists(outdir) ) :
    read_units(simdir, outnum)
    numlines = 0 
    nH_star = []
    T_star  = []
    for iproc in range(1, nproc+1) :
      iproc_char = str(iproc).zfill(5)
      starfile   = outdir + "/stars_" + outnum_char + ".out" + iproc_char
      if(os.path.exists(starfile)) :
        with open(starfile, 'r') as file:
          lines  = file.readlines()
          nline = 0 
          for line in lines :
            if(nline > 0) :
              parts = line.split()
              ndens = float(parts[10])*scale_nH
              temp = float(parts[14])
              nH_star.append(ndens)
              T_star.append(temp)
              #print(numlines, ndesn, T)
            nline = nline + 1  
            numlines = numlines +1
      else  :
        print(starfile, "doesnt exist")       
  else :
    print(outdir, "doesnt exist") 
  return [nH_star, T_star]             








def plot_sfrlog(simdir_arr, simtyp_arr, cbar_arr, overwrite) :

  fig1, ax1 = plt.subplots()
  fig2, ax2 = plt.subplots()
  fig3, ax3 = plt.subplots()

  isim  = 0 
  for simdir in simdir_arr :
    print("simdir :", simdir)
    mfblog_all = simdir  + "/sfrlog_all.txt" 
    if  (not  os.path.exists(mfblog_all) or overwrite ) :
      if(os.path.exists(mfblog_all) and overwrite) :
        command = "rm " + mfblog_all
        print('command :', command)
        os.system(command)

      command = "touch " + mfblog_all
      print('command :', command)
      os.system(command)

      nlogfile = 1 
      for ifile in range(1,5) :
        simlog = simdir + "/simulation_" + str(ifile) + ".log"
        mfblog = simdir + "/sfrlog_"     + str(ifile) + ".txt"
        if(os.path.exists(simlog)) :
          command = "grep -rH 'SF z=' " +  simlog + " > " + mfblog
          print('command :', command)
          print('writing :',  mfblog)
          os.system(command)
          command = "cat " + mfblog + " >> " + mfblog_all
          print('command :', command)
          os.system(command)
          nlogfile = nlogfile + 1
        else :
          simlog = simdir + "/simulation" + ".log"
          mfblog = simdir + "/sfrlog"     + ".txt"
          command = "grep -rH 'SF z=' " +  simlog + " > " + mfblog
          print('command :', command)
          print('writing :',  mfblog)
          os.system(command)
          command = "cat " + mfblog + " >> " + mfblog_all
          print('command :', command)
          os.system(command)
          nlogfile = nlogfile + 1
          break
    else : 
      print(mfblog_all, 'exists')

    with open(mfblog_all, 'r') as file:  
      lines  = file.readlines()
      Nlines = len(lines)
      z_arr  = []
      N_arr  = []
      nH_arr = []
      T_arr  = []
      Z_arr  = []
      dx_arr = []
      E_arr  = []
      for line in lines:
        #SF z=   7.61010  nH=   1.650 T=  -0.904 N=     1 eps_ff=   0.13134 alpha=   0.00000
        parts = line.split()
        if(float(parts[2]) > 8.4) :
          z_arr.append  (  float(parts[2])      )
          nH_arr.append (  float(parts[4])      )
          T_arr.append  (  float(parts[6])      )
      ax1.hist(z_arr , bins=256, density=True, color=cbar_arr[isim], label=simtyp_arr[isim], edgecolor=cbar_arr[isim], histtype='step')
      ax2.hist(nH_arr, bins=256, density=True, color=cbar_arr[isim], label=simtyp_arr[isim], edgecolor=cbar_arr[isim], histtype='step')
      ax3.hist(T_arr , bins=256, density=True, color=cbar_arr[isim], label=simtyp_arr[isim], edgecolor=cbar_arr[isim], histtype='step')
      isim = isim + 1
        #print("iline :", parts[0], parts[1] , parts[2], parts[4], parts[5], parts[6], parts[7], parts[8], parts[9], parts[10]  )

  imgname_1 = "hist_SFz.png"
  imgname_2 = "hist_SFnH.png"
  imgname_3 = "hist_SFT.png"
  
  ax1.set_xlabel("Redshift(z)")  
  ax1.set_ylabel("")
  ax1.set_title("Redshift of star formation")
  ax1.set_yscale('log')
  ax1.legend(loc='best')
  fig1.gca().invert_xaxis()
  fig1.tight_layout()

  ax2.set_xlabel(" log10 ( SF_nH[cc] )")  
  ax2.set_ylabel("")
  ax2.set_title("Number density of star forming cell")
  #ax2.set_yscale('log')
  ax2.legend(loc='best')
  fig2.tight_layout()

  ax3.set_xlabel("log10 (SF_T[K] ) ")  
  ax3.set_ylabel("")
  ax3.set_title("Temperature of star forming cell")
  #ax2.set_yscale('log')
  ax3.legend(loc='best')
  fig3.tight_layout()

  print("Saving : ", imgname_1)
  fig1.savefig(imgname_1)
  fig1.clf()

  print("Saving : ", imgname_2)
  fig2.savefig(imgname_2)
  fig2.clf()

  print("Saving : ", imgname_3)
  fig3.savefig(imgname_3)
  fig3.clf()
