import os 
import struct
import math 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import yt_funcs
#plt.rcParams['text.usetex'] = True
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
mP = 1.6726219e-24
X = 0.76
h = 0.703100835997471
Om = 0.276
G_cgs = 6.674e-8
H0 = 70.3/h # h(km/s)/Mpc

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

#code units 
Z         = -1.0
scale_M   = -1.0
scale_L   = -1.0
scale_T   = -1.0
scale_d   = -1.0
scale_nH  = -1.0
scale_vel = -1.0
def give_units(aexp, boxlength):
  global Z
  global scale_M
  global scale_L
  global scale_T
  global scale_d
  global scale_nH
  global scale_vel

  Lx          = boxlength # h^-1 mpc
  Ly          = Lx
  Lz          = Lx
  V           = Lx*Ly*Lz
  M_V         = rho_crit0*Om*V

  scale_L     = (aexp*boxlength/h)*unit_L
  scale_M     = (M_V/h)*unit_M
  scale_d     = (rho_crit0*Om*h*h/( aexp*aexp*aexp )) *unit_rho
  scale_T     = aexp*aexp/( h * 1e5) *unit_L
  scale_nH    = X*scale_d/mH
  scale_vel   = (scale_L/scale_T)*unit_vel

def read_units(simdir, outnum):
  global Z
  global scale_M
  global scale_L
  global scale_T
  global scale_d
  global scale_nH
  global scale_vel
  global boxlen

  outnumn_char = str(outnum).zfill(5)
  outdir = simdir+'/output_' + outnumn_char
 
  aexp    = None
  scale_L = None 
  scale_d = None
  scale_T = None 
  with open(outdir+'/info_'+outnumn_char+'.txt' , 'r') as file:
    for line in file:
      if 'aexp' in line:
        aexp = float(line.split()[2])
      if 'unit_l' in line:
        scale_L = float(line.split()[2])
      if 'unit_d' in line:
        scale_d = float(line.split()[2])
      if 'unit_t' in line:
        scale_T = float(line.split()[2])
        break
  if(aexp == None or scale_L == None or scale_d == None or scale_T == None ) :
    raise ValueError(' units not read ')
  Z           = 1.0/aexp - 1.0
  scale_nH    = X/mH * scale_d
  scale_M     = scale_d*scale_L*scale_L*scale_L
  scale_vel   = (scale_L/scale_T)*unit_vel
 
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
  outdir=simdir+"/sinkdata" 
  sinknum_char = str(sinknum).zfill(5)
  outfile=outdir+"/sink_"+sinknum_char+".csv"
  #print(outfile) 
  prop_arr = []
  with open(outfile, 'r') as sinkfile:
    for line in sinkfile:
      propval = float(line.split(',')[propnum])      
      prop_arr.append(propval) 
  return prop_arr


def give_massive_sinkid(simdir):

  sinkdir   = simdir+'/sinkdata/'
  csv_files = [file for file in os.listdir(sinkdir) if file.endswith('.csv')]
  Nfiles = len(csv_files)

  Nsink = Nfiles
  massive_id = 1
  msinkmax = 0
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
  plt.xlabel('Redshift')
  plt.ylabel('Msink[Msun]')
  plt.title('Most Massive sink')
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
  plt.legend()
  plt.gca().invert_xaxis()
  plt.tight_layout()
  imgname = imgname + ".png"
  print("Saving :" , imgname)
  plt.savefig(imgname)

def compare_sinkdmdt(simdir_arr, label_arr, sinkid_arr) :
  plt.xlabel('Redshift')
  plt.ylabel('dmdsink_dt[Msun/Myr]')
  plt.title('Sink accretion rate')
  plt.yscale('log')

  isink =0
  imgname = "dmdt_z" 
  for simdir in simdir_arr :
    sinkid      = sinkid_arr[isink]
    plabel      = label_arr[isink]
    sink_dmBHdt = sink_prop_return(simdir, sinkid, 9)
    aexp_arr    = sink_prop_return(simdir, sinkid, 10)

    Ntime        = len(aexp_arr)
    Z_arr        = []

    for itime in range(0, Ntime):
      aexp = aexp_arr[itime]
      z    = 1.0/aexp -1.0
      give_units(aexp)
      sink_dmBHdt[itime]  = (scale_M/Msun)*(Myr/scale_T)*sink_dmBHdt[itime]
      Z_arr.append(z)
    plt.plot(Z_arr, sink_dmBHdt,  label=plabel, color=cbar[isink] ,marker='o', linestyle='-')
    imgname = imgname + "_" + plabel
    isink   = isink + 1
  plt.legend()
  plt.gca().invert_xaxis()
  plt.tight_layout()
  imgname = imgname + ".png"
  plt.savefig(imgname)

def plot_sinksmasses(simdir, labl, Nsinks) :

  for isink in range(1, Nsinks + 1):
    sinkmass    = sink_prop_return(simdir, isink, 1)
    aexp_arr    = sink_prop_return(simdir, isink, 10)
  
    Ntime       = len(aexp_arr)
    Z_arr       = [] 

    for itime in range(0, Ntime):
      aexp = aexp_arr[itime]
      z  = 1.0/aexp -1.0
      give_units(aexp)
      sinkmass[itime]   = (scale_M/Msun) * sinkmass[itime]
      Z_arr.append(z)

    plt.plot(Z_arr, sinkmass)
  plt.xlabel('Redshift')
  plt.ylabel('Msink[Msun')
  plt.title('')
  plt.yscale('log')
  plt.gca().invert_xaxis()
  plt.tight_layout()
  imgname = "sinksmasses_z_" + labl  + ".png"
  plt.savefig(imgname)


def plot_sinkdmdt(simdir, labl, sinkid) :

  sink_dmBHdt = sink_prop_return(simdir, sinkid, 9)
  aexp_arr    = sink_prop_return(simdir, sinkid, 10)
  
  Ntime       = len(aexp_arr)
  Z_arr       = [] 

  for itime in range(0, Ntime):
    aexp = aexp_arr[itime]
    z  = 1.0/aexp -1.0
    give_units(aexp)
    sink_dmBHdt[itime]  = (scale_M/Msun)*(Myr/scale_T)*sink_dmBHdt[itime]
    Z_arr.append(z)

  plt.plot(Z_arr, sink_dmBHdt, label=labl, color=cbar[0] ,marker='o', linestyle='-')
  plt.xlabel('Redshift')
  plt.ylabel('dmdsink_dt[Msun/Myr]')
  plt.title('')
  plt.yscale('log')
  plt.legend()
  plt.gca().invert_xaxis()
  plt.tight_layout()
  imgname = "dmdt_z_" + labl  + ".png"
  plt.savefig(imgname)


def compare_sinkmasstot(simdir_arr, label_arr, sinkid_arr) :
  plt.xlabel('Redshift')
  plt.ylabel('Msink[Msun]')
  plt.title('Sinkmasses')
  plt.yscale('log')

  isink =0
  imgname = "sinkmasstot" 
  for simdir in simdir_arr :
    print(simdir)
    sinkid   = sinkid_arr[isink]
    plabel   = label_arr[isink]
    sinkmasstot = sink_prop_return(simdir, sinkid, 11)
    aexp_arr    = sink_prop_return(simdir, sinkid, 10)

    Ntime        = len(aexp_arr)
    Z_arr        = []

    for itime in range(0, Ntime):
      aexp = aexp_arr[itime]
      z    = 1.0/aexp -1.0
      give_units(aexp)
      sinkmasstot[itime]   = (scale_M/Msun) * sinkmasstot[itime]
      Z_arr.append(z)
    plt.plot(Z_arr, sinkmasstot,  label=plabel, color=cbar[isink] ,marker='o', linestyle='-')
    imgname = imgname + "_" + plabel
    isink   = isink + 1
  plt.legend()
  plt.gca().invert_xaxis()
  plt.tight_layout()
  imgname = imgname + ".png"
  print("Saving : ", imgname)
  plt.savefig(imgname)

#plot total star mass vs total sink mass w.r.t redshift
def star_vs_sink_masstot(simdir, labl, sinkid):
  aexp_arr     = sink_prop_return(simdir, sinkid, 10)
  sinkmasstot  = sink_prop_return(simdir, sinkid, 11)
  starmasstot  = sink_prop_return(simdir, sinkid, 12)
  Ntime        = len(aexp_arr)
  Z_arr        = [] 

  for itime in range(0, Ntime):
    aexp = aexp_arr[itime]
    z    = 1.0/aexp -1.0
    give_units(aexp)
    sinkmasstot[itime]   = (scale_M/Msun) * sinkmasstot[itime]
    starmasstot[itime]   = (scale_M/Msun) * starmasstot[itime]
    Z_arr.append(z)


  plt.plot(Z_arr, sinkmasstot,  label='M_sinktot', color=cbar[0], linestyle='-')
  plt.plot(Z_arr, starmasstot,  label='M_startot', color=cbar[1], linestyle='-')
  plt.gca().invert_xaxis()  
  plt.xlabel('Redshift')
  plt.yscale('log')
  plt.ylabel('Mass[Msun]')
  plt.legend()
  plt.title('')
  imgname = "sinkvsstar_masstot_"+ labl +".png"
  plt.savefig(imgname)


#give file with most number of lines
def give_longest_file(simdir) :
  max_lines    = 0
  max_lines_id = 0

  sinkdir   = simdir+'/sinkdata/'
  csv_files = [file for file in os.listdir(sinkdir) if file.endswith('.csv')]
  Nfiles = len(csv_files)

  for iout in range(1, Nfiles+1):     
    outnumn_char = str(iout).zfill(5)
    sinkfile = simdir+'/sinkdata/sink_' + outnumn_char + ".csv"
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

def plot_halodens(simdir, outnum, halonum, plotsinks, plotcircle, plotstars, overwrite_densbin, overwrite_sink, overwrite_star, boxlen_comov, dir):
  read_units(simdir, outnum)
  print("scale_M : ", scale_L)
  print("scale_L : ", scale_M)
  print("scale_T : ", scale_T)
  print("Z       : ", Z)

  outnum_char    = str(outnum).zfill(5)
  outdir        = simdir+'/output_' + outnum_char
  halofile       = outdir + '/info_'+ outnum_char + "/halodata.csv"

  if  (not os.path.exists(halofile) ) : 
    yt_funcs.write_halos(simdir, outnum)
  else :
    print(halofile, "already exists")  

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
  rvir = (halo_rvir/1000.0) /boxlen_comov

  xmin = max(halo_posx - 2.0*rvir,0)
  xmax = min(halo_posx + 2.0*rvir,1)
  ymin = max(halo_posy - 2.0*rvir,0)
  ymax = min(halo_posy + 2.0*rvir,1)
  zmin = max(halo_posz - 2.0*rvir,0)
  zmax = min(halo_posz + 2.0*rvir,1)
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

  print("xmin, xmax, ymin, ymax : ", xmin, xmax, ymin, ymax)
  print("xbox, ybox, zbox :", xbox, ybox, zbox)

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

  print("************************8")
  print("xmin, xmax, ymin, ymax : ", xmin, xmax, ymin, ymax)
  print("xbox, ybox, zbox :", xbox, ybox, zbox)
  print("rbox, rvir: ",  rbox, rvir)


  lbox = math.sqrt(xbox*xbox + ybox*ybox + zbox*zbox)
  

  #xmin = 0
  #xmax = 0.5
  #ymin = 0
  #ymax = 0.5
  #zmin = 0
  #zmax = 1
  #extent  = [xmin, xmax, ymin, ymax ]
  extent  = [(xmin-halo_posx)*boxlen_comov, (xmax-halo_posx)*boxlen_comov, (ymin-halo_posy)*boxlen_comov, (ymax-halo_posy)*boxlen_comov ]
  lmax = 13
  nx = int ( (2**lmax)*xbox - 1)
  ny = nx

  print("nx, ny : ", nx, ny)

  densbin = outdir + '/info_'+ outnum_char + "/halodens_" + str(halonum) + dir
  if ( (not os.path.exists(densbin) ) or (os.path.exists(densbin) and overwrite_densbin == True) ) :
    command = ("./amr2map -inp " + outdir + " -out " + densbin + " -dir " + dir + " -xmi " + str(xmin) + " -xma " + str(xmax) + " -ymi " + 
               str(ymin) + " -yma " + str(ymax) + " -zmi " + str(zmin) + " -zma " + str(zmax) + " -typ 1 " + "-lma " + str(lmax) + 
               " -nx " + str(nx) + " -ny " + str(ny))
    print("command : ", command)
    os.system(command)
  else :
    print(densbin, "already exists")  

  dens_arr  = read_amr2map_binary(densbin)

  if(Nx < 0 or Ny < 0 or Nz < 0):
    print("Nx, Ny, Nz", Nx, Ny, Nz)
    raise ValueError("Invalid grid size.")

  dens_arr = dens_arr*scale_d
  dens_reshape = dens_arr.reshape(Nx, Ny)
  dens_trans   = dens_reshape
  
  #plt.figure(figsize=(10, 8))

  if(plotstars == True) :
    mstar = []
    xstar = []
    ystar = []
    zstar = []
    nstar=0
    starfile = outdir + "/starpos.csv"
    if ( (not os.path.exists(starfile) ) or (os.path.exists(starfile) and overwrite_star == True) ) :  
      command = "./read_star -inp " + outdir + " -out " + starfile  
      os.system(command)
    else :
      print(starfile, "already exists")

    with open(starfile, 'r') as file:
      for line in file:
        mass = float(line.split(',')[1])
        posx = float(line.split(',')[2]) 
        posy = float(line.split(',')[3]) 
        posz = float(line.split(',')[4]) 
        if( posx>xmin and posx<xmax and posy>ymin and posy<ymax and posz>zmin and posz<zmax ) :
          star_posx = (posx - halo_posx)*boxlen_comov 
          star_posy = (posy - halo_posy)*boxlen_comov
          star_posz = (posz - halo_posz)*boxlen_comov
          mstar.append ( mass )
          xstar.append( star_posx )
          ystar.append( star_posy ) 
          zstar.append( star_posz )
          #print(nsink, sink_posx, sink_posy, sink_posz)
          nstar=nstar+1
    if(nstar > 0) :    
      if(dir == 'z') :  
        plt.scatter(xstar, ystar, marker='.', color="yellow" ,s=1)
        #plt.scatter(xstar, ystar, marker='.', c=mstar, cmap='viridis' ,s=1)
      elif (dir == 'x') :
        plt.scatter(ystar, zstar, marker='.', color="yellow" ,s=1)
        #plt.scatter(ystar, zstar, marker='.', c=mstar, cmap='viridis' ,s=1)  
      elif (dir == 'y') :
        plt.scatter(xstar, zstar, marker='.', color="yellow" ,s=1)
        #plt.scatter(xstar, zstar, marker='.', c=mstar, cmap='viridis',s=1)  

    else :
      print("No stars in the given region")

  if(plotsinks == True) :
    xsink = []
    ysink = []
    zsink = []
    nsink=0
    sinkfile = outdir + "/sinkpos.csv"
    if ( (not os.path.exists(sinkfile) ) or (os.path.exists(sinkfile) and overwrite_sink == True) ) :  
      command = "./read_sink -inp " + outdir + " -out " + sinkfile  
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
          #print(nsink, sink_posx, sink_posy, sink_posz)
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
  


  plt.imshow(dens_trans, cmap='viridis', norm=LogNorm(),extent = extent, origin='lower' ,label=f'Z = {Z:.2f}' )
  plt.text(0.95, 0.95, f'Z = {Z:.2f}', ha='right', va='top', color='black',fontweight='bold', fontsize=14, transform=plt.gca().transAxes)

  xtext = 1.5*rvir/lbox + 0.5
  ytext = 2.0*rvir/lbox + 0.5
  rgal = halo_rvir
  plt.text(xtext, ytext, f'{rgal:.2f} kpccm/h', ha='right', va='top', color='black',fontweight='bold', fontsize=14, transform=plt.gca().transAxes)

  if(plotcircle == True ):
    
    cen_x = 0
    cen_y = 0
    rad   = rvir*boxlen_comov 
    print("circle : ", cen_x, cen_y, rad)
    circle = plt.Circle((cen_x, cen_y), rad , edgecolor='red', facecolor='none', linewidth=2)
    plt.gca().add_patch(circle)
  
  cbar = plt.colorbar()
  cbar.set_label(r'density$ [ \textrm{gm}/\textrm{cm}^3 ] $ ')
  if(dir == 'z') :
    plt.xlabel(r' $\textrm{X} \left[  \textrm{Mpccm/h} \right]$ ')
    plt.ylabel(r' $\textrm{Y} \left[  \textrm{Mpccm/h} \right] $')
  elif (dir == 'x') :
    plt.xlabel(r' $\textrm{Y} \left[  \textrm{Mpccm/h} \right]$ ')
    plt.ylabel(r' $\textrm{Z} \left[  \textrm{Mpccm/h} \right] $')
  elif (dir == 'y') :
    plt.xlabel(r' $\textrm{X} \left[  \textrm{Mpccm/h} \right]$ ')
    plt.ylabel(r' $\textrm{Z} \left[  \textrm{Mpccm/h} \right] $')

  figname =  'dens_' + str(halonum) + dir + '.png'
  print("saving :", figname)  
  plt.savefig(figname)


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
          if(sink_dr <= rvir2)  :
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

