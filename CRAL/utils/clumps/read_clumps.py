import cfuncs

simdir              = "/u/sraghu/sim/SHUB/b128/smbhstar"
outnum              = 15
halonum             = 1
plotsinks           = True
plotstars           = True
plotcircle          = True
overwrite_densbin   = True
overwrite_sink      = False
overwrite_star      = False
boxlen_comov        = 10
dir                 = 'z'
#cfuncs.plot_halodens(simdir, outnum, halonum, plotsinks, plotcircle, plotstars, overwrite_densbin, overwrite_sink, overwrite_star, boxlen_comov, dir)


plotsinkvsstar = False
overwrite_halo = False
cfuncs.sink_star_halo_mass(simdir, outnum, overwrite_sink, overwrite_star, overwrite_halo, plotsinkvsstar ,boxlen_comov)

#simdir_1="/u/sraghu/sim/SHUB/b128/star"
#simdir_2="/u/sraghu/sim/SHUB/b128/smbh"
#simdir_3="/u/sraghu/sim/SHUB/b128/smbhstar"
#simdir_4="/u/sraghu/sim/SHUB/b128/smbhfeed"

'''
#compare total mass in sinks 
simdir_arr = [simdir_2, simdir_3, simdir_4]
simtyp_arr = ["sink_nostar", "sink_star", "sink_feed"]
sinkid_2   = cfuncs.give_longest_file(simdir_2)
sinkid_3   = cfuncs.give_longest_file(simdir_3)
sinkid_4   = cfuncs.give_longest_file(simdir_4) 
sinkid_arr = [sinkid_2, sinkid_3, sinkid_4]
cfuncs.compare_sinkmasstot(simdir_arr, simtyp_arr, sinkid_arr) 
'''


'''
#compare most massive sink mass
simdir_arr = [simdir_2, simdir_3, simdir_4]
simtyp_arr = ["sink_nostar", "sink_star", "sink_feed"]
sinkid_2   = cfuncs.give_massive_sinkid(simdir_2)
sinkid_3   = cfuncs.give_massive_sinkid(simdir_3)
sinkid_4   = cfuncs.give_massive_sinkid(simdir_4)
sinkid_arr = [sinkid_2, sinkid_3, sinkid_4]
cfuncs.compare_sinkmasses(simdir_arr, simtyp_arr, sinkid_arr)
'''

#plot total star vs sink mass

'''
#plot masses of Nsink number of sinks
simdir="/u/sraghu/sim/SHUB/b128/smbhstar"
Nsink=10
for isink in range (1, Nsink+1) :
  #sink_posx    = cfuncs.sink_prop_return(simdir, isink, 2)
  sinkmasses   = cfuncs.sink_prop_return(simdir, isink, 1)
  #aexp_arr     = cfuncs.sink_prop_return(simdir, isink, 22)
  aexp_arr     = cfuncs.sink_prop_return(simdir, isink, 10)
  Ntime        = len(aexp_arr)
  Z_arr        = [] 

  for itime in range(0, Ntime):
    aexp = aexp_arr[itime]
    z    = 1.0/aexp -1.0
    cfuncs.give_units(aexp)
    sinkmasses[itime]   = (cfuncs.scale_M/cfuncs.Msun) * sinkmasses[itime]
    #sink_posx[itime]   = (cfuncs.scale_L/cfuncs.Kpc) * sink_posx[itime]
    Z_arr.append(z)

  plt.plot(Z_arr, sinkmasses,  label='',  linestyle='-')
plt.gca().invert_xaxis()  
plt.xlabel('Redshift')
#plt.xscale('log')
plt.yscale('log')
plt.ylabel('Msink[Msun]')
#plt.ylabel('Xpos')
plt.title('')
plt.savefig("sinkmassall_smbhstar.png")
'''


'''
#plot total star mass vs total sink mass w.r.t redshift 
simdir="/u/sraghu/sim/SHUB/b128/smbhstar"
sinkid       = 1 
sinknum_char = str(sinkid).zfill(5)

aexp_arr      = cfuncs.sink_prop_return(simdir, sinkid, 10)
sinkmasstot   = cfuncs.sink_prop_return(simdir, sinkid, 11)
starmasstot   = cfuncs.sink_prop_return(simdir, sinkid, 12)
Ntime        = len(aexp_arr)
Z_arr        = [] 

for itime in range(0, Ntime):
  aexp = aexp_arr[itime]
  z    = 1.0/aexp -1.0
  cfuncs.give_units(aexp)
  sinkmasstot[itime]   = (cfuncs.scale_M/cfuncs.Msun) * sinkmasstot[itime]
  starmasstot[itime]   = (cfuncs.scale_M/cfuncs.Msun) * starmasstot[itime]
  Z_arr.append(z)


plt.plot(Z_arr, sinkmasstot,  label='M_sinktot',  linestyle='-')
plt.plot(Z_arr, starmasstot,  label='M_startot',  linestyle='-')
plt.gca().invert_xaxis()  
plt.xlabel('Redshift')
#plt.xscale('log')
plt.yscale('log')
plt.ylabel('Mass[Msun]')
#plt.ylabel('Xpos')
plt.legend()
plt.title('')
plt.savefig("sinkstar_masstot_128_croblisk.png")
'''

'''
NTask = 64
outnum=15
cfuncs.read_units(simdir, outnum)
print("")
print("------------------------------------------------")
print("Z :",         cfuncs.Z )
print("scale_L : ",  cfuncs.scale_L)
print("scale_d : ",  cfuncs.scale_d)
print("scale_T : ",  cfuncs.scale_T)
print("scale_nH : ", cfuncs.scale_nH)
print("scale_M : ", cfuncs.scale_M)
print("scale_M : ", "{:e}".format(cfuncs.scale_M*cfuncs.h/cfuncs.Msun), "Msun/h")
print("------------------------------------------------")
print("")

cfuncs.merge_clump_props(simdir, outnum, NTask )
cfuncs.merge_halo_props(simdir, outnum, NTask )

clump_masses = cfuncs.read_clump_prop(simdir, outnum, 10)
clump_denses = cfuncs.read_clump_prop(simdir, outnum, 9)
clump_indexs = cfuncs.read_clump_prop(simdir, outnum, 0)
clump_host   = cfuncs.read_clump_prop(simdir, outnum, 12)
halos_indexs = cfuncs.read_halo_prop(simdir, outnum, 0)
halos_indexs = cfuncs.read_halo_prop(simdir, outnum, 6)
sink_masses  = cfuncs.read_sink_prop(simdir, outnum, 1)

print("")
print("------------------------------------------------")
print("Number of clumps read : ", len(clump_masses))
print("Number of halos read : ", len(halos_indexs))
print("------------------------------------------------")
print("")


clump_masses = [cmass * cfuncs.scale_M/cfuncs.Msun for cmass in clump_masses]
Nclump = len(clump_masses)
Nhalo  = len(halos_indexs)
Nsink  = len(sink_masses)

Ncand =0 
Nhost=0
for iclump in range(0, Nclump) :
  for ihalo in range(0, Nhalo):
    if(  halos_indexs[ihalo] == clump_indexs[iclump]  ) :
      if( clump_masses[iclump] > cfuncs.mass_halo_AGN  and clump_masses[iclump] > cfuncs.mass_clump_AGN  and clump_denses[iclump] > cfuncs.n_star/cfuncs.scale_nH) :
        Ncand = Ncand + 1
      if(clump_host[iclump]==1) :
        Nhost = Nhost + 1
      break

print("Z :",         cfuncs.Z )
print("Number of candicate clump for sink : ",Ncand)        
print("Number of sinks : ",Nsink)
print("Number of hosts : ",Nhost)
print("------------------------------------------------")
print("")
'''
