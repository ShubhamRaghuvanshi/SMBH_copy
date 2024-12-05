import cfuncs
import matplotlib.pyplot as plt    

simdir_0='/u/sraghu/sim/SHUB/b256/sn_eq'
simdir_1='/ptmp/mpa/sraghu/tests_nosink/sn_eq'
simdir_2='/ptmp/mpa/sraghu/tests_nosink/sn_neq'
simdir_3='/ptmp/mpa/sraghu/tests_nosink/sn_neq_18'
simdir_4='/ptmp/mpa/sraghu/tests_nosink/sn_neq_lightstar'

simtyp_arr    = ['star']

simdir_arr  = [simdir_1]
cbar_arr    = ['k']

#normailsed histogram 
#check temp, maybe outputted before upper limit correction 
#  

#cfuncs.plot_mfblog(simdir_arr, simtyp_arr, cbar_arr, False)
#cfuncs.plot_sfrlog(simdir_arr, simtyp_arr, cbar_arr, False)

outnum_arr = [10]
isim =0

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()

for simdir in simdir_arr :
  outnum = outnum_arr[isim]
  starlog = cfuncs.read_starfiles(simdir, outnum ,False) 
  #starlog = cfuncs.read_sfrlog(simdir, False)
  nH_star = starlog[0] 
  T_star  = starlog[1]
  ax1.scatter(nH_star, T_star, s=5, label = simtyp_arr[isim])
  ax1.set_xlim(1e-7,1e4)
  ax1.set_ylim(1, 1e9)
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  # Add labels and title
  ax1.set_xlabel('nH[cc] of star forming cells')
  ax1.set_ylabel('T[K] of star forming cells')
  ax1.set_title('nH vs T for star forming cells')
  ax1.legend(loc='best')
  imgname = "scatter_star_" + simtyp_arr[isim] + ".png"
  print("Saving : ", imgname)
  fig1.savefig(imgname)
  fig1.clf()

  ax2.hist(nH_star, bins=256, density=True, color=cbar_arr[isim], label=simtyp_arr[isim], edgecolor=cbar_arr[isim], histtype='step')
  ax3.hist(T_star , bins=256, density=True, color=cbar_arr[isim], label=simtyp_arr[isim], edgecolor=cbar_arr[isim], histtype='step')
  isim = isim + 1
  
imgname_2 = "hist_SFnH.png"
imgname_3 = "hist_SFT.png"

ax2.set_xlabel(" log10 ( SF_nH[cc] )")  
ax2.set_ylabel("")
ax2.set_title("Number density of star forming cell")
ax2.set_xscale('log')
ax2.legend(loc='best')
fig2.tight_layout()

ax3.set_xlabel("log10 (SF_T[K] ) ")  
ax3.set_ylabel("")
ax3.set_title("Temperature of star forming cells")
ax2.set_xscale('log')
ax3.legend(loc='best')
fig3.tight_layout()

print("Saving : ", imgname_2)
fig2.savefig(imgname_2)
fig2.clf()

print("Saving : ", imgname_3)
fig3.savefig(imgname_3)
fig3.clf()


#halonum = 1
#outnum_arr = [8,9]
#cfuncs.plot_halohist(simdir_arr, simtyp_arr, outnum_arr, cbar_arr, halonum)
