import numpy as np
import math 

#constants
pi = 3.14159265

Kilo = 1e3
million = 1e6

AU_cgs = 1.495978707e13
Mpc = 3.0856e24 #cm
pc = 3.0856e18
Kpc = 3.0856e21

Msun = 1.988e33 #g
yr = 31536000.0 #sec
Myr = yr*million

G_cgs = 6.674e-8
c_cgs = 29979245800.0

mH = 1.6735575e-24
mP = 1.6726219e-24

kB    = 1.3806e-16

#units
unit_M = Msun
unit_L = Mpc
unit_T = 1.0
unit_vol = unit_L*unit_L*unit_L
unit_vel = unit_L/unit_T
unit_E = unit_M*unit_vel*unit_vel
unit_rho = unit_M/unit_vol

print('----------------Units in cgs-------------------------')
print("Unit_M  =", "{:e}".format(unit_M)   )
print("Unit_L  =", "{:e}".format(unit_L)   )
print("Unit_T  =", "{:e}".format(unit_T)   )
print("Unit_velocity  =", "{:e}".format(unit_vel)   )
print("Unit_E  =", "{:e}".format(unit_E)   )
print("Unit_rho  =", "{:e}".format(unit_rho)   )

#cosmological constants
h = 0.703100835997471
H0 = 70.3/h # h(km/s)/Mpc
Om = 0.276
a0 = 1
H0_sys = (H0*Kilo/Mpc)*unit_T*100.0

G = G_cgs*unit_M*unit_T*unit_T/( unit_L*unit_L*unit_L )
rho_crit0 = (3.0*H0_sys*H0_sys)/(8.0*pi*G)

print("G  = ", "{:e}".format(G) )
print("rho_crit0 :", "{:e}".format(rho_crit0) )
print('----------------------------------------------------')
print('')

print('----------------Comoving Units-------------------------')

# dimentional quantities
boxlength = 10.0 #h^-1 mpc
Lx = boxlength # h^-1 mpc
Ly = Lx
Lz = Lx
V = Lx*Ly*Lz

levelmin = 7
levelmax = 15
Nx = 2**levelmin
Ny = Nx
Nz = Nx
Ngrid = Nx*Ny*Nz
M_V = rho_crit0*Om*V
M_cell = M_V/Ngrid
dx_cell = Lx/Nx
V_cell = V/Ngrid 

#aexp = 0.0770442042993282
#Z = 1.0/aexp - 1.0
Z = 6.0
aexp = 1.0/(1.0 + Z)

h = 0.703000030517578
H0 = 100.0 # h(km/s)/Mpc
H0_sys = (H0*Kilo/Mpc)*unit_T*100.0
X = 0.7599999904632


print('Redshift = ',"{:e}".format(Z))
print('aexp = ',"{:e}".format(aexp))
print('Nx, Ny, Nz, Ngrid : ', Nx, Ny, Nz, Ngrid)
print('MTotGrid = ',"{:e}".format(M_V), "h^-1 Msun ")
print('MTotGrid = ',"{:e}".format(M_V*unit_M/h), "g ")
print('')

M_cell = M_V/Ngrid
dx_cell = Lx/Nx
V_cell = V/Ngrid 
print('M_cell at level ', levelmin, ' = ', "{:e}".format(M_cell) , "h^-1 Msun " )
print('dx_cell at level ', levelmin, ' = ', "{:e}".format(dx_cell) , "h^-1 Mpc     ", "{:e}".format(dx_cell*1000.0) , "h^-1 Kpc " )
print('V_cell at level ', levelmin, ' = ', "{:e}".format(V_cell) , "h^-3 Mpc^3 " )
print('')

fac                  = 2**(levelmax - levelmin)
fac3                 = fac*fac*fac
M_cell_levelmax      = M_cell/fac3
dx_cell_levelmax     = dx_cell/fac
V_cell_levelmax      = V_cell/fac3

print('fac, fac3 : ', fac, fac3)
print('M_cell at level ', levelmax, ' = ', "{:e}".format(M_cell_levelmax) , "h^-1 Msun " )
print('dx_cell at level ', levelmax, ' = ', "{:e}".format(dx_cell_levelmax) , "h^-1 Mpc     ", "{:e}".format(dx_cell_levelmax*1000.0) , "h^-1 Kpc " )
print('V_cell at level ', levelmax, ' = ', "{:e}".format(V_cell_levelmax) , "h^-3 Mpc^3 " )
print('----------------------------------------------------')
print('')

'''
print('----------------Code Units-------------------------')
scale_L = (aexp*boxlength/h)*unit_L
scale_V = scale_L*scale_L*scale_L
scale_M = (M_V/h)*unit_M
scale_d = (rho_crit0*Om*h*h/( aexp*aexp*aexp )) *unit_rho
scale_T = aexp*aexp/( H0_sys * 1e5) *unit_L
scale_nH = X*scale_d/mH
scale_vel = (scale_L/scale_T)*unit_vel

print("scale_M  =", "{:e}".format(scale_M)   )
print("scale_L  =", "{:e}".format(scale_L)   )
print("scale_T  =", "{:e}".format(scale_T)   )
print("scale_d  =", "{:e}".format(scale_d)   )
print("scale_vel  =", "{:e}".format(scale_vel)   )
print("scale_nH  =", "{:e}".format(scale_nH)   )
print('')

print('M_cell at level ', levelmin, ' = ', "{:e}".format(M_cell*unit_M/scale_M)  )
print('dx_cell at level ', levelmin, ' = ', "{:e}".format(dx_cell*unit_L/scale_L)  )
print('V_cell at level ', levelmin, ' = ', "{:e}".format(V_cell*unit_vol/scale_V) )
print('')

print('M_cell at level ', levelmax, ' = ', "{:e}".format(M_cell_levelmax*unit_M/scale_M)  )
print('dx_cell at level ', levelmax, ' = ', "{:e}".format(dx_cell_levelmax*unit_L/scale_L)  )
print('V_cell at level ', levelmax, ' = ', "{:e}".format(V_cell_levelmax*unit_vol/scale_V) )

print('----------------------------------------------------')
print('')
'''

print('----------------Sink and Resolution Numbers-------------------------')
J = 0.5
T = 10000.0
ktotwothirds = (5.0*kB/(G_cgs*mH)) * pow(3.0/(4.0*pi), 3)
ktoonethird  = math.sqrt(ktotwothirds)
dx_cell_levelmax_cgs = dx_cell_levelmax*Mpc/h
rho_Jeans_cgs = ktotwothirds*T/(dx_cell_levelmax_cgs*dx_cell_levelmax_cgs)
rho_crit_cgs  = J*J*rho_Jeans_cgs
R_Jeans_cgs   = ktoonethird * math.sqrt(T/rho_Jeans_cgs)
R_Jeans_cgs   = R_Jeans_cgs/pow(4.0*pi/3.0, 1.0/3.0)
n_Jeans_cgs   = rho_Jeans_cgs/mH


#M_Jeans = rho_Jeans_cgs*dx_cell_levelmax_cgs*dx_cell_levelmax_cgs*dx_cell_levelmax_cgs

M_Jeans = rho_Jeans_cgs* R_Jeans_cgs*R_Jeans_cgs*R_Jeans_cgs
M_Jeans = M_Jeans/Msun

M_Sink  = (1.0 - J*J)*rho_Jeans_cgs*pow(dx_cell_levelmax_cgs, 3.0)
M_Sink  = M_Sink/Msun

print("rho_Jeans_cgs", "{:e}".format(rho_Jeans_cgs)   )
print("R_Jeans_cgs=", "{:e}".format(R_Jeans_cgs)   )
print("n_Jeans_cgs =", "{:e}".format(n_Jeans_cgs)   )
print("M_Jeans=", "{:e}".format(M_Jeans), "Msun"   )
print("M_Sink (J=", J, ") =", "{:e}".format(M_Sink), "Msun"   )
print("R_gal_Kpc=", "{:e}".format(R_Jeans_cgs/Kpc)   )
#R_J_cgs = 15.0*kB*T/(4.0*pi*G_cgs)

rho_star_cgs = rho_Jeans_cgs/10.0
M_star  = rho_star_cgs*pow(dx_cell_levelmax_cgs, 3.0)
M_star  = M_star/Msun

print("M_star (J=", J, ") =", "{:e}".format(M_star), "Msun"   )

print('----------------Typical BH numbers-------------------------')
M_seed  = 3.7e6 
M_seed  = M_seed*Msun 
rfrac    = 0.5  
msq_rho_inv = (32.0/3.0) * (pi*pow(G_cgs, 3.0)/pow(c_cgs, 6.0))
rho_bh_cgs =  1.0/(M_seed*M_seed*msq_rho_inv)
n_bh_cgs   = rho_bh_cgs/mH
Rs_cgs  = 2.0*G_cgs*M_seed/(c_cgs*c_cgs)


print("Mseed : ", M_seed/Msun)
print("Rs : ", Rs_cgs/AU_cgs)
print("rho_bh_cgs : ", rho_bh_cgs)
print("n_bh_cgs : ", n_bh_cgs)


#n_sink_cgs  = 0.01
#rho_sink_cgs  = n_sink_cgs*mH


#print("d_sink_cgs  =", "{:e}".format(d_sink_cgs)   )
#print('V_cell at level ', levelmax, ' = ', "{:e}".format(V_cell_levelmax) , "h^-3 Mpc^3 " )
#min((d-d_thres/4.0)*dx_loc**3,Mseed*2d33/scale_m)
#at d=1 e-22 T=10000
























