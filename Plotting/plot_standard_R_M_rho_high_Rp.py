import numpy as np
import matplotlib
import matplotlib.image as image
matplotlib.rcParams['backend'] = "TkAgg"
import matplotlib.pyplot as plt
import glob


sil_fR=1
sil_fRQ=1
sil_GR=1
diff=0
profile = 0


###############
## Get data from files to arrays
###############


fR_Rp = glob.glob('Results/TOV_output_fR_*')[0][25:32] #35]
fR_Rp +=   '{' + glob.glob('Results/TOV_output_fR_*')[0][32:35] + '}'
print("Rp=", fR_Rp)

fRQ_Rp = glob.glob('Results/TOV_output_fRQ_*')[0][26:33] #35]
fRQ_Rp +=   '{' + glob.glob('Results/TOV_output_fRQ_*')[0][33:36] + '}'
print("fRQ,Rp=", fRQ_Rp)

fRQ_Rq = glob.glob('Results/TOV_output_fRQ_*')[0][39:46]
fRQ_Rq +=   '{' + glob.glob('Results/TOV_output_fRQ_*')[0][46:49] + '}'
print("fRQ,Rq=", fRQ_Rq)


#fRQ_Rq = glob.glob('Results/TOV_output_fRQ_*')[0][39:49]



with open(glob.glob('Results/TOV_output_fR_*')[0]) as f:
#with open('Results/TOV_output_fR_Rp_2.1*10^-10') as f:
    lines = f.readlines()
    x = [line.split()[2] for line in lines]
    y = [line.split()[1] for line in lines]
    z = [line.split()[0] for line in lines]

with open(glob.glob('Results/TOV_output_fRQ_*')[0]) as f:
    lines = f.readlines()
    x2 = [line.split()[2] for line in lines]
    y2 = [line.split()[1] for line in lines]
    z2 = [line.split()[0] for line in lines]


with open('Results/TOV_output_GR') as f:
    lines = f.readlines()
    x5 = [line.split()[2] for line in lines]
    y5 = [line.split()[1] for line in lines]
    z5 = [line.split()[0] for line in lines]







for i in range(len(x)):
    x[i]=np.log10(float(x[i])/100000.)
    z[i] = np.log10(float(z[i]))
for i in range(len(x5)):
    x5[i]=np.log10(float(x5[i])/100000.)
    z5[i]=np.log10(float(z5[i]))
for i in range(len(x2)):
    x2[i]=np.log10(float(x2[i])/100000.)
    z2[i]=np.log10(float(z2[i]))


#im = image.imread('/home/silvester/Desktop/Masterthesis_Grav/Code/grav_master_git/TOV_Baym.png')
#fig, ax = plt.subplots()
#ax.imshow(im, aspect='auto', extent=(5.2, 9.0, -0.1, 1.6), zorder=-1)





############################
## Plot fR and GR
#######################

fR_label = 'f(R)\n  $(R_p=$' + r'$%s$' %fR_Rp + r'$)$'
fRQ_label = 'f(R,Q) \n  $(R_p=$' + r'$%s$' %fRQ_Rp + '$)$\n  ' + r'$(R_q=$' + r'$%s$' %fRQ_Rq + '$)$'

# Plot R-M
if sil_fR==1:
    plt.plot(x,y, label=fR_label,color='blue')

if sil_fRQ==1:
    plt.plot(x2, y2, label=fRQ_label,color='orange')
if sil_GR==1:
    plt.plot(x5,y5, label='GR',color='black')
    # TEXT
    plt.xlabel('log R [km]')
    plt.ylabel('$M/M_0$')
    plt.title('fR/GR comparison for Mass and Radius for $R_p=2.1*10^{80}$')
    plt.legend()
    plt.show()

#Plot rho-M
if sil_fR==1:
    plt.plot(z,y, label='fR')
if sil_fRQ==1:
    plt.plot(z2, y2, label='fRQ')
if sil_GR==1:
    plt.plot(z5,y5, label='GR')
    # TEXT
    plt.xlabel(r'log $\rho_c$ [$g/cm^3$]')
    plt.ylabel('$M/M_0$')
    plt.title('fR - GR comparison for Mass and Density  for $R_p=2.1*10^{80}$')
    plt.legend()
    plt.show()







print("yes baby")




#######################
## Plot diff of GR and fR
########################
if diff==1:
    plt.plot(x,y6)
    # TEXT
    plt.xlabel('log R [km]')
    plt.ylabel('$M/M_0$')
    plt.title('Differences between GR and fR: $fR_M - GR_M$')
    plt.legend()
    plt.show()



######################
## Plot p and rho profile
######################
neutron_star = plt.Circle((0,0), radius = r_fR_7[-1], color='#e6e6e6')
plt.gca().add_artist(neutron_star)


if profile==1:
    plt.plot(r_fR_7, rho_fR_7)
    # TEXT
    plt.xlabel('r [km]')
    plt.ylabel('rho [$g/cm^3$]')
    plt.title('Density-Radius Profile')
    plt.legend()
    plt.axis([-r_fR_7[-1], r_fR_7[-1], -2*rho_fR_7[0], 2*rho_fR_7[0]])
    plt.show()

    plt.plot(r_fR_7, press_fR_7)
    # TEXT
    plt.xlabel('r [km]')
    plt.ylabel('press [$g/(cm*s^2)$]')
    plt.title('Press-Radius Profile')
    plt.legend()
    plt.show()


#########################
## Plot multiple star profiles in GR
###########################

neutron_star = plt.Circle((0,0), radius = r_fR_7[-1], color='#e6e6e6')
#plt.gca().add_artist(neutron_star)


if profile==1:
    labely = r'$\rho _c= 10$^' + str(rho_fR_7[0])[0] + str(rho_fR_7[0])[1] + str(rho_fR_7[0])[2] + ' ' + r'$\frac{g}{cm^3}$'
    plt.plot(r_fR_7, rho_fR_7, label=labely)
    labely = r'$\rho _c= 10$^' + str(rho_fR_41[0])[0] + str(rho_fR_41[0])[1] + str(rho_fR_41[0])[2] + str(rho_fR_41[0])[3] + ' ' + r'$\frac{g}{cm^3}$'
    plt.plot(r_fR_41, rho_fR_41, label=labely)
    #x=0
    #y=0
    #plt.plot(x, y, marker='o', ms=x * 5, mfc=(1., 0., 0., 0.5), mec='None')
    labely = r'$\rho _c= 10$^' + str(rho_fR_49[0])[0] + str(rho_fR_49[0])[1] + str(rho_fR_49[0])[2] + str(rho_fR_49[0])[3] + ' ' + r'$\frac{g}{cm^3}$'
    plt.plot(r_fR_49, rho_fR_49, label = labely)
    labely = r'$\rho _c= 10$^' + str(rho_fR_54[0])[0] + str(rho_fR_54[0])[1] + str(rho_fR_54[0])[2] + str(rho_fR_54[0])[3] + ' ' + r'$\frac{g}{cm^3}$'
    plt.plot(r_fR_54, rho_fR_54, label = labely)
    labely = r'$\rho _c= 10$^' + str(rho_fR_55[0])[0] + str(rho_fR_55[0])[1] + str(rho_fR_55[0])[2] + str(rho_fR_55[0])[3] + ' ' + r'$\frac{g}{cm^3}$'
    plt.plot(r_fR_55, rho_fR_55, label = labely)
    # TEXT
    plt.xlabel('r [km]')
    plt.ylabel(r'$\rho$ [$g/cm^3$]')
    plt.title('Density-Radius Profile for fR')
    plt.legend()
    plt.axis([-2*r_fR_7[-1], 2*r_fR_7[-1], -4, 20])
    #plt.yscale('log')

    plt.show()

    plt.plot(r_fR_7, press_fR_7)
    # TEXT
    plt.xlabel('r [km]')
    plt.ylabel('press [$g/(cm*s^2)$]')
    plt.title('Press-Radius Profile')
    plt.legend()
    plt.show()
