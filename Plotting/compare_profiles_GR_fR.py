import numpy as np
import matplotlib
import matplotlib.image as image
matplotlib.rcParams['backend'] = "TkAgg"
import matplotlib.pyplot as plt
from itertools import cycle


sil_fR=0
sil_GR=0
diff=0
profile = 1
lines = ["-","--","-.",":"]
linecycler = cycle(lines)


###############
## Get data from files to arrays
###############

###### STAR 1

with open('star2/TOV_output_GR_SLY') as f:
    lines = f.readlines()
    r_GR = [line.split()[2] for line in lines]
    rho_GR = [line.split()[1] for line in lines]
    press_GR = [line.split()[0] for line in lines]

with open('star2/TOV_output_fR_Rp_2.0*10^+9.3_SLY') as f:
    lines = f.readlines()
    r_fR = [line.split()[2] for line in lines]
    rho_fR = [line.split()[1] for line in lines]
    press_fR = [line.split()[0] for line in lines]


for i in range(len(r_GR)):
    r_GR[i]=float(r_GR[i])/100000.
    rho_GR[i] =   float(rho_GR[i])
    press_GR[i] =  float(press_GR[i])

for i in range(len(r_fR)):
    r_fR[i]=float(r_fR[i])/100000.
    rho_fR[i] =   float(rho_fR[i])
    press_fR[i] =  float(press_fR[i])

'''
################################
## Create new list for difference of GR and fR
###############################
delta_rho = list(range(len(rho_GR)))
for i in range(len(rho_GR)):
    delta_rho[i]=float(rho_GR[i])-float(rho_fR[i])
'''



plt.plot(r_GR, rho_GR, label='GR')#, marker='o', ms=450, mec='None', color='#e6e6e6')
plt.plot(r_fR, rho_fR,label='f(R)')
plt.xlabel(r'$r [km]$',fontsize=30)
plt.ylabel(r'$\rho [g/cm^3]$',fontsize=30)
plt.legend(fontsize=30, loc='lower left')
plt.show()

plt.plot(r_GR, press_GR, label='GR')#, marker='o', ms=450, mec='None', color='#e6e6e6')
plt.plot(r_fR, press_fR,label='f(R)')
plt.xlabel(r'$r [km]$',fontsize=30)
plt.ylabel(r'$P [g/cm s^2]$',fontsize=30)
plt.legend(fontsize=30, loc='lower left')
plt.show()





######## Plotting maximum masses

with open('../fR_SLY_Maxmass') as f:
    lines = f.readlines()
    alphas = [line.split()[0] for line in lines]
    maxmasses = [line.split()[1] for line in lines]
    maxmass_radi = [line.split()[2] for line in lines]



for i in range(len(alphas)):
    alphas[i]=float(alphas[i])/(10**9)
    maxmasses[i] =   float(maxmasses[i])
    maxmass_radi[i] =  float(maxmass_radi[i])/100000.

print('alphas= ',alphas)
print('amaxmasses= ',maxmasses)

plt.plot(alphas, maxmasses)#, marker='o', ms=450, mec='None', color='#e6e6e6')
plt.xlabel(r'$\alpha/\alpha_9$',fontsize=30)
plt.ylabel(r'$M_{max} [M_\odot]$',fontsize=30)
#plt.legend(fontsize=30, loc='lower left')
plt.show()


plt.plot(alphas, maxmass_radi)#, marker='o', ms=450, mec='None', color='#e6e6e6')
plt.xlabel(r'$\alpha/\alpha_9$',fontsize=30)
plt.ylabel(r'$ R(M_{max}) [km]$',fontsize=30)
#plt.legend(fontsize=30, loc='lower left')
plt.show()











'''
#Get Profiles for fR
with open('Profiles/p_rho_profile_Rp_2_2_10_ccount_7') as f:
    lines = f.readlines()
    press_fR_7 = [line.split()[0] for line in lines]
    rho_fR_7 = [line.split()[1] for line in lines]
    r_fR_7 = [line.split()[2] for line in lines]

with open('Profiles/p_rho_profile_Rp_2_2_10_ccount_41') as f:
    lines = f.readlines()
    press_fR_41 = [line.split()[0] for line in lines]
    rho_fR_41 = [line.split()[1] for line in lines]
    r_fR_41 = [line.split()[2] for line in lines]

with open('Profiles/p_rho_profile_Rp_2_2_10_ccount_49') as f:
    lines = f.readlines()
    press_fR_49 = [line.split()[0] for line in lines]
    rho_fR_49 = [line.split()[1] for line in lines]
    r_fR_49 = [line.split()[2] for line in lines]

with open('Profiles/p_rho_profile_Rp_2_2_10_ccount_54') as f:
    lines = f.readlines()
    press_fR_54 = [line.split()[0] for line in lines]
    rho_fR_54 = [line.split()[1] for line in lines]
    r_fR_54 = [line.split()[2] for line in lines]

with open('Profiles/p_rho_profile_Rp_2_2_10_ccount_55') as f:
    lines = f.readlines()
    press_fR_55 = [line.split()[0] for line in lines]
    rho_fR_55 = [line.split()[1] for line in lines]
    r_fR_55 = [line.split()[2] for line in lines]

with open('Profiles/p_rho_profile_Rp_2_2_10_ccount_56') as f:
    lines = f.readlines()
    press_fR_56 = [line.split()[0] for line in lines]
    rho_fR_56 = [line.split()[1] for line in lines]
    r_fR_56 = [line.split()[2] for line in lines]





#Get Profiles for GR
with open('Profiles/p_rho_profile_Rp_2_2_80_ccount_7') as f:
    lines = f.readlines()
    press_GR_7 = [line.split()[0] for line in lines]
    rho_GR_7 = [line.split()[1] for line in lines]
    r_GR_7 = [line.split()[2] for line in lines]

with open('Profiles/p_rho_profile_Rp_2_2_80_ccount_41') as f:
    lines = f.readlines()
    press_GR_41 = [line.split()[0] for line in lines]
    rho_GR_41 = [line.split()[1] for line in lines]
    r_GR_41 = [line.split()[2] for line in lines]

with open('Profiles/p_rho_profile_Rp_2_2_80_ccount_49') as f:
    lines = f.readlines()
    press_GR_49 = [line.split()[0] for line in lines]
    rho_GR_49 = [line.split()[1] for line in lines]
    r_GR_49 = [line.split()[2] for line in lines]

with open('Profiles/p_rho_profile_Rp_2_2_80_ccount_54') as f:
    lines = f.readlines()
    press_GR_54 = [line.split()[0] for line in lines]
    rho_GR_54 = [line.split()[1] for line in lines]
    r_GR_54 = [line.split()[2] for line in lines]

with open('Profiles/p_rho_profile_Rp_2_2_80_ccount_55') as f:
    lines = f.readlines()
    press_GR_55 = [line.split()[0] for line in lines]
    rho_GR_55 = [line.split()[1] for line in lines]
    r_GR_55 = [line.split()[2] for line in lines]

with open('Profiles/p_rho_profile_Rp_2_2_80_ccount_56') as f:
    lines = f.readlines()
    press_GR_56 = [line.split()[0] for line in lines]
    rho_GR_56 = [line.split()[1] for line in lines]
    r_GR_56 = [line.split()[2] for line in lines]

'''





















'''
##############################
## Convert from String to Float (and optional x to log scale)
############################

# Convert fR
for i in range(len(press_fR_7)):
    press_fR_7[i]=float(press_fR_7[i])
    rho_fR_7[i] =   float(rho_fR_7[i])
    r_fR_7[i] =  float(r_fR_7[i])/100000.

for i in range(len(press_fR_41)):
    press_fR_41[i]=float(press_fR_41[i])
    rho_fR_41[i] =  float(rho_fR_41[i])
    r_fR_41[i] =  float(r_fR_41[i])/100000.

for i in range(len(press_fR_49)):
    press_fR_49[i]=float(press_fR_49[i])
    rho_fR_49[i] = float(rho_fR_49[i])
    r_fR_49[i] =  float(r_fR_49[i])/100000.

for i in range(len(press_fR_54)):
    press_fR_54[i]=float(press_fR_54[i])
    #rho_fR_54[i] =  np.log10( float(rho_fR_54[i]) )
    rho_fR_54[i] =  float(rho_fR_54[i])
    r_fR_54[i] =  float(r_fR_54[i])/100000.

for i in range(len(press_fR_55)):
    press_fR_55[i]=float(press_fR_55[i])
    rho_fR_55[i] =  float(rho_fR_55[i])
    r_fR_55[i] =  float(r_fR_55[i])/100000.

for i in range(len(press_fR_56)):
    press_fR_56[i]=float(press_fR_56[i])
    rho_fR_56[i] =  float(rho_fR_56[i])
    r_fR_56[i] =  float(r_fR_56[i])/100000.



# Convert GR
for i in range(len(press_GR_7)):
    press_GR_7[i]=float(press_GR_7[i])
    rho_GR_7[i] =  float(rho_GR_7[i])
    r_GR_7[i] =  float(r_GR_7[i])/100000.

for i in range(len(press_GR_41)):
    press_GR_41[i]=float(press_GR_41[i])
    rho_GR_41[i] =  float(rho_GR_41[i])
    r_GR_41[i] =  float(r_GR_41[i])/100000.

for i in range(len(press_GR_49)):
    press_GR_49[i]=float(press_GR_49[i])
    rho_GR_49[i] = float(rho_GR_49[i])
    r_GR_49[i] =  float(r_GR_49[i])/100000.

for i in range(len(press_GR_54)):
    press_GR_54[i]=float(press_GR_54[i])
    #rho_GR_54[i] = np.log10(float(rho_GR_54[i]))
    rho_GR_54[i] = float(rho_GR_54[i])
    r_GR_54[i] =  float(r_GR_54[i])/100000.

for i in range(len(press_GR_55)):
    press_GR_55[i]=float(press_GR_55[i])
    rho_GR_55[i] =  float(rho_GR_55[i])
    r_GR_55[i] =  float(r_GR_55[i])/100000.

for i in range(len(press_GR_56)):
    press_GR_56[i]=float(press_GR_56[i])
    rho_GR_56[i] =  float(rho_GR_56[i])
    r_GR_56[i] =  float(r_GR_56[i])/100000.

'''




#im = image.imread('/home/silvester/Desktop/Masterthesis_Grav/Code/grav_master_git/TOV_Baym.png')
#fig, ax = plt.subplots()
#ax.imshow(im, aspect='auto', extent=(5.2, 9.0, -0.1, 1.6), zorder=-1)





############################
## Plot fR and GR
#######################
if sil_fR==1:
    plt.plot(x,y, label='Silvester fR')
if sil_GR==1:
    plt.plot(x5,y5, label='Silvester GR')
    # TEXT
    plt.xlabel('R [km]')
    plt.ylabel('$M/M_0$')
    plt.title('fR - GR comparison')
    plt.legend()
    plt.show()



#######################
## Plot diff of GR and fR
########################
if diff==1:
    plt.plot(x,y6)
    # TEXT
    plt.xlabel('R [km]')
    plt.ylabel('$M/M_0$')
    plt.title('Differences between GR and fR: $fR_M - GR_M$')
    plt.legend()
    plt.show()



######################
## Plot p and rho profile
######################
#neutron_star = plt.Circle((0,0), radius = r_GR_54[-1], color='#e6e6e6')
#plt.gca().add_artist(neutron_star)

'''
if profile==1:

    # Star 7
    # Rho
    x = 0
    y = 0
    plt.plot(x, y, marker='o', ms=450, mec='None', color='#e6e6e6')

    plt.plot(r_fR_7, rho_fR_7, label='fR')
    plt.plot(r_GR_7, rho_GR_7, label='GR')

    # TEXT
    labely = r'$\rho _c= 10$^' + str(np.log10(rho_fR_7[0]))[0] \
             + str(np.log10(rho_fR_7[0]))[1] \
             + str(np.log10(rho_fR_7[0]))[2] \
             + str(np.log10(rho_fR_7[0]))[3] \
             + str(np.log10(rho_fR_7[0]))[4] \
             + ' ' + r'$\frac{g}{cm^3}$'
    plt.xlabel('r [km]')
    plt.ylabel('rho [$g/cm^3$]')
    plt.title('GR/fR-Comparison of Density-Radius Profile1 ' + labely)
    plt.legend()
    plt.axis([-r_GR_7[-1], r_GR_7[-1], -2 * rho_GR_7[0], 2 * rho_GR_7[0]])
    plt.show()

    #Press
    plt.plot(x, y, marker='o', ms=450, mec='None', color='#e6e6e6')
    plt.plot(r_fR_7, press_fR_7, label='fR')
    plt.plot(r_GR_7, press_GR_7, label='GR')
    # TEXT
    plt.xlabel('r [km]')
    plt.ylabel('press [$g/(cm*s^2)$]')
    plt.title('Comparison of Press-Radius Profile of GR and fR' + labely)
    plt.legend()
    plt.axis([-r_GR_7[-1], r_GR_7[-1], -2 * press_GR_7[0], 2 * press_GR_7[0]])
    plt.show()


    # Star 49
    # Rho
    x=0
    y=0
    plt.plot(x, y, marker='o', ms= 450, mec='None', color='#e6e6e6')

    plt.plot(r_fR_49, rho_fR_49, label='fR')
    plt.plot(r_GR_49, rho_GR_49, label='GR')

    # TEXT
    labely = r'$\rho _c= 10$^' + str(np.log10( rho_fR_49[0]) )[0] \
                                + str(np.log10( rho_fR_49[0]) )[1] \
                                + str(np.log10(  rho_fR_49[0]) )[2] \
                                + str(np.log10(  rho_fR_49[0]) )[3] \
                                + str(np.log10(rho_fR_49[0]))[4] \
                                +' ' + r'$\frac{g}{cm^3}$'
    plt.xlabel('r [km]')
    plt.ylabel('rho [$g/cm^3$]')
    plt.title('GR/fR-Comparison of Density-Radius Profile '+ labely)
    plt.legend()
    plt.axis([-r_GR_49[-1], r_GR_49[-1], -2*rho_GR_49[0], 2*rho_GR_49[0]])
    plt.show()

    #Press
    plt.plot(x, y, marker='o', ms=450, mec='None', color='#e6e6e6')
    plt.plot(r_fR_49, press_fR_49, label='fR')
    plt.plot(r_GR_49, press_GR_49, label='GR')
    # TEXT
    plt.xlabel('r [km]')
    plt.ylabel('press [$g/(cm*s^2)$]')
    plt.title('Comparison of Press-Radius Profile of GR and fR' + labely)
    plt.legend()
    plt.axis([-r_GR_49[-1], r_GR_49[-1], -2 * press_GR_49[0], 2 * press_GR_49[0]])
    plt.show()

    # Star 54
    # Rho
    x=0
    y=0
    plt.plot(x, y, marker='o', ms=450, mec='None', color='#e6e6e6')

    plt.plot(r_fR_54, rho_fR_54, label='fR')
    plt.plot(r_GR_54, rho_GR_54, label='GR')

    # TEXT
    labely = r'$\rho _c= 10$^' + str(np.log10( rho_fR_54[0]) )[0] \
                                + str(np.log10( rho_fR_54[0]) )[1] \
                                + str(np.log10(  rho_fR_54[0]) )[2] \
                                + str(np.log10(  rho_fR_54[0]) )[3] \
                                + str(np.log10(rho_fR_54[0]))[4] \
                                +' ' + r'$\frac{g}{cm^3}$'
    plt.xlabel('r [km]')
    plt.ylabel('rho [$g/cm^3$]')
    plt.title('GR/fR-Comparison of Density-Radius Profile ' + labely)
    plt.legend()
    plt.axis([-r_GR_54[-1], r_GR_54[-1], -2*rho_GR_54[0], 2*rho_GR_54[0]])
    plt.show()

    #Press
    plt.plot(x, y, marker='o', ms=450, mec='None', color='#e6e6e6')
    plt.plot(r_fR_54, press_fR_54, label='fR')
    plt.plot(r_GR_54, press_GR_54, label='GR')
    # TEXT
    plt.xlabel('r [km]')
    plt.ylabel('press [$g/(cm*s^2)$]')
    plt.title('Comparison of Press-Radius Profile of GR and fR' + labely)
    plt.legend()
    plt.axis([-r_GR_54[-1], r_GR_54[-1], -2 * press_GR_54[0], 2 * press_GR_54[0]])
    plt.show()

# Star 55
    # Rho
    x=0
    y=0
    plt.plot(x, y, marker='o', ms=450, mec='None', color='#e6e6e6')

    plt.plot(r_fR_55, rho_fR_55, label='fR')
    plt.plot(r_GR_55, rho_GR_55, label='GR')

    # TEXT
    labely = r'$\rho _c= 10$^' + str(np.log10( rho_fR_55[0]) )[0] \
                                + str(np.log10( rho_fR_55[0]) )[1] \
                                + str(np.log10(  rho_fR_55[0]) )[2] \
                                + str(np.log10(  rho_fR_55[0]) )[3] \
                                + str(np.log10(rho_fR_55[0]))[4] \
                                +' ' + r'$\frac{g}{cm^3}$'
    plt.xlabel('r [km]')
    plt.ylabel('rho [$g/cm^3$]')
    plt.title('GR/fR-Comparison of Density-Radius Profile ' + labely)
    plt.legend()
    plt.axis([-r_GR_55[-1], r_GR_55[-1], -2*rho_GR_55[0], 2*rho_GR_55[0]])
    plt.show()

    #Press
    plt.plot(x, y, marker='o', ms=450, mec='None', color='#e6e6e6')
    plt.plot(r_fR_55, press_fR_55, label='fR')
    plt.plot(r_GR_55, press_GR_55, label='GR')
    # TEXT
    plt.xlabel('r [km]')
    plt.ylabel('press [$g/(cm*s^2)$]')
    plt.title('Comparison of Press-Radius Profile of GR and fR' + labely)
    plt.legend()
    plt.axis([-r_GR_55[-1], r_GR_55[-1], -2 * press_GR_55[0], 2 * press_GR_55[0]])
    plt.show()


    # Star 56
    # Rho
    x=0
    y=0
    plt.plot(x, y, marker='o', ms=450, mec='None', color='#e6e6e6')

    plt.plot(r_fR_56, rho_fR_56, label='fR')
    plt.plot(r_GR_56, rho_GR_56, label='GR')

    # TEXT
    labely = r'$\rho _c= 10$^' + str(np.log10( rho_fR_56[0]) )[0] \
                                + str(np.log10( rho_fR_56[0]) )[1] \
                                + str(np.log10(  rho_fR_56[0]) )[2] \
                                + str(np.log10(  rho_fR_56[0]) )[3] \
                                + str(np.log10(rho_fR_56[0]))[4] \
                                +' ' + r'$\frac{g}{cm^3}$'
    plt.xlabel('r [km]')
    plt.ylabel('rho [$g/cm^3$]')
    plt.title('GR/fR-Comparison of Density-Radius Profile ' + labely)
    plt.legend()
    plt.axis([-r_GR_56[-1], r_GR_56[-1], -2*rho_GR_56[0], 2*rho_GR_56[0]])
    plt.show()

    #Press
    plt.plot(x, y, marker='o', ms=450, mec='None', color='#e6e6e6')
    plt.plot(r_fR_56, press_fR_56, label='fR')
    plt.plot(r_GR_56, press_GR_56, label='GR')
    # TEXT
    plt.xlabel('r [km]')
    plt.ylabel('press [$g/(cm*s^2)$]')
    plt.title('Comparison of Press-Radius Profile of GR and fR' + labely)
    plt.legend()
    plt.axis([-r_GR_56[-1], r_GR_56[-1], -2 * press_GR_56[0], 2 * press_GR_56[0]])
    plt.show()


'''