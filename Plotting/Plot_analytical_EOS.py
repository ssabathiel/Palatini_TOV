import numpy as np
import matplotlib
import matplotlib.image as image
matplotlib.rcParams['backend'] = "TkAgg"
import matplotlib.pyplot as plt


rad_mass = 0

def plot_data(analyt,derivatives,second_derivatives,f,ax):
    with open(analyt) as f:
        lines = f.readlines()
        xi1 = [line.split()[rad_mass] for line in lines]
        zeta = [line.split()[1] for line in lines]

    with open(derivatives) as f:
        lines = f.readlines()
        xi2 = [line.split()[rad_mass] for line in lines]
        dzetadxi = [line.split()[1] for line in lines]

    with open(second_derivatives) as f:
        lines = f.readlines()
        xi3 = [line.split()[rad_mass] for line in lines]
        ddzetaddxi = [line.split()[1] for line in lines]



    for i in range(len(xi1)):
        xi1[i]=float(xi1[i])
        zeta[i]=float(zeta[i])
    for i in range(len(xi2)):
        xi2[i]=float(xi2[i])
        dzetadxi[i]=float(dzetadxi[i])
    for i in range(len(xi3)):
        xi3[i]=float(xi3[i])
        ddzetaddxi[i]=float(ddzetaddxi[i])

    labelly = ""
    if(analyt=='analyt_EOS_sly'):
        labelly = 'SLY'
    if(analyt=='analyt_EOS_fps'):
        labelly = 'FPS'
    if (analyt == 'analyt_EOS_ply'):
        labelly = 'PLY'

    ax[0].plot(xi1, zeta,label=labelly)
    ax[1].plot(xi2, dzetadxi)
    ax[2].plot(xi3, ddzetaddxi)

#im = image.imread('/home/silvester/Desktop/Masterthesis_Grav/Code/grav_master_git/TOV_Baym.png')
#fig, ax = plt.subplots()
#ax.imshow(im, aspect='auto', extent=(5.2, 9.0, -0.1, 1.6), zorder=-1)



f, ax = plt.subplots(3, sharex=True)
#f.suptitle('Analytical EOS and its derivatives')

plot_data('analyt_EOS_sly','derivatives_sly','2nd_derivatives_sly',f,ax)
plot_data('analyt_EOS_fps','derivatives_fps','2nd_derivatives_fps',f,ax)
plot_data('analyt_EOS_ply','derivatives_ply','2nd_derivatives_ply',f,ax)

ax[0].invert_xaxis()
ax[0].set_xlim([15, 10])
ax[0].set_ylim([24, 36])
ax[1].set_ylim([0, 4])
ax[2].set_ylim([-4, 7])
plt.xlabel(r'$\xi$',fontsize=20)
ax[0].set_ylabel(r'$\zeta$   ',rotation=0,fontsize=20)
ax[1].set_ylabel(r'$\frac{d\zeta}{d\xi}$   ',rotation=0,fontsize=20)
ax[2].set_ylabel(r'$\frac{d^2\zeta}{d\xi^2}$   ',rotation=0,fontsize=20)
ax[0].legend()



plt.show()


