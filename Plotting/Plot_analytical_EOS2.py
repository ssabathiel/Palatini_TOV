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



def plot_single_EOS(analyt,labelly):
    with open(analyt) as f:
        lines = f.readlines()
        rho = [line.split()[0] for line in lines]
        p = [line.split()[1] for line in lines]


    for i in range(len(rho)):
        rho[i]=np.log10( float(rho[i]) )
        p[i]=np.log10( float(p[i]) )


    plt.plot(rho, p,label=labelly)
    plt.xlabel(r'log$\rho(g/cm^3)$')
    plt.ylabel(r'log$p (dyn/cm^2)$')




def compare_analyt_spline(analyt,spline):
    with open(analyt) as f:
        lines = f.readlines()
        rho = [line.split()[0] for line in lines]
        p = [line.split()[1] for line in lines]

    with open(spline) as f:
        lines = f.readlines()
        rho_splined = [line.split()[0] for line in lines]
        p_splined = [line.split()[1] for line in lines]

    logg=1
    sparse_spline=1
    nr_points = 10

    if(logg==1):
        for i in range(len(rho)):
            rho[i]=np.log10( float(rho[i]) )
            p[i]=np.log10( float(p[i]) )
        for i in range(len(rho_splined)):
            rho_splined[i]=np.log10( float(rho_splined[i]) )
            p_splined[i]=np.log10( float(p_splined[i]) )

    else:
        for i in range(len(rho)):
            rho[i]= float(rho[i])
            p[i]= float(p[i])
        for i in range(len(rho_splined)):
            rho_splined[i]= float(rho_splined[i])
            p_splined[i]=float(p_splined[i])


    labelly = ""
    labelly2 = ""
    if(analyt=='rho_of_p_analytical_sly'):
        labelly = 'SLY Analytical'
        labelly2 = 'SLY Splined'
    if(analyt=='rho_of_p_analytical_fps'):
        labelly = 'FPS Analytical'
        labelly2 = 'FPS Splined'
    if (analyt == 'rho_of_p_analytical_ply'):
        labelly = 'PLY Analytical'
        labelly2 = 'PLY Splined'


    plt.plot(rho, p,label=labelly)
    plt.plot(rho_splined[::100], p_splined[::100],'*',label=labelly2)
    plt.xlabel(r'$\log\rho(g/cm^3)$')
    plt.ylabel(r'$\log p (dyn/cm^2)$')

#im = image.imread('/home/silvester/Desktop/Masterthesis_Grav/Code/grav_master_git/TOV_Baym.png')
#fig, ax = plt.subplots()
#ax.imshow(im, aspect='auto', extent=(5.2, 9.0, -0.1, 1.6), zorder=-1)



f, ax = plt.subplots(3, sharex=True)
#f2, ax2 = plt.subplots(3, sharex=True)
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




compare_analyt_spline("rho_of_p_analytical_sly","rho_of_p_splined_sly")





#plot_single_EOS("rho_of_p_analytical_sly","SLY")
#plot_single_EOS("rho_of_p_analytical_ply","PLY")
#plot_single_EOS("rho_of_p_analytical_fps","FPS")
plt.legend(loc='1')
plt.show()


E = np.exp(1.0)

def ddpsi_ddeps(eps):
    ddpsiddeps = (11.690904 * E ** (6.54 * (11.495 - eps))) / (1 + E ** (6.54 * (11.495 - eps))) ** 2 + \
      (13.50802 * E ** (4.3 * (14.08 - eps))) / (1 + E ** (4.3 * (14.08 - eps))) ** 2 -   \
      (4.959 * E ** (1.5 * (14.67 - eps))) / (1 + E ** (1.5 * (14.67 - eps))) ** 2 +     \
      (4.5 * E ** (3. * (14.67 - eps)) * (27.8 - 1.653 * eps)) / (1 + E ** (1.5 * (14.67 - eps))) ** 3 - \
      (2.25 * E ** (1.5 * (14.67 - eps)) * (27.8 - 1.653 * eps)) / (1 + E ** (1.5 * (14.67 - eps))) ** 2 +\
      (85.5432 * E ** (13.08 * (11.495 - eps)) * (19.105 + 0.8938 * eps)) /  \
       (1 + E ** (6.54 * (11.495 - eps))) ** 3 -   \
      (42.7716 * E ** (6.54 * (11.495 - eps)) * (19.105 + 0.8938 * eps)) /   \
       (1 + E ** (6.54 * (11.495 - eps))) ** 2 +   \
      (0.03555 * eps) / ((1 + E ** (6.48 * (-11.4971 + eps))) * (1 + 0.16326 * eps)) +   \
      (36.98 * E ** (8.6 * (14.08 - eps)) * (-22.775 + 1.5707 * eps)) / (1 + E ** (4.3 * (14.08 - eps))) ** 3 -\
      (18.49 * E ** (4.3 * (14.08 - eps)) * (-22.775 + 1.5707 * eps)) / (1 + E ** (4.3 * (14.08 - eps))) ** 2 -\
      (0.32652 * (6.121 + 0.017775 * eps ** 2)) /  \
       ((1 + E ** (6.48 * (-11.4971 + eps))) * (1 + 0.16326 * eps) ** 2) -    \
      (12.96 * E ** (6.48 * (-11.4971 + eps)) * (6.121 + 0.017775 * eps ** 2)) /   \
       ((1 + E ** (6.48 * (-11.4971 + eps))) ** 2 * (1 + 0.16326 * eps)) +    \
      (0.05330765519999999 * (6.22 + 6.121 * eps + 0.005925 * eps ** 3)) /  \
       ((1 + E ** (6.48 * (-11.4971 + eps))) * (1 + 0.16326 * eps) ** 3) +  \
      (2.1158495999999998 * E ** (6.48 * (-11.4971 + eps)) * (6.22 + 6.121 * eps + 0.005925 * eps ** 3)) /   \
       ((1 + E ** (6.48 * (-11.4971 + eps))) ** 2 * (1 + 0.16326 * eps) ** 2) +   \
      (83.98080000000002 * E ** (12.96 * (-11.4971 + eps)) * (6.22 + 6.121 * eps + 0.005925 * eps ** 3)) / \
       ((1 + E ** (6.48 * (-11.4971 + eps))) ** 3 * (1 + 0.16326 * eps)) -  \
      (41.99040000000001 * E ** (6.48 * (-11.4971 + eps)) * (6.22 + 6.121 * eps + 0.005925 * eps ** 3)) /   \
       ((1 + E ** (6.48 * (-11.4971 + eps))) ** 2 * (1 + 0.16326 * eps))   \


    return ddpsiddeps



eps=np.linspace(15,10,100)

#plt.plot(eps, ddpsi_ddeps(eps))
#plt.show()



x = linspa





