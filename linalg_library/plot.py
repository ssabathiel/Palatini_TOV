import numpy as np
import matplotlib
import matplotlib.image as image
matplotlib.rcParams['backend'] = "TkAgg"
import matplotlib.pyplot as plt


arnett=0
weily=0
baym=0
sil_GR=1
sil_fR=1

diff=1

rad_mass=2     #if turned on 2: rad_mass is on
rad_mass=2   #if turned on 0: rho

weily_rad_mass=0

if rad_mass==2:
    weily_rad_mass=1



with open('TOV_output_fR') as f:
    lines = f.readlines()
    x = [line.split()[rad_mass] for line in lines]
    y = [line.split()[1] for line in lines]

with open('TOV_Baym_output') as f:
    lines = f.readlines()
    x2 = [line.split()[rad_mass] for line in lines]
    y2 = [line.split()[1] for line in lines]

with open('TOV_Arnett_output') as f:
    lines = f.readlines()
    x3 = [line.split()[rad_mass] for line in lines]
    y3 = [line.split()[1] for line in lines]

with open('TOV_Weily_output') as f:
    lines = f.readlines()
    x4 = [line.split()[weily_rad_mass] for line in lines]
    y4 = [line.split()[2] for line in lines]

with open('TOV_output_GR') as f:
    lines = f.readlines()
    x5 = [line.split()[rad_mass] for line in lines]
    y5 = [line.split()[1] for line in lines]

with open('TOV_output_diff') as f:
    lines = f.readlines()
    x6 = [line.split()[1] for line in lines]
    y6 = [line.split()[4] for line in lines]



for i in range(len(x)):
    x[i]=np.log10(float(x[i]))
for i in range(len(x2)):
    x2[i]=np.log10(float(x2[i]))
for i in range(len(x3)):
    x3[i]=np.log10(float(x3[i]))
for i in range(len(x4)):
    x4[i]=np.log10(float(x4[i])*100000)
for i in range(len(x5)):
    x5[i]=np.log10(float(x5[i]))
for i in range(len(x6)):
    x6[i]=np.log10(float(x6[i]))
#print(type(x[0]))

im = image.imread('/home/silvester/Desktop/Masterthesis_Grav/Code/grav_master_git/TOV_Baym.png')
fig, ax = plt.subplots()
#ax.imshow(im, aspect='auto', extent=(5.2, 9.0, -0.1, 1.6), zorder=-1)

#plt.semilogx(x,y)

if sil_fR==1:
    plt.plot(x,y, label='Silvester fR')
if baym==1:
    plt.plot(x2,y2, label='Baym')
if arnett==1:
    plt.plot(x3,y3, label='Arnett')
if weily==1:
    plt.plot(x4,y4, label='Weilnboeck')
if sil_GR==1:
    plt.plot(x,y, label='Silvester GR')
legend = ax.legend(loc='upper center', shadow=True)
plt.show()


if diff==1:
    plt.plot(x6,y6, label='Differences')
plt.show()



