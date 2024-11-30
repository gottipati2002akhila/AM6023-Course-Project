import numpy as np
import matplotlib.pyplot as plt

def buckling_modes(k_p,k_m,s):
    term1=(k_m*np.sin(k_p)-k_p*np.sin(k_m))*(np.cos(k_p*s)-np.cos(k_m*s))
    term2=(np.cos(k_p)-np.cos(k_m))*(k_m*np.sin(k_p*s)-k_p*np.sin(k_m*s))
    return term1 -term2

def kplus(eta):
    f=np.sqrt(1+(eta/np.pi)**2)
    km=np.pi*(f-1)+0.5*np.arccos(1-((1-np.cos(2*np.pi*f))/f**2))
    kp=2*np.pi+np.pi*(f-1)-0.5*np.arccos(1-((1-np.cos(2*np.pi*f))/f**2))
    return kp,km

eta=[2,8.88,5.44,15.39,13.5,18.59]
s=np.linspace(0,1,100)
fig, axs = plt.subplots(6, 1, figsize=(8, 12), sharex=True)
font=20

for i in range(6):
    kp, km = kplus(eta[5 - i])  
    y_vals = buckling_modes(kp, km, s) 
    
    if i % 2 == 0:
        axs[i].plot(s, y_vals, color='r', label=r'$i = %d$' % (6 - i))
    else:
        axs[i].plot(s, -y_vals, color='r', label=r'$i = %d$' % (6 - i))
    axs[i].set_xlabel(r'$s$', fontsize=font)
    axs[i].set_ylabel(r'$\theta(s)$', fontsize=font)
    
    axs[i].set_yticks([])
    axs[i].set_xticks([])

    axs[i].legend(loc='lower right', fontsize=font-3,frameon=False)

plt.tight_layout()
plt.show()