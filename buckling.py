import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def buckling_equation(k_minus, k_plus):
    return 2 * k_plus * k_minus * (np.cos(k_plus) * np.cos(k_minus) - 1) + (k_plus**2 + k_minus**2) * np.sin(k_plus) * np.sin(k_minus)

def solve_for_k_minus(k_plus_initial, k_minus_initial_guess):
    solution = fsolve(buckling_equation, k_minus_initial_guess, args=(k_plus_initial), xtol=1e-6)
    return solution[0]

k_plus_range_hyperbolas = np.linspace(0.05, 8*np.pi, 500)
etas = [0, 2, 5.44, 7, 8.89, 10, 12.17, 15.39, 18.59]



font = 30
plt.figure(figsize=(10,12))

for eta in etas:
    k_minus_hyperbola = (eta**2) / k_plus_range_hyperbolas
    plt.plot(k_plus_range_hyperbolas, k_minus_hyperbola, 'black', alpha=0.7)

x_positions = [0.001, 0.123 * np.pi, 0.52 * np.pi, 0.9 * np.pi, 1.4 * np.pi, 
               1.9 * np.pi, 2.5 * np.pi, 3.85 * np.pi, 5.3 * np.pi]
y_positions = [1.2, 5 * np.pi, 5.2 * np.pi, 5.4 * np.pi, 5.6 * np.pi, 
               5.8 * np.pi, 6.5 * np.pi, 6.5 * np.pi, 7 * np.pi]

for i, eta in enumerate(etas):
    x_pos = x_positions[i]
    y_pos = y_positions[i]
    plt.text(x_pos, y_pos, f'{eta:.2f}', fontsize=font, rotation=270, verticalalignment='center')

tick_positions = np.arange(0, 9, 2) * np.pi  
tick_labels = [f'${i}\\pi$' if i > 0 else '0' for i in range(0, 9, 2)]  
plt.xticks(tick_positions, tick_labels, fontsize=font)
plt.yticks(tick_positions, tick_labels, fontsize=font)

plt.fill_between([3*np.pi, 4*np.pi], 0, np.pi, color='gray', alpha=0.5)
plt.fill_between([5*np.pi, 6*np.pi], 0, np.pi, color='gray', alpha=0.5)
plt.fill_between([7*np.pi, 8*np.pi], 0, np.pi, color='gray', alpha=0.5)
plt.fill_between([4*np.pi, 5*np.pi],  np.pi, 2*np.pi, color='gray', alpha=0.5)
plt.fill_between([6*np.pi, 7*np.pi], np.pi, 2*np.pi, color='gray', alpha=0.5)
plt.fill_between([5*np.pi, 6*np.pi], 2*np.pi, 3*np.pi, color='gray', alpha=0.5)
plt.fill_between([7*np.pi, 8*np.pi], 2*np.pi, 3*np.pi, color='gray', alpha=0.5)
plt.fill_between([6*np.pi, 7*np.pi], 3*np.pi, 4*np.pi, color='gray', alpha=0.5)
plt.fill_between([7*np.pi, 8*np.pi], 4*np.pi, 5*np.pi, color='gray', alpha=0.5)

def buckling_modes(k_plus, k_minus):
    term1 = 2 * k_plus * k_minus * (np.cos(k_plus) * np.cos(k_minus) - 1)
    term2 = (k_plus**2 + k_minus**2) * np.sin(k_plus) * np.sin(k_minus)
    return term1 + term2

k_plus_vals = np.linspace(0, 8*np.pi, 2000)  
k_minus_vals = np.linspace(0, 8*np.pi, 2000)
Kp, Km = np.meshgrid(k_plus_vals, k_minus_vals)
Z = buckling_modes(Kp, Km)
Z_masked = np.ma.masked_where(Kp <= Km, Z)

plt.scatter([2*np.pi], [0], color='red', label=r'Planar Elastica ($k_+ = 2\pi$, $k_- = 0$)', zorder=5)
plt.contour(Kp, Km, Z_masked, levels=[0], colors='blue', label='Antisymmetric modes')
plt.xlabel('$k_+$', fontsize=font)
plt.ylabel('$k_-$', fontsize=font)
plt.xlim([0, 8*np.pi])
plt.ylim([-0.5, 8*np.pi])

plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12), fontsize=font,frameon=False)

plt.grid(True)
plt.xticks(fontsize=font)  
plt.yticks(fontsize=font)  

plt.gca().spines['top'].set_linewidth(3)
plt.gca().spines['right'].set_linewidth(3)
plt.gca().spines['left'].set_linewidth(3)
plt.gca().spines['bottom'].set_linewidth(3)

plt.gca().spines['top'].set_color('black')
plt.gca().spines['right'].set_color('black')
plt.gca().spines['left'].set_color('black')
plt.gca().spines['bottom'].set_color('black')

plt.tight_layout()
plt.show()







