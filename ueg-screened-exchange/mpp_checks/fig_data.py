"""Data behind the plotted curves of main-text Figures 1, 2 and 4 (CSV files in
mpp_checks/figure_data/).  The figures are drawn from these files by the scripts in
mpp_figures/ (Fig. 4 evaluates its closed form directly).

  fig1_poles.csv          s = N(d-1) (N = 1..8), s = 2, s = d+1, s = 4 on 1.2 <= d <= 5
  fig1_coincidences.csv   the marked coincidences (open circles at d = 2, 3; filled dots at d = 5/3, 3/2)
  fig2_R2.csv             R_2(z), h_2(z), R_3(z) = R^(1)(z)/(4 pi), h_3(z) on 0 <= z <= 4
  fig2_R2_loglog.csv      R_2(z) on 1e-2 <= z <= 1e2 with 2/pi and 5/(12 pi z^4)
  fig4_delta.csv          delta_p/(q^2/2) on the cross section through the q axis, q = 3.2, inside the two balls
"""
import os, csv
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'figure_data'); os.makedirs(OUT, exist_ok=True)
def write(name, header, rows):
    with open(os.path.join(OUT, name), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(header); [w.writerow([f'{x:.12g}' for x in r]) for r in rows]
d = np.linspace(1.2, 5, 381)
write('fig1_poles.csv', ['d'] + [f's_N{N}' for N in range(1, 9)] + ['s_insertion', 's_rings', 's_q2_insertion'],
      [[x] + [N * (x - 1) for N in range(1, 9)] + [2, x + 1, 4] for x in d])
pts = []
for kern, s_of_d in (('insertion s=2', lambda x: 2.0), ('rings s=d+1', lambda x: x + 1), ('q^2 insertion s=4', lambda x: 4.0)):
    for N in range(1, 7):
        # N (d - 1) = s(d)
        if kern == 'rings s=d+1':
            if N == 1: continue
            x = (N + 1) / (N - 1)
        else:
            x = 1 + s_of_d(0) / N
        if 1.2 <= x <= 5: pts.append((kern, N, x, N * (x - 1)))
with open(os.path.join(OUT, 'fig1_coincidences.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['kernel_pole', 'N', 'd', 's', 'integer_d'])
    for kern, N, x, s in pts: w.writerow([kern, N, f'{x:.12g}', f'{s:.12g}', abs(x - round(x)) < 1e-12])
def h2(z): return 1 - z / np.sqrt(1 + z * z)
def h3(z): return np.where(z == 0, 1.0, 1 - z * np.arctan2(1, z))
def R2(z):
    s = np.sqrt(1 + z * z); r = s - z; a = z / s
    return 2 / np.pi * (h2(z) - z / s ** 3 + 2 * a * a * r * np.arctan(r) / (1 + r * r))
def R3(z):
    with np.errstate(divide='ignore', invalid='ignore'):
        v = 4 / (1 + z * z) * ((1 + 3 * z * z) - (2 + 3 * z * z) * z * np.arctan2(1, z)) / (4 * np.pi)
    return v
z = np.linspace(0, 4, 401)
write('fig2_R2.csv', ['z', 'R2', 'h2', 'R3', 'h3'], zip(z, R2(z), h2(z), R3(z), h3(z)))
zl = np.logspace(-2, 2, 401)
write('fig2_R2_loglog.csv', ['z', 'R2', 'static_2_over_pi', 'tail_5_over_12pi_z4'], zip(zl, R2(zl), np.full_like(zl, 2 / np.pi), 5 / (12 * np.pi * zl ** 4)))
q = 3.2; xs = np.linspace(-q - 1.2, 1.2, 161); ys = np.linspace(-1.2, 1.2, 97)
rows = []
for x in xs:
    for y in ys:
        inB = x * x + y * y < 1; inBq = (x + q) ** 2 + y * y < 1
        if inB or inBq:
            delta = x * q + q * q / 2     # p.q + q^2/2 with q along x
            rows.append((x, y, 1 if inB else -1, delta / (q * q / 2)))
write('fig4_delta.csv', ['p_parallel', 'p_perp', 'Delta_p', 'delta_over_q2_half'], rows)
print('wrote', sorted(os.listdir(OUT)))
