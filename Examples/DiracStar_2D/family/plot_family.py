"""Central amplitude vs ADM mass and binding energy for ground-state Dirac stars (data: dirac_family.dat).
Run: python3 plot_family.py [output.png]"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dirac_family.dat'))
A, M, N, Eb = d[:, 0], d[:, 1], d[:, 2], d[:, 3]

def extrema(y):
    out = []
    for i in range(1, len(A) - 1):
        if (y[i] - y[i - 1]) * (y[i + 1] - y[i]) < 0:
            c = np.polyfit(A[i - 1:i + 2], y[i - 1:i + 2], 2)
            a0 = -c[1] / (2 * c[0])
            out.append((a0, np.polyval(c, a0)))
    return out

tp = extrema(M)                                   # turning points (M and N extremal together)
A_max, M_max = tp[0]
i0 = np.where(np.diff(np.sign(Eb)))[0][0]
A_unb = A[i0] - Eb[i0] * (A[i0 + 1] - A[i0]) / (Eb[i0 + 1] - Eb[i0])

INK, INK2, GRID, SURF = '#1f1f1e', '#5f5e5a', '#e4e3df', '#fcfcfb'
STABLE, BLUE, ORANGE = '#1baf7a', '#2a78d6', '#eb6834'
plt.rcParams.update({'font.size': 10, 'text.color': INK, 'axes.labelcolor': INK})
fig, axs = plt.subplots(2, 1, figsize=(9, 8.2), sharex=True, facecolor=SURF, gridspec_kw=dict(hspace=0.08))
for ax in axs:
    ax.set_facecolor(SURF)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(True, color=GRID, lw=0.8, which='major')
    ax.set_axisbelow(True)
    ax.set_xscale('log')
    ax.axvspan(A[0] * 0.9, A_max, color=STABLE, alpha=0.12, lw=0)
    ax.axvline(A_max, color=STABLE, lw=1.4)
    for a, _ in tp[1:]:
        ax.axvline(a, color=INK2, lw=1, ls=':')

stab = A <= A_max
for ax, y in ((axs[0], M), (axs[1], Eb)):
    ax.plot(A[~stab], y[~stab], color=INK, lw=1.6)
    ax.plot(A[stab], y[stab], color=STABLE, lw=2.6)

# test models evolved with Spherical_BS_Evolver (1D): marker = what happened
tests = [(0.05, 'S1', 'stationary'), (0.15, 'S2', 'stationary'), (0.30, 'U1', 'collapse or migrate'),
         (0.40, 'U2', 'collapse or migrate'), (0.60, 'U3', 'collapses')]
for a, name, out in tests:
    m = np.interp(a, A, M); e = np.interp(a, A, Eb)
    col = STABLE if a < A_max else (ORANGE if 'migrate' in out else BLUE)
    for ax, y in ((axs[0], m), (axs[1], e)):
        ax.plot(a, y, 'o', ms=8, color=col, mec=SURF, mew=1.5, zorder=5)
    axs[0].annotate(name, (a, m), textcoords='offset points', xytext=(0, 9), ha='center', fontsize=9, color=INK)

ax = axs[0]
ax.set_ylabel(r'ADM mass  $M\mu/M_\mathrm{Pl}^2$')
ax.set_ylim(0.2, 0.78)
ax.annotate('maximum mass\n$A_c$ = %.4f, $M$ = %.4f' % (A_max, M_max), (A_max, M_max), textcoords='data',
            xytext=(0.25, 0.30), ha='left', fontsize=9, color=INK, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
ax.text(0.0045, 0.765, 'perturbatively stable branch', color='#0d7a54', fontsize=10, va='top')
ax.text(0.25, 0.215, '1 unstable radial mode', color=INK, fontsize=9, va='bottom')
ax.text(1.32, 0.215, '2 modes', color=INK2, fontsize=9, va='bottom')
ax.text(2.42, 0.215, '3', color=INK2, fontsize=9, va='bottom')
ax.set_title('Ground-state Dirac stars (two spinors, $\\mu = G = 1$): mass and binding energy vs central amplitude',
             loc='left', fontsize=11)

ax = axs[1]
ax.axhline(0, color=INK2, lw=1)
ax.set_ylabel(r'binding energy  $E_b = M - \mu N$')
ax.set_xlabel(r'central amplitude  $A_c = g(0)$   (log scale)')
ax.set_xlim(A[0] * 0.9, 3.0)
ax.annotate('$E_b = 0$ at $A_c$ = %.3f\n(beyond: unbound)' % A_unb, (A_unb, 0), textcoords='offset points',
            xytext=(12, -38), fontsize=9, color=INK, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
ax.annotate('most bound at the maximum mass\n$E_b$ = %.4f' % Eb.min(), (A[np.argmin(Eb)], Eb.min()),
            textcoords='offset points', xytext=(-150, 0), ha='right', va='center', fontsize=9, color=INK,
            arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
from matplotlib.lines import Line2D
handles = [Line2D([], [], color=STABLE, lw=2.6, label='stable branch ($A_c < %.4f$)' % A_max),
           Line2D([], [], color=INK, lw=1.6, label='unstable branch'),
           Line2D([], [], color=INK2, lw=1, ls=':', label='further turning points of $M$ and $N$'),
           Line2D([], [], marker='o', ls='', color=STABLE, ms=7, label='1D runs: stationary'),
           Line2D([], [], marker='o', ls='', color=ORANGE, ms=7, label='1D runs: collapse (a > 0) or migrate (a < 0)'),
           Line2D([], [], marker='o', ls='', color=BLUE, ms=7, label='1D runs: collapse (unbound)')]
ax.legend(handles=handles, frameon=False, fontsize=8.5, loc='upper left')
fig.savefig(sys.argv[1] if len(sys.argv) > 1 else 'dirac_family.png', dpi=150, facecolor=SURF, bbox_inches='tight')
print('A_max %.4f M_max %.5f A(E_b=0) %.4f, other turning points %s' % (A_max, M_max, A_unb, [(round(a, 3), round(m, 4)) for a, m in tp[1:]]))
