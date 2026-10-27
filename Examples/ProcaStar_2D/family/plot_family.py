"""Central amplitude vs ADM mass and binding energy for ground-state mini Proca stars (data: proca_family.dat,
regenerate with make_family.cpp). Run: python3 plot_family.py [output.png]"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

d = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'proca_family.dat'))
A, M, N, Eb, f0 = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 6]


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
f0_max = np.interp(A_max, A, f0)

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

# the default star of the 2D/3D examples, and the unstable models evolved with Spherical_BS_Evolver (1D, perturbation
# amplitude a = +-0.01, width 3): a > 0 formed an apparent horizon, a < 0 expanded off the unstable branch
tests = [(0.02, 'default', 'example'), (0.06, 'U1', '1D'), (0.1005, 'U2', '1D'), (0.1995, 'U3', '1D (unbound)')]
for a, name, kind in tests:
    m = np.interp(a, A, M); e = np.interp(a, A, Eb)
    col = STABLE if kind == 'example' else (ORANGE if a < A_unb else BLUE)
    for ax, y in ((axs[0], m), (axs[1], e)):
        ax.plot(a, y, 'o', ms=8, color=col, mec=SURF, mew=1.5, zorder=5)
    axs[0].annotate(name, (a, m), textcoords='offset points', xytext=(0, 9), ha='center', fontsize=9, color=INK)

ax = axs[0]
ax.set_ylabel(r'ADM mass  $M\mu/M_\mathrm{Pl}^2$')
ax.set_ylim(0.3, 1.12)
ax.annotate('maximum mass\n$A_c$ = %.4f ($f_0$ = %.4f), $M$ = %.4f' % (A_max, f0_max, M_max), (A_max, M_max),
            textcoords='data', xytext=(0.042, 0.44), ha='left', fontsize=9, color=INK,
            arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
ax.text(0.00105, 1.105, 'perturbatively stable branch', color='#0d7a54', fontsize=10, va='top')
ax.text(0.045, 0.31, '1 unstable radial mode', color=INK, fontsize=9, va='bottom')
ax.text(tp[1][0] * 1.06, 0.31, '2 modes', color=INK2, fontsize=9, va='bottom')
ax.set_title('Ground-state Proca stars ($\\mu = G = 1$): mass and binding energy vs central amplitude', loc='left',
             fontsize=11)

ax = axs[1]
ax.axhline(0, color=INK2, lw=1)
ax.set_ylim(-0.046, 0.092)
ax.set_ylabel(r'binding energy  $E_b = M - \mu N$')
ax.set_xlabel(r'central amplitude  $A_c = f(0)/\alpha(0) = |X_\phi(0)|$   (log scale)')
ax.set_xlim(A[0] * 0.9, 4.2)
ax.annotate('$E_b = 0$ at $A_c$ = %.3f\n(beyond: unbound)' % A_unb, (A_unb, 0), textcoords='offset points',
            xytext=(14, -40), fontsize=9, color=INK, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
ax.annotate('most bound (at the maximum mass): $E_b$ = %.4f' % Eb.min(), (A[np.argmin(Eb)], Eb.min()),
            textcoords='offset points', xytext=(12, -20), ha='left', va='center', fontsize=9, color=INK,
            arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
handles = [Line2D([], [], color=STABLE, lw=2.6, label='stable branch ($A_c < %.4f$)' % A_max),
           Line2D([], [], color=INK, lw=1.6, label='unstable branch'),
           Line2D([], [], color=INK2, lw=1, ls=':', label='second turning point of $M$ and $N$ (minimum)'),
           Line2D([], [], marker='o', ls='', color=STABLE, ms=7, label='default star of the 2D/3D examples'),
           Line2D([], [], marker='o', ls='', color=ORANGE, ms=7, label='1D runs: collapse (a > 0) or expand (a < 0)'),
           Line2D([], [], marker='o', ls='', color=BLUE, ms=7, label='1D run, unbound: collapse (a > 0) or disperse (a < 0)')]
ax.legend(handles=handles, frameon=False, fontsize=8.5, loc='upper left')
fig.savefig(sys.argv[1] if len(sys.argv) > 1 else 'proca_family.png', dpi=150, facecolor=SURF, bbox_inches='tight')
print('A_max %.4f (f0 %.5f) M_max %.5f N_max %.5f, A(E_b=0) %.4f, other turning points %s'
      % (A_max, f0_max, M_max, np.interp(A_max, A, N), A_unb, [(round(a, 3), round(m, 4)) for a, m in tp[1:]]))
