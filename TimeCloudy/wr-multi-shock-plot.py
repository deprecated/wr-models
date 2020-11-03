models=[["5.77", "62.3", "17", "3", "11", "100", "1.13e+05", "0.00401", "231"], ["6.32", "68.2", "18.3", "2.5", "9.3", "100", "1.34e+05", "0.00766", "408"], ["7.07", "76.3", "20.2", "2", "7.55", "100", "1.65e+05", "0.0167", "808"], ["7.56", "81.6", "21.5", "1.75", "6.65", "100", "1.87e+05", "0.0279", "1.27e+03"], ["8.16", "88.1", "23", "1.5", "5.74", "100", "2.17e+05", "0.0518", "2.2e+03"], ["8.94", "96.5", "25", "1.25", "4.82", "100", "2.59e+05", "0.114", "4.46e+03"], ["10", "108", "27.8", "1", "3.88", "100", "3.21e+05", "0.322", "1.13e+04"], ["11.5", "125", "31.8", "0.75", "2.93", "100", "4.25e+05", "1.19", "3.65e+04"], ["14.1", "153", "38.7", "0.5", "1.97", "100", "6.34e+05", "5.18", "1.31e+05"], ["18.3", "197", "49.7", "0.3", "1.19", "100", "1.05e+06", "23.4", "4.61e+05"]]
import os
import numpy as np
from scipy import interpolate, optimize, integrate
from astropy.table import Table
from matplotlib import pyplot as plt
import seaborn as sns
import sys
sys.path.append('/Users/will/Work/CLOUDY/claudia/')
from claudia import CloudyModel

CloudyModel.skipsaves.append('continuum')
CloudyModel.skipsaves.remove(".tim")


k = 1.3806503e-16                         # Boltzmann's constant [cgs]
mp = 1.67262158e-24                       # Proton rest mass [cgs]
yHe = 0.162                               # He/H abundance
xHeplus = 1.0                             # He is all singly ionized
xH = 1.0                                  # H is all ionized
mu = 1.0 + 3.0*yHe                        # Mean mass per nucleon
gamma = 5./3.                             # adiabatic index
yr = 3.15576e7                            # Year in seconds
pc = 3.085677582e18                       # Parsec in cm
km = 1e5

phiH = 11.0

# Set up graph for temperature and density
pltfile = 'wr-multi-shock-distance.pdf'
fig, (axtop, axbot) = plt.subplots(2, 1, sharex=True)

pltfile_em = pltfile.replace('distance', 'em-distance')
fig_em, axes_em = plt.subplots(9, 1, sharex=True)
fig_em.set_size_inches(10, 27)
ax6563, ax5007, ax4363, axO3Ha, axLcool, ax5007frac, axOcharge, axTagain, axNagain = axes_em

pltfile_em2 = pltfile.replace('distance', 'em2-distance')
fig_em2, (ax6563_2, ax5007_2) = plt.subplots(2, 1, sharex=True)

pltfile_emcum = pltfile.replace('distance', 'emcum')
fig_cum, [ax5007_cum, ax4363_cum, axT_cum] = plt.subplots(3, 1, sharex=True)

# Loop over all the shock velocities
colors = sns.dark_palette('orange', len(models[:-1]))
print(models)
for row, c in zip(models[:-1], colors):
    M0, u0, v1, n0, n1, N2, T1, dcool, tcool = [float(x) for x in row]
    model_id = 'wr-phi{:02.0f}-shock-v{:03.0f}'.format(phiH, u0)
    label = 'Vs = {:.0f} km/s'.format(u0)

    try:
        m = CloudyModel(model_id, niter=0)
    except:
        print('Failed to read', model_id)
        continue
    print(m)
    # Net cooling coefficient for all times
    NeNp = m.ovr.HII*m.ovr.hden*m.ovr.eden
    Lambda_full = (m.cool.Ctot_ergcm3s - m.cool.Htot_ergcm3s)/NeNp
    # index corresponding to initial post-shock state
    # Heuristic is that it is point where net cooling is highest
    istart = np.argmax(Lambda_full)
    # And corresponding T, which should be more or less T1
    Tstart = m.cool.Temp_K[istart]
    # Photoionization equilibrium T
    Teq = m.cool.Temp_K.min()
    print(istart, Teq, Tstart)
    # Now restrict to the post-shock zone
    T_grid = m.cool.Temp_K[istart:]
    Lambda_grid = Lambda_full[istart:]
    integrand_grid = T_grid**2 / Lambda_grid
    integral_grid = integrate.cumtrapz(integrand_grid, T_grid, initial=0.0)
    T = T_grid
    s = (2./3.)*(Lambda_grid[0]/Tstart**3)*(integral_grid[0] - integral_grid)

    # We need to recalculate tcool and dcool because the Lambda(T1) is
    # now very different - it is much higher because of the under-ionization
    Lambda1 = Lambda_grid[0]
    Pressure = (m.ovr.hden*(1.0 + yHe) + m.ovr.eden)*k*m.cool.Temp_K
    P1 = Pressure[istart]
    L1 = Lambda1*NeNp[istart]
    # Cooling time in seconds
    tcool = P1/((gamma - 1.)*L1)
    # Cooling distance in parsecs
    dcool = v1*km*tcool/pc

    x = np.hstack([[-0.05, 0.0], dcool*s]) 
    axtop.semilogy(x, np.hstack([[Teq, Teq], T]), color=c)
    den = n1*Tstart/T
    axbot.semilogy(x, np.hstack([[n0, n0], den]), label=label, color=c)

    # And plot the emissivities too
    Lcool = m.cool.Ctot_ergcm3s[istart:]*(den/n1)**2
    em5007 = (10**m.ems.O__3_500684A[istart:])*(den/n1)**2 
    em4363 = (10**m.ems.O__3_436321A[istart:])*(den/n1)**2 
    em6563 = (10**m.ems.H__1_656285A[istart:])*(den/n1)**2 
    Ostack = np.vstack([m.ovr["O"+j] for j in "123456"])
    O789 = 1.0 - Ostack.sum(axis=0)
    Ostack = np.vstack([m.ovr["O"+j] for j in "123456"] + [O789])
    Ocharge = np.sum(Ostack*np.arange(7)[:, None], axis=0)[istart:]
    istop = np.nanargmax(s[T > 1.05*Teq])
    ss = s/s[istop]

    # Fractional cumulative emissivity of [O III]
    dss = np.diff(ss, prepend=0.0)
    cumem = integrate.cumtrapz(em5007, s, initial=0.0)
    cumem /= cumem[istop]
    print(cumem[::10])

    ax5007.plot(ss, em5007, color=c)
    ax6563.plot(ss, em6563, color=c)
    ax5007_2.plot(ss, em5007, label=label, color=c)
    ax6563_2.plot(ss, em6563, color=c)
    ax4363.plot(ss, em4363/em5007, label=label, color=c)
    axO3Ha.plot(ss, em5007/em6563, color=c)
    axLcool.plot(ss, Lcool, color=c)
    ax5007frac.plot(ss, em5007/Lcool, color=c)
    axOcharge.plot(ss, Ocharge, color=c)
    axTagain.plot(ss, T, color=c)
    axNagain.plot(ss, den, color=c)

    ax5007_cum.plot(cumem, em5007, label=label, color=c)
    ax4363_cum.plot(cumem, em4363/em5007, color=c)
    axT_cum.plot(cumem, T, color=c)


axtop.set_ylim(9000, 1.5e6)
axbot.set_ylim(0.3, 200.0)
axbot.set_xlabel('Distance, pc')
axbot.set_ylabel('Density, pcc')
axtop.set_ylabel('Temperature, K')
axbot.set_xscale('symlog', linthreshx=1.e-5)
axtop.set_xscale('symlog', linthreshx=1.e-5)
axbot.legend(ncol=2, fontsize='x-small', loc='upper left')
fig.savefig(pltfile)

axes_em[-1].set_xlabel('Fraction of total cooling distance')
ax6563.set_ylabel('Hα 6563 emissivity')
ax4363.legend(ncol=2, fontsize='x-small', loc='lower left')
ax4363.set_ylabel('[O III] 4363/5007 ratio')
axO3Ha.set_ylabel('[O III] 5007/Hα ratio')
ax5007.set_ylabel('[O III] 5007 emissivity')
axLcool.set_ylabel('Total cooling, erg/cm³/s')
axTagain.set_ylabel('Temperature, K')
axNagain.set_ylabel('Total Hydrogen density, /cm³')
ax5007frac.set_ylabel('[O III] 5007 fraction of cooling')
axOcharge.set_ylabel('Mean charge of Oxygen')
for ax in axes_em:
    ax.set_xscale('linear')
    ax.set_yscale('log')
    ax.set_xlim(0.0, 1.2)
ax5007.set_ylim(3e-25, 1.5e-20)
axO3Ha.set_ylim(0.1, 150)
axOcharge.set_yscale('linear')
axOcharge.set_ylim(0.0, 8.0)
for ax in axLcool, ax4363, ax5007, ax6563, axNagain, axTagain, ax5007frac, axO3Ha:
    ax.set_yscale('linear')
    ax.set_ylim(0.0, None)
axO3Ha.set_ylim(0.0, 40.0)



fig_em.tight_layout()
fig_em.savefig(pltfile_em)


ax5007_2.set_ylim(0.0, None)
ax6563_2.set_ylim(0.0, None)
ax5007_2.set_xlabel('Fraction of total cooling distance')
ax5007_2.set_ylabel('[O III] 5007 emissivity')
ax6563_2.set_ylabel('Hα 6563 emissivity')
ax5007_2.set_xlim(0.0, 1.2)
ax5007_2.legend(ncol=2, fontsize='x-small', loc='upper left')
fig_em2.savefig(pltfile_em2)

ax5007_cum.set_xlim(0.0, 1.2)
ax5007_cum.set_ylim(0.0, None)
ax4363_cum.set_ylim(0.0, None)
axT_cum.set_ylim(0.0, 1e5)
fig_cum.tight_layout()
fig_cum.savefig(pltfile_emcum)
