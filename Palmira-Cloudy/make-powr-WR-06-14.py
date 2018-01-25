import numpy as np
import textwrap

modelnames = 'WR-06-14',
outfilename = 'powr-WR-06-14.ascii'

output = ''


SMALL_FLOAT = 1.e-30
wavelengths = None
for mn in modelnames:
    sedfile = 'WR-06-14/sed.txt'
    log_wavs, log_fluxes = np.loadtxt(sedfile, unpack=True)
    wavs = 10**log_wavs
    fluxes = 10**log_fluxes
    fluxes = np.maximum(SMALL_FLOAT, fluxes)
    if wavelengths is None:
        wavelengths = wavs
        s = ' '.join(['{:.6e}'.format(x) for x in wavelengths])
        output += textwrap.fill(s) + '\n'
    else:
        assert np.alltrue(wavs == wavelengths)
    s = ' '.join(['{:.6e}'.format(x) for x in fluxes])
    output += textwrap.fill(s) + '\n'

header = ''
header += '{magicnumber}'.format(magicnumber=20060612) + '\n'
header += '{ndim}'.format(ndim=1) + '\n'
header += '{npar}'.format(npar=1) + '\n'
header += '{label}'.format(label='age') + '\n'
header += '{nmod}'.format(nmod=len(modelnames)) + '\n'
header += '{nfreq}'.format(nfreq=len(wavelengths)) + '\n'
header += 'lambda' + '\n'
header += '1.0' + '\n'
header += 'F_lambda' + '\n'
header += '1.0' + '\n'
header += '1.0' + '\n'

with open(outfilename, 'w') as f:
    f.write(header+output)
