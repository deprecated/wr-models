import numpy as np
import textwrap

modelnames = 'WR-6-14', 'WR-10-16', 'WR-12-18'
outfilename = 'powr-models.ascii'

output = ''
output += '{magicnumber}'.format(magicnumber=20060612) + '\n'
output += '{ndim}'.format(ndim=1) + '\n'
output += '{npar}'.format(npar=1) + '\n'
output += '{label}'.format(label='age') + '\n'
output += '{nmod}'.format(nmod=len(modelnames)) + '\n'
output += '{nfreq}'.format(nfreq=852) + '\n'
output += 'lambda' + '\n'
output += '1.0' + '\n'
output += 'F_lambda' + '\n'
output += '1.0' + '\n'
output += '1.0 2.0 3.0' + '\n'


SMALL_FLOAT = 1.e-30
wavelengths = None
for mn in modelnames:
    sedfile = mn + '/sed.txt'
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

with open(outfilename, 'w') as f:
    f.write(output)
