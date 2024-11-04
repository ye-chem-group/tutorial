import sys
import numpy as np

from pyscf import lib

log = lib.logger.Logger(sys.stdout, 6)


# fit BM eos
def e_from_bm(v, e0, v0, b0, b0p):
    vr = (v0 / v)**(2./3.)
    e = e0 + (9./16.)*v0*b0 * (
        (vr-1)**3.*b0p + (vr-1)**2.*(6.-4.*vr)
    )
    return e
def fit_eos_bm(vs, es, ret_raw=False, ndata=None):
    """ Fit equation of state according to the BM formula

    Args:
        vs (numpy array):
            cell volumes in Ang^3
        es (numpy array):
            cell energies in Hartree
        ret_raw (bool; default: False):
            See Returns.

    Returns:
        Return fitted equilibrium cell energy (in Hartree), cell volume
        (in Ang^3) and bulk modulus (in GPa).
        If ret_raw is True, output of scipy curve_fit is also returned.
    """
    vs = np.asarray(vs)
    es = np.asarray(es)
    n = len(vs)

# keep ndata (if provided)
    if isinstance(ndata, int):
        nkeep = min(ndata, n)
    else:
        nkeep = n

# dim check
    if vs.ndim != 1 or es.ndim != 1 or len(es) != n:
        raise RuntimeError('Input arrays must be 1d of same length.')
    if nkeep < 4:
        raise RuntimeError('At least 4 data points must be provided '
                           '(got only %d)' % n)
# sort by volume
    keep = np.argsort(es)[:nkeep]
    vs = np.array(vs[keep])
    es = np.array(es[keep])
    order = np.argsort(vs)
    vs = np.array(vs[order])
    es = np.array(es[order])

# init guess
    order = np.argsort(es)
    imin = order[0]
    v0 = vs[imin]
    e0 = es[imin]
    if imin == 0:
        i0,i1,i2 = imin,imin+1,imin+2
    elif imin == nkeep-1:
        i0,i1,i2 = imin-2,imin-1,imin
    else:
        i0,i1,i2 = imin-1,imin,imin+1
    b0 = (es[i0] + es[i2] - es[i1]*2.) / (vs[i2]-vs[i1])**2. * v0
    b0p = 0.1
    p0 = [e0, v0, b0, b0p]

    from scipy.optimize import curve_fit
    p, varp = curve_fit(e_from_bm, vs, es, p0=p0)

    e = p[0]
    v = p[1]
    b = p[2] * 4359.74  # to GPa
    ret = (e, v, b)
    if ret_raw:
        ret += (p, varp)

    return ret


if __name__ == '__main__':
    try:
        finp = sys.argv[1]
    except:
        log.info('Usage: file_energy_vs_cellvol')
        sys.exit(1)

    vs, es = np.loadtxt(finp).T
    e0, v0, b0 = fit_eos_bm(vs, es)
    a0 = (v0*4)**(1./3.)    # cubic lattice constant
    log.info('Equilibrium energy= %.10f  lattice const= %.3f  bulk modulus= %.1f', e0, a0, b0)

    # plot
    param = fit_eos_bm(vs, es, ret_raw=True)[3]
    vs_fit = np.linspace(vs.min(), vs.max(), 100)
    es_fit = e_from_bm(vs_fit, *param)

    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=(3.2, 2.8))
    ax = fig.gca()

    ys = (es - e0) * 27.211399  # relative energy in eV
    ys_fit = (es_fit - e0) * 27.211399

    ax.plot(vs, ys, 'o', color='blue', label='data')
    ax.plot(vs_fit, ys_fit, '--', color='gray', label='fit', zorder=-1)

    ax.set_xlabel('Cell volume (Ang^3)')
    ax.set_ylabel('Relative energy (eV)')

    ax.legend(frameon=False)

    plt.tight_layout()
    plt.savefig('fit.pdf')
    plt.close()
