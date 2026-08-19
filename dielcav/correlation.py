
import matplotlib.pyplot as plt
import numpy as np
import time
import logging

from .utils import make_dit, legendre

def correlation_self(fname, Ndit=210, Li0_step=50, legendre_rank=1, save=True):
    """ 
    Compute the self dipole correlation function

    Parameters
    ----------
    fname : str
        .npz file containing the dipole vectors and positions
    Ndit : int
        nnumber of log spaced point
    Li0_step : int
        frame step
    legendre_rank : int, optional
        rank of the Legendre polynomial (default: 1)
    save : bool, optional
        save results into a .txt file (default: True)

    Returns
    -------
    t: np.array
        time
    c: np.array
        c(t) is the self dipole correlation function

    Example
    -------
    >>> t, c = dielcav.correlation_self(f'./out/dipoles_320K_1', Ndit=150, Li0_step=1, save=True)
    """

    with np.load(f'{fname}.npz') as data:
        t = data['t']
        dips = data['dipoles']

    imax, Nmol, Ndip_mol = dips.shape[:3]
    dips = dips.reshape(len(dips), Nmol*Ndip_mol, 3)
   
    dit = make_dit(imax-1, Ndit)
    Li0 = range(0, imax//2, Li0_step)

    res = np.empty((len(Li0), len(dit)))
    res[:] = np.nan

    for j, i0 in enumerate(Li0):
        dit0 = dit + i0
        dit0 = dit0[dit0<imax]
        corr = legendre(np.sum(dips[i0,:,:][np.newaxis,:,:]*dips[dit0,:,:], axis=-1), rank=legendre_rank)
        res[j, :len(corr)] = np.mean(corr, axis = 1)

    dt = np.around(t[dit]-t[0],4)
    c = np.nanmean(res, axis=0)

    if save:
        np.savetxt(f'{fname}_self_corr.txt',
               np.array([dt, c]).T)
    return dt, c


def correlation_total_box(fname, Ndit=210, Li0_step=50, save=True):
    """ 
    Compute the total dipole correlation function

    Parameters
    ----------
    fname : str
        .npz file containing the dipole vectors and positions
    Ndit : int
        number of log spaced point
    Li0_step : int
        frame step
    save : bool, optional
        save results into a .txt file (default: True)

    Returns
    -------
    t: np.array
        time
    c: np.array
        c(t) is the total dipole correlation function

    Example
    -------
    >>> t, c = dielcav.correlation_total_box(f'./out/dipoles_320K_1', Ndit=150, Li0_step=1, save=True)
    """
    with np.load(f'{fname}.npz') as data:
        t = data['t']
        dips = data['dipoles']

    imax, Nmol, Ndip_mol = dips.shape[:3]
    dips = dips.reshape(len(dips), Nmol*Ndip_mol, 3)
   
    dit = make_dit(imax-1, Ndit)
    Li0 = range(0, imax//2, Li0_step)

    M = np.sum(dips, axis=1)

    res = np.empty((len(Li0), len(dit)))
    res[:] = np.nan

    for j, i0 in enumerate(Li0):
        dit0 = dit + i0
        dit0 = dit0[dit0<imax]
        corr = np.sum(M[i0,:][np.newaxis,:]*M[dit0,:], axis=-1)
        res[j, :len(corr)] = corr
    dt = np.around(t[dit]-t[0],4)
    c = np.nanmean(res, axis=0)
    if save:
        np.savetxt(f'{fname}_self_total.txt',
        np.array([dt, c]).T)
    return dt, c

