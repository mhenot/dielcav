import numpy as np

def legendre(x, rank):
    """
    Legendre polynomial

    Parameters
    ----------
    x: float
        abscisse (dans [-1 ,1])
    rank: int
        polynomial rank

    Returns
    -------
    y: float
        P(x)
    """
    if rank == 1:
        return x
    elif rank == 2:
        return 1 / 2 * (3 * x**2 - 1)
    else:  
        raise NotImplementedError

def make_dit(imax, Ndit):
    """
    Makes an array of index for the spacing first linearly spaced between 0 and 10, then logarithmically spaced.

    Parameters
    ----------
    imax: int
        last value
    Ndit: int
        number of values

    Returns
    -------
    dit: np.array of int
        indexes
    """
    return np.unique(np.concatenate([np.arange(10), np.geomspace(10, imax - 1, Ndit).astype(int), [imax - 1]]))

def get_pos_dec(pos, cells, j):
    """
    Return all dipole positions shifted so that dipole j is centered

    Parameters
    ----------
    pos : np.ndarray((Nt, Ndipoles, 3)) or (Ndipoles, 3)
        dipoles positions
    cells : np.ndarray((Nt, 3+)) or np.ndarray((3+))
        simulation box size
    j : int
        index of the dipole to be centered

    Returns
    -------
    pos_dec : same as pos
        shifted dipoles positions
    """
    
    if len(cells.shape) == 1:
        dec = cells[:3]/2 - pos[j,:]
        pos_dec = pos + dec[np.newaxis,:]
        pos_dec[pos_dec[:,0]>cells[0, np.newaxis],0] -= cells[0]
        pos_dec[pos_dec[:,0]<0,0] += cells[0]
        pos_dec[pos_dec[:,1]>cells[0, np.newaxis],1] -= cells[1]
        pos_dec[pos_dec[:,1]<0,1] += cells[1]
        pos_dec[pos_dec[:,2]>cells[0, np.newaxis],2] -= cells[2]
        pos_dec[pos_dec[:,2]<0,2] += cells[2]
    else:
        dec = cells[:,:3]/2 - pos[:,j,:]
        pos_dec = pos + dec[:,np.newaxis,:]
        for i in range(pos_dec.shape[0]):
            pos_dec[i,pos_dec[i,:,0]>cells[i,0, np.newaxis],0] -= cells[i,0]
            pos_dec[i,pos_dec[i,:,0]<0,0] += cells[i,0]
            pos_dec[i,pos_dec[i,:,1]>cells[i,0, np.newaxis],1] -= cells[i,1]
            pos_dec[i,pos_dec[i,:,1]<0,1] += cells[i,1]
            pos_dec[i,pos_dec[i,:,2]>cells[i,0, np.newaxis],2] -= cells[i,2]
            pos_dec[i,pos_dec[i,:,2]<0,2] += cells[i,2]
    return pos_dec