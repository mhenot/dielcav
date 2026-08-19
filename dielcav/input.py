 
import MDAnalysis as mda

import numpy as np
import time
import logging


def compute_dipoles(top, traj, savefile, partial=None, dipoles_struct=None):
    """ 
    Compute the dipole vectors and the position of the center of charge for each time step.
    A .npz file is saved containing : simulation time, box size, dipoles vectors and positions.

    Parameters
    ----------
    top : str
        path of the topology file (for instance .tpr)
    traj : str
        path of the trajectory file (for instance .dcd)
    partial : list of int, optional
        [start, stop] to analyze only a subset of the trajectory.
        by default: all timesteps are considered
    dipoles_struct : list of list of str, optional
        structures of the dipoles to be considered. Molecules can be divided into subdipoles
        for instance [[at1, at2], [at3, at4, at5]] will compute for each molecules, two dipoles : one between at1 and at2 and one between at3, at4 and at5
        for instance [[at1, at2, at3]] will compute one dipole per molecule containing only at1, at2 and at3
        by default, all atoms of each molecule are considered

    Example
    -------
    >>> logging.basicConfig(level=logging.INFO)
    >>> path = '../examples/'
    >>> ref = '320K_1'
    >>> dielcav.compute_dipoles(f'{path}in/water.tpr', f'{path}out/trajectory_{ref}.dcd', savefile=f'{path}out/dipoles_{ref}')

    """


    # Loading the trajectory
    u = mda.Universe(top, traj)
    logging.info(f'{traj} : loaded')

    if dipoles_struct is None: # if the user did not give a structure for the dipoles, we will take the whole molecule
        # we first check that there is only one kind of molecule
        resnames = np.unique(u.residues.resnames)
        if len(resnames) != 1:
            logging.info('More than one type of residues was found')
            raise NotImplementedError

        first_residue = u.select_atoms(f'resname {resnames[0]}').residues[0]
        logging.info(f'One residue found: {resnames[0]}')
        for atom in first_residue.atoms:
            logging.info(f"{atom.name:6s}  type={atom.type:6s}  masse={atom.mass:7.3f}  charge={atom.charge:+.4f}")
        dipoles_struct = [[atom.name for atom in first_residue.atoms]]


    if partial is None:
        partial = []
    t, cell, dipoles, dipoles_pos = [], [], [], []
    k = 0
    for ts in u.trajectory:
        if len(partial)==0 or (k<partial[1] and k>partial[0]):
            t.append(u.trajectory.time)
            cell.append(u.dimensions)
            d_, d_pos_ = [], []
            for dipole_struct in dipoles_struct:
                pos_ch, ch = [], []
                for els in dipole_struct: 
                    sel = u.select_atoms('name {}'.format(els))
                    pos_ch.append(sel.positions * sel.charges[:, np.newaxis])
                    ch.append(sel[0].charge)
                pos_ch, ch = np.array(pos_ch), np.array(ch)
                bary_pos, bary_neg = np.sum(pos_ch[ch>0],axis=0) / np.sum(ch[ch>0]), np.sum(pos_ch[ch<0],axis=0)/np.sum(ch[ch<0])
                d_.append(np.sum(np.abs(ch))/2 * (bary_pos-bary_neg))
                d_pos_.append((bary_pos+bary_neg)/2)
            dipoles.append(np.transpose(np.array(d_), axes=(1,0,2))) # back in order : res, dip
            dipoles_pos.append(np.transpose(np.array(d_pos_), axes=(1,0,2)))
        elif k>partial[1]:
           break
        k += 1
        if k%1000==0:
            logging.info(f'ts {k}')

    np.savez(savefile, t=np.array(t), cell=np.array(cell), dipoles=np.array(dipoles), dipoles_pos=np.array(dipoles_pos))

    