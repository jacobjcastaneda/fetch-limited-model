"""
@Author: Jacob Castaneda
@Date: 8-16-21
"""

import numpy as np
from stompy.grid.unstructured_grid import UnTRIM08Grid
from scipy.interpolate import NearestNDInterpolator, LinearNDInterpolator
import matplotlib.pyplot as plt
import pdb
import time
import pickle

def array_NDInterp(xf, yf, xy=None, Z=None, type='Linear', interpolator=None):
    """
    xy: ndarray of x,y coordinates
    Z: scalar value at x,y
    xf, yf: array of new xf, yf values to be interpolated on (not a meshgrid)
    type: specify if Nearest Neighbor or Linear interpolation desired
    interpolator: if you define an interpolator above the scope of this function
    then pass it in to avoid re-doing work
    -----
    RETURN: interpolated values in ndarray
    """
    if interpolator:
        interp = interpolator
    else:
        if type != 'Linear':
            interp = NearestNDInterpolator(xy, Z)
        else:
            interp = LinearNDInterpolator(xy, Z)
            
    res = np.zeros(xf.shape[0])
    for i in range(xf.shape[0]):
        res[i] = interp(xf[i], yf[i])
    return res


        
def fetch(grd, d, thetaw, res, test=None):
    """
    grd: pass in stompy UnstructuredGrid Object
    d: bathymetry corresponding to grd cell-centers (made for SWAN)
    thetaw: dir- of wind source 0 rad in x+ dir- (+) CCW
    res: resolution for the 1D grid to determine fetch
    
    ----
    RETURN: numpy array of the fetch for cell centers
    """
    ctrs = grd.cells_center()
    L = np.max([np.max(ctrs[:, 0]) - np.min(ctrs[:, 0]),
                np.max(ctrs[:, 1]) - np.min(ctrs[:, 1])])
    F = np.zeros(ctrs.shape[0])
    Df = np.zeros(ctrs.shape[0])
    interp = LinearNDInterpolator(ctrs, d)
    print("Begin computing Fetch...")
    print("Test Mode...") if test else print("\n")
    start = time.time()
    for i in range(ctrs.shape[0]):
        if test:
            if i%test != 0:
                F[i] = -99
                Df[i] = -99
                continue
        xc, yc = ctrs[i, :]
        xe = xc + L*np.cos(thetaw)
        ye = yc + L*np.sin(thetaw)
        xl = np.linspace(xc, xe, res)
        yl = np.linspace(yc, ye, res)
        rl = np.sqrt( (xl-xc)**2 + (yl-yc)**2 )
        dl = array_NDInterp(xl, yl, interpolator=interp)
        dl[np.isnan(dl)] = 0

        imin = np.min(np.where(dl == 0))
        F[i] = rl[imin - 1]
        Df[i] = np.mean(d[:imin]) if len(d[:imin]) > 0 else 0
        # print updates about progress
        if i != 0:
            if i%(5*ctrs.shape[0] // 100) == 0:
                print('{:.2f}  % done on i = {} of {}. Step Time: {} s'.format(100*i/ctrs.shape[0], i, 
                                                                     ctrs.shape[0], time.time()-start))
         
    return F, Df

def read_node_depths(fp):
    """
    fp: pointer to file with depth information... assumes format of depths
    passed to SWAN for UNSTRUC grid

    ---
    RETURN: numpy array of the depths at the grid nodes
    """
    with open(fp, 'rt') as src:
        lines = src.readlines()

    for idx, elem in enumerate(lines):
        tmp = elem.rstrip()
        lines[idx] = float(tmp)
    return np.asarray(lines)


def visualize_grid(grd, vals, xlim, ylim):
    fig, ax = plt.subplots()
    grd.plot_cells(centers=True, values=vals)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)



"""
n=1000 works fairly well for resolution here where the max Length is ~O(100,000)
... result is about +/- 100 m on a grid with O(50m-100m) resolution... seems fair 
for fetch computation...
"""
if __name__ == "__main__":
    start = time.time()
    print("Running in Main...")
    # compute for wind source from West, with n=1000 grid resolution
    THETAW = np.pi
    RES = 1000
    
    fpgrd = 'input/sfei_v20_net_mod.grd'
    fpbathy = 'input/sfei_v20_net_mod_depth'

    grd = UnTRIM08Grid(fpgrd)
    ctrs = grd.cells_center()
    d = array_NDInterp(ctrs[:, 0], ctrs[:, 1], xy=grd.nodes['x'],
                       Z=read_node_depths(fpbathy), type='Nearest')
#    pdb.set_trace()
    F, Df = fetch(grd, d, THETAW, RES)
 #   pdb.set_trace()
    #print(F[:20], Df[:20])
    print('Elapsed Time: ', time.time() - start)
    mapp = {'Fetch': F, 'AveDepth': Df}
    with open('fetch.pkl', 'wb') as dst:
        pickle.dump(mapp, dst)
    #pdb.set_trace()
    print('hello')
