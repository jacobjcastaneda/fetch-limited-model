"""
@Author: Jacob Castaneda
@Date: 8-16-21
"""

import numpy as np
import numpy.ma as ma
from stompy.grid.unstructured_grid import UnTRIM08Grid
from scipy.interpolate import NearestNDInterpolator, LinearNDInterpolator
import matplotlib.pyplot as plt
from shapely import geometry
import pdb
import copy
import time
import pickle

def array_NDInterp(xf, yf, xy=None, Z=None, type='Linear', interpolator=None, mask=None):
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
    if mask is not None:
        res[~mask] = np.nan
        return res
    else:
        return res


    
    
    
    
def find_bdry(grd, poly, xl, yl, L, xc, yc):
    lnLst = [geometry.LineString(zip(xl, yl))]
    lnlngth = len([x for x in lnLst[-1].coords])
    while lnlngth > 3:
        mdpt = lnlngth // 2
        lnLst.append(geometry.LineString(lnLst[-1].coords[:mdpt]))
        if poly.contains(lnLst[-1]):
            break
        lnlngth = len([x for x in lnLst[-1].coords])
    result = np.asarray(lnLst[-1].coords)
    return result[:, 0], result[:, 1]
        
def distance(x1, y1, x2, y2):
    return np.sqrt( (x2-x1)**2 + (y2-y1)**2 )

def fetch(grd, d, thetaw, res, mask=None, plot=False):
    """
    grd: pass in stompy UnstructuredGrid Object
    d: bathymetry corresponding to grd cell-centers (made for SWAN)
    thetaw: dir- of wind source 0 rad in x+ dir- (+) CCW
    res: resolution for the 1D grid to determine fetch
    
    ----
    RETURN: numpy array of the fetch for cell centers
    """
    ctrs = grd.cells_center()
    L_init = np.max([np.max(ctrs[:, 0]) - np.min(ctrs[:, 0]),
                np.max(ctrs[:, 1]) - np.min(ctrs[:, 1])])
    F = np.zeros(ctrs.shape[0])
    Df = np.zeros(ctrs.shape[0])
    
    #test = []
    
    #interp = LinearNDInterpolator(ctrs, d)
    
    nodes = grd.nodes['x'][grd.boundary_cycle()]
    poly = geometry.Polygon(nodes)
    
    print("Begin computing Fetch...")
    start = time.time()
    #pdb.set_trace()
    for i in range(ctrs.shape[0]):
        L = copy.deepcopy(L_init)
        if mask is not None:
            if ~mask[i]:
                F[i] = -99
                Df[i] = -99
                continue
        xc, yc = ctrs[i, :]
        xe = xc + L*np.cos(thetaw)
        ye = yc + L*np.sin(thetaw)
        xl = np.linspace(xc, xe, res)
        yl = np.linspace(yc, ye, res)
        ###
        #LOOP = True
        #count = 0
        #while LOOP:
        if mask[i]:
            #pdb.set_trace()
            pass
        start = time.time()
        for i in range(10):
            if i == 3:
                pass
                #pdb.set_trace()
            # cuts line in half until within polygon boundary
            xl, yl = find_bdry(grd, poly, xl, yl, L, xc, yc)
            # new length is average of old L with length of line found within Polygon
            L = np.mean([L, distance(xc, yc, xl[-1], yl[-1])])
            xe = xc + L*np.cos(thetaw)
            ye = yc + L*np.sin(thetaw)
            xl = np.linspace(xc, xe, res)
            yl = np.linspace(yc, ye, res)
            ln = geometry.LineString(zip(xl, yl))
            if poly.contains(ln):
                while poly.contains(ln):
                    L = distance(xc, yc, xl[-1], yl[-1])*1.05
                    xe = xc + L*np.cos(thetaw)
                    ye = yc + L*np.sin(thetaw)
                    xl = np.linspace(xc, xe, res)
                    yl = np.linspace(yc, ye, res)
                    ln = geometry.LineString(zip(xl, yl))
                break
                
            if i >= 3 and res < 100000:
                # begin increasing resolution
                res =  res*10
            if plot and i > 10:
                fig, ax = plt.subplots()
                grd.plot_edges(ax=ax)
                ax.scatter([xc,xe], [yc,ye], c='r', zorder=10)
                ax.plot(xl, yl, 'k--')
                plt.show()
                
        print('Fetch for Cell time: {}'.format(time.time() -  start))
        
            
            
        
        
        
        
        #####
        #rl = np.sqrt( (xl-xc)**2 + (yl-yc)**2 )
        #dl = array_NDInterp(xl, yl, interpolator=interp)
        #dl[np.isnan(dl)] = 0
        
        #imin = np.min(np.where(dl == 0))
        F[i] = distance(xc, yc, xl[-1], yl[-1]) #rl[imin - 1]
        Df[i] = 0 #np.mean(d[:imin]) if len(d[:imin]) > 0 else 0
        
        #test.append((xl, yl, rl, dl, imin))
        if plot:
            fig, ax = plt.subplots()
            grd.plot_edges(ax=ax)
            ax.scatter([xc,xe], [yc,ye], c='r', zorder=10)
            ax.plot(xl, yl, 'k--')
            plt.show()
        # print updates about progress
        if i != 0:
            if i%(ctrs.shape[0] // 100) == 0:
                print('{:.2f}  % done on i = {} of {}. Step Time: {} s'.format(100*i/ctrs.shape[0], i, 
                                                                    ctrs.shape[0], time.time()-start))
    if mask is not None:
        return ma.array(F, mask=~mask), ma.array(Df, mask=~mask) # , test
    else:
        return F, Df  # , test

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

def get_mask(xy, xlim, ylim, lighten=(None, 1)):
    """
    xy: ndarray of the x,y coordinates for points passed in
    xlim: ndarray limits on x to filter out ~AOI
    ylim: ndarray limites on y to filter out ~AOI
    lighten: tuple(boolean, iterations) whether to lighten the point 
    density and by how many iterations. Each iteration halves point density
    -------
    RETURN: A mask of True/False for the points passed in (built to handle 
    cell-centers of Unstructured Grid)... False indicates point is undesired, 
    and True indicates point is in AOI
    NOTE: this is opposite of mask for numpy.ma, should therefore set mask=~mask
    for numpy masked array
    """
    TWO = 2
    x1, x2 = xlim
    y1, y2 = ylim
    res = []
    for i in range(xy.shape[0]):
        if (xy[i, 0] > x1 and xy[i, 0] < x2) and (xy[i, 1] > y1 and xy[i, 1] < y2):
            res.append(True)
        else:
            res.append(False)
    # if lighten != None, then make more sparse by factor given
    if lighten[0]:
        for i in range(lighten[1]):
            cnt = 0
            for j in range(len(res)):
                if res[j]:
                    if cnt%TWO == 0:
                        res[j] = False
                    cnt = cnt + 1
    return np.asarray(res)
       
    

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
    RES = 100
    
    fpgrd = 'input/sfei_v20_net_mod.grd'
    fpbathy = 'input/sfei_v20_net_mod_depth'

    grd = UnTRIM08Grid(fpgrd)
    ctrs = grd.cells_center()
    m = get_mask(ctrs, np.array([550000, 577000]), np.array([4155000, 4175000]), lighten=(True, 10))
    d = array_NDInterp(ctrs[:, 0], ctrs[:, 1], xy=grd.nodes['x'],
                       Z=read_node_depths(fpbathy), type='Nearest', mask=m)
    #pdb.set_trace()
    F, Df = fetch(grd, d, THETAW, RES, mask=m, plot=True)
 #   pdb.set_trace()
    #print(F[:20], Df[:20])
    print('Elapsed Time: ', time.time() - start)
    #mapp = {'Fetch': F, 'AveDepth': Df}
    #with open('fetch.pkl', 'wb') as dst:
    #    pickle.dump(mapp, dst)
    #pdb.set_trace()
    print('hello')
