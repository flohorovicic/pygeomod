"""Analysis and modification of structural data exported from GeoModeller

All structural data from an entire GeoModeller project can be exported into ASCII
files using the function in the GUI:

Export -> 3D Structural Data

This method generates files for defined geological parameters:
"Points" (i.e. formation contact points) and
"Foliations" (i.e. orientations/ potential field gradients). 

Exported parameters include all those defined in sections as well as 3D data points.

This package contains methods to check, visualise, and extract/modify parts of these
exported data sets, for example to import them into a different Geomodeller project.
"""

# import os, sys
import numpy as np
import numpy.random
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

class Struct3DPoints():
    """Class container for 3D structural points data sets"""
    
    def __init__(self, **kwds):
        """Structural points data set
        
        **Optional keywords**:
            - *filename* = string : filename of csv file with exported points to load
        """
        # store point information in purpose defined numpy record
        self.ptype = np.dtype([('x', np.float32),
                          ('y', np.float32),
                          ('z', np.float32),
                          ('formation', np.str_, 32)])

        if kwds.has_key("filename"):
            self.filename = kwds['filename']
            # read data
            self.parse()
            self.get_formation_names()
            self.get_range()
            
    def parse(self):
        """Parse filename and load data into numpy record
        
        The point information is stored in a purpose defined numpy record
        self.points
        """
        f = open(self.filename, "r")
        lines = f.readlines()
        self.header = lines[0]
        # determine position of elements in header (for extension to foliations, etc.)
        h_elem = np.array(self.header.rstrip().split(','))
        x_id = np.where(h_elem == 'X')[0]
        y_id = np.where(h_elem == 'Y')[0]
        z_id = np.where(h_elem == 'Z')[0]
        form_id = np.where(h_elem == 'formation')[0]
        # print x_id
        # create numpy array for points
        self.len = (len(lines)-1)
        self.points = np.ndarray(self.len, dtype = self.ptype)
        for i,line in enumerate(lines[1:]):
            l = line.rstrip().split(',')
            self.points[i]['x'] = float(l[x_id]) 
            self.points[i]['y'] = float(l[y_id]) 
            self.points[i]['z'] = float(l[z_id]) 
            self.points[i]['formation'] = l[form_id]
    
    def get_formation_names(self):
        """Get names of all formations that have a point in this data set
        and store in:
        
        self.formation_names
        """
#        self.formation_names = np.unique(self.formations)
        self.formation_names = np.unique(self.points[:]['formation'])
        
    def get_range(self):
        """Update min, max for all coordinate axes and store in
        self.xmin, self.xmax, ..."""
        self.xmin = np.min(self.points['x'])
        self.ymin = np.min(self.points['y'])
        self.zmin = np.min(self.points['z'])
        self.xmax = np.max(self.points['x'])
        self.ymax = np.max(self.points['y'])
        self.zmax = np.max(self.points['z'])
        
    def create_formation_subset(self, formation_names):
        """Create a subset (as another Struct3DPoints object) with specified formations only
        
        **Arguments**:
            - *formation_names* : list of formation names
            
        **Returns**:
            Struct3DPoints object with subset of points
        """
        # create new object
        # reference to own class type for consistency with Struct3DFoliations
        pts_subset = self.__class__()        
        
        # determine ids for all points of these formations:
        ids = np.ndarray((self.len), dtype='bool') 
        ids[:] = False
        if type(formation_names) == list:
            for formation in formation_names:
                ids[self.points['formation'] == formation] = True
        else:
            ids[self.points['formation'] == formation_names] = True
        
        # new length is identical to sum of ids bool array (all True elements)
        pts_subset.len = np.sum(ids)
         
        # extract points
        pts_subset.points = self.points[ids]
        
        # update range
        pts_subset.get_range()
        
        # update formation names
        pts_subset.get_formation_names()
        
        # get header from original
        pts_subset.header = self.header
        return pts_subset
    
    def remove_formations(self, formation_names):
        """Remove points for specified formations from the point set
        
        This function can be useful, for example, to remove one formation, perform
        a thinning operation, and then add it back in with the `combine_with` function.
        
        **Arguments**:
            - *formation_names* = list of formations to be removed (or a single string to
            remove only one formation)
        """
        # Note: implementation is very similar to create_formation_subset, only inverse
        # and changes in original point set!
        
        # determine ids for all points of these formations:
        ids = np.ndarray((self.len), dtype='bool') 
        ids[:] = True
        if type(formation_names) == list:
            for formation in formation_names:
                ids[self.points['formation'] == formation] = False
        else:
            ids[self.points['formation'] == formation_names] = False
 
        self.len = np.sum(ids)
         
        # extract points
        self.points = self.points[ids]
        
        # update range
        self.get_range()
        
        # update formation names
        self.get_formation_names()
        
    def rename_formations(self, rename_dict):
        """Rename formation according to assignments in dictionary
        
        Mapping in dictionary is of the form:
            old_name_1 : new_name_1, old_name_2 : new_name_2, ...
        """
        for k,v in rename_dict.items():
            
            print("Change name from %s to %s" % (k,v))
            for p in self.points:
                if p['formation'] == k: p['formation'] = v
         
        # update formation names
        self.get_formation_names()
        
    def extract_range(self, **kwds):
        """Extract subset for defined ranges
        
        Pass ranges as keywords: from_x, to_x, from_y, to_y, from_z, to_z
        All not defined ranges are simply kept as before
        
        **Returns**:
            pts_subset : Struct3DPoints data subset
        """
        from_x = kwds.get("from_x", self.xmin)
        from_y = kwds.get("from_y", self.ymin)
        from_z = kwds.get("from_z", self.zmin)
        to_x = kwds.get("to_x", self.xmax)
        to_y = kwds.get("to_y", self.ymax)
        to_z = kwds.get("to_z", self.zmax)
        
        # create new object
        # pts_subset = Struct3DPoints()
        pts_subset = self.__class__()

        # determine ids for points in range
        ids = np.ndarray((self.len), dtype='bool') 
        ids[:] = False
        ids[(self.points['x'] >= from_x) *
            (self.points['y'] >= from_y) *
            (self.points['z'] >= from_z) *
            (self.points['x'] <= to_x) *
            (self.points['y'] <= to_y) *
            (self.points['z'] <= to_z)] = True
        
        # new length is identical to sum of ids bool array (all True elements)
        pts_subset.len = np.sum(ids)
         
        # extract points
        pts_subset.points = self.points[ids]
        
        # update range
        pts_subset.get_range()
        
        # update formation names
        pts_subset.get_formation_names()
        
        # get header from original
        pts_subset.header = self.header
        
        return pts_subset
        
    
    def thin(self, nx, ny, nz, **kwds):
        """Thin data for one formations on grid with defined number of cells and store as subset
        
        **Arguments**:
            - *nx*, *ny*, *nz* = int : number of cells in each direction for thinning grid
        
        The thinning is performed on a raster and not 'formation-aware', 
        following this simple procedure:

        (1) Iterate through grid
        (2) If multiple points for formation in this cell: thin
        (3a) If thin: Select one point in cell at random and keep this one!
        (3b) else: if one point in raneg, keep it!
        
        Note: Thinning is performed for all formations, so make sure to create a subset
        for a single formation first! 
        
        **Returns**:
            pts_subset = Struct3DPoints : subset with thinned data for formation
        """
        # DEVNOTE: This would be an awesome function to parallelise! Should be quite simple!
        
        # first step: generate subset
        # pts_subset = self.create_formation_subset([formation])
         
        # create new pointset:
        # reference to own class type for consistency with Struct3DFoliations
        pts_subset = self.__class__()
         
        # determine cell boundaries of subset for thinning:
        delx = np.ones(nx) * (self.xmax - self.xmin) / nx
        bound_x = self.xmin + np.cumsum(delx)
        dely = np.ones(ny) * (self.ymax - self.ymin) / ny
        bound_y = self.ymin + np.cumsum(dely)
        delz = np.ones(nz) * (self.zmax - self.zmin) / nz
        bound_z = self.zmin + np.cumsum(delz)
        
        ids_to_keep = []
        
        for i in range(nx-1):
            for j in range(ny-1):
                for k in range(nz-1):
                    # determin number of points in this cell
                    ids = np.ndarray((self.len), dtype='bool') 
                    ids[:] = False
                    ids[(self.points['x'] > bound_x[i]) *
                        (self.points['y'] > bound_y[j]) *
                        (self.points['z'] > bound_z[k]) *
                        (self.points['x'] < bound_x[i+1]) *
                        (self.points['y'] < bound_y[j+1]) *
                        (self.points['z'] < bound_z[k+1])] = True
                    if np.sum(ids) > 1:
                        # Thinning required!
                        # keep random point
                        ids_to_keep.append(numpy.random.choice(np.where(ids)[0]))
                        
#                        pts_subset.points[nx * ny * i + ny * j + k] = self.points[id_to_keep]
                        # assign to new pointset:
                    elif np.sum(ids) == 1: 
                        # keep the one point, of course!
#                        pts_subset.points[nx * ny * i + ny * j + k] = self.points[ids[0]]
                        ids_to_keep.append(ids[0])
        
        # now get points for all those ids:
        # extract points
        pts_subset.points = self.points[np.array(ids_to_keep)]
        
        
        # update range
        pts_subset.get_range()
        
        # update length
        pts_subset.len = len(pts_subset.points)
        
        # update formation names
        pts_subset.get_formation_names()
        
        # get header from original
        pts_subset.header = self.header
                        
        return pts_subset
                        
                                    
    def combine_with(self, pts_set):
        """Combine this point set with another point set
        
        **Arguments**:
            - *pts_set* = Struct3DPoints : points set to combine
        """
        self.points = np.concatenate((self.points, pts_set.points))
        # update range and everything
        self.get_range()
        self.get_formation_names()
        self.len = len(self.points)
        
        
    def plot_plane(self, plane=('x','y'), **kwds):
        """Create 2-D plots for point distribution
        
        **Arguments**:
            - *plane* = tuple of plane axes directions, e.g. ('x','y') (default)
        
        **Optional Keywords**:
            - *ax* = matplotlib axis object: if provided, plot is attached to this axis
            - *formation_names* = list of formations : plot only points for specific formations
        """
        color = kwds.get("color", 'b')
        if kwds.has_key("ax"):
            # axis is provided, attach here
            ax = kwds['ax']
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        
        if kwds.has_key("formation_names"):
            pts_subset = self.create_formation_subset(kwds['formation_names'])
            ax.plot(pts_subset.points[:][plane[0]], pts_subset.points[:][plane[1]], '.', color = color)            
        else:
            ax.plot(self.points[:][plane[0]], self.points[:][plane[1]], '.', color = color)
        
    def plot_3D(self, **kwds):
        """Create a plot of points in 3-D
        
        **Optional keywords**:
            - *ax* = matplotlib axis object: if provided, plot is attached to this axis
            - *formation_names* = list of formations : plot only points for specific formations
        """
        if kwds.has_key("ax"):
            # axis is provided, attach here
            ax = kwds['ax']
        else:
            fig = plt.figure(figsize = (10,8))
            ax = fig.add_subplot(111, projection='3d')
        if kwds.has_key("formation_names"):
            # create a subset with speficied formations, only
            pts_subset = self.create_formation_subset(kwds['formation_names'])
            pts_subset.plot_3D(ax = ax)
        else:
            # plot all    
            ax.scatter(self.points['x'], self.points['y'], self.points['z'])
            
    def save(self, filename):
        """Save points set to file
        
        **Arguments**:
            - *filename* = string : name of new file
        """
        f = open(filename, 'w')
        f.write(self.header)
        for point in self.points:
            f.write("%.2f,%.2f,%.3f,%s\n" % (point['x'], point['y'], point['z'], point['formation']))
        f.close()
        
class Struct3DFoliations(Struct3DPoints):
    """Class container for foliations (i.e. orientations) exported from GeoModeller
    
    Mainly based on Struct3DPoints as must required functionality 
    for location of elements - some functions overwritten, e.g. save and parse to read orientation data,
    as well!
    
    However, further methods might be added or adapted in the future, for example:
    - downsampling according to (eigen)vector methods, e.g. the work from the Monash guys, etc.
    - ploting of orientations in 2-D and 3-D
    """
    
    def __init__(self, **kwds):
        """Structural points data set
        
        **Optional keywords**:
            - *filename* = string : filename of csv file with exported points to load
        """
        # store point information in purpose defined numpy record
        self.ftype = np.dtype([('x', np.float32),
                          ('y', np.float32),
                          ('z', np.float32),
                          ('azimuth', np.float32),
                          ('dip', np.float32),
                          ('polarity', np.int),
                          ('formation', np.str_, 32)])

        if kwds.has_key("filename"):
            self.filename = kwds['filename']
            # read data
            self.parse()
            self.get_formation_names()
            self.get_range()

    
    def parse(self):
        """Parse filename and load data into numpy record
        
        The point information is stored in a purpose defined numpy record
        self.points
        """
        f = open(self.filename, "r")
        lines = f.readlines()
        self.header = lines[0]
        # determine position of elements in header (for extension to foliations, etc.)
        h_elem = np.array(self.header.rstrip().split(','))
        x_id = np.where(h_elem == 'X')[0]
        y_id = np.where(h_elem == 'Y')[0]
        z_id = np.where(h_elem == 'Z')[0]
        azi_id = np.where(h_elem == 'azimuth')[0]
        dip_id = np.where(h_elem == 'dip')[0]
        pol_id = np.where(h_elem == 'polarity')[0]
        form_id = np.where(h_elem == 'formation')[0]
        # print x_id
        # create numpy array for points
        self.len = (len(lines)-1)
        self.points = np.ndarray(self.len, dtype = self.ftype)
        for i,line in enumerate(lines[1:]):
            l = line.rstrip().split(',')
            self.points[i]['x'] = float(l[x_id]) 
            self.points[i]['y'] = float(l[y_id]) 
            self.points[i]['z'] = float(l[z_id]) 
            self.points[i]['azimuth'] = float(l[azi_id]) 
            self.points[i]['dip'] = float(l[dip_id]) 
            self.points[i]['polarity'] = float(l[pol_id]) 
            self.points[i]['formation'] = l[form_id]

    
    def save(self, filename):
        """Save points set to file
        
        **Arguments**:
            - *filename* = string : name of new file
        """
        f = open(filename, 'w')
        f.write(self.header)
        for point in self.points:
            f.write("%.2f,%.2f,%.3f,%.3f,%.3f,%d,%s\n" % (point['x'], point['y'], point['z'], 
                                             point['azimuth'], point['dip'], point['polarity'],
                                             point['formation']))
        f.close()
     
    
    
        


if __name__ == '__main__':
    pass

