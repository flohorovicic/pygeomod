
Creating an irregular mesh from Geomodeller for SHEMAT simulations
==================================================================

Regular meshes can be exported directly from within the Geomodeller GUI.
However, in many cases, a more flexible solution is required, for
example to:

-  update a mesh automatically (without using the GUI), or
-  create an irregular mesh with refined regions

These steps can easily be performed with a set of Python scripts and C
programs that access Geomodellers funcionality through the API.

The main funcionality required here is combined in the Python package
``pygeomod``. Two main packages are required: ``geogrid.py`` is the most
recent development and contains a (relatively general) class definition
for rectangular grids in general, with the link to Geomodeller in
particular. The package ´geomodeller\_xml\_obj.py\` contains methods to
access and modify information stored in the Geomodeller xml
Project-files. This functionality can be used, for example, to change
geological input parameters (e.g. dip of a fault) directly from the
Python script.

--------------

A note on installation:

The most tricky part is to get the API properly installed, all libraries
linked, and compiled on a system. On esim39, the required library path
settings are defined in

``adjust_to_jni.sh``

Another important point (for now, should be fixed at some stage...) is
that the shared object ``libgeomod.so`` has to be located in the current
directory... time to write a proper make file, but to date that's the
stage the project is in.

--------------

We will first start here with an example for the generation of an
rectilinear refined mesh for a simulation with SHEMAT.

.. code:: python

    # first step: import standard libraries and set pylab for plotting functionalities
    %pylab inline
    import numpy as np
    import matplotlib.pyplot as plt
    import sys, os

.. parsed-literal::

    
    Welcome to pylab, a matplotlib-based Python environment [backend: module://IPython.zmq.pylab.backend_inline].
    For more information, type 'help(pylab)'.


.. code:: python

    # Add path to pygeomod and import module (note: this is only required because it can't be installed properly at the moment)
    sys.path.append(r'/home/jni/git/tmp/pygeomod_tmp')
    import geogrid

Creating a regular grid
-----------------------

The geogrid module contains a variety of methods to generate grids. In
combinbation wtih Geomodeller, the easiest thing to do is to create a
regular mesh from a Geomodeller project:

.. code:: python

    # Define path to geomodeller model file:
    geomodel = r'/home/jni/git/tmp/geomuce/gemuce_tmp/examples/simple_three_layer/simple_three_layer.xml'
.. code:: python

    reload(geogrid) # only required for development stage - can be removed afterwards 
    # Now: define a GeoGrid object:
    G1 = geogrid.GeoGrid()
    # and set the boundaries/ model extent according to the Geomodeller model:
    G1.get_dimensions_from_geomodeller_xml_project(geomodel)
.. code:: python

    # and create a regular grid for a defined number of cells in each direction:
    nx = 25
    ny = 2
    nz = 25
    G1.define_regular_grid(nx, ny, nz)
.. code:: python

    # ...and, finally, update the grid properties on the base of the Geomodeller model:
    G1.update_from_geomodeller_project(geomodel)
The grid is stored in the object variable ``G1.grid`` as a numpy array.

.. code:: python

    type(G1.grid)



.. parsed-literal::

    numpy.ndarray



So the grid can directly be used to create slices, plots, further
caluclations, etc. However, a lot of functionality is alread implemented
in the geogrid package. For example, slice plots through the model can
simply be generated with:

.. code:: python

    G1.plot_section('y', colorbar=False, cmap='RdBu') # more plotting options possible, generally following the logic of matplotlibs imshow function



.. image:: Geomodeller-Export_files/Geomodeller-Export_11_0.png


It is also possible to export the model directly to VTK - however, this
requires an installation of the pyevtk package which is not installed on
esim for now:

.. code:: python

    G1.export_to_vtk()

::


    ---------------------------------------------------------------------------
    ImportError                               Traceback (most recent call last)

    <ipython-input-71-972ad06a1420> in <module>()
    ----> 1 G1.export_to_vtk()
    

    /home/jni/git/tmp/pygeomod_tmp/geogrid.py in export_to_vtk(self, vtk_filename, real_coords, **kwds)
        327         grid = kwds.get("grid", self.grid)
        328         var_name = kwds.get("var_name", "Geology")
    --> 329         from evtk.hl import gridToVTK
        330         # define coordinates
        331         x = np.zeros(self.nx + 1)


    ImportError: No module named evtk.hl


Rectilinear grids
-----------------

Creating a rectilinear grid requires only that the cell spacings are
explicitly defined. Everything else is exactly the same as before. Note
that it is (at the moment) your responsibility to assing proper spacings
- if you go beyond the bounds of the Geomodel, the function will not
crash, but return the standard Geomodeller "out" value (usually the
number of stratigraphic units + 1).

One way to create meshes in the correct range is, of course, to use the
extent of the Geomodel, determined with the function:

.. code:: python

    reload(geogrid) # only required for development stage - can be removed afterwards 
    # Now: define a GeoGrid object:
    G1 = geogrid.GeoGrid()
    # and set the boundaries/ model extent according to the Geomodeller model:
    G1.get_dimensions_from_geomodeller_xml_project(geomodel)
.. code:: python

    # The extent of the Geomodeller model can be obtained with:
    G1.xmin, G1.xmax



.. parsed-literal::

    (0, 1000)



.. code:: python

    # and the extent with:
    G1.extent_x



.. parsed-literal::

    1000



Let's be a bit fancy and create the horizontal (x,y) grid with a core
region of high refinement and increasing mesh sizes towards the
boundary. First, we define the geometry:

.. code:: python

    core_region = 100 # m
    # define cell width in core region:
    cell_width_core = 25 # m
    del_core = np.ones(int(core_region / cell_width_core)) * cell_width_core
    # and the number of cells in the boundary regions (the innermost cell has the size of the core cells):
    n_boundary = 10
    # now determine the boundary width on both sides of the core region:
    width_boundary_x = (G1.extent_x - core_region) / 2. 
    width_boundary_y = (G1.extent_y - core_region) / 2.
A little helper function in the ``geogrid`` package can be used to
determine an optimal cell increase factor for the boundary cells for a
given width an a number of cells, and a fixed inner cell width which we
take as the width of the core cells for a neat transition:

.. code:: python

    dx_boundary = geogrid.optimial_cell_increase(cell_width_core, n_boundary, width_boundary_x)
    dy_boundary = geogrid.optimial_cell_increase(cell_width_core, n_boundary, width_boundary_y)
We now simply combine the boundary and core cells for the complete
discretisation array:

.. code:: python

    delx = np.concatenate((dx_boundary[::-1], del_core, dx_boundary)) # first array reversed from large to small cells
    dely = np.concatenate((dy_boundary[::-1], del_core, dy_boundary))
A plot of the grid:

.. code:: python

    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(111)
    for dx in np.cumsum(delx):
        ax.axvline(dx, color = 'k')
    for dy in np.cumsum(dely):
        ax.axhline(dy, color = 'k')
    
    ax.set_xlim((0,sum(delx)))
    ax.set_ylim((0,sum(dely)))



.. parsed-literal::

    (0, 999.99999999999864)




.. image:: Geomodeller-Export_files/Geomodeller-Export_25_1.png


In z-direction we will create a regular mesh:

.. code:: python

    nz = 20
    delz = np.ones(nz) * G1.extent_z / nz
Ok, back to the geogrid package: we now assign the cell discretisation
arrays to the geogrid object and populate the grid with geology ids
determined from the Geomodeller model:

.. code:: python

    G1.define_irregular_grid(delx, dely, delz)
    G1.update_from_geomodeller_project(geomodel)
.. code:: python

    G1.grid



.. parsed-literal::

    array([[[ 1.,  1.,  1., ...,  1.,  1.,  1.],
            [ 1.,  1.,  1., ...,  1.,  1.,  1.],
            [ 1.,  1.,  1., ...,  1.,  1.,  1.],
            ..., 
            [ 1.,  1.,  1., ...,  1.,  1.,  1.],
            [ 1.,  1.,  1., ...,  1.,  1.,  1.],
            [ 1.,  1.,  1., ...,  1.,  1.,  1.]],
    
           [[ 1.,  1.,  1., ...,  1.,  1.,  1.],
            [ 1.,  1.,  1., ...,  1.,  1.,  1.],
            [ 1.,  1.,  1., ...,  1.,  1.,  1.],
            ..., 
            [ 1.,  1.,  1., ...,  1.,  1.,  1.],
            [ 1.,  1.,  1., ...,  1.,  1.,  1.],
            [ 1.,  1.,  1., ...,  1.,  1.,  1.]],
    
           [[ 1.,  1.,  1., ...,  1.,  1.,  1.],
            [ 1.,  1.,  1., ...,  1.,  1.,  1.],
            [ 1.,  1.,  1., ...,  1.,  1.,  1.],
            ..., 
            [ 1.,  1.,  1., ...,  1.,  1.,  1.],
            [ 1.,  1.,  1., ...,  1.,  1.,  1.],
            [ 1.,  1.,  1., ...,  1.,  1.,  1.]],
    
           ..., 
           [[ 1.,  1.,  1., ...,  3.,  3.,  3.],
            [ 1.,  1.,  1., ...,  3.,  3.,  3.],
            [ 1.,  1.,  1., ...,  3.,  3.,  3.],
            ..., 
            [ 1.,  1.,  1., ...,  3.,  3.,  3.],
            [ 1.,  1.,  1., ...,  3.,  3.,  3.],
            [ 1.,  1.,  1., ...,  3.,  3.,  3.]],
    
           [[ 1.,  1.,  1., ...,  3.,  3.,  3.],
            [ 1.,  1.,  1., ...,  3.,  3.,  3.],
            [ 1.,  1.,  1., ...,  3.,  3.,  3.],
            ..., 
            [ 1.,  1.,  1., ...,  3.,  3.,  3.],
            [ 1.,  1.,  1., ...,  3.,  3.,  3.],
            [ 1.,  1.,  1., ...,  3.,  3.,  3.]],
    
           [[ 1.,  1.,  1., ...,  3.,  3.,  3.],
            [ 1.,  1.,  1., ...,  3.,  3.,  3.],
            [ 1.,  1.,  1., ...,  3.,  3.,  3.],
            ..., 
            [ 1.,  1.,  1., ...,  3.,  3.,  3.],
            [ 1.,  1.,  1., ...,  3.,  3.,  3.],
            [ 1.,  1.,  1., ...,  3.,  3.,  3.]]])



The simple plotting functions don't work for irregular/ rectilinear
grids at to date (as imshow can only plot regular grids). Export to VTK
would work, in principle.

What we can do, however, is create a SHEMAT nml file (for the old SHEMAT
version) directly from the grid:

.. code:: python

    sys.path.append(r'/home/jni/git/tmp/PySHEMAT/PySHEMAT-master')
    import PySHEMAT
.. code:: python

    S1 = PySHEMAT.Shemat_file(from_geogrid = G1, nml_filename = 'updated_model.nml')

.. parsed-literal::

    create empty file

