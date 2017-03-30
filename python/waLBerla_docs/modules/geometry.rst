***************
Geometry module
***************


.. py:class:: TriangleMesh

   Corresponds to C++ class :doxylink:class:`walberla::geometry::TriangleMesh`

   .. py:attribute:: numTriangles
   
      Number of triangles.
   
   .. py:attribute:: numVertices
   
      Number of vertices.
   
   .. py:attribute:: numVertexNormals
   
      Number of vertex normals.
   
   .. py:method:: getAABB()
   
      Returns the axis aligned bounding box of the mesh.
      
   .. py:method:: volume()

      Volume of the Mesh.

   .. py:method:: scale(factor)

      Scales the complete mesh by the given factor.

   .. py:method:: scaleXYZ( factors )
   
      Scales the mesh by different factors in x,y,z direction.

      :param factors: tuple or list with 3 entries corresponding to x,y,z factors

   .. py:method:: exchangeAxes( xAxisId, yAxisId, zAxisId )

      Permutes the coordinate order of each vertex. e.g. ``m.exchangeAxes(0,2,1)`` exchanges z and y axis.


   .. py:method:: removeDuplicateVertices( tolerance = 1e-4)
      
      Merges vertices with a distance smaller than tolerance


   .. py:method:: merge( other, offset=(0,0,0 ) )
      
      Merges another mesh into the current mesh. During the merging all vertices of the
      other mesh are shifted by the given offset.


   .. py:method:: save( filename )
      
      Saves the mesh to a file. The mesh format is deduced using the filename extension.
      Supported formats are: obj,pov,off and vtp.

   .. py:staticmethod:: load( filename, broadcast=True )

      Loads mesh from a file. The mesh format is deduced using the filename extension.
      Supported formats are obj, pov and off.
   
      :param broadcast: If True the mesh is read on the root system only and broadcasted
                        to all other processes using MPI to reduce file system load.
