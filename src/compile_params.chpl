/* compile_params.chpl — compile-time parameters for Chaa.
 *
 * These are `config param`s: edit this file, or override on the chpl
 * command line without touching the source, e.g.
 *     chpl ... -sNG=3 -shdf5Enabled=true
 * (CMake exposes the matching -DCHAA_* options.)
 */
module CompileParams {

  /* number of ghost-cell layers per active dimension (2 is required by
     the piecewise-linear reconstruction; raise for wider stencils) */
  config param NG = 2;

  /* compile in the HDF5 writer (links against -lhdf5) */
  config param hdf5Enabled = false;
}
