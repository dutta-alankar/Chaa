/* compile_params.chpl — compile-time parameters for Chaa.
 *
 * These are `config param`s: edit this file, or override on the chpl
 * command line without touching the source, e.g.
 *     chpl ... -sNG=3 -shdf5Enabled=true
 * (CMake exposes the matching -DCHAA_* options.)
 */
module CompileParams {

  /* number of ghost-cell layers per active dimension: 2 suffices for
     constant/linear/limo3 reconstruction, 3 is needed for ppm (default) */
  config param NG = 3;

  /* compile in the HDF5 writer (links against -lhdf5) */
  config param hdf5Enabled = false;

  /* number of passive tracer (scalar colour) fields advected with the
     flow; they ride in the state vector after the five hydro slots */
  config param NSCAL = 1;
}
