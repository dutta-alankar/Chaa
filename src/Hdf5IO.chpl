/* Hdf5IO.chpl — minimal HDF5 writer through direct C bindings.
 *
 * Only a handful of H5 core functions are needed, so they are declared
 * as extern procs (no header required) and the code links with -lhdf5.
 * Compile with  -shdf5Enabled=true  to activate; when disabled the
 * param-folded branches keep the binary free of any HDF5 dependency.
 */
module Hdf5IO {
  use Params, Grid;
  use CTypes;
  import CompileParams.hdf5Enabled;   // compile_params.chpl

  type hid = int(64);

  // prototypes live in src/chaa_h5.h (passed on the chpl command line)
  extern proc H5open(): c_int;
  extern proc H5Fcreate(name: c_ptrConst(c_char), flags: c_uint,
                        fcpl: hid, fapl: hid): hid;
  extern proc H5Fclose(f: hid): c_int;
  extern proc H5Screate_simple(rank: c_int, dims: c_ptrConst(uint(64)),
                               maxdims: c_ptrConst(uint(64))): hid;
  extern proc H5Sclose(s: hid): c_int;
  extern proc H5Dcreate2(loc: hid, name: c_ptrConst(c_char), dtype: hid,
                         space: hid, lcpl: hid, dcpl: hid, dapl: hid): hid;
  extern proc H5Dwrite(dset: hid, memtype: hid, memspace: hid,
                       filespace: hid, xfer: hid,
                       buf: c_ptrConst(void)): c_int;
  extern proc H5Dclose(d: hid): c_int;
  extern var H5T_NATIVE_DOUBLE_g: hid;

  param H5F_ACC_TRUNC: c_uint = 2;
  param H5P_DEFAULT: hid = 0;
  param H5S_ALL: hid = 0;

  proc h5PutArray(fid: hid, name: string, rank: int,
                  d0: int, d1: int, d2: int, ref buf: [] real) {
    if hdf5Enabled {
      var dims: [0..2] uint(64);
      dims[0] = d0: uint(64);
      dims[1] = d1: uint(64);
      dims[2] = d2: uint(64);
      var nilDims: c_ptrConst(uint(64));
      const sid = H5Screate_simple(rank: c_int,
                                   c_ptrTo(dims): c_ptrConst(uint(64)),
                                   nilDims);
      const did = H5Dcreate2(fid, name.c_str(), H5T_NATIVE_DOUBLE_g, sid,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if did < 0 then halt("HDF5: failed to create dataset " + name);
      H5Dwrite(did, H5T_NATIVE_DOUBLE_g, H5S_ALL, H5S_ALL, H5P_DEFAULT,
               c_ptrTo(buf): c_ptrConst(void));
      H5Dclose(did);
      H5Sclose(sid);
    }
  }

  /* write one dump: fields (C order, x1 fastest), cell-centre and node
     coordinates, mapped node positions for curvilinear meshes, time */
  proc dumpHdf5(path: string, time: real, ref locV) {
    if !hdf5Enabled {
      halt("hdf5 output requested but the binary was compiled without " +
           "HDF5 support; rebuild with HDF5=1 (-shdf5Enabled=true -lhdf5)");
    } else {
      H5open();
      const fid = H5Fcreate(path.c_str(), H5F_ACC_TRUNC,
                            H5P_DEFAULT, H5P_DEFAULT);
      if fid < 0 then halt("HDF5: cannot create " + path);

      // ---- fields ----
      const n = nx1*nx2*nx3;
      var buf: [0..#n] real;
      const names = ("rho", "vx1", "vx2", "vx3", "prs");
      for param c in 0..NVAR-1 {
        forall (i, j, k) in {1..nx1, 1..nx2, 1..nx3} with (ref buf) do
          buf[((k-1)*nx2 + (j-1))*nx1 + (i-1)] = locV[i, j, k](c);
        select ndim {
          when 1 do h5PutArray(fid, names(c), 1, nx1, 1, 1, buf);
          when 2 do h5PutArray(fid, names(c), 2, nx2, nx1, 1, buf);
          otherwise do h5PutArray(fid, names(c), 3, nx3, nx2, nx1, buf);
        }
      }

      // ---- cell-centre coordinates (1D per axis) ----
      {
        var c1: [0..#nx1] real = [i in 0..#nx1] x1c(i+1);
        h5PutArray(fid, "cc_x1", 1, nx1, 1, 1, c1);
        var c2: [0..#nx2] real = [j in 0..#nx2] x2c(j+1);
        h5PutArray(fid, "cc_x2", 1, nx2, 1, 1, c2);
        var c3: [0..#nx3] real = [kk in 0..#nx3] x3c(kk+1);
        h5PutArray(fid, "cc_x3", 1, nx3, 1, 1, c3);
      }

      // ---- node coordinates ----
      {
        var f1: [0..#(nx1+1)] real = [i in 0..#(nx1+1)] x1f(i+1);
        h5PutArray(fid, "node_x1", 1, nx1+1, 1, 1, f1);
        var f2: [0..#(nx2+1)] real = [j in 0..#(nx2+1)] x2f(j+1);
        h5PutArray(fid, "node_x2", 1, nx2+1, 1, 1, f2);
        var f3: [0..#(nx3+1)] real = [kk in 0..#(nx3+1)] x3f(kk+1);
        h5PutArray(fid, "node_x3", 1, nx3+1, 1, 1, f3);
      }

      // ---- mapped node positions for curvilinear meshes (2D/3D) ----
      if geom != Geom.cartesian && ndim >= 2 {
        const m1 = nx1 + 1,
              m2 = nx2 + 1,
              m3 = if act3 then nx3 + 1 else 1;
        var nbx, nby, nbz: [0..#(m1*m2*m3)] real;
        forall (i, j, k) in {1..m1, 1..m2, 1..m3}
            with (ref nbx, ref nby, ref nbz) {
          const p = nodePos(i, j, k);
          const idx = ((k-1)*m2 + (j-1))*m1 + (i-1);
          nbx[idx] = p(0); nby[idx] = p(1); nbz[idx] = p(2);
        }
        if ndim == 2 {
          h5PutArray(fid, "nodes_x", 2, m2, m1, 1, nbx);
          h5PutArray(fid, "nodes_y", 2, m2, m1, 1, nby);
        } else {
          h5PutArray(fid, "nodes_x", 3, m3, m2, m1, nbx);
          h5PutArray(fid, "nodes_y", 3, m3, m2, m1, nby);
          h5PutArray(fid, "nodes_z", 3, m3, m2, m1, nbz);
        }
      }

      // ---- time stamp ----
      {
        var tbuf: [0..0] real = time;
        h5PutArray(fid, "time", 1, 1, 1, 1, tbuf);
      }

      H5Fclose(fid);
    }
  }
}
