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

  /* write one dump over the cell ranges (r1,r2,r3) — the full interior
     on a single locale, or this locale's block in parallel piece
     output.  Layout: fields in C order with x1 fastest, cell-centre and
     node coordinates, mapped node positions for curvilinear meshes,
     time stamp. */
  proc dumpHdf5(path: string, time: real,
                r1: range, r2: range, r3: range) {
    if !hdf5Enabled {
      halt("hdf5 output requested but the binary was compiled without " +
           "HDF5 support; rebuild with -DCHAA_HDF5=ON");
    } else {
      use State;
      H5open();
      const fid = H5Fcreate(path.c_str(), H5F_ACC_TRUNC,
                            H5P_DEFAULT, H5P_DEFAULT);
      if fid < 0 then halt("HDF5: cannot create " + path);

      const c1 = r1.size, c2 = r2.size, c3 = r3.size;

      // ---- fields ----
      var buf: [0..#(c1*c2*c3)] real;
      for param c in 0..NTOT-1 {
        forall (i, j, k) in {r1, r2, r3} with (ref buf) do
          buf[((k-r3.low)*c2 + (j-r2.low))*c1 + (i-r1.low)] = V[i, j, k](c);
        select ndim {
          when 1 do h5PutArray(fid, fieldName(c), 1, c1, 1, 1, buf);
          when 2 do h5PutArray(fid, fieldName(c), 2, c2, c1, 1, buf);
          otherwise do h5PutArray(fid, fieldName(c), 3, c3, c2, c1, buf);
        }
      }

      // ---- cell-centre coordinates (1D per axis) ----
      {
        var cc1: [0..#c1] real = [i in 0..#c1] x1c(r1.low + i);
        h5PutArray(fid, "cc_x1", 1, c1, 1, 1, cc1);
        var cc2: [0..#c2] real = [j in 0..#c2] x2c(r2.low + j);
        h5PutArray(fid, "cc_x2", 1, c2, 1, 1, cc2);
        var cc3: [0..#c3] real = [kk in 0..#c3] x3c(r3.low + kk);
        h5PutArray(fid, "cc_x3", 1, c3, 1, 1, cc3);
      }

      // ---- node coordinates ----
      {
        var f1: [0..#(c1+1)] real = [i in 0..#(c1+1)] x1f(r1.low + i);
        h5PutArray(fid, "node_x1", 1, c1+1, 1, 1, f1);
        var f2: [0..#(c2+1)] real = [j in 0..#(c2+1)] x2f(r2.low + j);
        h5PutArray(fid, "node_x2", 1, c2+1, 1, 1, f2);
        var f3: [0..#(c3+1)] real = [kk in 0..#(c3+1)] x3f(r3.low + kk);
        h5PutArray(fid, "node_x3", 1, c3+1, 1, 1, f3);
      }

      // ---- mapped node positions for curvilinear meshes (2D/3D) ----
      if geom != Geom.cartesian && ndim >= 2 {
        const m1 = c1 + 1,
              m2 = c2 + 1,
              m3 = if act3 then c3 + 1 else 1;
        var nbx, nby, nbz: [0..#(m1*m2*m3)] real;
        forall (i, j, k) in {1..m1, 1..m2, 1..m3}
            with (ref nbx, ref nby, ref nbz) {
          const p = nodePos(r1.low + i - 1, r2.low + j - 1, r3.low + k - 1);
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
