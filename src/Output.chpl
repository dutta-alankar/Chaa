/* Output.chpl — simulation output:
 *   txt   column ASCII            (1D only)
 *   vtk   legacy ASCII VTK        (rectilinear for Cartesian, structured
 *                                  grid mapped to physical space otherwise)
 *   hdf5  HDF5 + XDMF (.xmf) companion for ParaView / VisIt (2D/3D)
 *
 * Single locale: one file per dump.
 * Multiple locales: parallel piece output — every locale writes its own
 * block of the domain concurrently (independent files, safe with both
 * serial and parallel libhdf5), and locale 0 writes an XDMF *spatial
 * collection* that stitches the pieces back together for ParaView/VisIt.
 * VTK pieces are written per locale the same way.
 */
module Output {
  use Params, Grid, State;
  use Hdf5IO, Particles;
  use IO, Math;
  import CompileParams;
  import Gpu;

  const fmtList = outFormats.split(",");

  proc hasFmt(s: string): bool {
    for f in fmtList do
      if f.strip() == s then return true;
    return false;
  }

  const doTxt = hasFmt("txt"),
        doVtk = hasFmt("vtk"),
        doH5  = hasFmt("hdf5");

  proc pad4(n: int): string {
    var s = n:string;
    while s.size < 4 do s = "0" + s;
    return s;
  }

  proc outPath(ext: string, num: int): string do
    return outDir + "/" + problem + "." + pad4(num) + "." + ext;

  proc writeOutputs(num: int, time: real) {
    // GPU build: the devices own the interior — refresh the host copy
    if CompileParams.gpuEnabled then Gpu.gpuDownV();
    if numLocales == 1 {
      if doTxt && ndim == 1 then try! writeTxt(num, time);
      if doVtk then try! writeVtk(num, time, 1..nx1, 1..nx2, 1..nx3,
                                  outPath("vtk", num));
      if doH5 {
        dumpHdf5(outPath("h5", num), time, 1..nx1, 1..nx2, 1..nx3);
        if ndim >= 2 then try! writeXmf(num, time);
      }
    } else {
      /* parallel piece output: one writer per locale, concurrently */
      coforall loc in Locales do on loc {
        const mine = DInt.localSubdomain();
        if mine.size > 0 {
          if doH5 then
            dumpHdf5(piecePath("h5", num, here.id), time,
                     mine.dim(0), mine.dim(1), mine.dim(2));
          if doVtk then
            try! writeVtk(num, time, mine.dim(0), mine.dim(1),
                          mine.dim(2), piecePath("vtk", num, here.id));
        }
      }
      if doH5 && ndim >= 2 then try! writeXmfCollection(num, time);
      if doTxt && ndim == 1 then try! writeTxt(num, time);
    }
    if nParticles > 0 then
      try! writeParticles(outDir + "/" + problem + ".particles." +
                          pad4(num) + ".txt");
    writeln("output ", num, " written at t = ", time);
  }

  proc piecePath(ext: string, num: int, lid: int): string do
    return outDir + "/" + problem + "." + pad4(num) + ".loc" + lid:string +
           "." + ext;

  /* ------------------------------- txt ------------------------------ */
  proc writeTxt(num: int, time: real) throws {
    var f = open(outPath("txt", num), ioMode.cw);
    var w = f.writer(locking=false);
    w.writeln("# Chaa  problem=", problem, "  geometry=", geometry,
              "  time=", time);
    w.write("# x1  rho  vx1  vx2  vx3  prs");
    for c in ISC..NTOT-1 do w.write("  ", fieldName(c));
    w.writeln();
    for i in 1..nx1 {
      const wv = V[i, 1, 1];
      w.writef("%.12er", x1c(i));
      for param c in 0..NTOT-1 do w.writef(" %.12er", wv(c));
      w.writeln();
    }
    w.close();
  }

  /* ------------------------------- vtk ------------------------------ */
  proc writeVtk(num: int, time: real, r1: range, r2: range, r3: range,
                path: string) throws {
    const n1 = r1.size + 1,
          n2 = if act2 then r2.size + 1 else 1,
          n3 = if act3 then r3.size + 1 else 1;

    var f = open(path, ioMode.cw);
    var w = f.writer(locking=false);
    w.writeln("# vtk DataFile Version 3.0");
    w.writeln("Chaa problem=", problem, " geometry=", geometry,
              " time=", time);
    w.writeln("ASCII");

    if geom == Geom.cartesian {
      w.writeln("DATASET RECTILINEAR_GRID");
      w.writeln("DIMENSIONS ", n1, " ", n2, " ", n3);
      w.writeln("X_COORDINATES ", n1, " double");
      for i in r1.low..r1.high+1 do w.writef("%.9er\n", x1f(i));
      w.writeln("Y_COORDINATES ", n2, " double");
      if act2 then for j in r2.low..r2.high+1 do w.writef("%.9er\n", x2f(j));
              else w.writeln("0.0");
      w.writeln("Z_COORDINATES ", n3, " double");
      if act3 then for k in r3.low..r3.high+1 do w.writef("%.9er\n", x3f(k));
              else w.writeln("0.0");
    } else {
      w.writeln("DATASET STRUCTURED_GRID");
      w.writeln("DIMENSIONS ", n1, " ", n2, " ", n3);
      w.writeln("POINTS ", n1*n2*n3, " double");
      for k in r3.low..(if act3 then r3.high+1 else r3.high) do
        for j in r2.low..(if act2 then r2.high+1 else r2.high) do
          for i in r1.low..r1.high+1 {
            const p = nodePos(i, j, k);
            w.writef("%.9er %.9er %.9er\n", p(0), p(1), p(2));
          }
    }

    w.writeln("CELL_DATA ", r1.size*r2.size*r3.size);
    for param c in 0..NTOT-1 {
      w.writeln("SCALARS ", fieldName(c), " double 1");
      w.writeln("LOOKUP_TABLE default");
      for k in r3 do
        for j in r2 do
          for i in r1 do
            w.writef("%.9er\n", V[i, j, k](c));
    }
    w.close();
  }

  /* ------------------------------- xmf ------------------------------ */

  /* emit one <Grid> block describing fields over (r1,r2,r3) stored in
     h5name (which already contains matching coordinate datasets) */
  proc xmfGrid(w, gname: string, h5name: string,
               r1: range, r2: range, r3: range) throws {
    const c1 = r1.size, c2 = r2.size, c3 = r3.size;
    const cellDims = if ndim == 2 then c2:string + " " + c1:string
                     else c3:string + " " + c2:string + " " + c1:string;
    const nodeDims = if ndim == 2
                     then (c2+1):string + " " + (c1+1):string
                     else (c3+1):string + " " + (c2+1):string + " "
                          + (c1+1):string;

    w.writeln("  <Grid Name=\"", gname, "\" GridType=\"Uniform\">");
    if geom == Geom.cartesian {
      const topo = if ndim == 2 then "2DRectMesh" else "3DRectMesh";
      const geo  = if ndim == 2 then "VXVY" else "VXVYVZ";
      w.writeln("   <Topology TopologyType=\"", topo,
                "\" Dimensions=\"", nodeDims, "\"/>");
      w.writeln("   <Geometry GeometryType=\"", geo, "\">");
      w.writeln("    <DataItem Dimensions=\"", c1+1,
                "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">",
                h5name, ":/node_x1</DataItem>");
      w.writeln("    <DataItem Dimensions=\"", c2+1,
                "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">",
                h5name, ":/node_x2</DataItem>");
      if ndim == 3 then
        w.writeln("    <DataItem Dimensions=\"", c3+1,
                  "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">",
                  h5name, ":/node_x3</DataItem>");
      w.writeln("   </Geometry>");
    } else {
      const topo = if ndim == 2 then "2DSMesh" else "3DSMesh";
      const geo  = if ndim == 2 then "X_Y" else "X_Y_Z";
      w.writeln("   <Topology TopologyType=\"", topo,
                "\" Dimensions=\"", nodeDims, "\"/>");
      w.writeln("   <Geometry GeometryType=\"", geo, "\">");
      const ncomp = if ndim == 2 then 2 else 3;
      const comps = ("nodes_x", "nodes_y", "nodes_z");
      for param c in 0..2 {
        if c < ncomp then
          w.writeln("    <DataItem Dimensions=\"", nodeDims,
                    "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">",
                    h5name, ":/", comps(c), "</DataItem>");
      }
      w.writeln("   </Geometry>");
    }

    for param c in 0..NTOT-1 {
      w.writeln("   <Attribute Name=\"", fieldName(c),
                "\" AttributeType=\"Scalar\" Center=\"Cell\">");
      w.writeln("    <DataItem Dimensions=\"", cellDims,
                "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">",
                h5name, ":/", fieldName(c), "</DataItem>");
      w.writeln("   </Attribute>");
    }
    w.writeln("  </Grid>");
  }

  proc writeXmf(num: int, time: real) throws {
    const h5name = problem + "." + pad4(num) + ".h5";
    var f = open(outPath("xmf", num), ioMode.cw);
    var w = f.writer(locking=false);
    w.writeln("<?xml version=\"1.0\" ?>");
    w.writeln("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>");
    w.writeln("<Xdmf Version=\"2.0\">");
    w.writeln(" <Domain>");
    w.writeln("  <Grid Name=\"chaa\" GridType=\"Collection\"",
              " CollectionType=\"Spatial\">");
    w.writeln("   <Time Value=\"", time, "\"/>");
    xmfGrid(w, "mesh", h5name, 1..nx1, 1..nx2, 1..nx3);
    w.writeln("  </Grid>");
    w.writeln(" </Domain>");
    w.writeln("</Xdmf>");
    w.close();
  }

  /* master file stitching the per-locale pieces */
  proc writeXmfCollection(num: int, time: real) throws {
    var f = open(outPath("xmf", num), ioMode.cw);
    var w = f.writer(locking=false);
    w.writeln("<?xml version=\"1.0\" ?>");
    w.writeln("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>");
    w.writeln("<Xdmf Version=\"2.0\">");
    w.writeln(" <Domain>");
    w.writeln("  <Grid Name=\"chaa\" GridType=\"Collection\"",
              " CollectionType=\"Spatial\">");
    w.writeln("   <Time Value=\"", time, "\"/>");
    for loc in Locales {
      const sub = DInt.localSubdomain(loc);
      if sub.size > 0 {
        const h5name = problem + "." + pad4(num) + ".loc" +
                       loc.id:string + ".h5";
        xmfGrid(w, "piece" + loc.id:string, h5name,
                sub.dim(0), sub.dim(1), sub.dim(2));
      }
    }
    w.writeln("  </Grid>");
    w.writeln(" </Domain>");
    w.writeln("</Xdmf>");
    w.close();
  }
}
