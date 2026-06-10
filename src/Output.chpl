/* Output.chpl — simulation output:
 *   txt   column ASCII            (1D only)
 *   vtk   legacy ASCII VTK        (rectilinear for Cartesian, structured
 *                                  grid mapped to physical space otherwise)
 *   hdf5  HDF5 + XDMF (.xmf) companion for ParaView / VisIt (2D/3D)
 *
 * Fields are gathered to locale 0 and written serially; output volumes
 * in the test suite are small.
 */
module Output {
  use Params, Grid, State;
  use Hdf5IO;
  use IO, Math;

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
    // gather primitives once
    var locV: [1..nx1, 1..nx2, 1..nx3] StateVec;
    forall idx in DInt with (ref locV) do locV[idx] = V[idx];

    if doTxt && ndim == 1 then try! writeTxt(num, time, locV);
    if doVtk               then try! writeVtk(num, time, locV);
    if doH5 {
      dumpHdf5(outPath("h5", num), time, locV);
      if ndim >= 2 then try! writeXmf(num, time);
    }
    writeln("output ", num, " written at t = ", time);
  }

  /* ------------------------------- txt ------------------------------ */
  proc writeTxt(num: int, time: real, ref locV) throws {
    var f = open(outPath("txt", num), ioMode.cw);
    var w = f.writer(locking=false);
    w.writeln("# Chaa  problem=", problem, "  geometry=", geometry,
              "  time=", time);
    w.writeln("# x1  rho  vx1  vx2  vx3  prs");
    for i in 1..nx1 {
      const wv = locV[i, 1, 1];
      w.writef("%.12er %.12er %.12er %.12er %.12er %.12er\n",
               x1c(i), wv(IRHO), wv(IVX1), wv(IVX2), wv(IVX3), wv(IPRS));
    }
    w.close();
  }

  /* ------------------------------- vtk ------------------------------ */
  proc writeVtk(num: int, time: real, ref locV) throws {
    const n1 = nx1 + 1,
          n2 = if act2 then nx2 + 1 else 1,
          n3 = if act3 then nx3 + 1 else 1;

    var f = open(outPath("vtk", num), ioMode.cw);
    var w = f.writer(locking=false);
    w.writeln("# vtk DataFile Version 3.0");
    w.writeln("Chaa problem=", problem, " geometry=", geometry,
              " time=", time);
    w.writeln("ASCII");

    if geom == Geom.cartesian {
      w.writeln("DATASET RECTILINEAR_GRID");
      w.writeln("DIMENSIONS ", n1, " ", n2, " ", n3);
      w.writeln("X_COORDINATES ", n1, " double");
      for i in 1..n1 do w.writef("%.9er\n", x1f(i));
      w.writeln("Y_COORDINATES ", n2, " double");
      if act2 then for j in 1..n2 do w.writef("%.9er\n", x2f(j));
              else w.writeln("0.0");
      w.writeln("Z_COORDINATES ", n3, " double");
      if act3 then for k in 1..n3 do w.writef("%.9er\n", x3f(k));
              else w.writeln("0.0");
    } else {
      w.writeln("DATASET STRUCTURED_GRID");
      w.writeln("DIMENSIONS ", n1, " ", n2, " ", n3);
      w.writeln("POINTS ", n1*n2*n3, " double");
      for k in 1..n3 do
        for j in 1..n2 do
          for i in 1..n1 {
            const p = nodePos(i, j, k);
            w.writef("%.9er %.9er %.9er\n", p(0), p(1), p(2));
          }
    }

    w.writeln("CELL_DATA ", nx1*nx2*nx3);
    const names = ("rho", "vx1", "vx2", "vx3", "prs");
    for param c in 0..NVAR-1 {
      w.writeln("SCALARS ", names(c), " double 1");
      w.writeln("LOOKUP_TABLE default");
      for k in 1..nx3 do
        for j in 1..nx2 do
          for i in 1..nx1 do
            w.writef("%.9er\n", locV[i, j, k](c));
    }
    w.close();
  }

  /* ------------------------------- xmf ------------------------------ */
  proc writeXmf(num: int, time: real) throws {
    const h5name = problem + "." + pad4(num) + ".h5";
    var f = open(outPath("xmf", num), ioMode.cw);
    var w = f.writer(locking=false);

    const cellDims = if ndim == 2 then nx2:string + " " + nx1:string
                     else nx3:string + " " + nx2:string + " " + nx1:string;
    const nodeDims = if ndim == 2
                     then (nx2+1):string + " " + (nx1+1):string
                     else (nx3+1):string + " " + (nx2+1):string + " "
                          + (nx1+1):string;

    w.writeln("<?xml version=\"1.0\" ?>");
    w.writeln("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>");
    w.writeln("<Xdmf Version=\"2.0\">");
    w.writeln(" <Domain>");
    w.writeln("  <Grid Name=\"mesh\" GridType=\"Uniform\">");
    w.writeln("   <Time Value=\"", time, "\"/>");

    if geom == Geom.cartesian {
      const topo = if ndim == 2 then "2DRectMesh" else "3DRectMesh";
      const geo  = if ndim == 2 then "VXVY" else "VXVYVZ";
      w.writeln("   <Topology TopologyType=\"", topo,
                "\" Dimensions=\"", nodeDims, "\"/>");
      w.writeln("   <Geometry GeometryType=\"", geo, "\">");
      w.writeln("    <DataItem Dimensions=\"", nx1+1,
                "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">",
                h5name, ":/node_x1</DataItem>");
      w.writeln("    <DataItem Dimensions=\"", nx2+1,
                "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">",
                h5name, ":/node_x2</DataItem>");
      if ndim == 3 then
        w.writeln("    <DataItem Dimensions=\"", nx3+1,
                  "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">",
                  h5name, ":/node_x3</DataItem>");
      w.writeln("   </Geometry>");
    } else {
      const topo = if ndim == 2 then "2DSMesh" else "3DSMesh";
      const geo  = if ndim == 2 then "X_Y" else "X_Y_Z";
      w.writeln("   <Topology TopologyType=\"", topo,
                "\" Dimensions=\"", nodeDims, "\"/>");
      w.writeln("   <Geometry GeometryType=\"", geo, "\">");
      const comps = if ndim == 2 then ("nodes_x", "nodes_y", "")
                                 else ("nodes_x", "nodes_y", "nodes_z");
      for param c in 0..2 {
        if comps(c) != "" then
          w.writeln("    <DataItem Dimensions=\"", nodeDims,
                    "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">",
                    h5name, ":/", comps(c), "</DataItem>");
      }
      w.writeln("   </Geometry>");
    }

    const names = ("rho", "vx1", "vx2", "vx3", "prs");
    for param c in 0..NVAR-1 {
      w.writeln("   <Attribute Name=\"", names(c),
                "\" AttributeType=\"Scalar\" Center=\"Cell\">");
      w.writeln("    <DataItem Dimensions=\"", cellDims,
                "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">",
                h5name, ":/", names(c), "</DataItem>");
      w.writeln("   </Attribute>");
    }

    w.writeln("  </Grid>");
    w.writeln(" </Domain>");
    w.writeln("</Xdmf>");
    w.close();
  }
}
