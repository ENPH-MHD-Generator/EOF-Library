import gmsh
import sys

def main(out_msh: str):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("elmer_em_3d_quasi2d")

    # Geometry
    L = 1.0    # x
    H = 0.1    # y
    W = 0.1    # z

    lc_yz = 0.01
    Nx = 1

    # ---- YZ rectangle at x=0 ----
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc_yz)
    p2 = gmsh.model.geo.addPoint(0, H, 0, lc_yz)
    p3 = gmsh.model.geo.addPoint(0, H, W, lc_yz)
    p4 = gmsh.model.geo.addPoint(0, 0, W, lc_yz)

    l1 = gmsh.model.geo.addLine(p1, p2)  # y
    l2 = gmsh.model.geo.addLine(p2, p3)  # z
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    s0 = gmsh.model.geo.addPlaneSurface([cl])

    gmsh.model.geo.synchronize()

    # ---- Extrude along X ----
    ext = gmsh.model.geo.extrude(
        [(2, s0)], L, 0, 0,
        numElements=[Nx],
        recombine=False
    )
    gmsh.model.geo.synchronize()

    # Extract volume
    vol = [tag for dim, tag in ext if dim == 3][0]

    # Extract surfaces by location
    surfaces = gmsh.model.getEntities(2)

    inlet  = []
    outlet = []
    y0     = []
    yH     = []
    z0     = []
    zW     = []

    for dim, s in surfaces:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, s)
        x = 0.5 * (xmin + xmax)
        y = 0.5 * (ymin + ymax)
        z = 0.5 * (zmin + zmax)

        if abs(x - 0.0) < 1e-8:
            inlet.append(s)
        elif abs(x - L) < 1e-8:
            outlet.append(s)
        elif abs(y - 0.0) < 1e-8:
            y0.append(s)
        elif abs(y - H) < 1e-8:
            yH.append(s)
        elif abs(z - 0.0) < 1e-8:
            z0.append(s)
        elif abs(z - W) < 1e-8:
            zW.append(s)


    # ---- Physical groups (DETERMINISTIC IDS) ----

    # Volume
    gmsh.model.addPhysicalGroup(3, [vol], tag=1)
    gmsh.model.setPhysicalName(3, 1, "fluid")

    # Faraday electrodes
    gmsh.model.addPhysicalGroup(2, y0, tag=10)
    gmsh.model.setPhysicalName(2, 10, "FaradayMinusY")

    gmsh.model.addPhysicalGroup(2, yH, tag=11)
    gmsh.model.setPhysicalName(2, 11, "FaradayPlusY")

    # Inlet / outlet (Hall electrodes later if desired)
    gmsh.model.addPhysicalGroup(2, inlet, tag=20)
    gmsh.model.setPhysicalName(2, 20, "InletX")

    gmsh.model.addPhysicalGroup(2, outlet, tag=21)
    gmsh.model.setPhysicalName(2, 21, "OutletX")

    # Insulating side walls
    gmsh.model.addPhysicalGroup(2, z0 + zW, tag=30)
    gmsh.model.setPhysicalName(2, 30, "SideWallsZ")

    # ---- Mesh ----
    gmsh.model.mesh.generate(3)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(out_msh)

    gmsh.finalize()


if __name__ == "__main__":
    main(sys.argv[1])
