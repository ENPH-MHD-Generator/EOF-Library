import gmsh
import sys
import math


# Note, I'm using Python 3.11.4 - F

# Usage:
#   python3 stl_to_elmer_msh.py MHD-Channel-x-axis.stl channel.msh
#
# Notes:
# - Works best when STL is a CLOSED, MANIFOLD surface (watertight).
# - Tags inlet/outlet by min/max X (assuming channel axis is X).
# - Everything else becomes "walls".

def main(stl_path: str, out_msh: str):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)

    gmsh.model.add("stl_model")

    # --- Import STL (discrete surface mesh)
    gmsh.merge(stl_path)

    # --- Classify the discrete triangles into surfaces, then create CAD geometry from them
    # Angle in radians. Increase if STL is noisy / faceted.
    angle = 40.0 * math.pi / 180.0
    force_parametrizable_patches = True
    include_boundary = True

    gmsh.model.mesh.classifySurfaces(
        angle,
        include_boundary,
        force_parametrizable_patches,
        angle
    )

    gmsh.model.mesh.createGeometry()

    # --- Collect all surface entities created from STL
    surfaces = gmsh.model.getEntities(2)
    if not surfaces:
        raise RuntimeError("No surfaces found after createGeometry(). Is STL empty/bad?")

    # --- Create a volume from all surfaces (watertight requirement!)
    surface_tags = [s[1] for s in surfaces]
    sl = gmsh.model.geo.addSurfaceLoop(surface_tags)
    vol = gmsh.model.geo.addVolume([sl])
    gmsh.model.geo.synchronize()

    # --- Determine bounding box of the entire model
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)

    # --- Identify inlet/outlet surfaces by X-min and X-max planes (tolerance-based)
    tol = 1e-6 * max(1.0, (xmax - xmin))
    inlet_surfs = gmsh.model.getEntitiesInBoundingBox(xmin - tol, ymin - 1e9, zmin - 1e9,
                                                      xmin + tol, ymax + 1e9, zmax + 1e9, 2)
    outlet_surfs = gmsh.model.getEntitiesInBoundingBox(xmax - tol, ymin - 1e9, zmin - 1e9,
                                                       xmax + tol, ymax + 1e9, zmax + 1e9, 2)

    inlet_tags = set([t for (d, t) in inlet_surfs])
    outlet_tags = set([t for (d, t) in outlet_surfs])

    # Walls = all other surfaces
    all_tags = set(surface_tags)
    wall_tags = list(all_tags - inlet_tags - outlet_tags)

    # --- Define Physical groups (these are what ElmerGrid and gmshToFoam rely on)
    # Choose stable IDs (or let gmsh assign). Names help OpenFOAM patch naming.
    phys_vol = gmsh.model.addPhysicalGroup(3, [vol], tag=1)
    gmsh.model.setPhysicalName(3, phys_vol, "fluid")

    if inlet_tags:
        pg = gmsh.model.addPhysicalGroup(2, list(inlet_tags), tag=101)
        gmsh.model.setPhysicalName(2, pg, "inlet")
    if outlet_tags:
        pg = gmsh.model.addPhysicalGroup(2, list(outlet_tags), tag=102)
        gmsh.model.setPhysicalName(2, pg, "outlet")
    if wall_tags:
        pg = gmsh.model.addPhysicalGroup(2, wall_tags, tag=103)
        gmsh.model.setPhysicalName(2, pg, "walls")

    # --- Meshing controls (tune these)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.02)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.1)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)  # 10 = HXT (if available), else choose 1/4/7...

    # --- Generate 3D mesh
    gmsh.model.mesh.generate(3)

    # --- Write msh v2.2 ASCII (ElmerGrid 14 and gmshToFoam are typically happiest with msh2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(out_msh)

    gmsh.finalize()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 stl_to_elmer_msh.py input.stl output.msh")
        sys.exit(2)
    main(sys.argv[1], sys.argv[2])
