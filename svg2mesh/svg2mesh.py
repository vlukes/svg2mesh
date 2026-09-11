import sys
import os
import argparse
from ast import literal_eval
import numpy as nm
import gmsh
import meshio
from svgelements import SVG


def set_periodic(model, dx, dy, dz=0, dim=2):
    """
    Set mesh periodic.
    
    Inputs
    ======
    model: gmsh.model
        GMSH model
    dx, dy, dz: floats
        Domain dimensions
    """
    def set_periodic(master, slave, dx, dy, dz):
        for m, s in zip(master, slave):
            gmsh.model.mesh.setPeriodic(
                dim - 1, [s], [m],
                [1, 0, 0, dx,
                0, 1, 0, dy,
                0, 0, 1, dz,
                0, 0, 0, 1]
            )

    tol = 1e-6
    bnd = nm.array([(tag,) + tuple(gmsh.model.occ.getCenterOfMass(d, tag))
                    for d, tag in model.getEntities(dim=dim-1)])

    for k in range(dim):
        cmin, cmax = nm.min(bnd[:, k + 1]), nm.max(bnd[:, k + 1])
        midxs = nm.where(nm.abs(bnd[:, k + 1] - cmin) < tol)[0]
        sidxs = nm.where(nm.abs(bnd[:, k + 1] - cmax) < tol)[0]

        dv = nm.array([0] * 3)
        dd = cmax - cmin
        dv[k] = dd

        pairs = []
        for m in midxs:
            dist = nm.linalg.norm(bnd[sidxs, 1:] - dv - bnd[m, 1:], axis=1)
            idx = nm.where(dist < (dd * tol))[0]
            pairs.append((m, sidxs[idx[0]]))
        
        pairs = nm.array(pairs)
        args = tuple(bnd[pairs, 0].T) + tuple(dv)
        set_periodic(*args)


def set_mesh_size(model, pgroups, msize, dim=2):
    """
    Set mesh size individually to each physical group.
    
    Inputs
    ======
    model: gmsh.model
        GMSH model
    pgroups: list
        Physical group IDs
    msize: list of floats
        Mesh size for each physical group
    """
    mesh = model.mesh
    fields = []
    for pg, ms in zip(pgroups, msize):
        entities = model.getEntitiesForPhysicalGroup(dim, pg)
        field = mesh.field.add("Constant")
        mesh.field.setNumbers(field, "SurfacesList", entities)
        mesh.field.setNumber(field, "VIn", ms)
        fields.append(field)

    field = mesh.field.add("Min")
    mesh.field.setNumbers(field, "FieldsList", fields)
    mesh.field.setAsBackgroundMesh(field)


def gmsh_init():
    "Initiate GMSH module."
    gmsh.initialize()
    gmsh.logger.start()
    gmsh.model.add('boolean')

    return gmsh, gmsh.model, gmsh.model.occ


def get_xy(p, dec=6):
    return (float(nm.round(p.x, dec)), float(nm.round(p.y, dec)))


def quadratic_to_cubic(p0, c1, p1):
    "Convert SVG quadratic bezier into equivalent cubic bezier."

    c1 = (p0[0] + 2.0/3.0*(c1[0] - p0[0]),
          p0[1] + 2.0/3.0*(c1[1] - p0[1]))
    c2 = (p1[0] + 2.0/3.0*(c1[0] - p1[0]),
          p1[1] + 2.0/3.0*(c1[1] - p1[1]))

    return c1, c2


def rotate_obj(obj, angle, occ, h, cx=0, cy=0):
    "Rotate OCC object around Z axis."
    if nm.abs(angle) > 1e-3:
        occ.rotate([obj], cx, h - cy, 0, 0, 0, 1, angle)


def get_group_elements(elements, occ, cargs, merge=False, oflag=''):
    out = []
    for el in elements:
        if isinstance(el, list):
            out += get_group_elements(el, occ, cargs, True, '    ')
        else:
            shape = el.__class__.__name__
            print(f'{oflag}{shape}: {el}')
            if shape == 'Rect':
                y = cargs[4] - el.y - el.height
                out.append((2, occ.addRectangle(el.x, y, 0,
                                                el.width, el.height)))
                # rotate_obj(out[-1], -el.rotation, occ, cargs[4])
                #            el.transform[4], el.transform[5])

            elif shape == 'Ellipse':
                out.append((2, occ.addDisk(el.cx, cargs[4] - el.cy, 0,
                                           el.rx, el.ry)))
                rotate_obj(out[-1], -el.rotation, occ, cargs[4])                           

            elif shape == 'Path':
                point_keys = {get_xy(p) for p in el.as_points()}
                points = {p: occ.addPoint(p[0], cargs[4] - p[1], 0)
                          for p in point_keys}
                
                loops, lines = [], []
                for seg in el:
                    seg_name = seg.__class__.__name__
                    if seg_name == 'Move':
                        if len(lines) > 0:
                            if  el.id.startswith('text'):
                                loops.append(occ.addCurveLoop(lines))
                                out.append((2, occ.addPlaneSurface(loops)))
                            else:
                                loops.append(occ.addCurveLoop(lines))
                            lines = []

                    elif seg_name in ['Line', 'Close']:
                        p1, p2 = get_xy(seg.start), get_xy(seg.end)
                        if p1 != p2:
                            lines.append(occ.addLine(points[p1], points[p2]))

                    elif seg_name == 'QuadraticBezier':
                        p1, c12, p2 = (get_xy(seg.start), get_xy(seg.control),
                                       get_xy(seg.end))
                        c1, c2 = quadratic_to_cubic(p1, c12, p2)
                        points[c1] = occ.addPoint(c1[0], cargs[4] - c1[1], 0)
                        points[c2] = occ.addPoint(c2[0], cargs[4] - c2[1], 0)
                        lines.append(occ.addBezier([points[p1], points[c1],
                                                    points[c2], points[p2]]))
                    
                    elif seg_name == 'CubicBezier':
                        p1, p2 = get_xy(seg.start), get_xy(seg.end)
                        c1, c2 = get_xy(seg.control1), get_xy(seg.control2)
                        lines.append(occ.addBezier([points[p1], points[c1],
                                                    points[c2], points[p2]]))
                    
                    elif seg_name == 'Arc':
                        apoints = [get_xy(seg.point(t))
                                   for t in nm.linspace(0, 1, 30)]
                        spoints = []
                        for apt in apoints:
                            if apt not in points:
                                y = cargs[4] - apt[1]
                                points[apt] = occ.addPoint(apt[0], y, 0)
                            spoints.append(points[apt])
                        
                        lines.append(occ.addBSpline(spoints))

                    else:
                        raise ValueError(f'Unknown Path element "{seg_name}"!')

                loops.append(occ.addCurveLoop(lines))
                out.append((2, occ.addPlaneSurface(loops)))

            else:
                raise ValueError(f'Unknown SVG element "{shape}"!')


    # merge grouped objects
    if merge and len(out) > 1:
        out, _ = occ.fuse(out[:1], out[1:])

    # crop to view box
    if cargs and len(out) > 0:
        cbox1 = (2, occ.addRectangle(*cargs))
        cbox2 = (2, occ.addRectangle(-cargs[3], -cargs[4], cargs[2],
                                     3*cargs[3], 3*cargs[4]))
        cbox, _ = occ.cut([cbox2], [cbox1])
        out, _ = occ.cut(out, cbox)

    return out


def save_png(filename_out, filename_png):
    """
    Save mesh screenshot as PNG file.

    Inputs
    ======
    filename_out: str
        Output filename
    filename_png: str
        PNG output filename
    """
    import pyvista as pv

    mesh = meshio.read(filename_out)

    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(mesh, show_edges=True, scalars='mat_id', show_scalar_bar=False)

    plotter.view_xy()
    plotter.camera.zoom(1.4)
    plotter.image_scale = 2
    plotter.show(screenshot=filename_png)


def sum_dict(d):
    return sum(d, [])


def gen_mesh_from_svg(filename_svg, filename_out=None,
                      mesh_size=None, unit_cell=False,
                      periodic=False, export_png=False):
    """
    Generate FE mesh from geometry defined by a SVG file.

    Inputs
    ======
    filename_svg: str
        SVG filename
    filename_out: str
        Output filename
    mesh_size: float
        Size of mesh elements
    unit_cell: bool
        If True, generate periodic unit cell
    periodic: bool
        If True, generate periodic mesh
    export_png: bool
        If True, save mesh screenshot (PNG)
    """
    gmsh, model, occ = gmsh_init()

    # set ppi=25.4 to keep units in mm
    svg = SVG.parse(filename_svg, ppi=25.4)

    vbox = svg.viewbox
    cargs = (vbox.x, vbox.y, 0, vbox.width, vbox.height)

    layers = []
    for el in svg:
        print(f'group {el.id}:')
        layers.append(get_group_elements(el, occ, cargs, oflag='  '))

    if len(layers) > 1:
        _, ovv = occ.fragment(layers[0], sum_dict(layers[1:]))
    else:
        ovv = layers

    model.occ.synchronize()

    mids = sum_dict([[ik] * len(k) for ik, k in enumerate(layers)])
    tags_1 = set()
    pgs = [[] for _ in range(len(layers))]
    for mid in nm.arange(len(ovv))[::-1]:
        ids = [k for _, k in ovv[mid] if k not in tags_1]
        tags_1 = tags_1 | set(ids)
        pgs[mids[mid]].append(ids)

    pgs1 = []
    for mid, pg in enumerate(pgs):
        pgs1.append(model.addPhysicalGroup(2, sum_dict(pg), mid + 1))

    if mesh_size is None:
        mesh_size = vbox.width / 10

    if isinstance(mesh_size, (int, float)):
        ms = [mesh_size] * len(pgs1)
    else:
        ms = mesh_size

    if periodic or unit_cell:
        set_periodic(model, vbox.width, vbox.height, 0)

    set_mesh_size(model, pgs1, ms)
    gmsh.option.setNumber("Mesh.Algorithm", 2)

    model.mesh.generate(2)

    filename_base = os.path.splitext(filename_svg)[0]
    filename_msh = f'{filename_base}.msh'

    gmsh.write(filename_msh)
    gmsh.finalize()

    mesh = meshio.read(filename_msh)
    mesh.cell_data = {'mat_id': mesh.cell_data['gmsh:physical']}
    mesh.point_data = {}
    mesh.cell_sets = None

    if unit_cell:
        points = mesh.points
        bbox = [nm.min(points, axis=0), nm.max(points, axis=0)]
        points -= bbox[0]
        dd = bbox[1] - bbox[0]
        for k, d in enumerate(dd[:2]):
            points[:, k] /= d
            points[nm.abs(points[:, k]) < 1e-6, k] = 0.

    print(mesh)

    if filename_out is None:
        filename_out = f'{filename_base}.vtk'

    mesh.write(filename_out, binary=False)   
    print(f'output file: {filename_out}')

    if export_png:
        filename_png = f'{filename_base}.png'
        save_png(filename_out, filename_png)
        print(f'screenshot file: {filename_png}')


def parse_args():
    parser = argparse.ArgumentParser(
        prog='SVG2MESH',
        description='SVG2MESH - Generate FE mesh using SVG defined geometry.')

    parser.add_argument('filename_svg')
    parser.add_argument('-o', '--output', action='store',
                        dest='filename_out', default=None)
    parser.add_argument('-u', '--unit-cell', action='store_true',
                        dest='unit_cell', default=False)
    parser.add_argument('-p', '--periodic', action='store_true',
                        dest='periodic', default=False)
    parser.add_argument('-s', '--mesh-size', action='store',
                        dest='mesh_size', default=None)
    parser.add_argument('-e', '--export-png', action='store_true',
                        dest='export_png', default=False)


    return parser.parse_args()


def main():
    args = parse_args()

    mesh_size = args.mesh_size
    if mesh_size is not None:
        mesh_size = literal_eval(mesh_size)

    gen_mesh_from_svg(args.filename_svg, filename_out=args.filename_out,
                      mesh_size=mesh_size, unit_cell=args.unit_cell,
                      periodic=args.periodic, export_png=args.export_png)


if __name__ == '__main__':
    main()
