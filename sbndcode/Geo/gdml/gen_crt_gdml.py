'''Generate the GDML for the SBND CRT geometry.

We define a global XML DOM which the function for each volume will append to.
Volume-building functions are called hierarchically, so if you call module(),
it will construct the module and all the parts that make it up, so you end
up with complete GDML for one module.

Each physical volume has a corresponding unique logical volume, as required
by LArG4 to keep track of energy depositions. The solids, however, can safely
be referenced many times, and so are stored only once (using a hash keyed on
the the linear dimensions).

The default is that the "unit cell" is a module (made of a cover and one
monolithic chunk of scintillator, where downstream code will figure out which
strip should be hit). To build the individual strips as part of the geometry,
set do_strips to True (this will require changes in the aforementioned
downstream code).

The output of this code is a file "crt.gdml" which contains the GDML snippets
to paste into the full SBND geometry.

:todo: Add bottom tagger (which is not simply a grid)
:todo: Make volume names based on something repeatable, not UUIDs

A. Mastbaum <mastbaum@uchicago.edu>, 2016/10/27
'''

import xml.etree.cElementTree as ET
from xml.dom import minidom
import uuid
from hashlib import md5

# Options
PAD = 0.4  # Padding between strips and module (Al thickness)
do_strips = False  # Build individual scintillator strips

##########

solids = ET.Element('solids')
structure = ET.Element('structure')
solids_store = {}

def strip(x, y, z=1.0):
    '''Build one scintillator strip.'''
    x = str(x)
    y = str(y - 2 * PAD)
    z = str(z)
    sname = 'CRTStrip_' + md5('_'.join((xx, yy, zz))).hexdigest()[:8]

    if not sname in solids_store:
        s = ET.SubElement(solids, 'box', name=sname, lunit="cm", x=x, y=y, z=z)
        solids_store[sname] = s
    else:
        s = solids_store[name]

    name = 'CRTStrip_' + uuid.uuid4().hex[:8]
    vname = 'volAuxDetSensitive' + name
    v = ET.SubElement(structure, 'volume', name=vname)
    ET.SubElement(v, 'materialref', ref='Polystyrene')
    ET.SubElement(v, 'solidref', ref=sname)

    return s, v


def strip_array(x, y, z=1.0, nx=16):
    '''Build an edge-to-edge array of scintillator strips.'''
    xx = str(x * nx)
    yy = str(y - 2 * PAD)
    zz = str(z)
    sname = 'CRTStripArray_' + md5('_'.join((xx, yy, zz))).hexdigest()[:8]

    if not sname in solids_store:
        s = ET.SubElement(solids, 'box', name=sname, lunit="cm", x=xx, y=yy, z=zz)
        solids_store[sname] = s
    else:
        s = solids_store[sname]

    if do_strips:
        for i in range(nx):
            es, ev = strip(x, y, z)

            pv = ET.SubElement(v, 'physvol')
            ET.SubElement(pv, 'volumeref', ref=ev.attrib['name'])

            dx = 1.0 * (2*i - nx + 1) / 2 * x
            posname = 'pos' + ev.attrib['name']
            ET.SubElement(pv, 'position', name=posname,
                          unit="cm", x=str(dx), y='0', z='0')

    name = 'CRTStripArray_' + uuid.uuid4().hex[:8]
    vname = 'vol' + '' if do_strips else 'AuxDetSensitive' + name
    v = ET.SubElement(structure, 'volume', name=vname)
    ET.SubElement(v, 'materialref', ref='Polystyrene')
    ET.SubElement(v, 'solidref', ref=sname)

    return s, v

def module(x, y, z=1.0, nx=16):
    '''Build a module: a strip array and an aluminum cover.'''
    # Module
    xx = str(x * nx + 2 * PAD)
    yy = str(y)
    zz = str(z + 2 * PAD)
    sname = 'CRTModule_' + md5('_'.join((xx, yy, zz))).hexdigest()[:8]

    if not sname in solids_store:
        s = ET.SubElement(solids, 'box', name=sname, lunit="cm", x=xx, y=yy, z=zz)
        solids_store[sname] = s
    else:
        s = solids_store[sname]

    # Strip array
    es, ev = strip_array(x, y, z)

    # Cover (module - strip array)
    cover_sname = 'CRTModuleCover_' + md5('_'.join((xx, yy, zz))).hexdigest()[:8]

    if not cover_sname in solids_store:
        #s = ET.SubElement(solids, 'box', name=cover_sname, lunit="cm", x=xx, y=yy, z=zz)
        scover = ET.SubElement(solids, 'subtraction', name=cover_sname)
        ET.SubElement(scover, 'first', ref=sname)
        ET.SubElement(scover, 'second', ref=es.attrib['name'])
        solids_store[cover_sname] = scover
    else:
        scover = solids_store[cover_sname]

    # Volumes
    # Cover
    cover_name = 'CRTModuleCover_' + uuid.uuid4().hex[:8]
    vcname = 'vol' + cover_name
    vc = ET.SubElement(structure, 'volume', name=vcname)
    ET.SubElement(vc, 'materialref', ref='ALUMINUM_Al')
    ET.SubElement(vc, 'solidref', ref=cover_sname)

    # Module
    name = 'CRTModule_' + uuid.uuid4().hex[:8]
    vname = 'vol' + name
    v = ET.SubElement(structure, 'volume', name=vname)
    ET.SubElement(v, 'materialref', ref='ALUMINUM_Al')
    ET.SubElement(v, 'solidref', ref=sname)

    pvc = ET.SubElement(v, 'physvol')
    poscname = 'pos' + cover_name
    ET.SubElement(pvc, 'volumeref', ref=vcname)
    ET.SubElement(pvc, 'position', name=poscname, unit="cm", x='0', y='0', z='0')

    pvs = ET.SubElement(v, 'physvol')
    possname = 'pos' + ev.attrib['name']
    ET.SubElement(pvs, 'volumeref', ref=ev.attrib['name'])
    ET.SubElement(pvs, 'position', name=possname, unit="cm", x='0', y='0', z='0')

    return s, v


def plane(x, y, nx, ny, z=1.0, nxs=16):
    '''Build a plane: an Nx by Ny grid of modules.'''
    xx = str((x * nxs + 2 * PAD) * nx)
    yy = str(y * ny)
    zz = str(z + 2 * PAD)
    sname = 'CRTPlane_' + md5('_'.join((xx, yy, zz))).hexdigest()[:8]

    if not sname in solids_store:
        s = ET.SubElement(solids, 'box', name=sname, lunit="cm", x=xx, y=yy, z=zz)
        solids_store[sname] = s
    else:
        s = solids_store[sname]

    modules = []
    for i in range(nx):
        for j in range(ny):
            modules.append(module(x, y, z, nxs))
 
    name = 'CRTPlane_' + uuid.uuid4().hex[:8]
    vname = 'vol' + name
    v = ET.SubElement(structure, 'volume', name=vname)
    ET.SubElement(v, 'materialref', ref='ALUMINUM_Al')
    ET.SubElement(v, 'solidref', ref=sname)

    for i in range(nx):
        for j in range(ny):
            es, ev = modules[j + i*ny]

            pv = ET.SubElement(v, 'physvol')
            ET.SubElement(pv, 'volumeref', ref=ev.attrib['name'])

            dx = 1.0 * (2*i - nx + 1) / 2 * (x * nxs + 2 * PAD)
            dy = 1.0 * (2*j - ny + 1) / 2 * y

            posname = 'pos' + ev.attrib['name']
            ET.SubElement(pv, 'position', name=posname,
                          unit="cm", x=str(dx), y=str(dy), z='0')

    return s, v


def tagger(x1, y1, nx1, ny1, x2, y2, nx2, ny2, z=1.0, nxs=16):
    '''Build a tagger: a stack of two perpendicular planes.'''
    xx = str((x1 * nxs + 2 * PAD) * nx1)
    yy = str(y1 * ny1)
    zz = str(2.0 * (z + 2 * PAD))
    sname = 'Tagger_' + md5('_'.join((xx, yy, zz))).hexdigest()[:8]

    if not sname in solids_store:
        s = ET.SubElement(solids, 'box', name=sname, lunit="cm", x=xx, y=yy, z=zz)
        solids_store[sname] = s
    else:
        s = solids_store[sname]

    s1, v1 = plane(x1, y1, nx1, ny1)
    s2, v2 = plane(x2, y2, nx2, ny2)
 
    name = 'Tagger_' + uuid.uuid4().hex[:8]
    vname = 'volAuxDet' + name
    v = ET.SubElement(structure, 'volume', name=vname)
    ET.SubElement(v, 'materialref', ref='ALUMINUM_Al')
    ET.SubElement(v, 'solidref', ref=sname)

    dz = (z + 2 * PAD) / 2

    # Plane 1: 0 degrees
    pv = ET.SubElement(v, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=v1.attrib['name'])
    posname = 'pos' + v1.attrib['name']
    ET.SubElement(pv, 'position', name=posname,
                  unit="cm", x='0', y='0', z=str(-dz))
    
    # Plane 2: 90 degrees
    pv = ET.SubElement(v, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=v2.attrib['name'])
    posname = 'pos' + v2.attrib['name']
    ET.SubElement(pv, 'position', name=posname,
                  unit="cm", x='0', y='0', z=str(dz))
    posname = 'rot' + v2.attrib['name']
    ET.SubElement(pv, 'rotation', name=posname,
                  unit="deg", x='0', y='0', z='90')

    return s, v


# Build the taggers
tths, tthv = tagger(11.2, 450.0, 5, 2, 11.2, 450.0, 5, 2)  # Top H
ttls, ttlv = tagger(11.2, 450.0, 5, 2, 11.2, 450.0, 5, 2)  # Top L
tsls, tslv = tagger(11.2, 360.0, 5, 2, 11.2, 450.0, 4, 2)  # Side L
tsrs, tsrv = tagger(11.2, 360.0, 5, 2, 11.2, 450.0, 4, 2)  # Side R
tffs, tffv = tagger(11.2, 360.0, 4, 2, 11.2, 360.0, 4, 2)  # Face F
tfbs, tfbv = tagger(11.2, 360.0, 4, 2, 11.2, 360.0, 4, 2)  # Face B
# FIXME: Add bottom tagger, which isn't simply a grid.

# Generate GDML for the world volume, for testing
#ws = ET.SubElement(solids, 'box', name='World', lunit="cm", x='10000', y='10000', z='10000')
#w = ET.SubElement(structure, 'volume', name='volWorld')
#ET.SubElement(w, 'materialref', ref='Air')
#ET.SubElement(w, 'solidref', ref='World')
#pv = ET.SubElement(w, 'physvol')
#ET.SubElement(pv, 'volumeref', ref=tffv.attrib['name'])
#ET.SubElement(pv, 'position', name='posA', unit="cm", x='0', y='0', z='0')
#setup = ET.Element('setup', name='Default', version='1.0')
#ET.SubElement(setup, 'world', ref='volWorld')

with open('crt.gdml', 'w') as f:
    f.write(minidom.parseString(ET.tostring(solids)).toprettyxml(indent=' '))
    f.write(minidom.parseString(ET.tostring(structure)).toprettyxml(indent=' '))

