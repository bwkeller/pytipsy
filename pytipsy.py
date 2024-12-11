import re
import struct
import numpy as np
import os

#STANDARD = False
# < means little endian (tipsy bindary, intel)
#STANDARD = True
# > means big endian (tipsy std, SUN)

def checkarray(filename, VERBOSE=True):
    """checkarray Checks tipsy array files detecting the format: 
    big endian, little endian, number 
    but does not read the data

    Usage: 
          checktipsy(filename, VERBOSE=False)

    Input parameters: 
    filename  filename string
    VERBOSE  print messages (optional)
    Return values:
    (file,header,endianswap)
    """
    try:
        f = open(filename, 'rb')
    except:
        print("array ERROR: Can't open file",filename)
        return (None,None,None)

    fs = os.fstat(f.fileno()).st_size
    #Read in the header
    n, = struct.unpack("<i", f.read(4))
    endianswap = False
    #Confirm Endianness
    if (fs != 4+4*n):
        endianswap = True
        f.seek(0)
        nswap, = struct.unpack(">i", f.read(4))
        if (fs != 4+4*nswap):
        f.close()
        if (VERBOSE):
            print("RTIPSY ERROR: Header (native: %d std: %d) and file size n (%d) inconsistent" % (n,nswap,(fs-4)//4) )
        return (None,None,None)
        n=nswap

    return(f,n,endianswap)

def rarray(filename, STANDARD_list=None, INTEGER=False, VERBOSE=False ):
    """rarray reads tipsy array files detecting the format: 
    big endian, little endian, number and returns the data
    
    Usage: 
    rarray(filename, [STANDARD_list,] VERBOSE=False)
    
    Input parameters: 
    filename  filename string
    STANDARD_list  optional list where you can save endianswap True/False
    INTEGER  assume integer data
    VERBOSE  print messages (optional)
    Return values:
    read array or 1 for fail
    """

    f,n,endianswap = checkarray(filename, VERBOSE=VERBOSE)
    if (f == None): return 1

    if (type(STANDARD_list) == list):
    STANDARD_list.insert(0,endianswap)

    readformat = '>%s' % (n) if endianswap else '<%s' % (n) 
    readformat += 'i' if INTEGER else 'f' 

    data = np.array(struct.unpack(readformat, f.read(n*4)))
    if (len(data)!=n):
    print("Failed to read all data",len(data),n)
    elif (VERBOSE):
    print("Succesfully read all data",len(data),n)
    f.close()
    return data

def warray(filename, data, STANDARD=True, VERBOSE=False):
    """warray writes tipsy array files using the given format
    default (STANDARD  big endian)
    
    Usage: 
    warray(filename, data, [STANDARD=True/False,] [VERBOSE=True/False])
    
    Input parameters: 
    filename  filename string
    data is a numpy array of floats or ints
    STANDARD=True  means write standard big endian
    VERBOSE=True   reports success
    Return values:
    0 success  1 fail
    """

    try:
    f = open(filename, 'wb')
    except:
    print("warray ERROR: Can't open file")
    return 1

    if (STANDARD):                
    f.write(struct.pack(">i", len(data)))
    writeformat = '>'
    else:                
    f.write(struct.pack("<i", len(data)))
    writeformat = '<'

    if (len(data)>0):
    writeformat += '%s' % (len(data)) 
    writeformat += 'i' if 'int' in str(type(data[0])) else 'f'

    f.write(struct.pack(writeformat, *data))
    if (VERBOSE):
        print("Wrote all data",len(data))

    f.close()
    return 0

def subset( a, i ):
   b={}
   for key, value in a.items():
      b[ key ] = a[ key ][i]
   return b
            
def checktipsy(filename, VERBOSE=False):
    """checktipsy Checks tipsy files detecting the format: 
    big endian, little endian, padded (standard) or non-padded header 
    but does not read the data

    Usage: 
          checktipsy(filename, VERBOSE=False)

    Input parameters: 
    filename  filename string
    VERBOSE  print messages (optional)
    Return values:
    (file,header,endianswap)
    """

    try:
        f = open(filename, 'rb')
    except:
        print("TIPSY ERROR: Can't open file",filename)
        return (None,None,None)
    fs = os.fstat(f.fileno()).st_size
    #Read in the header
    t, n, ndim, ng, nd, ns = struct.unpack("<diiiii", f.read(28))
    endianswap = False
    #Check Endianness
    if (ndim < 1 or ndim > 3):
        endianswap = True
        f.seek(0)
        t, n, ndim, ng, nd, ns = struct.unpack(">diiiii", f.read(28))
        if VERBOSE:
            print("SWAP_ENDIAN")
    if VERBOSE:
        print("Header: time,n,ngas,ndark,nstar: ", t, n, ng, nd, ns)
    #Catch for 4 byte padding
    if (fs == 32+48*ng+36*nd+44*ns):
        f.read(4)
    #File is borked if this is true
    elif (fs != 28+48*ng+36*nd+44*ns):
        print("TIPSY ERROR: Header and file size inconsistent")
        print("Estimates: Header bytes:  28 or 32 (either is OK)")
        print("     ngas: ",ng," bytes:",48*ng)
        print("    ndark: ",nd," bytes:",36*nd)
        print("    nstar: ",ns," bytes:",44*ns)
        print("Actual File bytes:",fs,"  not one of:",28+48*ng+36*nd+44*ns,32+48*ng+36*nd+44*ns)
        f.close()
        return (None,None,None)

    return(f,(t, n, ndim, ng, nd, ns), endianswap)

def rtipsy(filename, STANDARD_list=None, VERBOSE=False):
    """rtipsy Reads tipsy files detecting the format: 
    big endian, little endian, padded (standard) or non-padded header 

    Usage: 
          rtipsy(filename, VERBOSE=False)

    Input parameters: 
    filename  filename string
    VERBOSE  print messages (optional)
    Return values:
    (header,g,d,s)
    header    tipsy header struct
    g,d,s     gas, dark and star structures
    Please read rtipsy.py for the structure definitions

    Example: 
    h,g,d,s = rtipsy('/home/wadsley/usr5/mihos/mihos.std')
    print, h['ndark']
    plt.plot(d['x'], d['y'], 'k,')"""

    f,header,endianswap = checktipsy(filename, VERBOSE=VERBOSE)
    t, n, ndim, ng, nd, ns = header
    if (type(STANDARD_list) == list):
        STANDARD_list.insert(0,endianswap)

    catg = {'mass':np.zeros(ng), 'pos':np.zeros((ng,3)),
            'vel':np.zeros((ng,3)), 'dens':np.zeros(ng),
            'tempg':np.zeros(ng), 'h':np.zeros(ng), 'zmetal':np.zeros(ng),
            'phi':np.zeros(ng)}
    catd = {'mass':np.zeros(nd), 'pos':np.zeros((nd,3)),
            'vel':np.zeros((nd,3)),          'eps':np.zeros(nd),
            'phi':np.zeros(nd)}
    cats = {'mass':np.zeros(ns), 'pos':np.zeros((ns,3)),
            'vel':np.zeros((ns,3)),          'metals':np.zeros(ns),
            'tform':np.zeros(ns), 'eps':np.zeros(ns), 'phi':np.zeros(ns)}
    for cat in ['g','d','s']:
        j = 0
        for qty in ['x','y','z']:
            locals()['cat'+cat][qty] = locals()['cat'+cat]['pos'][:,j]
            locals()['cat'+cat]['v'+qty] = locals()['cat'+cat]['vel'][:,j]
            j += 1

    if (ng > 0):
        for i in range(ng):
            if endianswap:
                mass, x, y, z, vx, vy, vz, dens, tempg, h, zmetal, phi = struct.unpack(">ffffffffffff", f.read(48))
            else:
                mass, x, y, z, vx, vy, vz, dens, tempg, h, zmetal, phi = struct.unpack("<ffffffffffff", f.read(48))
            catg['mass'][i] = mass
            catg['x'][i] = x
            catg['y'][i] = y
            catg['z'][i] = z
            catg['vx'][i] = vx
            catg['vy'][i] = vy
            catg['vz'][i] = vz
            catg['dens'][i] = dens
            catg['tempg'][i] = tempg
            catg['h'][i] = h
            catg['zmetal'][i] = zmetal
            catg['phi'][i] = phi
    if (nd > 0):
        for i in range(nd):
            if endianswap:
                mass, x, y, z, vx, vy, vz, eps, phi = struct.unpack(">fffffffff", f.read(36))
            else:
                mass, x, y, z, vx, vy, vz, eps, phi = struct.unpack("<fffffffff", f.read(36))
            catd['mass'][i] = mass
            catd['x'][i] = x
            catd['y'][i] = y
            catd['z'][i] = z
            catd['vx'][i] = vx
            catd['vy'][i] = vy
            catd['vz'][i] = vz
            catd['eps'][i] = eps
            catd['phi'][i] = phi
    if (ns > 0):
        for i in range(ns):
            if endianswap:
                mass, x, y, z, vx, vy, vz, metals, tform, eps, phi = struct.unpack(">fffffffffff", f.read(44))
            else:
                mass, x, y, z, vx, vy, vz, metals, tform, eps, phi = struct.unpack("<fffffffffff", f.read(44))
            cats['mass'][i] = mass
            cats['x'][i] = x
            cats['y'][i] = y
            cats['z'][i] = z
            cats['vx'][i] = vx
            cats['vy'][i] = vy
            cats['vz'][i] = vz
            cats['metals'][i] = metals
            cats['tform'][i] = tform
            cats['eps'][i] = eps
            cats['phi'][i] = phi

            
    return (header,catg,catd,cats)

def wtipsy(filename, header, catg, catd, cats, STANDARD=True, VERBOSE=False):
    """wtipsy  Write tipsy files in selected format
    big endian, little endian, padded (standard) or non-padded header 

    Usage: 
         wtipsy(filename, header, g, d, s, STANDARD=True, VERBOSE=False)

    Input parameters: 
    filename  filename string
    header    tipsy header struct
    g,d,s     gas, dark and star structures
    STANDARD True for standard big endian
    VERBOSE  print messages (optional)
    Return values:
      0 success  1 fail
    Please read pytipsy.py for the structure definitions
    """

    try:
        f = open(filename, 'wb')
    except:
        print("WTIPSY ERROR: Can't open file")
        return 1

    endian='>' if STANDARD else '<'

    f.write(struct.pack(endian+"diiiii", header['time'], header['n'], header['ndim'], header['ngas'], header['ndark'], header['nstar']))
    if (STANDARD):
        f.write(struct.pack("xxxx"))
        if (VERBOSE):
        print("STANDARD Write. Header: ",header)
    elif (VERBOSE):
        print("Native Write. Header: ",header)

    for i in range(header['ngas']):
        f.write(struct.pack(endian+"ffffffffffff", catg['mass'][i], catg['x'][i], catg['y'][i], catg['z'][i], catg['vx'][i], catg['vy'][i], 
            catg['vz'][i], catg['dens'][i], catg['tempg'][i], catg['h'][i], catg['zmetal'][i], catg['phi'][i]))
    for i in range(header['ndark']):
        f.write(struct.pack(endian+"fffffffff", catd['mass'][i], catd['x'][i], catd['y'][i], catd['z'][i], catd['vx'][i], catd['vy'][i], 
            catd['vz'][i], catd['eps'][i], catd['phi'][i]))
    for i in range(header['nstar']):
        f.write(struct.pack(endian+"fffffffffff", cats['mass'][i], cats['x'][i], cats['y'][i], cats['z'][i], cats['vx'][i], cats['vy'][i], 
            cats['vz'][i], cats['metals'][i], cats['tform'][i], cats['eps'][i], cats['phi'][i]))
    f.close()
    return 0

class gaslog(dict):
    def __init__(self, fname):
        self.rawdata = np.genfromtxt(fname, comments='#', dtype=None, names=['dTime', 'z', 'E', 'T', 'U', 'Eth', 'Lx', 'Ly', 'Lz',
            'WallTime', 'dwMax', 'dIMax', 'dEMax', 'dMultiEff'])
        for name in self.rawdata.dtype.names:
            self[name] = self.rawdata[name]
        self.units = {'erg': float(re.findall('dErgPerGmUnit:\s*[0-9,.,e,+,-]*', open(fname).read())[0].split()[1]) * \
                float(re.findall('dMsolUnit:\s*[0-9,.,e,+,-]*', open(fname).read())[0].split()[1]) * 1.9891e33, 'yr': \
        float(re.findall('dSecUnit:\s*[0-9,.,e,+,-]*', open(fname).read())[0].split()[1]) / 3.1557e7}
        self['dTime'] *= self.units['yr']
        self['E'] *= self.units['erg']
        self['T'] *= self.units['erg']
        self['U'] *= self.units['erg']
        self['Eth'] *= self.units['erg']

class starlog(dict):
    def __init__(self, fname):
        try:
            f = open(fname, 'rb')
        except:
            print("Cannot open starlog!")
            return 1
        # Calculate number of entries
        n_sf = int((f.seek(0, os.SEEK_END)-4)/96)
        self['iOrdStar'] = np.zeros(n_sf, dtype=np.int64)
        self['iOrdGas'] = np.zeros(n_sf, dtype=np.int64)
        self['timeForm'] = np.zeros(n_sf)
        self['xForm'] = np.zeros(n_sf)
        self['yForm'] = np.zeros(n_sf)
        self['zForm'] = np.zeros(n_sf)
        self['vxForm'] = np.zeros(n_sf)
        self['vyForm'] = np.zeros(n_sf)
        self['vzForm'] = np.zeros(n_sf)
        self['massForm'] = np.zeros(n_sf)
        self['rhoForm'] = np.zeros(n_sf)
        self['TForm'] = np.zeros(n_sf)
        f.seek(4) # Remove 4-byte pad
        for i in range(n_sf):
            iOrdStar, iOrdGas, timeForm, xForm, yForm, zForm, vxForm, vyForm, vzForm, massForm, \
            rhoForm, Tform = struct.unpack('>qqdddddddddd', f.read(96))
            self['iOrdStar'][i] = iOrdStar
            self['iOrdGas'][i] = iOrdGas
            self['timeForm'][i] = timeForm
            self['xForm'][i] = xForm
            self['yForm'][i] = yForm
            self['zForm'][i] = zForm
            self['vxForm'][i] = vxForm
            self['vyForm'][i] = vyForm
            self['vzForm'][i] = vzForm
            self['massForm'][i] = massForm
            self['rhoForm'][i] = rhoForm
            self['TForm'][i] = Tform
