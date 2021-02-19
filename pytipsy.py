import re
import struct
import numpy as np

def rtipsy(filename, VERBOSE=False):
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
	try:
		f = open(filename, 'rb')
	except:
		print("RTIPSY ERROR: Can't open file")
		return 1
	fs = len(f.read())
	f.seek(0)
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
		print("Read time,n,ngas,ndark,nstar: ", t, n, ng, nd, ns)
	#Catch for 4 byte padding
	if (fs == 32+48*ng+36*nd+44*ns):
		f.read(4)
	#File is borked if this is true
	elif (fs != 28+48*ng+36*nd+44*ns):
		print("RTIPSY ERROR: Header and file size inconsistent")
		print("Estimates: Header bytes:  28 or 32 (either is OK)")
		print("     ngas: ",ng," bytes:",48*ng)
		print("    ndark: ",nd," bytes:",36*nd)
		print("    nstar: ",ns," bytes:",44*ns)
		print("Actual File bytes:",fs,"  not one of:",28+48*ng+36*nd+44*ns,32+48*ng+36*nd+44*ns)
		f.close()
		return 1

	catg = {'mass':np.zeros(ng), 'pos':np.zeros((ng,3)), 'vel':np.zeros((ng,3)), 'dens':np.zeros(ng),
			'tempg':np.zeros(ng), 'h':np.zeros(ng), 'zmetal':np.zeros(ng),	'phi':np.zeros(ng)}
	catd = {'mass':np.zeros(nd), 'pos':np.zeros((nd,3)), 'vel':np.zeros((nd,3)),
			'eps':np.zeros(nd), 'phi':np.zeros(nd)}
	cats = {'mass':np.zeros(ns), 'pos':np.zeros((ns,3)), 'vel':np.zeros((ns,3)),
			'metals':np.zeros(ns), 'tform':np.zeros(ns), 'eps':np.zeros(ns), 'phi':np.zeros(ns)}
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
	header = {'time':t, 'n':n, 'ndim':ndim, 'ngas':ng, 'ndark':nd, 'nstar':ns}
	return (header,catg,catd,cats)

def wtipsy(filename, header, catg, catd, cats, STANDARD=True):
	try:
		f = open(filename, 'wb')
	except:
		print("WTIPSY ERROR: Can't open file")
		return 1
	f.write(struct.pack(">diiiii", header['time'], header['n'], header['ndim'], header['ngas'], header['ndark'], header['nstar']))
	if STANDARD:
		f.write(struct.pack("xxxx"))
	for i in range(header['ngas']):
		f.write(struct.pack(">ffffffffffff", catg['mass'][i], catg['x'][i], catg['y'][i], catg['z'][i], catg['vx'][i], catg['vy'][i], 
			catg['vz'][i], catg['dens'][i], catg['tempg'][i], catg['h'][i], catg['zmetal'][i], catg['phi'][i]))
	for i in range(header['ndark']):
		f.write(struct.pack(">fffffffff", catd['mass'][i], catd['x'][i], catd['y'][i], catd['z'][i], catd['vx'][i], catd['vy'][i], 
			catd['vz'][i], catd['eps'][i], catd['phi'][i]))
	for i in range(header['nstar']):
		f.write(struct.pack(">fffffffffff", cats['mass'][i], cats['x'][i], cats['y'][i], cats['z'][i], cats['vx'][i], cats['vy'][i], 
			cats['vz'][i], cats['metals'][i], cats['tform'][i], cats['eps'][i], cats['phi'][i]))
	f.close()

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
