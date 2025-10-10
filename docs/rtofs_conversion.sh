# !/bin/bash
# This example script is used to convert RTOFS model output files to netcdf format using HYCOM-Tools.
# Assuming the bathymetry and grid files indicated below.

### Get bathymetry and grid files
# Regional Depth
wget https://data.hycom.org/datasets/GLBb0.08/expt_53.X/topo/depth_GLBb0.08_09m11.a -O /data/regional.depth.a

# Note: this could be downloaded from https://data.hycom.org/datasets/GLBb0.08/expt_53.X/topo/depth_GLBb0.08_09m11.b
cat > /data/regional.depth.b << EOF
bathymetery from 30-secnd GEBCO_08 20091120 global dataset, new 11x11 average
i/jdm = 4500 3298; plon,plat range = 74.12003 434.12000   -78.64000 89.97772
Filled single-width inlets and (B&C grid) enclosed seas except the Black Sea.
Clipped to range:    5.00 m to 11000.00 m. Compatible with CICE.
Topography mods to Faroe Bank Channel. Fixed Yenisei River. depth_11 coastline.
min,max depth =     5.00000  7199.97998
EOF

# Regional Grid
wget https://data.hycom.org/datasets/GLBb0.08/expt_53.X/topo/regional.grid.a -O /data/regional.grid.a

# Note: this could be downloaded from https://data.hycom.org/datasets/GLBb0.08/expt_53.X/topo/regional.grid.b
cat > /data/regional.grid.b << EOF
 4500    'idm   ' = longitudinal array size
 3298    'jdm   ' = latitudinal  array size
   12    'mapflg' = map flag (0=mercator,10=panam,12=ulon-panam)
plon:  min,max =        74.12003      434.12000
plat:  min,max =       -78.64000       89.97772
qlon:  min,max =        74.12000      434.12000
qlat:  min,max =       -78.65600       90.00000
ulon:  min,max =        74.12000      434.11993
ulat:  min,max =       -78.64000       89.98425
vlon:  min,max =        74.12006      434.12000
vlat:  min,max =       -78.65600       89.98425
pang:  min,max =        -2.67803        2.67803
pscx:  min,max =      1751.89807     8898.68848
pscy:  min,max =         0.00000     8895.61328
qscx:  min,max =         0.00000     8897.49316
qscy:  min,max =         0.00000     8895.61328
uscx:  min,max =         2.30988     8897.49512
uscy:  min,max =         0.00000     8895.61328
vscx:  min,max =      1749.46313     8898.68555
vscy:  min,max =         0.00000     8895.61328
cori:  min,max =   -0.0001429933   0.0001458425
pasp:  min,max =         0.49239       99.00000
EOF

### Get RTOFS Model Output Files

# RTOFS Binary [.ab] Files From Amazon S3 storage.
setenv RTOFS_DATE 20251008 # YYYYMMDD
# ~ 5GB download. ~12GB unzipped.
aws s3 cp s3://noaa-nws-rtofs-pds/rtofs.${RTOFS_DATE}/rtofs_glo.t00z.f06.archv.a.tgz /data/
tar -xzf /data/rtofs_glo.t00z.f06.archv.a.tgz
aws s3 cp s3://noaa-nws-rtofs-pds/rtofs.${RTOFS_DATE}/rtofs_glo.t00z.f06.archv.b /data/

### Setup conversion
set echo
setenv R GLBb0.08
#
# --- optional title and institution.
#
setenv CDF_TITLE        "HYCOM ${R}"
setenv CDF_INST         "Naval Research Laboratory"
    
# Set Output file names
# # Note: CDF030 captures any var configured in archv2ncdf3z.IN  with the '30' I/O unit... and so forth.
setenv CDF030  /data/rtofs_out.nc 

  cat > /data/archv2ncdf3z.IN << 'E-o-D'
rtofs_glo.t00z.f06.archv.a
netCDF
 000    'iexpt ' = experiment number x10 (000=from archive file)
   0    'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
 4500 'idm   ' = longitudinal array size
 3298 'jdm   ' = latitudinal  array size
  41  'kdm   ' = number of layers
  34.0  'thbase' = reference density (sigma units)
   0    'smooth' = smooth the layered fields (0=F,1=T)
   1    'iorign' = i-origin of plotted subregion
   1    'jorign' = j-origin of plotted subregion
   0    'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
   0    'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)
   3    'itype ' = interpolation type (0=sample,1=linear)
  33    'kz    ' = number of depths to sample
   0.0  'z     ' = sample depth  1
  10.0  'z     ' = sample depth  2
  20.0  'z     ' = sample depth  3
  30.0  'z     ' = sample depth  4
  50.0  'z     ' = sample depth  5
  75.0  'z     ' = sample depth  6
 100.0  'z     ' = sample depth  7
 125.0  'z     ' = sample depth  8
 150.0  'z     ' = sample depth  9
 200.0  'z     ' = sample depth 10
 250.0  'z     ' = sample depth 11
 300.0  'z     ' = sample depth 12
 400.0  'z     ' = sample depth 13
 500.0  'z     ' = sample depth 14
 600.0  'z     ' = sample depth 15
 700.0  'z     ' = sample depth 16
 800.0  'z     ' = sample depth 17
 900.0  'z     ' = sample depth 18
1000.0  'z     ' = sample depth 19
1100.0  'z     ' = sample depth 20
1200.0  'z     ' = sample depth 21
1300.0  'z     ' = sample depth 22
1400.0  'z     ' = sample depth 23
1500.0  'z     ' = sample depth 24
1750.0  'z     ' = sample depth 25
2000.0  'z     ' = sample depth 26
2500.0  'z     ' = sample depth 27
3000.0  'z     ' = sample depth 28
3500.0  'z     ' = sample depth 29
4000.0  'z     ' = sample depth 30
4500.0  'z     ' = sample depth 31
5000.0  'z     ' = sample depth 32
5500.0  'z     ' = sample depth 33
   30    'botio ' = bathymetry  I/O unit (0 no I/O)
   30    'mltio ' = mix.l.thk.  I/O unit (0 no I/O)
   1.0  'tempml' = temperature jump across mixed-layer (degC,  0 no I/O)
   0.05 'densml' =     density jump across mixed-layer (kg/m3, 0 no I/O)
   30    'infio ' = intf. depth I/O unit (0 no I/O, <0 label with layer #)
   0    'wviio ' = intf. veloc I/O unit (0 no I/O)
  0    'wvlio ' = w-velocity  I/O unit (0 no I/O)
  30    'uvlio ' = u-velocity  I/O unit (0 no I/O)
  30    'vvlio ' = v-velocity  I/O unit (0 no I/O)
   0    'splio ' = speed       I/O unit (0 no I/O)
  30    'temio ' = temperature I/O unit (0 no I/O)
  30    'salio ' = salinity    I/O unit (0 no I/O)
  30    'tthio ' = density     I/O unit (0 no I/O)
  0    'keio  ' = kinetic egy I/O unit (0 no I/O)
E-o-D
    
# Run Conversion
archv2ncdf3z <  /data/archv2ncdf3z.IN