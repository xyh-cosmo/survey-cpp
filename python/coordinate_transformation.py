from astropy import units as u
from astropy.coordinates import SkyCoord,Galactic,GeocentricTrueEcliptic

print('> Testing transformation from Galactic to Ecliptic coordinate frame')
c_gal = SkyCoord(l=180*u.degree,b=0*u.degree,frame='galactic')
c_epl_g = c_gal.geocentrictrueecliptic
c_epl_b = c_gal.barycentrictrueecliptic

print(c_gal)
print(c_epl_g.lon.deg,c_epl_g.lat.deg)
print(c_epl_b.lon.deg,c_epl_b.lat.deg)


print('> Testing inverse transformation, from Eclptic to Galactic')
lon,lat=c_epl_g.lon.deg,c_epl_g.lat.deg
c_epl_g = GeocentricTrueEcliptic(lon=lon*u.degree,lat=lat*u.degree)
c_gal = SkyCoord(c_epl_g).galactic  # this is fine
# c_gal = Galactic(c_epl_g)  # this does not work!

print(c_epl_g.lon.deg,c_epl_g.lat.deg)
print(c_gal)
