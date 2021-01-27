import sys
from astral import LocationInfo
import datetime
from calendar import monthrange
from astral.sun import sun
from astral.geocoder import database, lookup
from geopy.geocoders import Nominatim
from timezonefinder import TimezoneFinder

incity = sys.argv[1]
geolocator = Nominatim(user_agent="Google Maps")
location = geolocator.geocode(incity)
real_lat = location.latitude
real_long = location.longitude
obj = TimezoneFinder() 
mytmz = obj.timezone_at(lng=real_long, lat=real_lat)
loc = LocationInfo(name='NA', region='NA', timezone=mytmz,latitude=real_lat, longitude=real_long)
# non-leap year for reference
ryear = 2018
for month in range(1, 13): # Month is always 1..12
    for day in range(1, monthrange(ryear, month)[1] + 1):
        s = sun(loc.observer, date=datetime.date(ryear, month, day), tzinfo=mytmz)
        suntime = s["sunrise"].hour * 60 + s["sunrise"].minute
        mydate = str(datetime.date(ryear, month, day)) 
        outline = [mydate,incity,suntime]
        print(*outline,sep=",")
