import sys
import numpy as np
import calendar
from calendar import monthrange
from datetime import datetime
from geopy.geocoders import Nominatim
def Sort(sub_li): 
    # reverse = None (Sorts in Ascending order) 
    # key is set to sort using second element of  
    # sublist lambda has been used 
    sub_li.sort(key = lambda x: x[2]) 
    return sub_li 
def get_best_neighbor(humidity_df,tl_lat,tl_long):
    scan = []
    for i in range(-25,26,1):
        for m in range(-25,26,1):
            if humidity_df[800,tl_lat+i,tl_long+m] != -1.0:
                scan.append([tl_lat+i,tl_long+m,abs(i)+abs(m)])
    sscan = Sort(scan)
    return sscan[0] 

def get_daily_mean_humidity(humidity_df,myyear,mymonth,myday,tl_lat,tl_long):
    # Get the day number in the dataframe
    x = datetime(myyear, mymonth, myday)
    # get number in the year
    day_of_year = x.timetuple().tm_yday
    # add days in previous years
    start_year = 2014
    increment_days = 0
    
    if myyear == start_year:
        tl_day = day_of_year - 1
    else:
        for i in range(0,myyear - start_year):
            cur_year = start_year + i
             
    # Feb 29 is the extra day
            if calendar.isleap(cur_year):
               days = 366
            else:
               days = 365
            increment_days = increment_days + days
        tl_day = increment_days
    return humidity_df[tl_day,tl_lat,tl_long]


# Generate the humidity dataframe using loop below
#import netCDF4 as nc
#months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
#humidity_df = np.empty((0, 360, 720))
#for yyyy in range(2014, 2019):
#    for mm in months:
#        ds = nc.Dataset(f'/local1/home/mo/Rob_COVID/ERA5_humidity/Qair_WFDE5_CRU_{yyyy}{mm}_v1.0.nc')
#        Qair = np.flip(ds['Qair'][:], 1)
#        daily_Qair = np.mean(Qair.reshape((-1, 24, 360, 720)), axis=1)
#        humidity_df = np.concatenate((humidity_df, daily_Qair.filled(fill_value=-1)))

# dataframe is 3-dimensional array:
# 1826 days from 2014 to 2019
# 360 grids lat
# 720 grids long


humidity_df = np.load("daily_hum_14_18_final.npy")
#indate = "2020-01-01"
#incity = "London"
inyear = int(sys.argv[1])
incity = sys.argv[2]

geolocator = Nominatim(user_agent="Google Maps")
location = geolocator.geocode(incity)
real_lat = location.latitude
real_long = location.longitude

tl_lat = round(abs(180 - real_lat*2))
tl_long = round(abs(real_long*2 + 360))
if humidity_df[800,tl_lat,tl_long] == -1.0:
    nearest = get_best_neighbor(humidity_df,tl_lat,tl_long)
    tl_lat = nearest[0]
    tl_long = nearest[1]


#print("Real",real_lat,real_long)
#print("Grid",tl_lat,tl_long)

years = [2014,2015,2016,2017,2018]
#mydate = indate.split("-")
#myyear = int(mydate[0])
#mymonth = int(mydate[1])
#myday = int(mydate[2])

for month in range(1, 13): # Month is always 1..12
    for day in range(1, monthrange(inyear, month)[1] + 1):
        if calendar.isleap(inyear) and month == 2 and day == 29:
            humidity = get_daily_mean_humidity(humidity_df,2016,month,day,tl_lat,tl_long)
        else:
            humlist = []
            for y in years:
                yhumidity = get_daily_mean_humidity(humidity_df,y,month,day,tl_lat,tl_long) 
                humlist.append(yhumidity)
                
            humidity = sum(humlist) / len(humlist)
        humidity = round(humidity,10)
        omonth = month
        oday = day
        if month < 10:
            omonth = "0" + str(month)
        if day < 10:
            oday = "0" + str(day)
        indate = str(inyear) + "-" + str(omonth) + "-" + str(oday)
        print(*[indate,incity,humidity],sep=",")


