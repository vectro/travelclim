#!env/bin/python

import datetime

import georasters 
from shapely.geometry import Point
import geopandas
import matplotlib
import matplotlib.pyplot
import numpy
import pandas

def open_file(month, attribute):
    # WorldClim 2.0 data from here:
    # http://worldclim.org/version2
    fname = 'wc2.0_10m_{}_{:02d}.tif'.format(attribute, month)
    g = georasters.from_file(fname)
    df = g.to_pandas()
    index = ['longitude', 'latitude']

    df = df.drop(['row', 'col'], axis='columns')
    df.columns = [attribute] + index

    if attribute in ('tmin', 'tmax', 'tavg'):
        minval = -273
    elif attribute in ('vapr', 'wind'):
        minval = 0
    else:
        assert False
    df[attribute] = df[attribute].where(df[attribute] > minval)

    df = df.set_index(index)

    return df[attribute]

def dew_point(vapr):
    # Based on formula here, this is the inverse.
    #https://www.weather.gov/media/epz/wxcalc/vaporPressure.pdf

    vapr = vapr * 10 # Convert to millibars / hectopascals
    l = numpy.log(0.16367 * vapr)
    dp = 237.3 * l / (17.26939 - l)

    dp = c_to_f(dp)
    
    dp.name = 'dp'
    return dp

def c_to_f(temp):
    # Convert temp to F.
    # https://www.weather.gov/media/epz/wxcalc/tempConvert.pdf
    return temp * 9 / 5. + 32
    

def wind_chill(temp, wind):
    temp = c_to_f(temp)

    # Convert wind to mph.
    # https://www.weather.gov/media/epz/wxcalc/windConversion.pdf
    wind = wind * 2.23694

    # From formula here
    # https://www.weather.gov/media/epz/wxcalc/windChill.pdf
    we = wind ** 0.16
    wc = 35.74 + (0.6215 * temp) - (35.75 * we) + (0.4275 * temp * we)

    wc.name = 'wc'
    return wc

def rh(temp, vapr_actual):
    # From formula here
    # https://www.weather.gov/media/epz/wxcalc/vaporPressure.pdf

    vapr_sat = 6.11 * (10 ** ((7.5 * temp) / (237.3 + temp)))
    vapr_actual = vapr_actual * 10 # Convert to millibars / hectopascals
    return vapr_actual / vapr_sat

def heat_index(temp, rh):
    # From formula here
    # https://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml
    temp = c_to_f(temp)
    hi1 = (-42.379 + 2.04901523*temp + 10.14333127*rh - .22475541*temp*rh
           - .00683783*temp*temp - .05481717*rh*rh + .00122874*temp*temp*rh
           + .00085282*temp*rh*rh - .00000199*temp*temp*rh*rh)

    adj1_where = rh.lt(.13) & temp.between(80, 112)
    adj1 = hi1 - ((13-rh)/4)*numpy.sqrt((17-numpy.abs(temp-95.))/17,
                                        where=adj1_where)
    adj2 = hi1 + ((rh-85)/10) * ((87-temp)/5)

    hi1 = adj1.where(adj1_where, hi1)
    hi1 = adj2.where(rh.gt(.85) & temp.between(80, 87), hi1)

    mhi = (hi1 + temp) / 2
    hi2 = 0.5 * (temp + 61.0 + ((temp-68.0)*1.2) + (rh*0.094))
    hi3 = hi1.where(mhi.ge(80), hi2)
    hi = hi3.where(temp.ge(70), temp)
    hi.name = 'hi'
    return hi

def draw(month):
    tmin = open_file(month, 'tmin')
    tmax = open_file(month, 'tmax')
    vapr = open_file(month, 'vapr')
    wind = open_file(month, 'wind')

    wc = wind_chill(tmin, wind)
    hi = heat_index(tmax, rh(tmax, vapr))

    ok_hi = hi.lt(80)
    ok_hi[hi.isnull()] = numpy.nan
    ok_hi.name = 'hi'
    ok_wc = wc.gt(45)
    ok_wc[wc.isnull()] = numpy.nan
    ok_wc.name = 'wc'

    gsdf = pandas.concat([ok_hi, ok_wc], axis=1)
    #gsdf = gsdf.sample(1000)
    #gsdf = gsdf.sort_index().loc(axis=0)[slice(49.,51.)]

    geom = [Point(lon, lat) for lon, lat in gsdf.index]
    gsdf = geopandas.GeoDataFrame(gsdf.reset_index(drop=True),
                                  geometry=geom,
                                  crs={'init': 'espg:4326'})

    # Natural Earth country borders from here:
    # https://www.naturalearthdata.com/downloads/50m-cultural-vectors/50m-admin-0-countries-2/
    coastline = geopandas.read_file('ne_50m_admin_0_countries.shp')

    fig, ax = matplotlib.pyplot.subplots()
    ax.set_aspect('equal')
    for dataset, color, label in [
        ((gsdf['hi'] == False) & (gsdf['wc'] == False), "Black", None),
        ((gsdf['hi'] != False) & (gsdf['wc'] == False), "Blue", "Too Cold"),
        ((gsdf['hi'] == False) & (gsdf['wc'] != False), "Red", "Too Hot"),
        ((gsdf['hi'] == True) & (gsdf['wc'] == True), "Green", "Just Right")
        ]:
        dataset = gsdf[dataset]['geometry']
        if not dataset.empty:
            dataset.plot(ax=ax, color=color, label=label, markersize=1)
    transparent = (0, 0, 0, 0)
    coastline.plot(ax=ax, color=transparent, edgecolor='black', linewidth=1)
    ax.legend(loc=3, framealpha=1.0)
    ax.set_title(datetime.date(1, month, 1).strftime("%B"))
    return fig

def main():
    for month in range(1, 13):
        print("Starting ", month)
        fig = draw(month)
        fig.set_size_inches(9.5, 4.75)
        matplotlib.pyplot.savefig("{:02d}.png".format(month), dpi=300)
        print("Wrote ", month)

def interactive_main(month):
    fig = draw(month)
    matplotlib.pyplot.show()

main()
#interactive_main(7)
