"""
Append a feature to the geojson file
"""

# Libraries
import json
import unicodedata
import pandas as pd

# Functions
def strip_accents(string_to_strip):
    """
    This function takes a string that potentially has a diacritic or accent (i.e.,
    a unicode string) and converts it to ASCII text.
    """

    return ''.join(c for c in unicodedata.normalize('NFD', string_to_strip)
                  if unicodedata.category(c) != 'Mn')

def append_geojson(country, fitted):
    """
    Add a feature from the Stan output files to a geojson

    Inputs
    country: one of ['Argentina', 'Brazil', 'country_geojsonguay']
    fitted: a pd.Dataframe with the fitted values from the Stan model

    Returns
    country_geojson: a geojson formatted object with the appended fitted feature
   """

    fitted['AL1'] = fitted['AL1'].str.lower()

    # Read in geojson
    country_geojson = json.load(open("../data/geojson/" + country + ".geojson"))

    names = [country_geojson['features'][i]['properties']['name'].lower() \
            for i in range(len(country_geojson['features']))]

    # Append feature to geojson
    for i in range(len(country_geojson['features'])):

        name = names[i]

        if name == 'buenos aires':
            name = 'caba'

        index = fitted['AL1'] == name

        if any(index):

            print(name)
            country_geojson['features'][i]['properties']['r'] = fitted.loc[index]['.mean'].iloc[0]
            country_geojson['features'][i]['properties']['clean_name'] = name

    return country_geojson

# Globals
countries = ['argentina', 'brazil', 'paraguay']

# Model output
SA_fit = pd.read_csv('../data/BA5_fit.csv')

geojson = {}

brazil = append_geojson('brazil', SA_fit)
argentina = append_geojson('argentina', SA_fit)
paraguay = append_geojson('paraguay', SA_fit)

with open('brazil.geojson', 'w') as f:
    json.dump(brazil, f)
with open('argentina.geojson', 'w') as f:
    json.dump(argentina, f)
with open('paraguay.geojson', 'w') as f:
    json.dump(paraguay, f)
