""" FIXME """
try:
    import cPickle as pickle
except:
    import pickle

""" load test data """
class Dict2Struc(object):
    """ all the variables from a dict in a "matlab-like-structure" """
    def __init__(self, adict):
        self.__dict__.update(adict)

data = pickle.load( open('gsw_cv.pkl','rb') )
gsw_cv = Dict2Struc(data) # then type dat.<tab> to navigate through your variables
SP = gsw_cv.SP_chck_cast
p = gsw_cv.p_chck_cast
lon = gsw_cv.long_chck_cast
lat  = gsw_cv.lat_chck_cast
SA_from_SP(SP, p, lon, lat)