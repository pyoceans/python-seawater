def test_print(STP, method, comp_value=None):
    """
    Run a function test mimicking the original logic. This is done to allow for a direct comparison of the result from the Matlab to the python package.
    """

    if comp_value is None:
        comp_value = method

    # test for floating differences with: computed - check_value >= defined_precision
    exec( "unequal = (gsw_cv." +comp_value+ " - STP." +method+ "() ) >= gsw_cv." +comp_value+ "_ca")

    width = 23
    if unequal.any():
        print "%s: Failed" % method.rjust(width)
    else:
        # test if check value is identical to computed value
        if eval( "( gsw_cv." +comp_value+ "[~np.isnan(gsw_cv."+comp_value+")] == STP." +method+ "()[~np.isnan(STP."+method+"())] ).all()" ):
            print "%s: Passed" % method.rjust(width)
        else:
            # test for differences in case their aren't equal. This is an attempt to place all tests together (i.e. term25 and small float differences that will appear)
            exec("nmax = ( gsw_cv."+comp_value+" - STP."+method+"() )[~np.isnan(gsw_cv."+comp_value+")].max()")
            exec("nmin = ( gsw_cv."+comp_value+" - STP."+method+"() )[~np.isnan(gsw_cv."+comp_value+")].min()")
            print "%s: Passed, but small diff ranging from: %s to %s" % ( method.rjust(width), nmax, nmin)

#load test data
data = pickle.load( open('gsw_cv.pkl','rb') )
gsw_cv = Dict2Struc(data) # then type dat.<tab> to navigate through your variables

""" z_from_p """
z_from_p = gsw.z_from_p(gsw_cv.p_chck_cast, gsw_cv.lat_chck_cast)
test_print("z_from_p")

""" p_from_z """ #NOTE: show diff not present in the original
p_from_z = gsw.p_from_z( z_from_p, gsw_cv.lat_chck_cast )
test_print("p_from_z")

""" grav """
grav = gsw.grav(gsw_cv.lat_chck_cast, gsw_cv.p_chck_cast )
test_print("grav")

""" distance """ #NOTE: show diff not present in the original
distance = gsw.distance(gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast, gsw_cv.p_chck_cast)
test_print("distance")

""" CT_from_pt """
CT_from_pt = gsw.CT_from_pt(SA_chck_cast, pt_chck_cast)
test_print("CT_from_pt")

""" pt_from_CT """
pt_from_CT = gsw.pt_from_CT(SA_chck_cast, CT_chck_cast)
test_print("pt_from_CT", "pt")

""" t_from_CT """
t_from_CT =  gsw.t_from_CT(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
test_print("t_from_CT", "t_chck_cast")

""" pt_from_t """
pt_from_t = gsw.pt_from_t(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast, gsw_cv.pr)
test_print("pt_from_t")

""" pt0_from_t """
pt0_from_t = gsw.pt0_from_t(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("pt0_from_t", "pt0")

""" pt_from_entropy (t_from_entropy) """ #NOTE: show small diff not present in the original
pt_from_entropy =  gsw.t_from_entropy(SA_chck_cast, entropy_chck_cast, 'pt')
test_print("pt_from_entropy")

""" CT_from_entropy (t_from_entropy) """ #NOTE: show diff not present in the original
CT_from_entropy =  gsw.t_from_entropy(SA_chck_cast, entropy_chck_cast, 'CT')
test_print("CT_from_entropy")

""" entropy_from_t (entropy_from_pt) """
entropy_from_pt =  gsw.entropy_from_t(SA_chck_cast, pt_chck_cast, 'pt')
test_print("entropy_from_pt")

""" entropy_from_t (entropy_from_CT) """
entropy_from_CT =  gsw.entropy_from_t(SA_chck_cast, CT_chck_cast, 'CT')
test_print("entropy_from_CT")

""" SA_from_SP """ #TODO: test the second output (in_ocean)
SA_from_SP  = gsw.SA_from_SP(gsw_cv.SP_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0]
test_print("SA_from_SP")

""" SP_from_SA """ #TODO: test the second output (in_ocean)
SA_chck_cast = sio.loadmat("derived_prop.mat", squeeze_me=True)['SA_chck_cast']
SP_from_SA = gsw.SP_from_SA(SA_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0]
test_print("SP_from_SA") #NOTE: the original make the comparison with SP_chck_cast, why?

""" SA_from_rho """
rho_chck_cast = sio.loadmat("derived_prop.mat", squeeze_me=True)['rho']
SA_from_rho = gsw.SA_from_rho(rho_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("SA_from_rho")

""" SA_from_Sstar """
Sstar_from_SA = gsw.Sstar_from_SA(SA_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0]
test_print("Sstar_from_SA")

""" Sstar_from_SA """
#FIXME: why SA_chck_cast instead of Sstar?
SA_from_Sstar = gsw.SA_from_Sstar(SA_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0]
test_print("SA_from_Sstar")

""" SP_from_Sstar """
#FIXME: original run SP_from_SA instead of SP_from_Sstar
SP_from_Sstar = gsw.SP_from_SA(SA_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0]
test_print("SP_from_Sstar")

""" Sstar_from_SP """
Sstar_from_SP = gsw.Sstar_from_SP(gsw_cv.SP_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0]
test_print("Sstar_from_SP")

""" cndr_from_SP """ #NOTE: show diff not present in the original
cndr = gsw.cndr_from_SP(gsw_cv.SP_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("cndr")

""" SP_from_cndr """
cndr = sio.loadmat("derived_prop.mat", squeeze_me=True)['cndr']
SP_from_cndr = gsw.SP_from_cndr(cndr, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("SP_from_cndr")

""" SA_Sstar_from_SP """
SA_SA_Sstar_from_SP, Sstar_SA_Sstar_from_SP = gsw.SA_Sstar_from_SP(gsw_cv.SP_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0:2]
test_print("SA_SA_Sstar_from_SP")
test_print("Sstar_SA_Sstar_from_SP")