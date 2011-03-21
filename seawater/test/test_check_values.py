# -*- coding: utf-8 -*-



import unittest
import numpy as np
import numpy.testing as npt
import seawater.gibbs as gsw
import seawater.library as gswl
import seawater.gibbs25 as gsw25

# Standard values for arguments from 
# http://www.teos-10.org/pubs/gsw/html/gsw_contents.html

# Salinities
SP    = [34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600]
SA    = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
Sstar = [34.7115, 34.8912, 35.0247, 34.8436, 34.7291, 34.7197]

# Temperatures
t     = [28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036]
pt    = [28.7832, 28.4209, 22.7850, 10.2305,  6.8292,  4.3245]
CT    = [28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236]
t48   = [29,      28,      23,      10,       7,       4]
t68   = [29,      28,      23,      10,       7,       4]

# Other
p       = [  10,       50,      125,      250,      600,    1000]
z       = [  10,       50,      125,      250,      600,    1000]
rho     = [1021.839, 1022.262, 1024.426, 1027.792, 1029.839, 1032.002]
entropy = [ 400.3892, 395.4378, 319.8668, 146.7910,  98.6471,  62.7919]
lon = 188
lat = 4

# Check values named out_check comes from the page
# http://www.teos-10.org/pubs/gsw/html/gsw_contents.html
#
# If check values above are missing or given with too
# low precision, check values named out_octave are
# calculated independently by octave using
# the GSW toolbox version 2.0 available from
# http://www.teos-10.org/pubs/gsw/html/gsw_contents.html

class Test_standard(unittest.TestCase):
#class Test_standard(npt.TestCase):

    # -----------------------------------------------------
    # Absolute Salinity (SA) and Preformed Salinity (Sstar)
    # -----------------------------------------------------
    
    def test_SA_from_SP(self):
        """Absolute Salinity from Practical Salinity"""
        out_check = np.array((34.711779712893147,
                              34.891523721622328,
                              35.025547737221643,
                              34.847230081708567,
                              34.736629598766619,
                              34.732361864645107))
        out = gsw.SA_from_SP(SP, p, lon, lat)
        npt.assert_array_equal(out, out_check)

    def test_Sstar_from_SP(self):
        """Preformed Salinity from Practical Salinity"""
        out_check = np.array((34.711553202053111,
                              34.891161009146472,
                              35.024649258886718,
                              34.843592772087710,
                              34.729033602488833,
                              34.719676382802788))
        out = gsw.Sstar_from_SP(SP, p, lon,lat)
        npt.assert_array_equal(out, out_check)

    #def test_SA_Sstar_from_SP 
        """Absolute Salinity & Preformed Salinity from Practical Salinity"""
        

    # -----------------------------
    # Conservative Temperature (CT)
    # -----------------------------

    
    def test_CT_from_t(self):
        """Conservative Temperature from in-situ temperature"""
        out = gsw.CT_from_t(SA, t, p)
        out_check = np.array((28.809919826700281,
                              28.439227816091140,
                              22.786176893078498,
                              10.226189266620782,
                               6.827213633479988,
                               4.323575748610455))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)

        
    # -----------------------------------------------------------------------
    # Other conversions between temperatures, salinities, pressure and height
    # -----------------------------------------------------------------------
    
    def test_t_from_CT(self):           # 
        """in-situ temperature from Conservative Temperature"""
        out = gsw.t_from_CT(SA,CT,p)
        out_check = np.array((28.785580227725703,
                              28.432872246163946,
                              22.810323087627076,
                              10.260010752788906,
                               6.886286301029376,
                               4.403624452383043))
        npt.assert_array_equal(out, out_check)
    
    def test_pt_from_t(self):
        """potential temperature"""
        pr = 0
        out = gsw.pt_from_t(SA, t, p, pr)
        out_check = np.array((28.783196819670632,
                              28.420983342398962,
                              22.784930399117108,
                              10.230523661095731,
                               6.829230224409661,
                               4.324510571845719))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)
    
    def test_CT_from_pt(self):
        """Conservative Temperature from potential temperature"""
        out = gsw.CT_from_pt(SA, pt)
        out_check = np.array((28.809923015982083,
                              28.439144260767169,
                              22.786246608464264,
                              10.226165605435785,
                               6.827183417643142,
                               4.323565182322069))     
        npt.assert_array_equal(out, out_check)    

    def test_pot_enthalpy_from_pt(self):
        """potential enthalpy from potential temperature"""
        out = gsw.pot_enthalpy_from_pt(SA,pt)
        out_check = np.array((115005.4085345822,
                              113525.3087024591,
                               90959.6876993543,
                               40821.5028045380,
                               27253.2147222681,
                               17259.1013118296))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=10)
    
    def test_pt0_from_t(self):
        """potential temperature with a reference pressure of zero dbar"""
        out = gsw.pt0_from_t(SA, t, p)
        out_check = np.array((28.783196819670632,
                              28.420983342398962,
                              22.784930399117108,
                              10.230523661095731,
                               6.829230224409661,
                               4.324510571845719))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)
        
    def test_pt_from_CT(self):
        """potential temperature from Conservative Temperature"""
        out = gsw.pt_from_CT(SA, CT)
        out_check = np.array((28.783177048624573,
                              28.420955597191984,
                              22.784953468087107,
                              10.230534394434429,
                               6.829216587061605,
                               4.324534835990236))
        npt.assert_array_equal(out, out_check)        
    def test_SP_from_SA(self):
        """Practical Salinity from Absolute Salinity"""
        out = gsw.SP_from_SA(SA, p, lon, lat)
        out_check = np.array((34.548720191893423,
                              34.727476389710311,
                              34.860552017493589,
                              34.680970059473843,
                              34.567970540149211,
                              34.560037956374323))
        npt.assert_array_equal(out, out_check)        

    def test_Sstar_from_SA(self):
        """Preformed Salinity from Absolute Salinity"""
        out = gsw.Sstar_from_SA(SA, p, lon, lat)
        out_check = np.array((34.711573489159953,
                              34.891137287524145,
                              35.024701521665072,
                              34.843562690379144,
                              34.729004003722217,
                              34.719714518157680))
        npt.assert_array_equal(out, out_check)
        
    def test_SA_from_Sstar(self):
        """Absolute Salinity from Preformed Salinity"""
        out = gsw.SA_from_Sstar(Sstar, p, lon, lat)
        out_check = np.array([34.711726510840045,
                             34.891562712475853,
                             35.025598478334928,
                             34.847237309620859,
                             34.736695996277788,
                             34.732385481842321])
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=10)

    def test_SP_from_Sstar(self):
        """Practical Salinity from Preformed Salinity"""
        out = gsw.SP_from_Sstar(Sstar, p, lon, lat)
        out_check = np.array([34.548647047639960,
                              34.727538807857847,
                              34.860550502970142,
                              34.681007193989544,
                              34.568066085887892,
                              34.560023506354682])
        npt.assert_array_equal(out, out_check)    

    def test_z_from_p(self):
        """height from pressure"""
        out = gsw.z_from_p(p, lat)
        out_check = np.array((  -9.9446,  -49.7182, -124.2728,
                              -248.4704, -595.8262, -992.0932))
        out_octave = [  -9.94460074378793,
                       -49.71817464749508,
                      -124.27282749612912,
                      -248.47044828162981,
                      -595.82618013969216,
                      -992.09317479591675]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check,  decimal=4)
        npt.assert_array_almost_equal(out, out_octave, decimal=14)
    
    def test_p_from_z(self):
        """pressure from height"""
        out = gsw.p_from_z(z,lat)
        out_check =  [  -10.055217939908,
                        -50.271175099730,
                       -125.654885697343,
                       -251.232845042860,
                       -602.440507522266,
                      -1003.076098067266]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=12)

    def test_t90_from_t48(self):
        """ITS-90 temperature from IPTS-48 temperature"""
        out = gsw.t90_from_t48(t48)
        out_check = np.array((28.9840, 27.9844, 22.9867,
                               9.9936,  6.9955,  3.9974))
        out_octave = np.array((28.98398424378149,
                               27.98441334079821,
                               22.98669079420939,
                                9.99364152603375,
                                6.99545669039431,
                                3.99735103575142))
        #print np.max(abs(out-out_octave))
        npt.assert_array_almost_equal(out, out_check,  decimal=4)
        npt.assert_array_almost_equal(out, out_octave, decimal=14)

    def test_t90_from_t68(self):
        """ITS-90 temperature from IPTS-68 temperature"""
        out = gsw.t90_from_t68(t68)
        out_check = np.array((28.9930, 27.9933, 22.9945,
                               9.9976,  6.9983,  3.9990))
        out_octave = np.array((28.99304166999919,
                               27.99328161241301,
                               22.99448132448212,
                                9.99760057586179,
                                6.99832040310325,
                                3.99904023034472))
        #print np.max(abs(out-out_octave))
        npt.assert_array_almost_equal(out, out_check,  decimal=4)
        npt.assert_array_almost_equal(out, out_octave, decimal=14)

    # -----------------------------------------------------------------
    # density and enthalpy, based on the 25-term expression for density
    # -----------------------------------------------------------------

    def test_rho_CT25(self):
        """in-situ density from 25 term expression"""
        out = gsw25.rho_CT25(SA, CT, p)
        out_check = np.array((1021.839541498949,
                              1022.261845123599,
                              1024.426245497197,
                              1027.792152543827,
                              1029.838876568916,
                              1032.002445797635))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=12)

    # Individual terms are tested
    #def test_rho_alpha_beta_CT25(self):
        """in-situ density, expansison and contraction coefficients"""

    def test_alpha_CT25(self):
        """thermal expansison coefficient by 25 term equation"""
        out = gsw25.alpha_CT25(SA, t, p)
        out_check = 1.0e-3 * np.array((0.324356899200044,
                                       0.322411511462120,
                                       0.281313314947018,
                                       0.173211327216534,
                                       0.146078398825498,
                                       0.129148950914940))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=6)

    def test_beta_CT25(self):
        """saline contraction coefficient by 25 term equation"""
        out = gsw25.beta_CT25(SA, t, p)
        out_check = 1.0e-3 * np.array((0.717340839421976,
                                       0.717514858800028,
                                       0.726299524800543,
                                       0.750622769257820,
                                       0.755000880950473,
                                       0.756865290355974))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=6)

    def test_specvol_CT25(self):
        """specific volume by 25 term equation"""
        out = gsw25.specvol_CT25(SA, t, p)
        out_check = 1.0e-3 * np.array((0.978627229998447,
                                       0.978222952143042,
                                       0.976156169753986,
                                       0.972959364911436,
                                       0.971025684456263,
                                       0.968989951595608))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=7)
        
    def test_specvol_anom_CT25(self):
        """specific volume anomaly by 25 term equation"""
        out = gsw25.specvol_anom_CT25(SA, t, p)
        out_check = 1.0e-5 * np.array((0.600921244780790,
                                       0.578505932118849,
                                       0.405538759163434,
                                       0.141863283479346,
                                       0.104119730440851,
                                       0.076279577097848))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=7)
        
    def test_enthalpy_CT25(self):
        """enthalpy by 25 term formulation"""
        out = gsw25.enthalpy_CT25(SA, CT, p)
        out_check = [115103.1813925039,
                     114014.6929487063,
                      92180.0148514069,
                      43255.3660538842,
                      33087.1524935877,
                      26970.6799089967]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=10)

    def test_enthalpy_diff_CT25(self):
        """difference of enthalpy between two pressures"""

        p_shallow =  [ 10,  50, 125, 250, 600, 1000]
        p_deep    =  [110, 150, 225, 350, 700, 1100]
        out = gsw25.enthalpy_diff_CT25(SA, CT, p_shallow, p_deep)
        out2 = ( gsw25.enthalpy_CT25(SA, CT, p_deep) -
                 gsw25.enthalpy_CT25(SA, CT, p_shallow) )
        out3 = ( gsw.enthalpy_CT25(SA, CT, p_deep) -
                 gsw.enthalpy_CT25(SA, CT, p_shallow) )
        out_check = [978.4262827561797,
                     978.0221873372307,
                     975.9531141438936,
                     972.7477253978700,
                     970.8128690308879,
                     968.7770260315075]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=13)
        npt.assert_array_almost_equal(out2, out_check, decimal=8)
        

    # --------------------------------------------------------------------
    # water column properties, based on the 25-term expression for density
    # --------------------------------------------------------------------
    
    #def test_Nsquared_CT25(self):
        """buoyancy (Brunt-Vaisala) frequency squared (N^2)"""
        
    #def test_Turner_Rsubrho_CT25(self):
        """Turner angle & Rsubrho"""
        
    #def test_IPV_vs_fNsquared_ratio_CT25(self):
        """potential density gradient ratio"""

    # ------------------------------------------------------------------------
    # geostrophic streamfunctions, based on the 25-term expression for density
    # ------------------------------------------------------------------------
    
    #def test_geo_strf_dyn_height(self):
        """dynamic height anomaly"""
        
    #def test_geo_strf_dyn_height_pc(self):
        """dynamic height anomaly for piecewise constant profiles"""

    #def test_geo_strf_McD_Klocker(self):
        """McDougall-Klocker geostrophic streamfunction"""
        
    #def test_geof_str_McD_Klocker_pc(self):
        """McDougall-Klocker for piecewise constant profiles"""
        
    #def test_geo_strf_Montgomery(self):
        """Montgomery geostrophic streamfunction"""
        
    #def test_geo_strf_Cunningham(self):
        """Cunningham geostrophic streamfunction"""
    #def test_geostrophic_velocity(self):
        """geostrophic velocity"""

    # ------------------------------------------------------------------
    # neutral and non-linear properties, based on the 25-term expression
    # ------------------------------------------------------------------
    
    def test_cabbeling_CT25(self):
        """cabbeling coefficient"""
        out = gsw25.cabbeling_CT25(SA, CT, p)
        # Buggy values, computed with buggy matlab routine
        out_check = 1.0e-4 * np.array((0.071255687860750,
                                       0.071469675168983,
                                       0.077959301008454,
                                       0.098539728631524,
                                       0.103657650078995,
                                       0.106537287638893))
        # Computed by bug-fixed octave
        out_octave = np.array((8.62401756238290e-06,
                               8.64247376958721e-06,
                               9.22258732263748e-06,
                               1.09466959666811e-05,
                               1.13365671614248e-05,
                               1.15414936408258e-05))

        #print np.max(abs(out-out_check))
        #print np.max(abs(out-out_octave))
        npt.assert_array_almost_equal(out, out_octave, decimal=14)

        
    def test_thermobaric_CT25(self):
        """thermobaric coefficient"""
        out = gsw25.thermobaric_CT25(SA, CT, p)
        out_check = 1.0e-11 * np.array((0.143712914501030,
                                        0.144525225115917,
                                        0.163502412605500,
                                        0.226474876009325,
                                        0.247208470958881,
                                        0.262712952286233))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=20)

    def test_isopycnal_slope_ratio_CT25(self):
        """ratio of the slopes of isopycnals on the SA-CT diagram"""
        out = gsw25.isopycnal_slope_ratio_CT25(SA, CT, p, 0)
        out_check = [1.000443300473299,
                     1.002247196919216,
                     1.007323843523325,
                     1.033841543679809,
                     1.113488165501534,
                     1.257603349708546]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=14)
        

    #def test_isopycnal_vs_ntp_CT_ratio_CT25(self):
        """ratio of the gradient of Conservative Temperature"""

    #def test_ntp_pt_vs_CT_ratio_CT25(self):
        """ratio of gradients of potential temperature"""

    # -----------------------------------------------------
    # basic thermodynamic properties in terms of (SA, t, p)
    # -----------------------------------------------------
    
    def test_rho(self):
        """in-situ density"""
        out = gsw.rho(SA, t, p)
        out_check = [1021.840173185531,
                     1022.262689926782,
                     1024.427715941676,
                     1027.790201811623,
                     1029.837714725961,
                     1032.002404116447]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=12)

    def test_pot_rho(self):
        """potential density"""
        pr = 0
        out = gsw.pot_rho(SA, t, p, pr)
        out_check = [1021.798145811089,
                     1022.052484416980,
                     1023.893583651958,
                     1026.667621124443,
                     1027.107230868492,
                     1027.409631264134]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=12)


    def test_specvol(self):
        """specific volume"""
        out = gsw.specvol(SA, t, p)
        out_check = np.array((0.000978626625025472,
                              0.000978222143734527,
                              0.000976154768597586,
                              0.000972961211575438,
                              0.000971026779948624,
                              0.000968989990731808))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=18)

    def test_specvol_anom(self):
        """specific volume anomaly"""
        out = gsw.specvol_anom(SA,t,p)
        out_check = 1.0e-5 * np.array((0.601044462887525,
                                       0.578602431524256,
                                       0.405564998680916,
                                       0.142198661956343,
                                       0.104351837470899,
                                       0.076396485029176))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=20)

    def test_alpha_wrt_CT(self):
        """thermal expansion coefficient wrt Conservative Temperature"""
        out = gsw.alpha_wrt_CT(SA, t, p)
        out_check = 1.0e-3 * np.array((0.324707942226035,
                                       0.322724175251768,
                                       0.281178901333399,
                                       0.173138222139544,
                                       0.146269210200123,
                                       0.129426849037398))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=18)

    def test_alpha_wrt_pt(self):
        """thermal expansion coefficient wrt potential temperature"""

        out = gsw.alpha_wrt_pt(SA, t, p)
        out_check = 1.0e-3 * np.array((0.325621974596088,
                                       0.323548679831175,
                                       0.281641477318055,
                                       0.173138875637364,
                                       0.146227722831487,
                                       0.129358812895628))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=18)
        
    def test_alpha_wrt_t(self):
        """thermal expansion coefficient wrt in-situ temperature"""
        out = gsw.alpha_wrt_t(SA, t, p)
        out_check = 1.0e-3 * np.array((0.325601747227247,
                                       0.323448083851267,
                                       0.281413883319329,
                                       0.172825692975230,
                                       0.145569941503599,
                                       0.128362986933288))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=18)
        
    def test_beta_const_CT(self):
        """saline contraction coefficient at constant Conservative Temperature"""
        out = gsw.beta_const_CT(SA, t, p)
        out_check = 1.0e-3 * np.array((0.717486774409917,
                                       0.717647567079959,
                                       0.726219047211703,
                                       0.750509250912251,
                                       0.755062962544804,
                                       0.757065895554659))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=18)
        
    def test_beta_const_pt(self):
        """saline contraction coefficient at constant potential temperature"""
        out = gsw.beta_const_pt(SA, t, p)
        out_check = 1.0e-3 * np.array((0.731118367454125,
                                       0.731059415831298,
                                       0.735986955829449,
                                       0.753748863626444,
                                       0.757121845539448,
                                       0.758434161456822))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=18)

    def test_beta_const_t(self):
        """saline contraction coefficient at constant in-situ temperature"""
        out = gsw.beta_const_t(SA, t, p)
        out_check = 1.0e-3 * np.array((0.731120837010429,
                                       0.731071779078011,
                                       0.736019128913071,
                                       0.753810501711847,
                                       0.757259405338257,
                                       0.758649268096996))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=18)

    def test_entropy(self):
        """entropy"""
        out = gsw.entropy(SA, t, p)
        out_check = [400.3894252787245,
                     395.4381784340642,
                     319.8664981986740,
                     146.7908815899072,
                      98.6473408657975,
                      62.7915087346090]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=13)

    def test_internal_energy(self):
        """internal energy"""
        out = gsw.internal_energy(SA,t,p)
        out_check = 1.0e5 * np.array((1.149062384730930,
                                      1.134265741706209,
                                      0.908608185884215,
                                      0.407243400571882,
                                      0.271626660018474,
                                      0.171825052266728))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=10)
        
    def test_enthalpy(self):
        """enthalpy"""
        out = gsw.enthalpy(SA, t, p)
        out_check = 1.0e5 * np.array((1.151032604783763,
                                      1.140148036012021,
                                      0.921799209310966,
                                      0.432553283808897,
                                      0.330872159700175,
                                      0.269705880448018))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=10)

    def test_cp(self):
        """isobaric heat capacity"""
        out = gsw.cp(SA, t, p)
        out_check = [4002.888003958537,
                     4000.980283927373,
                     3995.546468894633,
                     3985.076769021370,
                     3973.593843482723,
                     3960.184084786622]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=12)

    def test_isochoric_heat_cap(self):
        """isochoric heat capacity"""
        out = gsw.isochoric_heat_cap(SA, t, p)
        out_check = [3928.137087019246,
                     3927.273816327282,
                     3941.364185254816,
                     3966.261261456055,
                     3960.509032220448,
                     3950.139013423850]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=12)
        
    def test_chem_potential_relative(self):
        """relative chemical potential"""
        out = gsw.chem_potential_relative(SA, t, p)
        out_check = np.array((79.425448103352224,
                              79.259892142259332,
                              74.691548585676486,
                              65.640637188627224,
                              61.226856561881860,
                              57.212985572773491))
        npt.assert_array_equal(out, out_check)        
    def test_chem_potential_water(self):
        """chemical potential of water in seawater"""
        out = gsw.chem_potential_water(SA, t, p)
        out_check = np.array((-8545.561146284534,
                              -8008.085548342105,
                              -5103.980139874876,
                               -634.067782745442,
                               3335.566803473286,
                               7555.434445971858))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=12)
        
    def test_chem_potential_salt(self):
        """chemical potential of salt in seawater"""
        out  = gsw.chem_potential_salt(SA, t, p)
        out_check = np.array((-8466.135698181180,
                              -7928.825656199846,
                              -5029.288591289200,
                               -568.427145556815,
                               3396.793660035168,
                               7612.647431544631))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=12)
        
    def test_Helmholtz_energy(self):
        """Helmholtz energy"""
        out = gsw.Helmholtz_energy(SA, t, p)
        #out = gsw.helmholtz_energy(SA, t, p)
        out_check = np.array((-5985.582882093846,
                              -5830.818452241629,
                              -3806.966178407540,
                               -877.663694207387,
                               -462.170339049288,
                               -245.504072049327))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=12)

    def test_sound_speed(self):
        """sound speed"""
        out = gsw.sound_speed(SA, t, p)
        out_check = np.array((1542.615803587414,
                              1542.703534065789,
                              1530.844979136360,
                              1494.409996920661,
                              1487.377102518027,
                              1483.934609078705))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=12)

    def test_kappa(self):
        """isentropic compressibility"""
        out = gsw.kappa(SA, t, p)
        out_check = 1.0e-9 * np.array((0.411245799180373,
                                       0.411029072229334,
                                       0.416539558054756,
                                       0.435668337689072,
                                       0.438923693006423,
                                       0.440037575765429))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=24)

    def test_kappa_const_t(self):
        """isothermal compressibility"""
        out = gsw.kappa_const_t(SA, t, p)
        out_check = 1.0e-9 * np.array((0.419071646368281,
                                       0.418743202287955,
                                       0.422265764368337,
                                       0.437735100406753,
                                       0.440373818138017,
                                       0.441156576599537))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=24)

    def test_adiabatic_lapse_rate(self):
        """adiabatic lapse rate"""
        out = gsw.adiabatic_lapse_rate(SA, t, p)
        out_check = 1.0e-7 * np.array((0.240350282348014,
                                       0.238496699896003,
                                       0.203479879742929,
                                       0.119586543071262,
                                       0.099617071808651,
                                       0.087174726986510))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=22)
        
    def test_molality(self):
        """molality of seawater"""
        out = gsw.molality(SA)
        out_check = np.array((1.145084756351415,
                              1.151227076302664,
                              1.155812234675510,
                              1.149712646992389,
                              1.145932308325875,
                              1.145788768234685))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)

    def test_ionic_strength(self):
        """ionic strength of seawater"""
        out = gsw.ionic_strength(SA)
        out_check = np.array((0.712981183609951,
                              0.716805667801764,
                              0.719660593278309,
                              0.715862716115311,
                              0.713508907524334,
                              0.713419533018609))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)
        
    def test_osmotic_coefficient(self):
        """osmotic coefficient of seawater"""
        out = gsw.osmotic_coefficient(SA, t, p)
        out_check = np.array((0.902847177948536,
                              0.902986238880573,
                              0.902388659002806,
                              0.898809266717742,
                              0.898010536533592,
                              0.897679115176392))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)

    # Individual functions are tested below
    #def test_temps_maxdensity(self):
    
    def test_t_maxdensity(self):
        """in-situ temperature of maximum density of seawater"""
        out = gsw.t_maxdensity(SA, p)
        out_check =  [-3.725008913916380,
                      -3.853714292562061,
                      -4.051433177642653,
                      -4.295542513200469,
                      -5.074746618871896,
                      -6.011701971231014]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)

    def test_pt_maxdensity(self):
        """potential temperature of maximum density of seawater"""
        out = gsw.pt_maxdensity(SA, p)
        out_check = [-3.725007817485227,
                     -3.853686822871752,
                     -4.051260804552126,
                     -4.294848434851046,
                     -5.070673022022400,
                     -6.000138120066936]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=32)

    def test_CT_maxdensity(self):
        """Conservative Temperature of maximum density of seawater"""
        out = gsw.CT_maxdensity(SA, p)
        out_check = [-3.720779524452287,
                     -3.849111898670119,
                     -4.046292804691113,
                     -4.290063909787930,
                     -5.065854290972595,
                     -5.995131932081709]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)

        
 
    # ----------------------------------------------------
    # basic thermodynamic properties in terms of CT and pt
    # ----------------------------------------------------

    def test_rho_CT(self):
        """in-situ density from CT"""
        out = gsw.rho_CT(SA, CT, p)
        out_check = np.array((1021.840179764021,
                              1022.262699103554,
                              1024.427709285783,
                              1027.790199901620,
                              1029.837716779620,
                              1032.002400877215))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=12)

    # Component functions are tested individually
    #def test_rho_alpha_beta_CT(self):
        """in-situ density, expansion & contraction coefficients from CT"""

    def test_alpha_CT(self):
        "thermal expansion coefficient from CT"""    
        out = gsw.alpha_CT(SA, CT, p)
        out_check = 1.0e-3 * np.array((0.324707799687110,
                                       0.322723974596565,
                                       0.281179083137791,
                                       0.173138326981408,
                                       0.146269069937061,
                                       0.129427106726016))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=18)

    def test_beta_CT(self):
        """saline contraction coefficient from CT"""
        out = gsw.beta_CT(SA, CT, p)
        out_check = 1.0e-3 * np.array((0.717486805409467,
                                       0.717647610706394,
                                       0.726219007035780,
                                       0.750509225553259,
                                       0.755062996958255,
                                       0.757065832066916))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=18)
                                      
    def test_specvol_CT(self):
        """specific volume from CT"""
        out = gsw.specvol_CT(SA, CT, p)
        out_check = 1.0e-3 * np.array((0.978626618725186,
                                       0.978222134953103,
                                       0.976154774939840,
                                       0.972961213383549,
                                       0.971026778012244,
                                       0.968989993773259))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=18)
        

    #def test_specvol_anom_CT(self):
        """specific volume anomaly from CT"""

    def test_sigma0_CT(self):
        """sigma0 in terms of SA & CT with reference pressure of 0 dbar"""
        out = gsw.sigma0_CT(SA, CT)
        out_check = np.array((21.798152387272808,
                              22.052493577393079,
                              23.893577037223736,
                              26.667619272292086,
                              27.107232725897575,
                              27.409628653483423))
        #print np.max(abs(out-out_check))
        npt.assert_array_equal(out, out_check)

    def test_sigma1_CT(self):
        """sigma1 in terms of SA & CT with reference pressure of 1000 dbar"""
        out = gsw.sigma1_CT(SA, CT)
        out_check = np.array((25.955551924297879,
                              26.213079338959005,
                              28.125625321852112,
                              31.120452618127047,
                              31.637685737359789,
                              32.002400877214541))
        #print np.max(abs(out-out_check))
        npt.assert_array_equal(out, out_check)

    def test_sigma2_CT(self):
        """sigma2 in terms of SA & CT with reference pressure of 2000 dbar"""
        out = gsw.sigma2_CT(SA, CT)
        out_check = np.array((30.023190232854176,
                              30.283832597402807,
                              32.265717400800895,
                              35.474625433103029,
                              36.067249427311708,
                              36.492547357788226))
        #print np.max(abs(out-out_check))
        npt.assert_array_equal(out, out_check)

    def test_sigma3_CT(self):
        """sigma3 in terms of SA & CT with reference pressure of 3000 dbar"""
        out = gsw.sigma3_CT(SA, CT)
        out_check = np.array((34.003709330818083,
                              34.267404341412885,
                              36.316609174494715,
                              39.732408885730138,
                              40.397870256006627,
                              40.881719119051240))
        #print np.max(abs(out-out_check))
        npt.assert_array_equal(out, out_check)

    def test_sigma4_CT(self):
        """sigma4 in terms of SA & CT with reference pressure of 4000 dbar"""
        out = gsw.sigma4_CT(SA, CT)
        out_check = np.array((37.899615932915594,
                              38.166310327096426,
                              40.280942547164841,
                              43.896123296674659,
                              44.631612266276306,
                              45.171743792926918))
        #print np.max(abs(out-out_check))
        npt.assert_array_equal(out, out_check)

    def test_enthalpy_CT(self):
        """enthalpy from CT"""
        out = gsw.enthalpy_CT(SA, CT, p)
        out_check = 1.0e+5 * np.array((1.151031813321767,
                                       1.140146925586514,
                                       0.921800131787836,
                                       0.432553712315790,
                                       0.330871615358722,
                                       0.269706848807403))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=10)
        

    # Not implemented yet
    def rest_enthalpy_diff_CT(self):
        """difference of enthalpy from CT between two pressures"""
        p_shallow = np.array(p)
        p_deep    = p_shallow + 100.0
        out = gsw.enthalpy_diff_CT(SA, CT, p_shallow, p_deep)
        out_check = np.array((978.4255933006788,
                              978.0212984387035,
                              975.9516776633390,
                              972.7494888619550,
                              970.8138965635516,
                              968.7770207718350))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=10)

    def test_entropy_from_pt(self):
        """entropy from potential temperature"""
        pt = [28.7832, 28.4210, 22.7850, 10.2305, 6.8292, 4.3245]
        out = gsw.entropy_from_pt(SA, pt)
        out_check = np.array([400.38946744, 395.43839949, 319.86743859,
                              146.79054828,  98.64691006,  62.79135672])
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=8)

    def test_entropy_from_CT(self):
        """entropy from Conservative Temperature"""
        out = gsw.entropy_from_CT(SA, CT)
        out_check = np.array([400.3892,  395.4378, 319.8668,
                              146.7910,   98.6471,  62.7919])
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=4)

    def test_pt_from_entropy(self):
        """potential temperature from entropy"""
        out = gsw.pt_from_entropy(SA, entropy)
        out_check = np.array((28.7832, 28.4210, 22.7850,
                              10.2305,  6.8292,  4.3245))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=4)

    def test_CT_from_entropy(self):
        """Conservative Temperature from entropy"""
        out = gsw.CT_from_entropy(SA, entropy)
        out_check = np.array((28.8099, 28.4392, 22.7862,
                              10.2262,  6.8272,  4.3236))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=5)

    # -------------------------------------------
    # derivatives of enthalpy, entropy, CT and pt
    # -------------------------------------------

    # Separate derivatives tested below
    #def test_CT_first_derivatives(self):
        """first derivatives of Conservative Temperature"""

    def test_CT_derivative_SA(self):
        """derivative of Conservative Temperature wrt absolute salinity"""
        out = gsw.CT_derivative_SA(SA, pt)
        out_check = np.array((-0.041981092877806,
                              -0.041558140199508,
                              -0.034739209004865,
                              -0.018711103772892,
                              -0.014075941811725,
                              -0.010571716552295))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)
    
    def test_CT_derivative_pt(self):
        """derivative of Conservative Temperature wrt potential temperature"""
        out = gsw.CT_derivative_pt(SA, pt)
        out_check = np.array((1.002814937296636,
                              1.002554817053239,
                              1.001645140295163,
                              1.000003771100520,
                              0.999716359504731,
                              0.999474326580093))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)
    
    #def test_CT_second_derivatives(self):
        """second derivatives of Conservative Temperature"""

    # Not implemented 2011-03-19
    def rest_CT_derivative_SA_SA(self):
        """second derivative of Cons. Temperature wrt Abs. Sallinity"""
        out = gsw.CT_derivative_SA_SA(SA, pt)
        out_check = 1.0e-3 * np.array((-0.060718502073595,
                                       -0.062065324397404,
                                       -0.084017055351272,
                                       -0.148436050118396,
                                       -0.171270386500246,
                                       -0.189920754897514))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)

    # Not implemented 2011-03-19
    def rest_CT_derivative_SA_pt(self):
        """mixed second derivative of Conservative Temperature"""
        out = gsw.CT_derivative_CT_SA(SA, pt)
        out_check = [-0.001197415000868,
                     -0.001198309530139,
                     -0.001226523296082,
                     -0.001335896286480,
                     -0.001380492698572,
                     -0.001417751669135]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)

    # Not implemented 2011-03-19
    def rest_CT_derivative_pt_pt(self):
        """second derivative of Cons. Temperature wrt pot. temperature"""
        out = gsw.CT_derivative_CT_SA(SA, pt)
        out_check = 1.0e-3 * np.array((0.123012754427146,
                                       0.124662008871271,
                                       0.140829458783443,
                                       0.140646803448166,
                                       0.113684095615077,
                                       0.082286843477998))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)


    #def test_enthalpy_first_derivatives(self):
        """first derivatives of enthalpy"""
    # Not implemented 2011-03-19
    def rest_enthalpy_derivative_SA(self):
        """derivative of enthalpy wrt Absolute Salinity"""
        out = gsw.enthalpy_derivative_SA(SA, CT, p)
        out_check = [-0.070221017857321,
                     -0.351155742717708,
                     -0.887081852469166,
                     -1.829811186103484,
                     -4.424369449732652,
                     -7.407229217110896]
        print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)

    # Not implemented 2011-03-19
    def rest_enthalpy_derivative_CT(self):
        """derivative of enthalpy wrt Conservative Temperature"""
        out = gsw.enthalpy_derivative_SA(SA, CT, p)
        out_check = [3991.899729625529,
                     3992.025696743769,
                     3992.210168035032,
                     3992.283178760239,
                     3992.681641668387,
                     3993.005774345293]
        print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)

    # Not implemented 2011-03-19
    def rest_enthalpy_derivative_p(self):
        """derivative of enthalpy wrt pressure"""
        out = gsw.enthalpy_derivative_SA(SA, CT, p)
        out_check = 1.0e-3 * np.array((0.978626618725186,
                                       0.978222134953103,
                                       0.976154774939840,
                                       0.972961213383549,
                                       0.971026778012244,
                                       0.968989993773259))
        print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)


    #def test_enthalpy_second_derivatives(self):
        """second derivatives of enthalpy"""



    # Tested individually below
    #def test_entropy_first_derivatives(self):
        """first derivatives of entropy"""

    def test_entropy_derivative_CT(self):
        """derivative of entropy wrt to Conservative Temperature"""
        out = gsw.entropy_derivative_CT(SA, CT)
        out_check = np.array((13.221031210083824,
                              13.236911191313675,
                              13.489004628681361,
                              14.086599016583795,
                              14.257729576432077,
                              14.386429945649411))
        npt.assert_array_equal(out, out_check)

    def test_entropy_derivative_SA(self):
        """derivative of entropy wrt to Absolute Salinity"""
        out = gsw.entropy_derivative_SA(SA, CT)
        out_check = np.array((-0.263286800711655,
                              -0.263977276574528,
                              -0.255367497912925,
                              -0.238066586439561,
                              -0.234438260606436,
                              -0.232820684341694))
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)

    #Tested indivdually below
    #def test_entropy_second_derivatives(self):
        """second derivatives of entropy"""

    def test_entropy_derivative_SA_SA(self):
        out = gsw.entropy_derivative_SA_SA(SA, CT)
        out_check = [-0.007627718929669,
                     -0.007591969960708,
                     -0.007528186784540,
                     -0.007455177590576,
                     -0.007441108287466,
                     -0.007414368396280]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)

    def test_entropy_derivative_SA_CT(self):
        out = gsw.entropy_derivative_SA_CT(SA, CT)
        out_check = [-0.001833104216751,
                     -0.001819473824306,
                     -0.001580843823414,
                     -0.000930111408561,
                     -0.000717011215195,
                     -0.000548410546830]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)

    def test_entropy_derivative_CT_CT(self):
        out = gsw.entropy_derivative_CT_CT(SA, CT)
        out_check = [-0.043665023731109,
                     -0.043781336189326,
                     -0.045506114440888,
                     -0.049708939454018,
                     -0.050938690879443,
                     -0.051875017843472]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)
        
    #def test_pt_first_derivatives(self):
        """first derivatives of potential temperature"""

    # Not implemented 2011-03-19
    def rest_pt_derivative_SA(self):
        """derivative of pot. temperature wrt Abs. Salinity"""
        out = gsw.pt_derivative_SA(SA, CT)
        out_check = [0.041863223165431,
                     0.041452303483011,
                     0.034682095247246,
                     0.018711079068408,
                     0.014079958329844,
                     0.010577326129948]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)

    # Not implemented 2011-03-10    
    def rest_pt_derivative_SA(self):
        """derivative of pot. temperature wrt. Cons. Temperature"""
        out = gsw.pt_derivative_CT(SA, CT)
        out_check = [0.997192967140242,
                     0.997451686508335,
                     0.998357568277750,
                     0.999996224076267,
                     1.000283719083268,
                     1.000525947028218]
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=15)
        

    #def test_pt_second_derivatives(self):
        """second derivatives of potential temperature"""

    # -----------------------
    # Planet Earth properties
    # -----------------------

    # Not implemented yet
    def rest_f(self):
        """Coriolis parameter (f)"""
        lat = [-90, -60, -30, 0]
        out = gsw.f(lat)
        out_check = np.array([-0.000145842300000000,
                              -0.000126303136746351,
                              -0.000072921150000000,
                               0.0])
        #print np.max(abs(out-out_check))
        npt.assert_array_almost_equal(out, out_check, decimal=8)

    def test_grav(self):
        """Gravitational acceleration"""
        lat = [-90, -60, -30, 0]
        p = 0
        out0 = gsw.grav(lat)
        out1 = gsw.grav(lat, p)
        out_check = np.array((9.832186205884799,
                              9.819178859991149,
                              9.793249257048750,
                              9.780327000000000))
        npt.assert_array_equal(out0, out_check)
        npt.assert_array_equal(out1, out_check)
        
    # ERROR in check values
    def test_distance(self):
        """spherical earth distance at a given pressure"""
        lon = [159, 220]
        lat = [-35,  35]
        out = gsw.distance(lon, lat)
        out_check  = 1.050955616404136e+007
        out_octave = 10030974.6529160
        #print abs(out-out_octave)
        npt.assert_almost_equal(out, out_octave, decimal=8)
        
        lon = [159,  220]
        lat = [-35,   35]
        p   = [200, 1000]
        out = gsw.distance(lon, lat, p)
        out_check  = 1.050922853324439e+007
        out_octave = 10030661.6387801
        #print abs(out-out_octave)
        npt.assert_almost_equal(out, out_octave, decimal=8)

    def test_SA_from_SP_Baltic(self):
        """Absolute Salinity in the Baltic Sea"""
        SP = [6.5683, 6.6719, 6.8108, 7.2629, 7.4825, 10.2796]
        lon, lat = 20, 59

        out = gswl._SA_from_SP_Baltic(SP,lon,lat)
        out_check = np.array([6.6699,  6.7738,
                              6.9130,  7.3661,
                              7.5862, 10.3895])
        out_octave = np.array(( 6.66994543234286,
                                6.77377643074286,
                                6.91298613805714,
                                7.36609419188571,
                                7.58618383714286,
                               10.38952057097143))
        npt.assert_array_almost_equal(out, out_check, 4)
        npt.assert_array_almost_equal(out, out_octave, 14)
        
    def test_SP_from_SA_Baltic(self):
        """Practical Salinity in the Baltic Sea"""
        SA = [6.6699, 6.7738, 6.9130, 7.3661, 7.5862, 10.3895]
        lon, lat = 20, 59
        out = gswl._SP_from_SA_Baltic(SA, lon, lat)
        out_check  = np.array([6.5683,  6.6719,
                               6.8108,  7.2629,
                               7.4825, 10.2796])
        out_octave = np.array(( 6.56825466873292,
                                6.67192351682135,
                                6.81081383110345,
                                7.26290579519266,
                                7.48251612689877,
                               10.27957947479392))
        npt.assert_array_almost_equal(out, out_check, decimal=4)
        npt.assert_array_almost_equal(out, out_octave, decimal=14)
        
    # No check values on the web site
    #def test_infunnel(self):
        """'oceanographic funnel' check for the 25-term equation"""

    # No check values on the web site
    #def test_entropy_part(self):
        """entropy minus the terms that are a function of only SA"""

    # No check values on the web site
    #def test_entropy_part_zerop(self):
        """entropy_part evaluated at 0 dbar"""
    
    # No check values on the web site
    #def test_interp_McD_Klocker(self):
        """linearly interpolates the reference cast"""

    # No check values on the web site
    #def test_interp_SA_CT(self):
        """linearly interpolates (SA,CT,p) to the desired p"""

    # No check values on the web site
    #def test_gibbs_pt0_pt0(self):
        """gibbs(0,2,0,SA,t,0)"""

    # No check values on the web site
    #def test_specvol_SSO_0_CT25(self):
        """spec_vol_CT25(35.16504,0,p)"""

    # No check values on the web site
    #def test_enthalpy_SSO_0_CT25(self):
        """enthalpy_CT25(35.16504,0,p)"""

# -----------------------------------------------

if __name__ == '__main__':
    # A more verbose output
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_standard)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
    #npt.test()

