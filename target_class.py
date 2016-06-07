####
# (DEPRECATED - Now I just call txt table files for target/arm info) originally, I've been just looking at the txt tables for information on the targets and arms and manually writing if statements in scripts on a target by target basis. Very inefficient and time consuming. So I've decided to just make a target class that contains all the information of all target for every cluster and observation block so that I can call it when I need to retrieve specific information on a target instead of writing if statements in every script.
####

class target:
    def __init__(self, cluster, OB,  name, arm, origarm, RA, DEC, redshift, redshiftquality, Zmag, Irac_chan_1_mag, z-irac_colour, z-irac_colour_error, Radius, MIPS_flux, D4000, M_L, log_10_Mstellar, EW_OII, SFR_OII, SFR_OII_uncorrected, EW_HDL):
        self.cluster = cluster
        self.OB = OB
        self.name = name
        self.arm = arm
        self.origarm = origarm
        self.RA = RA
        self.DEC = DEC
        self.rs = redshift
        self.rsq = redshiftquality
        self.Zmag = Zmag
        self.Irac = Irac_chan_1_mag
        self.Z_I = z-irac_colour
        self.Z_I_err = z-irac_colour_error
        self.R = Radius
        self.MIPS = MIPS_flux
        self.D4000 = D4000
        self.M_L = M_L
        self.log_10_Mstellar = log_10_Mstellar
        self.EW_OII = EW_OII
        self.SFR_OII = SFR_OII
        self.SFR_OII_uncorr = SFR_OII_uncorrected
        self.EW_HDL = EW_HDL        

##CL0034, OB1

CL0034_YJ_OB1_target_4 = target("CL0034", 1, 4, 1, 3, 8.66301, -43.0885, 0.854, 2.0, 21.71, 17.98, 3.97, 0.91, 1.2749, -1.0, 1.34837, 0.9515, 10.11589, 14.1294, 7.6854, 7.8249, -99.0)
CL0034_YJ_OB1_target_8 = target("CL0034", 1, 8, 3, 5, 8.66185, -43.1057, 0.868, 2.0, 22.71, 19.27, 3.62, 0.83, 0.8008, -1.0, 1.32452, 0.9526, 9.61778, 46.5875, 4.1742, 9.2399, -99.0)
#CL0034_YJ_OB1_target_3 = target("CL0034, 1, 
