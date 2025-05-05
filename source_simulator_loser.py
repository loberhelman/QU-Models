"""
There are three classes here the first two are parents of the final class. The two parent classes are called 
Simulated_Fg_Component and Simulated_Src_Component and simulate different foreground and source models respectively. 
The final class is called Simulated_Components and is a child of the two classes and will work to combine foreground 
and source models to create a complete line of sight. 

WHAT ABOUT THE BEAM MODELS?


"""

import numpy as np


class Simulated_Fg_Component:

    def __init__(self, lambda_sqr_array, RM_radm2, delta_rm_2, sigmaRM_radm2, psi0_deg):
        #Sets the source parameters
        #FIXME WHAT ARE WE DOING ABOUT RM? RM Center RM 

        self.lambda_sqr_array = lambda_sqr_array
        self.RM_radm2 = RM_radm2
        self.delta_rm_2 = delta_rm_2
        self.sigmaRM_radm2 = sigmaRM_radm2
        self.psi0_deg = psi0_deg
        self.quarr_fg = 0


    def set_a_linear_gradient_fg(self):
        """
        Linear Gradient
        Basically where delta RM 
        Sokoloff et al 1998
        """
        
        quarr = 1.0 * np.exp(-2.0 * self.RM_radm2 * self.lambda_sqr_array  - 2 * (self.delta_rm_2 * self.lambda_sqr_array)**2)
        self.quarr_fg = quarr

    def set_a_extrenal_dispersion_fg(self):
        """
        Sokoloff et al. (1998) Eq B3
        O'Sullivan et al. (2012) Eq 11
        """

        quarr = 1.0 * np.exp(-2.0 * self.sigmaRM_radm2**2.0 * self.lambda_sqr_array**2.0)
        self.quarr_fg = quarr

    def get_quarr(self):
        return self.quarr_fg



class Simulated_Src_Component:

    def __init__(self, lambda_sqr_array, frac_pol, RM_radm2, delta_rm_2, sigmaRM_radm2, psi0_deg):
        #Sets the source parameters
        #FIXME WHAT ARE WE DOING ABOUT RM? RM Center RM 

        self.lambda_sqr_array = lambda_sqr_array
        self.frac_pol = frac_pol
        self.RM_radm2 = RM_radm2
        self.delta_rm_2 = delta_rm_2
        self.sigmaRM_radm2 = sigmaRM_radm2
        self.psi0_deg = psi0_deg
        self.quarr_src = 0

    def set_thin_screen(self):
        """
    
        Simple Faraday thin source
    
        Ref:
        Sokoloff et al. (1998) Eq 2
        O'Sullivan et al. (2012) Eq 8
        Ma et al. (2019a) Eq 10
    
        """

        quarr = np.exp(2j * (np.radians(self.psi0_deg) + self.RM_radm2 * self.lambda_sqr_array))
        self.quarr_src = quarr


    def set_a_burn_slab_src(self):
        """
        Burn Slab
        Single Faraday component with differential Faraday rotation
        Burn (1966) Eq 18; with N >> (H_r/2H_z^0)^2
        Sokoloff et al. (1998) Eq 3
        O'Sullivan et al. (2012) Eq 9
        Ma et al. (2019a) Eq 12
        """

        quarr = np.exp( 2j * (np.radians(self.psi0_deg) + (0.5*self.sigmaRM_radm2 +RM) * self.lambda_sqr_array))
             * np.sin(self.sigmaRM_radm2 * self.lambda_sqr_array) / (self.sigmaRM_radm2 * self.lambda_sqr_array))
        

    def set_a_turbulent_slab_src(self):
        """
        Turbulent Slab
        Sokloff et al 1998 eqn 34
        delta_rm = 0
        sigRM_radm2 = nonzero
        """
        S = (2. * self.lambda_sqr_array**2 * self.sigmaRM_radm2**2 - 
             2j * self.lambda_sqr_array * self.delta_rm_2)
        quArr = (np.exp( 2j * (np.radians(self.psi0_deg) +
                                  self.RM_radm2 * self.lambda_sqr_array)) * 
            (1 - np.exp(-1.*S)) / S)
       
        self.quarr_src = quarr

    def set_a_mixed_slab_src(self):
        """
        Mixed Slab
        Sokloff et al 1998 eqn 34
        """
        S = (2. * self.lambda_sqr_array**2*self.sigmaRM_radm2**2)
        quArr = (np.exp( 2j * (np.radians(self.psi0_deg) +
                                  self.RM_radm2 * self.lambda_sqr_array)) * 
                                  (1 - np.exp(-1.*S)) / S)
       
        self.quarr_src = quarr
    
    def get_quarr(self):

        return self.quarr_src


class Simulated_LoS(Simulated_Src_Component, Simulated_Fg_Component):
    """
    Class provides the frame work to combine operators

    ### FIXME ADD THE TWO COMPONENT
    ### FIXME ADD IN THE TRACKING OF THE RANDOM SEEDS FOR NOISE GENERATION
    ### FIXME CONTEMPLATE ADDING NOISE - SIGMA OF NOISE DIST

    The theory behind it is that there are emission and propagation operators
    The emission operators need to be added and the propagation operators can be multiplied
    """
    def __init__(self, RM_screen, RM_src, sigma_RM_2):
        self.quarr_LoS = 0
        self.RM_screen_radm2 = RM_screen
        self.RM_src_radm2 = RM_src
        self.sigma_RM_2 = sigma_RM_2


    ##Dreaming big here
    def create_emission(self, src_model1,src_model2):
        #essentitally the emission operator
        self.quarr_LoS = src_model1 + src_model2


    def propagate_emission(self, src_model, fg_model):
        #propagation operator
        self.quarr_LoS = src_model * fg_model


    #But maybe too big and need to have set models

    def create_a_faraday_simple_LOS(self):
        """
        Faraday Simple: Uniform Emission Slab and the Uniform Foreground
        """

        self.quarrLos = self.frac_pol * np.exp(2j *(self.psi0_deg + self.RM_screen_radm2 * self.lambda_sqr_array))

    def create_a_burnslab_unf_fg_LOS(self):
        """
        Part 1 of the Uniform Burn Slab combined with the Uniform Foreground
        """

        self.quarrLos = self.frac_pol * ((np.sin(self.RM_src_radm2 * self.lambda_sqr_array))/(self.RM_src_radm2 * self.lambda_sqr_array)) * np.exp(2j * (self.psi0_deg + (1/2*self.RM_src_radm2 + self.RM_screen_radm2)*self.lambda_sqr_array))


    def create_a_unif_emiss_screen_lin_fg_LOS(self):
        """
        Uniform emission screen with a linear foreground screen convolved with a beam
        Inputs:
            -fractional polarization
            -polarization angle
            -RM_screen - RM at the center of the beam
            -sigma_RM_2 = the change in RM across the beam (sigma_b Delta RM/ D)
        """

        self.quarrLos = self.frac_pol * np.exp(2j * (self.psi0_deg + self.RM_screen_radm2 * self.lambda_sqr_array)) * np.exp(-2 * self.lambda_sqr_array**(2) * self.sigma_RM_2)
    
    def create_a_burn_slab_lin_fg_LOS(self):
        """
        Uniform emission screen with a linear foreground screen convolved with a beam
        Inputs:
            -fractional polarization
            -polarization angle
            -RM_screen - RM at the center of the beam
            -RM_src - RM of the bg source
            -sigma_RM_2 = the change in RM across the beam (sigma_b Delta RM/ D)
        """

        self.quarrLos = self.frac_pol * ((np.sin(self.RM_src_radm2 * self.lambda_sqr_array))/(self.RM_src_radm2 * self.lambda_sqr_array)) * np.exp(2j * (self.psi0_deg + 0.5*self.RM_src_radm2 * self.lambda_sqr_array)) * np.exp(-2 * self.lambda_sqr_array**(2) * self.sigma_RM_2)

    def create_a_turbulent_slab_unif_fg_LOS(self):
        """
        Uniform emission screen with a linear foreground screen convolved with a beam
        Inputs:
            -fractional polarization
            -polarization angle
            -RM_screen - RM at the center of the beam
            -sigma_RM_2 = the change in RM across the beam (sigma_b Delta RM/ D)
        """

        self.quarrLos = self.frac_pol * (1 - np.exp(-2*self.lambda_sqr_array^2 * self.sigma_RM_2))/(2*self.lambda_sqr_array^2 * self.sigma_RM_2) * (self.psi0_deg + self.RM_screen_radm2 * self.lambda_sqr_array))

    def get_quarr(self):
        return self.quarr_LoS


        

