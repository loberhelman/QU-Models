"""
There are three classes here the first two are parents of the final class. The two parent classes are called 
Simulated_Fg_Component and Simulated_Src_Component and simulate different foreground and source models respectively. 
The final class is called Simulated_Components and is a child of the two classes and will work to combine foreground 
and source models to create a complete line of sight. 

WHAT ABOUT THE BEAM MODELS?


"""

import numpy as np


class Simulated_LoS:
    """
    Class provides the frame work to combine operators


    The theory behind it is that there are emission and propagation operators
    The emission operators need to be added and the propagation operators can be multiplied
    """

    #model_list = {'faraday': self.create_a_faraday_simple_LOS}
               #'create_a_burnslab_unf_fg_LOS',
               #'create_a_unif_emiss_screen_lin_fg_LOS', 'create_a_burn_slab_lin_fg_LOS', 'create_a_turbulent_slab_ran_fg_LOS']


    def __init__(self, RM_screen, RM_src, sigma_RM_2, psi0, lambda_sqr_array, frac_pol_emn, sigma_RM_FG_2, frac_pol_second, psi0_2_deg, RM_screen_2_radm2, SNR, model_name):
        self.model_list = {'faraday_simple': self._create_a_faraday_simple_LOS,
                           'unif_screen':self._create_a_unif_emis_screen_LOS,
                           'burnslab': self._create_a_unifrom_slab_LOS,
                           'burnslab_unif_fg': self._create_a_burnslab_unf_fg_LOS,
                           'unif_screen_lin_fg': self._create_a_unif_emiss_screen_lin_fg_LOS,
                           'burnslab_lin_fg':self._create_a_burn_slab_lin_fg_LOS,
                           'turbslab_ran_fg':self._create_a_mixed_slab_ran_fg_LOS,
                           'mixedslab_unif_fg':self._create_a_mixed_slab_unf_fg_LOS,
                           'two_component':self._create_two_component_LOS}
        

        self.quarr_LoS = 0
        self.RM_screen_radm2 = RM_screen
        self.RM_src_radm2 = RM_src
        self.sigma_RM_2 = sigma_RM_2
        self.psi0_deg = psi0
        self.lambda_sqr_array = lambda_sqr_array
        self.frac_pol = frac_pol_emn
        self.frac_pol_2_comp = frac_pol_second
        self.sigma_RM_FG_2 = sigma_RM_FG_2
        self.psi0_2_deg = psi0_2_deg
        self.RM_screen_2_radm2 = RM_screen_2_radm2
        self.SNR = SNR
        self.model_n = model_name
    

        if model_name not in self.model_list:
            err = f'model name :{model_name} unknown, choose one of:\n{list(self.model_list.keys())}'
            raise TypeError(err)
        
        self.model = self.model_list[model_name]()


    #But maybe too big and need to have set models

    def _create_a_unif_emis_screen_LOS(self):
        """
        Faraday Simple: Uniform Emission Slab and the Uniform Foreground
        Inputs:
            -fractional polarization
            -polarization angle
            -RM_screen
        """

        self.quarr_LoS = self.frac_pol * np.exp(2j * np.radians(self.psi0_deg)*self.lambda_sqr_array)
        self.sigma_RM_2 = 0
        self.sigma_RM_FG_2 = 0
        self.RM_src_radm2 = 0
        self.RM_screen_radm2 = 0
        self.frac_pol_2_comp = 0
        self.psi0_2_deg=0
        self.RM_screen_2_radm2 = 0


    def _create_a_unifrom_slab_LOS(self):
        """
        Faraday Simple: Uniform Emission Slab and the Uniform Foreground
        Inputs:
            -fractional polarization
            -polarization angle
            -RM_screen
        """

        self.quarr_LoS = self.frac_pol * np.exp(2j * (np.radians(self.psi0_deg) + (0.5*self.RM_src_radm2 *self.lambda_sqr_array)) * ((np.sin(self.RM_src_radm2 * self.lambda_sqr_array))/(self.RM_src_radm2 * self.lambda_sqr_array)))
        self.sigma_RM_2 = 0
        self.sigma_RM_FG_2 = 0
        self.RM_screen_radm2 =0
        self.frac_pol_2_comp=0
        self.psi0_2_deg=0
        self.RM_screen_2_radm2=0

    def _create_a_faraday_simple_LOS(self):
        """
        Faraday Simple: Uniform Emission Slab and the Uniform Foreground
        Inputs:
            -fractional polarization
            -polarization angle
            -RM_screen
        """

        self.quarr_LoS = self.frac_pol * np.exp(2j * (np.radians(self.psi0_deg) + self.RM_screen_radm2 * self.lambda_sqr_array))
        self.sigma_RM_2 = 0
        self.sigma_RM_FG_2 = 0
        self.RM_src_radm2 = 0
        self.frac_pol_2_comp=0
        self.psi0_2_deg=0
        self.RM_screen_2_radm2=0


       
    def _create_a_burnslab_unf_fg_LOS(self):
        """
        Part 1 of the Uniform Burn Slab combined with the Uniform Foreground
        Inputs:
            -fractional polarization
            -polarization angle
            -RM_screen
            -RM_source
        """
        print('creating slab')
        self.quarr_LoS = self.frac_pol * np.exp(2j * (np.radians(self.psi0_deg) + (0.5*self.RM_src_radm2 + self.RM_screen_radm2)*self.lambda_sqr_array)) * ((np.sin(self.RM_src_radm2 * self.lambda_sqr_array))/(self.RM_src_radm2 * self.lambda_sqr_array)) 
        self.sigma_RM_2 = 0
        self.sigma_RM_FG_2 = 0
        self.frac_pol_2_comp=0
        self.psi0_2_deg=0
        self.RM_screen_2_radm2=0

    def _create_a_unif_emiss_screen_lin_fg_LOS(self):
        """
        Uniform emission screen with a linear foreground screen convolved with a beam
        Inputs:
            -fractional polarization
            -polarization angle
            -RM_screen - RM at the center of the beam
            -sigma_RM_2 = the change in RM across the beam (sigma_b Delta RM/ D)
        """
        
        self.quarr_LoS = self.frac_pol * np.exp(2j * (np.radians(self.psi0_deg) + self.RM_screen_radm2 * self.lambda_sqr_array)) * np.exp(-2 * self.lambda_sqr_array**(2) * self.sigma_RM_2**2.)
        self.sigma_RM_FG_2 = 0
        self.RM_src_radm2 = 0
        self.frac_pol_2_comp=0
        self.psi0_2_deg=0
        self.RM_screen_2_radm2=0

    def _create_a_burn_slab_lin_fg_LOS(self):
        """
        Uniform emission screen with a linear foreground screen convolved with a beam
        Inputs:
            -fractional polarization
            -polarization angle
            -RM_screen - RM at the center of the beam
            -RM_src - RM of the bg source
            -sigma_RM_2 = the change in RM across the beam (sigma_b Delta RM/ D)
        """

        #self.quarr_LoS =  ((np.sin(self.RM_src_radm2 * self.lambda_sqr_array))/(self.RM_src_radm2 * self.lambda_sqr_array)) * np.exp(2j * self.lambda_sqr_array * self.RM_screen_radm2 - 2 * self.lambda_sqr_array**(2) * self.sigma_RM_2**2.))
        self.quarr_LoS = (self.frac_pol * np.exp(2j * (np.radians(self.psi0_deg) + (0.5*self.RM_src_radm2 * self.lambda_sqr_array))) * ((np.sin(self.RM_src_radm2 * self.lambda_sqr_array))/(self.RM_src_radm2 * self.lambda_sqr_array)) * np.exp(2j * self.RM_screen_radm2 * self.lambda_sqr_array - 2 * self.lambda_sqr_array**2.0 * self.sigma_RM_FG_2**2.0))
        self.sigma_RM_FG_2 = 0
        self.psi0_2_deg=0
        self.frac_pol_2_comp=0
        self.RM_screen_2_radm2=0

    def _create_a_mixed_slab_ran_fg_LOS(self):
        """
        Turbulent Slab with a linear RM gradient foreground screen convolved with a beam
        Inputs:
            -fractional polarization
            -polarization angle
            -RM_screen - RM at the center of the beam
            -sigma_RM_2 = the change in RM across the beam (sigma_b Delta RM/ D)
        """
        
        para_S = (2. * self.lambda_sqr_array**2 * self.sigma_RM_2**2 -  2j * self.lambda_sqr_array * self.RM_src_radm2)

        #self.quarr_LoS  = self.frac_pol * np.exp( 2j * (np.radians(self.psi0_deg) + self.RM_screen_radm2*self.lambda_sqr_array)) * ((1 - np.exp(-1.*para_S)) / para_S) * np.exp(-2.0 * self.sigma_RM_FG_2 ** 2.0 * self.lambda_sqr_array**2.0)
        self.quarr_LoS  = self.frac_pol * np.exp( 2j * np.radians(self.psi0_deg)) * ((1 - np.exp(-1.*para_S)) / para_S) * np.exp(-2.0 * self.sigma_RM_FG_2 ** 2.0 * self.lambda_sqr_array**2.0)

        #self.RM_src_radm2 = 0
        self.frac_pol_2_comp=0
        self.psi0_2_deg=0
        self.RM_screen_2_radm2=0
        self.RM_screen_radm2 = 0

    def _create_a_mixed_slab_unf_fg_LOS(self):
        """
        Turbulent Slab with a uniform foreground screen convolved with a beam
        Inputs:
            -fractional polarization
            -polarization angle
            -RM_screen - RM at the center of the beam
            -sigma_RM_2 = the change in RM across the beam (sigma_b Delta RM/ D)
        """
    

        para_S = (2. * self.lambda_sqr_array**2 * self.sigma_RM_2**2 - 2j * self.lambda_sqr_array * self.RM_src_radm2)
        self.quarr_LoS = self.frac_pol * np.exp( 2j * (np.radians(self.psi0_deg) + self.RM_screen_radm2 * self.lambda_sqr_array)) * ((1 - np.exp(-1.*para_S)) / para_S) 
        self.frac_pol_2_comp=0
        self.psi0_2_deg = 0
        self.sigma_RM_FG_2=0
        self.RM_screen_2_radm2=0


    def _create_two_component_LOS(self):

        self.quarr_LoS = (self.frac_pol * np.exp(2j * (np.radians(self.psi0_deg) + self.RM_screen_radm2* self.lambda_sqr_array))) + (self.frac_pol_2_comp * np.exp(2j * (np.radians(self.psi0_2_deg) + self.RM_screen_2_radm2* self.lambda_sqr_array)))
        self.sigma_RM_FG_2 = 0

    def get_quarr(self):
     
        return self.quarr_LoS
    
    def get_model_parameters(self):
        return self.frac_pol, self.psi0_deg, self.psi0_2_deg, self.RM_screen_radm2, self.RM_screen_2_radm2, self.RM_src_radm2, self.sigma_RM_2, self.sigma_RM_FG_2, self.frac_pol_2_comp, self.SNR, self.model_n


        

