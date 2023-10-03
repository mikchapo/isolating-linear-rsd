# f_max HOD model
# v0.1.0, 2021-03-03 - Started code using snippets from halotools tutorial (https://halotools.readthedocs.io/en/stable/quickstart_and_tutorials/tutorials/model_building/composing_models/hod_modeling/hod_modeling_tutorial1.html)

# Imports
from halotools.custom_exceptions import HalotoolsError
from halotools.empirical_models import OccupationComponent, Zheng07Cens, Zheng07Sats
from scipy.special import erf
import numpy as np
import warnings

baseline_model_params = {'omegam': 0.30609356,
                         'ombh2': 0.04868216126407447,
                         'sigma8': 0.80294489,
                         'h': 0.6805418299999999,
                         'ns': 0.95655848,
                         'N_eff': 3.046,
                         'w': -1.0,
                         'logM_sat': 14.312161,
                         'alpha': 0.54126815,
                         'logM_cut': 13.810948,
                         'sigma_logM': 1.081985,
                         'v_bc': 0.56097753,
                         'v_bs': 0.7405123,
                         'c_vir': 1.6794195,
                         'f': 0.68929181}

reid_logM_min = 13.031
prim_haloprop_key = 'halo_mvir'
default_luminosity_threshold = -20


class FMaxAemulusCens(OccupationComponent):
    def __init__(self,
        prim_haloprop_key=prim_haloprop_key,
        threshold=default_luminosity_threshold,
        logM_min=reid_logM_min,
        sigma_logM=baseline_model_params['sigma_logM'],
        fmax=1.,
        **kwargs):
        upper_occupation_bound = 1.0

        # Call the super class constructor, which binds all the
        # arguments to the instance.
        super(FMaxAemulusCens, self).__init__(gal_type='centrals',
            upper_occupation_bound=upper_occupation_bound,
            prim_haloprop_key=prim_haloprop_key,
            threshold=threshold,
            logM_min=logM_min,
            sigma_logM=sigma_logM,
            fmax=fmax,
            **kwargs)

        self.param_dict = {'logM_min': logM_min,
            'sigma_logM': sigma_logM,
            'fmax': fmax}

    def mean_occupation(self, **kwargs):
        if 'table' in list(kwargs.keys()):
            mass = kwargs['table'][self.prim_haloprop_key]
        elif 'prim_haloprop' in list(kwargs.keys()):
            mass = np.atleast_1d(kwargs['prim_haloprop'])
        else:
            msg = ("\nYou must pass either a ``table`` or ``prim_haloprop`` argument \n"
                "to the ``mean_occupation`` function of the ``Zheng07Cens`` class.\n")
            raise HalotoolsError(msg)

        logM = np.log10(mass)
        mean_ncen = self.param_dict['fmax']*0.5*(1.0 + erf(
            (logM - self.param_dict['logM_min']) / self.param_dict['sigma_logM']))

        return mean_ncen


class FMaxAemulusSats(OccupationComponent):
    def __init__(self,
            prim_haloprop_key=prim_haloprop_key,
            threshold=default_luminosity_threshold,
            modulate_with_cenocc=True,
            cenocc_model=None,
            logM_min=reid_logM_min,
            sigma_logM=baseline_model_params['sigma_logM'],
            fmax=1.,
            logM_sat=baseline_model_params['logM_sat'],
            alpha=baseline_model_params['alpha'],
            logM_cut=baseline_model_params['logM_cut'],
            **kwargs):

        upper_occupation_bound = float("inf")

        super(FMaxAemulusSats, self).__init__(
            gal_type='satellites', threshold=threshold,
            upper_occupation_bound=upper_occupation_bound,
            prim_haloprop_key=prim_haloprop_key,
            **kwargs)

        self.param_dict = {'logM_min': logM_min,
            'sigma_logM': sigma_logM,
            'fmax': fmax,
            'logM_sat': logM_sat,
            'alpha': alpha,
            'logM_cut': logM_cut}

        if cenocc_model is None:
            cenocc_model = FMaxAemulusCens(
                prim_haloprop_key=prim_haloprop_key, threshold=threshold)
        else:
            if modulate_with_cenocc is False:
                msg = ("You chose to input a ``cenocc_model``, but you set the \n"
                    "``modulate_with_cenocc`` keyword to False, so your "
                    "``cenocc_model`` will have no impact on the model's behavior.\n"
                    "Be sure this is what you intend before proceeding.\n"
                    "Refer to the Zheng et al. (2007) composite model tutorial for details.\n")
                warnings.warn(msg)

        self.modulate_with_cenocc = modulate_with_cenocc
        if self.modulate_with_cenocc:
            try:
                assert isinstance(cenocc_model, OccupationComponent)
            except AssertionError:
                msg = ("The input ``cenocc_model`` must be an instance of \n"
                    "``OccupationComponent`` or one of its sub-classes.\n")
                raise HalotoolsError(msg)

            self.central_occupation_model = cenocc_model

            # self.param_dict.update(self.central_occupation_model.param_dict)


    def mean_occupation(self, **kwargs):
        if self.modulate_with_cenocc:
            for key, value in list(self.param_dict.items()):
                if key in self.central_occupation_model.param_dict:
                    if key == "fmax":
                        self.central_occupation_model.param_dict[key] = 1.
                    else:
                        self.central_occupation_model.param_dict[key] = value

        # Retrieve the array storing the mass-like variable
        if 'table' in list(kwargs.keys()):
            mass = kwargs['table'][self.prim_haloprop_key]
        elif 'prim_haloprop' in list(kwargs.keys()):
            mass = np.atleast_1d(kwargs['prim_haloprop'])
        else:
            msg = ("\nYou must pass either a ``table`` or ``prim_haloprop`` argument \n"
                "to the ``mean_occupation`` function of the ``Zheng07Sats`` class.\n")
            raise HalotoolsError(msg)

        M_sat = 10.**self.param_dict['logM_sat']
        M_cut = 10.**self.param_dict['logM_cut']

        # Call to np.where raises a harmless RuntimeWarning exception if
        # there are entries of input logM for which mean_nsat = 0
        # Evaluating mean_nsat using the catch_warnings context manager
        # suppresses this warning
        mean_nsat = np.zeros_like(mass)

        idx_nonzero = np.where(mass - M_cut > 0)[0]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)

            mean_nsat[idx_nonzero] = (mass[idx_nonzero]/M_sat)**self.param_dict['alpha']*np.exp(-M_cut/mass[idx_nonzero])

        # If a central occupation model was passed to the constructor,
        # multiply mean_nsat by an overall factor of mean_ncen
        if self.modulate_with_cenocc:
            # compatible with AB models
            mean_ncen = getattr(self.central_occupation_model, "baseline_mean_occupation", self.central_occupation_model.mean_occupation)(**kwargs)
            mean_nsat *= mean_ncen

        return mean_nsat