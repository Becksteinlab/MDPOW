import pytest

import numpy as np
from numpy.testing import assert_array_almost_equal

from . import CONFIGURATIONS

from mdpow import config

class TestAlteredConfig(object):
    params_altered = {
        'DEFAULT':
            {
            'qscripts':'custom.sh'
            },
        'setup':
            {
            'name': 'custom_name',
            'molecule': 'some_molecule_ident',
            'itp': 'some_molecules_itp',
            'structure': 'some_molecules_structure',
            'watermodel': 'spce',
            'maxwarn': 2,
            'distance': None,          # default (not in this input file)
            'boxtype': 'dodecahedron', # default (not in this input file)
            'gromacsoutput': True,
            },
        'energy_minimize':
            {
            'mdp': 'custom_emin.mdp'
            },
        'MD_relaxed':
            {
            'qscript': 'MD_relaxed.sge',
            'runtime': 10,
            'runlocal': False,
            'mdp': 'MD_relaxed_NPT_opls.mdp'
            },
        'MD_NPT':
            {
            'qscript': 'MD_NPT.sge',
            'runtime': 10000,
            'runlocal': True,
            'mdp': 'MD_NPT_opls.mdp',
            },
        'FEP':
            {
            'method': 'TI',
            'qscript': 'FEP.sge',
            'runtime': 1000,
            'runlocal': True,
            'maxwarn': 3,
            'mdp': 'fep_custom_opls.mdp'
            },
        'FEP_schedule_Coulomb':
            {
            'name': 'Coul',
            'description': 'transition_1',
            'label': 'coulomb',
            'couple_lambda0': 'vdw',
            'couple_lambda1': 'vdw-q',
            'sc_alpha': 0.2,
            'sc_power': 2,
            'sc_sigma': 0.6,
            'lambdas': np.array([ 0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1. ])
            },
        'FEP_schedule_VDW':
            {
            'name': 'VANDERWAALS',
            'description': 'transition_2',
            'label': 'vanderwaals',
            'couple_lambda0': 'none',
            'couple_lambda1': 'vdw',
            'sc_alpha': 0,
            'sc_power': 3,
            'sc_sigma': 0.1,
            'lambdas': np.array([ 0.0 , 0.25 , 0.50 , 0.75 , 1 ])
            },
        'mdrun':
            {
            'stepout': 12000,
            'verbose': False,
            'nice': 12,
            'maxthreads': 1
            }
    }

    @pytest.fixture
    def cfg(self):
        return config.get_configuration(str(CONFIGURATIONS / 'altered_runinput.yml'))

    def _test_section(self, cfg, section):
        section_dict = self.params_altered[section]
        for k in section_dict.keys():
            if k == 'lambdas':
                parsed = np.array([float(x.strip()) for x in cfg.get(section,k).split(",")])
                assert_array_almost_equal(parsed, section_dict[k],
                                          err_msg="mismatch in lambdas")
            else:
                assert cfg.get(section,k) == section_dict[k], \
                    "mismatch in {}:{}".format(section,k)

    def test_DEFAULT(self, cfg):
        return self._test_section(cfg, "DEFAULT")

    def test_setup(self, cfg):
        return self._test_section(cfg, "setup")

    def test_energy_minimize(self, cfg):
        return self._test_section(cfg, "energy_minimize")

    def test_MD_relaxed(self, cfg):
        return self._test_section(cfg, "MD_relaxed")

    def test_MD_NPT(self, cfg):
        return self._test_section(cfg, "MD_NPT")

    def test_FEP(self, cfg):
        return self._test_section(cfg, "FEP")

    def test_FEP_schedule_Coulomb(self, cfg):
        return self._test_section(cfg, "FEP_schedule_Coulomb")

    def test_FEP_schedule_VDW(self, cfg):
        return self._test_section(cfg, "FEP_schedule_VDW")

    def test_mdrun(self, cfg):
        return self._test_section(cfg, "mdrun")
