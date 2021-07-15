import os.path

import pytest
import gromacs.cbook
from mdpow.config import topfiles
from mdpow.forcefields import *

# currently supported
WATERMODELS = ('tip4p', 'tip3p', 'tip5p', 'spc', 'spce', 'm24', 'tip4pd')
SOLVENTMODELS = ('water', 'cyclohexane', 'octanol')


class TestIncludedForcefiels(object):
    @staticmethod
    def test_default_forcefield():
        assert DEFAULT_FORCEFIELD == "OPLS-AA"

    @staticmethod
    def test_oplsaa_itp():
        assert "ffoplsaa.itp" in topfiles
        assert topfiles["ffoplsaa.itp"].endswith(
            os.path.join('mdpow', 'top', 'ffoplsaa.itp'))

    @staticmethod
    def test_oplsaa_ff():
        assert "oplsaa.ff" in topfiles
        assert topfiles["oplsaa.ff"].endswith(
            os.path.join('mdpow', 'top', 'oplsaa.ff'))

class TestIncludedSolvents(object):
    solvents = {
        'tip4p': {
            'tip4p.itp': os.path.join('mdpow', 'top', 'oplsaa.ff', 'tip4p.itp'),
            'tip4p.gro': os.path.join('mdpow', 'top', 'tip4p.gro')
        },
        'octanol': {
            '1oct.gro': os.path.join('mdpow', 'top', '1oct.gro'),
            '1oct.itp': os.path.join('mdpow', 'top', 'oplsaa.ff', '1oct.itp'),
        },
        'cyclohexane': {
            '1cyclo.gro': os.path.join('mdpow', 'top', '1cyclo.gro'),
            '1cyclo.itp': os.path.join('mdpow', 'top', 'oplsaa.ff', '1cyclo.itp')
        },
    }

    # using nosetest-style test generators.. .feel free to rewrite for
    # py.test but I found the docs for parameterizing tests in py.test
    # too complicated

    def _test_solvent(self, name):
        solvent = self.solvents[name]
        def _assert_has_filename(filename):
            assert filename in topfiles
        def _assert_correct_path(filename, path):
            assert topfiles[filename].endswith(path)

        for filename, path in solvent:
            yield _assert_has_filename, filename
            yield _assert_correct_path, filename, path

    def test_tip4p(self):
        self._test_solvent('tip4p')

    def test_octanol(self):
        self._test_solvent('octanol')

    def test_cyclohexane(self):
        self._test_solvent('cyclohexane')

class TestWatermodels(object):
    watermodels = WATERMODELS

    @staticmethod
    def test_default_water_model():
        assert DEFAULT_WATER_MODEL == "tip4p"

    def test_watermodelsdat(self):
        included_watermodels = open(topfiles['watermodels.dat']).read()
        for line, ref in zip(self._simple_line_parser(GMX_WATERMODELS_DAT),
                             self._simple_line_parser(included_watermodels)):
            assert line.strip() == ref.strip()

    def test_gromacs_water_models(self):
        models = GROMACS_WATER_MODELS
        def has_identifier(identifier):
            assert identifier in models
        def itp_in_top(identifier):
            model = models[identifier]
            assert model.itp in topfiles
        def coordinates_in_top(identifier):
            model = models[identifier]
            assert model.coordinates in topfiles

        for identifier in self.watermodels:
            yield has_identifier, identifier
            yield itp_in_top, identifier
            yield coordinates_in_top, identifier

    @staticmethod
    def _simple_line_parser(string):
        for line in string.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            yield line

    @staticmethod
    def test_get_water_model():
        model = DEFAULT_WATER_MODEL
        assert (get_water_model(model) is
                GROMACS_WATER_MODELS[model])

    @staticmethod
    def test_get_water_model_ValueError():
        with pytest.raises(ValueError):
            get_water_model("The Jabberwock is an imaginary beast.")


class TestSolventModels(object):
    watermodels = WATERMODELS
    solventmodels = [model for model in SOLVENTMODELS if model != "water"]

    @staticmethod
    def test_get_solvent_default_water():
        model = "water"
        defaultmodel = DEFAULT_WATER_MODEL
        assert (get_solvent_model(model) is
                GROMACS_WATER_MODELS[defaultmodel])

    @staticmethod
    def test_get_solvent_model_ValueError():
        with pytest.raises(ValueError):
            get_solvent_model("The Jabberwock is an imaginary beast.")

    @staticmethod
    def test_get_solvent_cyclohexane():
        model = 'cyclohexane'
        forcefield = 'OPLS-AA'
        assert (get_solvent_model(model) is
                GROMACS_SOLVENT_MODELS[forcefield][model])

    @pytest.mark.parametrize("forcefield", ['OPLS-AA', 'CHARMM', 'AMBER'])
    @staticmethod
    def test_get_solvent_octanol(forcefield):
        model = 'octanol'
        assert (get_solvent_model(model, forcefield=forcefield) is
                GROMACS_SOLVENT_MODELS[forcefield][model])

    @pytest.mark.parametrize("forcefield", ['OPLS-AA', 'CHARMM', 'AMBER'])
    @staticmethod
    def test_get_solvent_wetoctanol(forcefield):
        model = 'wetoctanol'
        assert (get_solvent_model(model, forcefield=forcefield) is
                GROMACS_SOLVENT_MODELS[forcefield][model])

    @staticmethod
    def test_get_solvent_identifier_default_is_water():
        assert (get_solvent_identifier('water') is
                DEFAULT_WATER_MODEL)

    def test_get_solvent_identifier_water(self):
        def _assert_model(model):
            assert get_solvent_identifier('water', model=model) is model

        for model in self.watermodels:
            yield _assert_model, model

    def test_get_solvent_identifier_solvents(self):
        def _assert_model(solvent, model):
            assert get_solvent_identifier(solvent, model=model) is solvent

        for solvent in self.solventmodels:
            yield _assert_model, solvent, None

        # make sure that model is ignored
        for solvent in self.solventmodels:
            yield _assert_model, solvent, "Jabberwock model"

    @staticmethod
    def test_get_solvent_identifier_None():
        assert get_solvent_identifier('water', model="foobar") is None
        assert get_solvent_identifier('benzene') is None
