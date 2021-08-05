import os.path

import pytest

from .. import config
from .. import forcefields

# currently supported
WATERMODELS = ('tip4p', 'tip3p', 'tip5p', 'spc', 'spce', 'm24', 'tip4pd')
SOLVENTMODELS = ('water', 'cyclohexane', 'octanol')

class TestIncludedForcefiels(object):
    @staticmethod
    def test_default_forcefield():
        assert forcefields.DEFAULT_FORCEFIELD == "OPLS-AA"

    @staticmethod
    def test_oplsaa_itp():
        assert "ffoplsaa.itp" in config.topfiles
        assert config.topfiles["ffoplsaa.itp"].endswith(
            os.path.join('mdpow', 'top', 'ffoplsaa.itp'))

    @staticmethod
    def test_oplsaa_ff():
        assert "oplsaa.ff" in config.topfiles
        assert config.topfiles["oplsaa.ff"].endswith(
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

    @pytest.mark.parametrize("solvent_name", ["tip4p", "octanol", "cyclohexane"])
    def test_solvent(self, solvent_name):
        solvent = self.solvents[solvent_name]
        for filename, path in solvent.items():
            assert filename in config.topfiles
            assert config.topfiles[filename].endswith(path)


class TestWatermodels(object):
    @staticmethod
    def test_default_water_model():
        assert forcefields.DEFAULT_WATER_MODEL == "tip4p"

    def test_watermodelsdat(self):
        included_watermodels = open(config.topfiles['watermodels.dat']).read()
        for line, ref in zip(self._simple_line_parser(forcefields.GMX_WATERMODELS_DAT),
                             self._simple_line_parser(included_watermodels)):
            assert line.strip() == ref.strip()

    @pytest.mark.parametrize('identifier', WATERMODELS)
    def test_gromacs_water_models(self, identifier):
        models = forcefields.GROMACS_WATER_MODELS

        assert identifier in models

        model = models[identifier]
        assert model.itp in config.topfiles
        assert model.coordinates in config.topfiles


    @staticmethod
    def _simple_line_parser(string):
        for line in string.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            yield line

    @staticmethod
    def test_get_water_model():
        model = forcefields.DEFAULT_WATER_MODEL
        assert (forcefields.get_water_model(model) is
                forcefields.GROMACS_WATER_MODELS[model])

    @staticmethod
    def test_get_water_model_ValueError():
        with pytest.raises(ValueError):
            forcefields.get_water_model("The Jabberwock is an imaginary beast.")

class TestSolventModels(object):
    @staticmethod
    def test_get_solvent_default_water():
        model = "water"
        defaultmodel = forcefields.DEFAULT_WATER_MODEL
        assert (forcefields.get_solvent_model(model) is
                forcefields.GROMACS_WATER_MODELS[defaultmodel])

    @staticmethod
    def test_get_solvent_model_ValueError():
        with pytest.raises(ValueError):
            forcefields.get_solvent_model("The Jabberwock is an imaginary beast.")

    @staticmethod
    def test_get_solvent_cyclohexane():
        model = 'cyclohexane'
        forcefield = 'OPLS-AA'
        assert (forcefields.get_solvent_model(model) is
                forcefields.GROMACS_SOLVENT_MODELS[forcefield][model])

    @pytest.mark.parametrize("forcefield", ['OPLS-AA', 'CHARMM', 'AMBER'])
    def test_get_solvent_octanol(self, forcefield):
        model = 'octanol'
        assert (forcefields.get_solvent_model(model, forcefield=forcefield) is
                forcefields.GROMACS_SOLVENT_MODELS[forcefield][model])

    @pytest.mark.parametrize("forcefield", ['OPLS-AA', 'CHARMM', 'AMBER'])
    def test_get_solvent_wetoctanol(self, forcefield):
        model = 'wetoctanol'
        assert (forcefields.get_solvent_model(model, forcefield=forcefield) is
                forcefields.GROMACS_SOLVENT_MODELS[forcefield][model])

    @staticmethod
    def test_get_solvent_identifier_default_is_water():
        assert (forcefields.get_solvent_identifier('water') is
                forcefields.DEFAULT_WATER_MODEL)

    @pytest.mark.parametrize("model", WATERMODELS)
    def test_get_solvent_identifier_water(self, model):
        assert forcefields.get_solvent_identifier('water', model=model) is model

    @pytest.mark.parametrize('solvent',
                             [model for model in SOLVENTMODELS if model != "water"])
    @pytest.mark.parametrize('model', [None, "Jabberwock model"])
    def test_get_solvent_identifier_solvents(self, solvent, model):
        # The model="Jabberwock model" checks that "model" is properly ignored.
        assert forcefields.get_solvent_identifier(solvent, model=model) is solvent

    @staticmethod
    def test_get_solvent_identifier_None():
        assert forcefields.get_solvent_identifier('water', model="foobar") is None
        assert forcefields.get_solvent_identifier('benzene') is None
