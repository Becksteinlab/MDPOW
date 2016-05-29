import pytest

import mdpow.config
import mdpow.forcefields

# currently supported
WATERMODELS = ('tip4p', 'tip3p', 'tip5p', 'spc', 'spce')
SOLVENTMODELS = ('water', 'cyclohexane', 'octanol')

class TestIncludedForcefiels(object):
    def test_default_forcefield(self):
        assert mdpow.forcefields.DEFAULT_FORCEFIELD == "OPLS-AA"

    def test_oplsaa_itp(self):
        assert "ffoplsaa.itp" in mdpow.config.topfiles
        assert mdpow.config.topfiles["ffoplsaa.itp"].endswith('mdpow/top/ffoplsaa.itp')

    def test_oplsaa_ff(self):
        assert "oplsaa.ff" in mdpow.config.topfiles
        assert mdpow.config.topfiles["oplsaa.ff"].endswith('mdpow/top/oplsaa.ff')

class TestIncludedSolvents(object):
    solvents = {
        'tip4p': {
            'tip4p.itp': 'mdpow/top/oplsaa.ff/tip4p.itp',
            'tip4p.gro': 'mdpow/top/tip4p.gro'
        },
        'octanol': {
            '1oct.gro': 'mdpow/top/1oct.gro',
            '1oct.itp': 'mdpow/top/oplsaa.ff/1oct.itp',
        },
        'cyclohexane': {
            '1cyclo.gro': 'mdpow/top/1cyclo.gro',
            '1cyclo.itp': 'mdpow/top/oplsaa.ff/1cyclo.itp'
        },
    }

    # using nosetest-style test generators.. .feel free to rewrite for
    # py.test but I found the docs for parameterizing tests in py.test
    # too complicated

    def _test_solvent(self, name):
        solvent = self.solvents[name]
        def _assert_has_filename(filename):
            assert filename in mdpow.config.topfiles
        def _assert_correct_path(filename, path):
            assert mdpow.config.topfiles[filename].endswith(path)

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

    def test_default_water_model(self):
        assert mdpow.forcefields.DEFAULT_WATER_MODEL == "tip4p"

    def test_watermodelsdat(self):
        included_watermodels = open(mdpow.config.topfiles['watermodels.dat']).read()
        for line, ref in zip(self._simple_line_parser(mdpow.forcefields.GMX_WATERMODELS_DAT),
                             self._simple_line_parser(included_watermodels)):
            assert line.strip() == ref.strip()

    def test_gromacs_water_models(self):
        models = mdpow.forcefields.GROMACS_WATER_MODELS
        def has_identifier(identifier):
            assert identifier in models
        def itp_in_top(identifier):
            model = models[identifier]
            assert model.itp in mdpow.config.topfiles
        def coordinates_in_top(identifier):
            model = models[identifier]
            assert model.coordinates in mdpow.config.topfiles

        for identifier in self.watermodels:
            yield has_identifier, identifier
            yield itp_in_top, identifier
            yield coordinates_in_top, identifier

    def _simple_line_parser(self, string):
        for line in string.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            yield line

    def test_get_water_model(self):
        model = mdpow.forcefields.DEFAULT_WATER_MODEL
        assert mdpow.forcefields.get_water_model(model) is mdpow.forcefields.GROMACS_WATER_MODELS[model]

    def test_get_water_model_ValueError(self):
        with pytest.raises(ValueError):
            mdpow.forcefields.get_water_model("The Jabberwock is an imaginary beast.")


class TestSolventModels(object):
    watermodels = WATERMODELS
    solventmodels = [model for model in SOLVENTMODELS if model != "water"]

    def test_get_solvent_default_water(self):
        model = "water"
        defaultmodel = mdpow.forcefields.DEFAULT_WATER_MODEL
        assert mdpow.forcefields.get_solvent_model(model) is mdpow.forcefields.GROMACS_WATER_MODELS[defaultmodel]

    def test_get_solvent_model_ValueError(self):
        with pytest.raises(ValueError):
            mdpow.forcefields.get_solvent_model("The Jabberwock is an imaginary beast.")

    def test_get_solvent_cyclohexane(self):
        model = 'cyclohexane'
        assert mdpow.forcefields.get_solvent_model(model) is mdpow.forcefields.GROMACS_SOLVENT_MODELS[model]

    def test_get_solvent_octanol(self):
        model = 'octanol'
        assert mdpow.forcefields.get_solvent_model(model) is mdpow.forcefields.GROMACS_SOLVENT_MODELS[model]

    def test_get_solvent(self):
        def _assert_water_model(model):
            assert mdpow.forcefields.get_solvent_model(model) is mdpow.forcefields.GROMACS_WATER_MODELS[model]

        def _assert_solvent_model(model):
            assert mdpow.forcefields.get_solvent_model(model) is mdpow.forcefields.GROMACS_SOLVENT_MODELS[model]

        for model in mdpow.forcefields.GROMACS_WATER_MODELS:
            yield _assert_water_model, model

        for model in mdpow.forcefields.GROMACS_SOLVENT_MODELS:
            yield _assert_solvent_model, model

    def test_get_solvent_identifier_default_is_water(self):
        assert mdpow.forcefields.get_solvent_identifier('water') is mdpow.forcefields.DEFAULT_WATER_MODEL

    def test_get_solvent_identifier_water(self):
        def _assert_model(model):
            assert mdpow.forcefields.get_solvent_identifier('water', model=model) is model

        for model in self.watermodels:
            yield _assert_model, model

    def test_get_solvent_identifier_solvents(self):
        def _assert_model(solvent, model):
            assert mdpow.forcefields.get_solvent_identifier(solvent, model=model) is solvent

        for solvent in self.solventmodels:
            yield _assert_model, solvent, None

        # make sure that model is ignored
        for solvent in self.solventmodels:
            yield _assert_model, solvent, "Jabberwock model"

    def test_get_solvent_identifier_None(self):
        assert mdpow.forcefields.get_solvent_identifier('water', model="foobar") is None
        assert mdpow.forcefields.get_solvent_identifier('benzene') is None
