import mdpow.equil
import tempfile
import os
import shutil
import mdpow.run
import gromacs.environment
gromacs.environment.flags['capture_output'] = True


class TestEquilibration(object):
    sim_filename = ".simulation"
    sims = {"water": mdpow.equil.WaterSimulation,  "octanol":mdpow.equil.OctanolSimulation, "cyclohexane":mdpow.equil.CyclohexaneSimulation}
    #sims = {"water": mdpow.equil.WaterSimulation}
    
    def setup(self):
        self.basedir = os.getcwd()
        self.directory_name = tempfile.mkdtemp()
        print("basedir:  {}".format(self.basedir))
        print("Copying required files into:  {}".format(self.directory_name))
        files = ["benzene.itp","benzene.pdb"]
        for f in files:
            shutil.copy("mdpow/tests/benzene/{}".format(f),"{0}/{1}".format(self.directory_name,f))
            print("\tbenzene/{0} --> {1}/{0}".format(f,self.directory_name))
        os.chdir(self.directory_name)
        print("Changed directory:  {}".format(os.getcwd()))
        for k in self.sims.keys():
            W = self.sims[k](molecule="BNZ")
            W.topology(itp="benzene.itp")
            W.solvate(struct="benzene.pdb")
            W.energy_minimize()
            W.MD_relaxed()
            mdrun = mdpow.run.MDrunnerSimple(dirname="/".join([self.directory_name,"Equilibrium",k,"MD_relaxed"]),deffnm="md")
            simulation_done = mdrun.run_check()
            W.MD(runtime=5,qscript=['local.sh'])
            W.save(k + self.sim_filename)
        
    def teardown(self):
        os.chdir(self.basedir)
        shutil.rmtree(self.directory_name)
        
    def test_simulation_file(self):
        for k in self.sims.keys():
            assert os.path.exists("{0}/{1}".format(self.directory_name,k + self.sim_filename))
    
    def test_equilibrium_solvent_dir(self):
        contents = ["em","MD_relaxed","solvation","top"]
        for k in self.sims.keys():
            solvent_dir = "{0}/{1}/{2}/".format(self.directory_name, "Equilibrium", k)
            assert os.path.exists(solvent_dir)
            for c in contents:
                assert os.path.exists(solvent_dir+c)
