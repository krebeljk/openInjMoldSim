import shutil # copy files, folders
import os # file paths
import fileinput # write to file
import re # regular expressions
import subprocess as supr # run external proc
import numpy as np
import hashlib

class Case:
    """
    Simulation instance
    """

    # static data members
    source = 'case_00_org'
    flnm_inp_dir = 'org0'
    flnm_inp = 'T'
    cmds = [ './Allrun' ]

    def __init__(self
        , dest = 'case1'
        , par1 = 0.477e6
        ):

        self.dest = dest
        self.par1 = par1
        self.par_dict = { 'PY_PL_HLD_01' : self.par1}
        self.org_inp = os.path.join(self.dest, self.flnm_inp_dir, self.flnm_inp)

    def make_case(self):
        shutil.copytree(self.source, self.dest)

    def insert_param(self):
        for line in fileinput.input(self.org_inp, inplace=1): # modify file
            for key in self.par_dict:
                line = re.sub(key ,str(self.par_dict[key]), line.rstrip())
            print(line)

    def run(self): # command exec
        for cmd in self.cmds:
            pop = supr.check_call( [cmd] , shell=True , cwd=self.dest)

    def clean(self):
        shutil.rmtree(self.dest)

    def conduct(self):
        self.make_case()
        self.insert_param()
        self.run()

def design_of_experiment():

    # parameters to use
    params = np.arange(120.0, 350.0, 30.0) # degC

    # loop through params
    for i in range(len(params)):
        dest = 'case_' + str(i+1).zfill(2)
        par1 = params[i] + 273.15 # convert to kelvin
        mojCase = Case(dest, par1)

        mojCase.conduct()
        # mojCase.clean()

def main():
    design_of_experiment()


if __name__ == '__main__':
    main()
