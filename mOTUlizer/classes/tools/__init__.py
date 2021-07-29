from abc import ABC, abstractmethod
import sys, shutil
from mOTUlizer.errors import *
from subprocess import call

class Wrapper(ABC):

    @abstractmethod
    def get_command(self):
        pass

    @abstractmethod
    def parse_output(self):
        pass

    def __init__(self):
        tests = {exe : shutil.which(exe) for exe in self.__executables__}
        if not all(list(tests.values())):
            raise CantFindExeError(f"Some of the executables necessary for {type(self)} ain't ere.\nYou are missing {[k for k,v in tests if not v]}")

    def get_requs(self):
        return self.__executables__

    def run_command(self):
        try :
            call(self.get_command(), shell = True)
        except :
            raise CantRunError(f"Some error happened while executing an {type(self)} object.\nThe command it was trying to run was {self.get_command()}")
