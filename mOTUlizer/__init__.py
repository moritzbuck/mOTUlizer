import subprocess
from importlib.util import find_spec
import os

cwd = os.getcwd()
module_path = find_spec("mOTUlizer").submodule_search_locations[0]

os.chdir(module_path)

try :
    label = subprocess.check_output(["git", "describe", "--tags"]).strip().decode()
except :
    label = "0.2.0"

os.chdir(cwd)

__version__ =  label
