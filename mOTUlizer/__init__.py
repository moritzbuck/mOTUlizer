import subprocess
from importlib.util import find_spec
import os

cwd = os.getcwd()
module_path = find_spec("mOTUlizer").submodule_search_locations[0]

os.chdir(module_path)

if os.path.exists(".git"):
    label = subprocess.check_output(["git", "describe", "--tags"]).strip().decode()
else:
    label = "0.4.0"

os.chdir(cwd)

__version__ =  label

_quiet_ = True
_temp_folder_ = "/tmp/"

def set_quiet(value : bool):
    _quiet_ = value

def get_quiet() -> bool:
    return _quiet_

_threads_ = 24

def set_threads(value : int):
    _threads_ = value

def get_threads() -> int:
    return _threads_
