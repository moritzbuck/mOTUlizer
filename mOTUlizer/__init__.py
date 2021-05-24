import subprocess
try :
    label = subprocess.check_output(["git", "describe", "--tags"]).strip().decode()
except :
    label = "0.1.4"

__version__ =  label
