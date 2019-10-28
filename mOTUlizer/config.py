from tqdm import tqdm

THREADS = 4
VERBOSE = False
DB_FOLDER = "/home/moritz/dbs/"

def vtqdm(args):
    if VERBOSE:
        return tqdm(args)
    else :
        return args
