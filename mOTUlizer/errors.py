class FileError(Exception):
    """Something is wrong with a file..."""
    pass

class CantAminoAcidsError(FileError):
    """Raised when amino acid sequences cannot be inputed"""
    pass

class CantGFFError(FileError):
    """Raised when GFF is missing"""
    pass

class CantCDSError(FileError):
    """Raised when CDS-fasta is missing"""
    pass

class CantGeneClusterError(Exception):
    """Raised when something makes creating GeneClusters problematic"""
    pass

class CantLoadGeneClusterError(CantGeneClusterError):
    """Raised when something makes creating GeneClusters problematic"""
    pass

class CantExeError(Exception):
    """A program is missing"""
    pass

class CantFindExeError(CantExeError):
    """A program is missing"""
    pass

class CantRunError(CantExeError):
    """A program is somehow not working"""
    pass

class CantGenesError(FileError):
    """Raised when nucleotide file is missing"""
    pass

class CantParseError(FileError):
    """Raised when something can't be parsed"""
    pass

class CantNucleotideError(FileError):
    """Raised when nucleotide file is missing"""
    pass

class GenomeIdError(Exception):
    """Raised when your genome IDs are somehow wrong"""
    pass

class CantMethodError(Exception):
    """Raised when a computing method you want don't exist"""
    pass

class GenomeMismatchError(Exception):
    """Raised when your somehow genome-sets that should be the same are not"""
    pass

class CoreNotComputedError(Exception):
    """Core hasn't been computed yet error."""
    pass

class DataBaseError(Exception):
    """Some Generic DataBaseError."""
    pass

class DataBaseNotInitialisedError(DataBaseError):
    """Db not initialised"""
    pass

class DataBaseBadValuesError(DataBaseError):
    """Some values passed to the database don't work"""
    pass
