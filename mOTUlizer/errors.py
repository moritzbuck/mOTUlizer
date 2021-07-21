class CantAminoAcidsError(Exception):
    """Raised when amino acid sequences cannot be inputed"""
    pass

class CantGFFError(Exception):
    """Raised when GFF is missing"""
    pass

class CantNucleotideError(Exception):
    """Raised when nucleotide file is missing"""
    pass

class GenomeIdError(Exception):
    """Raised when your genome IDs are somehow wrong"""
    pass

class FileError(Exception):
    """Something is wrong with a file..."""
    pass
