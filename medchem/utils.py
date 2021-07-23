import os
import medchem


def get_data(file=None):
    """Return the folder that contains the package specific data"""

    path = os.path.join(medchem.__path__[0], "data/")
    if file is not None:
        path = os.path.join(path, file)
    return path
