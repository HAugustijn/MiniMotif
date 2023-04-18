"""  """

import os
import re
import shutil


def movetodir(outdir, dirname, pattern):
    """ moves files matching a pattern into new directory
    parameters
    ----------
    outdir
        string, the path to output directory
    dirname
        string, name of the new direcory
    pattern
        string, regex

    """
    try:
        os.mkdir(os.path.join(outdir, dirname))
    except:
        pass
    # Move files into new directory
    for f in os.listdir(outdir):
        if re.search(pattern, f):
            try:
                shutil.move(os.path.join(outdir, f), os.path.join(outdir, dirname))
            except:
                pass