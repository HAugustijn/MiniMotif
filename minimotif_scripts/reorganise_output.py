"""  """

import os
import re
import shutil


def movetodir(outdir, dirname, patterns):
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
        for pattern in patterns:
            if re.search(pattern, f):
                if "results.tsv" in f:
                    pass
                else:
                    try:
                        shutil.move(os.path.join(outdir, f), os.path.join(outdir, dirname))
                    except:
                        pass
