import os
os.environ['SPARK_HOME'] = os.path.join(os.path.expanduser('~'),'spark')

import subprocess, re
from IPython.display import HTML

import numpy as np
from scipy.sparse import csr_matrix 

import imnet

strings = imnet.random_strings.generate_random_sequences(5000)
