#RRT implementation for reorientation with collision-free
from openravepy import *
from pylab import *

import time
import string
import numpy as np
import copy
import random
import os

import l