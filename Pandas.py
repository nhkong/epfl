#!/usr/bin/python

import pandas as pd

import math
df = pd.read_table('15.1.tmp', sep="\t", low_memory=False)


df.to_csv('16.tmp', index=False, sep="\t")
