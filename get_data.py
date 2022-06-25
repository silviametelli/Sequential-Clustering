#! /usr/bin/env python
import sys


old_key = ""
for line in sys.stdin:
    [key, val] = line.strip().split(",")
    if key != old_key:
        if old_key != "":
            print(old_key + "\t" + " ".join(str(c) for c in sorted(comp)))
        comp = [int(val)]
        old_key = key
    else:
        comp += [int(val)]
if old_key != "":
    print(old_key + "\t" + " ".join(str(c) for c in sorted(comp)))


            
