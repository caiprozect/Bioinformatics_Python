import numpy as np 

fname = "sequence"

ascii_seq = np.fromfile(fname, dtype=np.uint8, count=-1, sep="")

a = (ascii_seq==65).sum()
c = (ascii_seq==67).sum()
g = (ascii_seq==71).sum()
t = (ascii_seq==84).sum()

print("# of A: ", a)
print("# of C: ", c)
print("# of G: ", g)
print("# of T: ", t)