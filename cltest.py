import pyopencl as cl 
import numpy as np 
import sys
from time import time

fname = "./chroms/chr1.fa"
csize = 1024*1024*64

f = open(fname, "rb")
data = f.read(csize)
f.close()
h_seq = np.frombuffer(data, dtype=np.uint8)
h_seq = h_seq.astype(dtype=np.int32)
h_mapped = np.empty(csize).astype(np.int32)

rtime = time()

kernelsource = '''
__kernel void freqtab(
	const int N,
	__global int* seq,
	__global int* mapped_seq
)
{
	int i = get_global_id(0);
	if ( i < N ) {
		mapped_seq[i] = seq[i] * 10 + 1;
	}
}
'''

context = cl.create_some_context()
queue = cl.CommandQueue(context)
program = cl.Program(context, kernelsource).build()
freqtab = program.freqtab
freqtab.set_scalar_arg_dtypes([np.int32, None, None, None])

d_seq = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf = h_seq)
d_mapped = cl.Buffer(context, cl.mem_flags.WRITE_ONLY, h_mapped.nbytes)

N = csize
globalsize = (N,)

freqtab(queue, globalsize, None, N, d_seq, d_freq, localmem)

cl.enqueue_copy(queue, h_freq, d_freq)

print(np.sum(h_freq.reshape((int(csize/N), N)), axis=0)[:4])

rtime = time() - rtime
print(rtime, " seconds")
