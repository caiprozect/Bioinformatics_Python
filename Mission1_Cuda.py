import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule

import numpy as np 
from itertools import product
from functools import reduce
from operator import mul

from time import time

def kMerCount(file, nK):
	fname = file
	K = nK
	f = open(fname, 'rb')
	f.readline()
	data = f.read().upper().splitlines()
	f.close()
	data = ("".encode('utf-8')).join(data)
	h_seq = np.frombuffer(data, dtype=np.uint8)
	h_seq = h_seq.astype(np.int32)

	mod = SourceModule("""
		__global__ void mapToNumb(
			const int N,
			const int M,
			const int K,
			const int numbKmer,
			int* seq,
			int* numb_seq,
			int* mono_freqs
		)
		{
			int gid = blockDim.x * blockIdx.x + threadIdx.x;
			int idx = gid * M;
			int i, letter;

			if(idx < N*M) {
			for(i=0; i < M; i++) {
				letter = seq[idx+i];
				switch(letter) {
					case 65:
						numb_seq[idx+i] = 0;
						atomicAdd(&mono_freqs[0], 1);
						break;
					case 67:
						numb_seq[idx+i] = 1;
						atomicAdd(&mono_freqs[1], 1);
						break;
					case 71:
						numb_seq[idx+i] = 2;
						atomicAdd(&mono_freqs[2], 1);
						break;
					case 84:
						numb_seq[idx+i] = 3;
						atomicAdd(&mono_freqs[3], 1);
						break;
					case 78:
						numb_seq[idx+i] = -1;
						break;
					default:
						numb_seq[idx+i] = -2;
						break;
				}
			}
		}
		}
		__global__ void freqTab(
		const int N,
		const int M,
		const int nK,
		const int numbKmer,
		int* numb_seq,
		int* freqs,
		int* checksum
	) {
		int gid = blockDim.x * blockIdx.x + threadIdx.x;
		int idx = gid * M;
		int i, numb;
		int k, loc_idx, ptn_idx;
		int dgt;
		int kmin;
		for(i=0; i < M; i++) {
			ptn_idx = 0;
			loc_idx = idx + i;
			kmin = 0;
			if(loc_idx < (N*M - nK + 1)) {
				for(k=0; k < nK; k++) {
					numb = numb_seq[loc_idx + k];
					switch(numb) {
						case (-1):
							atomicAdd(&checksum[1], 1);
							break;
						case (-2):
							atomicAdd(&checksum[0], 1);
							break;
						default:
							dgt = (int)(powf(4, (float)(nK-1-k)));
							ptn_idx += dgt * numb;
							break;
					}
					if(numb < kmin) {
						kmin = numb;
					} 
				}
				if(kmin >= 0) {
					atomicAdd(&freqs[ptn_idx], 1);
				}
			}
		}
	}
		""")

	gridSize = 1024
	blockSize = cuda.device_attribute.MAX_THREADS_PER_BLOCK

	seqLen = h_seq.size
	q, r = divmod(seqLen, gridSize*blockSize)
	q = q + 1
	h_seq = np.concatenate((h_seq, np.repeat(78, gridSize*blockSize-r).astype(np.int32)))

	h_numb_seq = np.zeros(h_seq.size).astype(np.int32)
	
	N = gridSize * blockSize
	M = q
	numbKmer = 4**K
	h_mono_freqs = np.zeros(4).astype(np.int32)

	mapToNumb = mod.get_function("mapToNumb")
	mapToNumb(np.int32(N), np.int32(M), np.int32(K), np.int32(numbKmer), cuda.In(h_seq), cuda.Out(h_numb_seq), cuda.Out(h_mono_freqs), block=(blockSize, 1, 1), grid=(gridSize, 1))

	h_freqs = np.zeros(numbKmer).astype(np.int32)
	h_checksum = np.zeros(2).astype(np.int32)

	freqTab = mod.get_function("freqTab")
	freqTab(np.int32(N), np.int32(M), np.int32(K), np.int32(numbKmer), cuda.In(h_numb_seq), cuda.Out(h_freqs), cuda.Out(h_checksum), block=(blockSize, 1, 1), grid=(gridSize, 1))

	assert(h_checksum[0] == 0), "File has unknown nucleotide character"
	print("Counting {} has been done".format(file))

	return h_mono_freqs.astype(np.int64), h_freqs.astype(np.int64)

def processOrganChromDict(organChromDict, nK, outfile):
	lNucls = ["A", "C", "G", "T"]
	lKmer = list(product(range(4), repeat=nK))
	with open(outfile, 'w') as f:
		seqNameList = organChromDict.keys()
		for seqName in seqNameList:
			chromNumbs = organChromDict[seqName]
			chromFileList = ["../data/{}.chroms/chr{}.fa".format(seqName, numb) for numb in chromNumbs]
			bothCnts = [kMerCount(file, nK) for file in chromFileList]
			monoCnts = np.sum(np.array([each[0] for each in bothCnts]), axis=0)
			monoFreqs = monoCnts / monoCnts.sum()
			polyCnts = np.sum(np.array([each[1] for each in bothCnts]), axis=0)
			polyFreqs = polyCnts / polyCnts.sum()
			polyExps = [reduce(mul, map((lambda x: monoFreqs[x]), kmer), 1) for kmer in lKmer]
			flds = "{:16}{:>16}{:>16}{:>16}\n"
			mono_entries = "{:16}{:16d}{:16.5f}\n"
			entries = "{:16}{:16d}{:16.5f}{:16.5f}\n"
			f.write("{} chromosomes statistics\n".format(seqName))
			f.write(flds.format("K-mer", "Counts", "Frequencies", "Expectations"))
			for i in range(4):
				f.write(mono_entries.format(lNucls[i], monoCnts[i], monoFreqs[i]))
			for i in range(4**nK):
				kmer = lKmer[i]
				kmer = reduce((lambda x,y: x+lNucls[y]), kmer, "")
				f.write(entries.format(kmer, polyCnts[i], polyFreqs[i], polyExps[i]))
			f.write("\n")

def main():
	dictChromCat = {'hg38': ([str(i) for i in range(1,23)]+['X','Y']), 'galGal3': ([str(i) for i in range(1,29)]+['32','W','Z']),
					'dm3': ['2R', '2L', '3R', '3L', '4', 'X'], 'ce10': ['I', 'II', 'III', 'IV', 'V', 'X']}
	K = 2
	outFile = "Mission1_Cuda.txt"
	processOrganChromDict(dictChromCat, K, outFile)

if __name__=="__main__":
	rtime = time()
	main()
	rtime = time() - rtime
	print("Run time: {} seconds".format(rtime))