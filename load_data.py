import numpy as np
import itertools
import subprocess
import pysam
import gzip
import os
import pdb

MIN_MAP_QUAL = 10
MAX_VAL = 65535
MIN_READ_LEN = 1
MAX_READ_LEN = 100000

class ZipFile():
    """File class to rapidly handle some file operations on 
    a gzip file. To greatly increase speed, methods in this class 
    use the system's `zcat` function to read each line of the 
    file, rather than the `gzip` module.

    Arguments
        filename : string
        name of the gzip file to parse.
    """

    def __init__(self, filename):

        if os.path.isfile(filename):
            pipe = subprocess.Popen(["zcat", filename], stdout=subprocess.PIPE)
            self.handle = pipe.stdout
            # remove header
            header = self.handle.next()
        else:
            raise IOError

    def _readline(self):

        for line in self.handle:
            yield line.strip().split('\t')

    def read(self, batch=None):
        """Reads in the lines of the file, either in batches
        or as a whole.

        Optional arguments
            batch : int
            read in `batch` number of lines at a time.

        """

        if batch is None:
            # read the whole file
            locations = [line for line in self._readline()]
        else:
            # read a chunk of the file
            locations = [loc for loc in itertools.islice(self._readline(), 0, batch)]

        for loc in locations:
            loc[1] = int(loc[1])
            loc[2] = int(loc[2])

        return locations

    def close(self): 
      
        pass

class BamFile():
    """File class to handle some file operations on 
    a bam file.

    Arguments
        filename : string
        name of the bam file to parse.

        protocol : string
        DNase_seq / ATAC_seq

    """

    def __init__(self, filename, protocol):

    	self._handle = pysam.Samfile(filename, "rb")
        self._protocol = protocol
        self._ref_lengths = dict([(r,l) \
            for r,l in zip(self._handle.references,self._handle.lengths)])

    def get_read_counts(self, locations, width=200):
        """Get the number of sequencing reads mapped to
        each base along a window centered at each of
        several motif instances.

        Arguments:
            locations : list
            each entry of the list is a list that specifies 
            information for a motif instance

            width : int
            length of the genomic window around the motif
            instance.

        """

        counts = []

        if self._protocol=='DNase_seq':

            for location in locations:

                chrom = location[0]
                strand = location[3]
                if strand=='+':
                    center = location[1]
                else:
                    center = location[2]
                left = max([1,center-width/2])
                right = min([center+width/2, self._ref_lengths[chrom]])

                # fetch all reads overlapping this genomic location
                sam_iter = self._handle.fetch(reference=chrom, start=left, end=right)
                forward = np.zeros((width,), dtype=np.float64)
                reverse = np.zeros((width,), dtype=np.float64)

                for read in sam_iter:

                    # skip read if unmapped
                    if read.is_unmapped:
                        continue

                    # skip read, if mapping quality is low
                    if read.mapq < MIN_MAP_QUAL:
                        continue

                    start = read.pos
                    end = start + read.alen - 1

                    # skip read, if site of cleavage is outside window
                    if (read.is_reverse and end >= right) or (not read.is_reverse and start < left):
                        continue

                    if read.is_reverse:
                        reverse[end-left] += 1
                    else:
                        forward[start-left] += 1

                # flip fwd and rev strand read counts, 
                # if the motif is on the opposite strand.
                if strand=='+':
                    count = [forward, reverse]
                else:
                    count = [reverse[::-1], forward[::-1]]

                count = np.hstack(count)

                # cap the read count at any location
                count[count>MAX_VAL] = MAX_VAL
                counts.append(count.astype(np.float64))

        else:

            for l,location in enumerate(locations):

                chrom = location[0]
                strand = location[3]
                if strand=='+':
                    center = location[1]
                else:
                    center = location[2]
                left = max([1,center-width/2])
                right = min([center+width/2, self._ref_lengths[chrom]])

                # fetch all reads overlapping this genomic location
                sam_iter = [read for read in self._handle.fetch(reference=chrom, start=left, end=right)]
                count = np.zeros((width,), dtype=np.float64)

                for read in sam_iter:

                    # discard anomalies
                    if not (read.is_read1 or read.is_read2):
                        continue

                    # require that (both) paired-reads are uniquely mapped
                    if read.is_unmapped or read.mate_is_unmapped:
                        continue

                    # paired reads mapped to same strand...
                    if read.is_reverse == read.mate_is_reverse:
                        continue

                    # read has poor mapping quality
                    if read.mapq < MIN_MAP_QUAL:
                        continue

                    # discard fragments that are too long or too short
                    isize = np.abs(read.isize)
                    if isize < MIN_READ_LEN or isize > MAX_READ_LEN:
                        continue

                    # pysam pos starts at 0, not 1
                    # site of adapter insertions are 9bp apart
                    # an offset of +4 / -5 gives approximate site of transposition
                    if read.is_reverse:
                        if read.pos<=read.mpos:
                            tpos = read.pos + 4
                        else:
                            tpos = read.mpos + isize - 5
                        # add read, if site of transposition is within window
                        if left<=tpos<right:
                            count[tpos-left] += 1
                    else:
                        if read.pos>=read.mpos:
                            tpos = read.mpos + isize - 5
                        else:
                            tpos = read.pos + 4
                        # add read, if site of transposition is within window
                        if left<=tpos<right:
                            count[tpos-left] += 1

                # flip read counts,
                # if the motif is on the opposite strand.
                if strand=='-':
                    count = count[::-1]

                counts.append(count)

        counts = np.array(counts).astype(np.float64)

        # cap the read count at all locations
        counts[counts>MAX_VAL] = MAX_VAL

        return counts

    def close(self):

        self._handle.close()
