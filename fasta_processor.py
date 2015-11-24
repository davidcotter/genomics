__author__ = 'DCotter'
import re
class FastProcessor:
    """Load the specified FASTA file for processing for processing
    """

    fasta_records = {}
    current_record = None
    start_codons = ["ATG"]
    stop_codons = ["TAA", "TAG", "TGA"]

    def read_fasta_file(self, filename):
        """ Load the specified FASTA file for processing
        :param filename: The filename of the FASTA file - can be a filename or a reletive or absolute path
        """
        with open(filename) as f:
           for line in [line.rstrip('\n') for line in open(filename)]:
               if line.startswith(">"):
                   name = re.findall('([^\s]+)' , line)[0][1:]
                   self.fasta_records[name] = ''
                   self.current_record = name
               else:
                    self.fasta_records[self.current_record] += line

    def record_count(self):
        """ Return the count of the number of records in the file
        :returns the record count
        """
        return len(self.fasta_records.keys())

    def record_lengths(self):
        """ Return the count of the number of records in the file
        :returns the record count
        """
        return list(self.map_lengths().values())

    def map_lengths(self):
        """ Map the name/value to name/length of the value:
        e.g.: A dict of {'a': 'aaa', 'b': 'bb'} would map to {'a': 3, 'b': 2}
        :returns a dict of name/length
        """
        return {k: len(v) for k, v in self.fasta_records.items()}

    def max_len_sequences(self):
        """ Return the names of the sequences with the max len
        e.g.: A dict of {'a': 'aaa', 'b': 'bb', 'c': 'ccc'} would return ['a', 'c'] since both a and c have the max len
        of 3
        :returns a list of names
        """
        return self.keys_for_aggregate_val(self.map_lengths(), max)

    def min_len_sequences(self):
        """ Return the names of the sequences with the min len
        e.g.: A dict of {'a': 'aaa', 'b': 'bb', 'c': 'ccc'} would return ['b'] since b has the min len of 2
        :returns a list of names
        """
        return self.keys_for_aggregate_val(self.map_lengths(), min)

    def keys_for_aggregate_val(self, dict, aggregate_fn):
        """ Return the names of the dict items satisfying aggregate_fn
        e.g.: A dict of {'a': 'aaa', 'b': 'bb', 'c': 'ccc'} would return ['a', 'c'] if 'max' was passed in as
        thr aggregate_fn since both a and c have the max len of 3
        :returns a list of names
        """
        _, the_value = aggregate_fn(dict.items(), key=lambda x:x[1])
        return [key for (key,value) in dict.items() if value == the_value]


    def open_reading_frames(self, sequence, frame_offset = 0):
        """ Return a list of open reading frames from a sequence - i.e. a sequence starting with a start codon
        and ending with a stop codon. Only take groups of 3 bps into account at a time and start at the specified
        frame_offset (0, 1 or 2 - although other values are allowed)
        e.g.:   "AATGGGGTAATTTATGTTTTAAGGG" would return ['ATGGGGTAA', 'ATGTTTTAA'] for frame offset 1
        thr aggregate_fn since both a and c have the max len of 3
        :returns a list open reading frames
        """
        orfs = []
        if not sequence:
            return orfs
        start_index = -1
        for i in range(frame_offset, len(sequence), 3):
            if sequence[i:i + 3] in self.start_codons:
               start_index = i
            if sequence[i:i + 3] in self.stop_codons and start_index > -1:
                orfs.append(sequence[start_index: i + 3])
                start_index = -1

        return orfs

    def repeats(self, sequence, length):
        """ Returns a dict of all sequences of the specified length with their repeat count. A repeat count of 2 means the sequence appears twice in total
        We do not return repeat counts of 1 as that is noisy and assumed for all sequences. Overlapping repeats are included
        e.g. searchig for repeats of length 2 in 'aaaa' yields {'aa': 3}
        :returns a dict of sequences with their repeat counts
        """
        repeats = {}
        if not sequence:
            return repeats
        for i in range(0, len(sequence) - (length - 1)):
            substring = sequence[i:i+length]
            if not substring in repeats:
                count = 1
                for j in range(i + 1, len(sequence) - (length - 1)):
                    if sequence[j:j+length] == substring:
                        count += 1
                if count > 1:
                    repeats[substring] = count
        return repeats





