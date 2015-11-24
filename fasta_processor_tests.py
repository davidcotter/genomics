import unittest
import fasta_processor

# Here's our "unit tests".
class TestMyFastaProcesssor(unittest.TestCase):

    the_processor = None
    def setUp(self):
        self.the_processor = fasta_processor.FastProcessor()
        self.the_processor.read_fasta_file("dna_sample.fasta")


    def test_read_in_fasta_file(self):
        self.failUnless(self.the_processor.fasta_records['gi|142022655|gb|EQ086233.1|43'].startswith('TCGG'))

    def test_count_records(self):
        self.failUnlessEqual(self.the_processor.record_count(), 25)

    def test_record_lengths(self):
        sorted_lengths = sorted(self.the_processor.record_lengths())
        self.failUnlessEqual(len(sorted_lengths), 25)
        self.failUnlessEqual(sorted_lengths[0], 512)

    def test_longest_sequences(self):
        self.failUnlessEqual(self.the_processor.max_len_sequences(), ['gi|142022655|gb|EQ086233.1|323'])

    def test_shortest_sequences(self):
        self.failUnlessEqual(self.the_processor.min_len_sequences(), ['gi|142022655|gb|EQ086233.1|521'])

    def test_open_reading_frames(self):
        self.failUnlessEqual(self.the_processor.open_reading_frames(None, 2), [])
        self.failUnlessEqual(self.the_processor.open_reading_frames("", 2), [])
        self.failUnlessEqual(self.the_processor.open_reading_frames("AAAAAAATGAAAAAAAAATAAGGGG", 0), ['ATGAAAAAAAAATAA'])
        self.failUnlessEqual(self.the_processor.open_reading_frames("AAAAATGAAAAAAAAATAAGGGG", 1), ['ATGAAAAAAAAATAA'])
        self.failUnlessEqual(self.the_processor.open_reading_frames("AAAAAATGAAAAAAAAATAAGGGG", 2), ['ATGAAAAAAAAATAA'])
        self.failUnlessEqual(self.the_processor.open_reading_frames("AAAAAAAAAAAAAA", 2), [])
        self.failUnlessEqual(self.the_processor.open_reading_frames("GGGATGAAAAAATAGTTTTTTATGGGGGGGGGGTAAAAAA"), ['ATGAAAAAATAG', 'ATGGGGGGGGGGTAA'])
        self.failUnlessEqual(self.the_processor.open_reading_frames("GGGATGAAAAAATAGTTTTTTATGGGGGGGGGGTAAAAAA", 1), [])
        self.failUnlessEqual(self.the_processor.open_reading_frames("GGGATGAAAAAATAGTTTTTTATGGGGGGGGGGTAAAAAA", 2), [])
        self.failUnlessEqual(self.the_processor.open_reading_frames("TAATAGTGA"), [])

    def test_repeats(self):
        self.failUnlessEqual(self.the_processor.repeats("aabbbaaaa", 2), {'aa': 4, 'bb': 2})

    def test_repeats_longest_sequence(self):
        repeats = self.the_processor.repeats(self.the_processor.fasta_records['gi|142022655|gb|EQ086233.1|521'], 6)
        self.failUnlessEqual(self.the_processor.keys_for_aggregate_val(repeats, max), ['CCGGGT'])

def main():
    unittest.main()

if __name__ == '__main__':
    main()