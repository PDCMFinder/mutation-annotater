from unittest import TestCase
import annotationMerger
import filter
import pandas as ps


class TestFilter(TestCase):

    def test_GivenChrPosKeyWithoutMatchingKeys_WhenChrPosKeyIsused_Then_ReturnChrPosKey(self):

        chromo = '1'
        seq_start = '10000'
        refAllele = 'A'
        altAllele = 'T'

        expectedChrPosKey = "{0}:{1}".format(chromo,seq_start)

        rowDic = {'chromosome': chromo, 'seq_start_position' : seq_start, 'ref_allele' : refAllele, 'alt_allele' : altAllele }
        row = ps.Series(rowDic)

        actualChrPos = annotationMerger.formatChrPosKey(row)

        self.assertEqual(expectedChrPosKey,actualChrPos)

    def test_GivenChrPosKeyWithInsertionWithMatchingFirstNucleotides_WhenChrPosKeyIsused_Then_ReturnShiftedAndHyphonatedResult(self):
        chromo = '1'
        seq_start = '10000'
        seq_shift = '10001'
        refAllele = 'A'
        altAllele = 'ATT'

        expectedChrPosKey = "{0}:{1}-{2}".format(chromo, seq_start,seq_shift)

        rowDic = {'chromosome': chromo, 'seq_start_position': seq_start, 'ref_allele': refAllele,
                  'alt_allele': altAllele}
        row = ps.Series(rowDic)

        actualChrPos = annotationMerger.formatChrPosKey(row)

        self.assertEqual(expectedChrPosKey, actualChrPos)

    def test_GivenChrPosKeyWithDeletionWithMatchingFirstNucleotides_WhenChrPosKeyIsused_Then_ReturnIncrementposBy1(self):
        chromo = '1'
        seq_start = '10000'
        seq_shift = '10001'
        refAllele = 'AT'
        altAllele = 'A'

        expectedChrPosKey = "{0}:{1}".format(chromo, seq_shift)

        rowDic = {'chromosome': chromo, 'seq_start_position': seq_start, 'ref_allele': refAllele,
                  'alt_allele': altAllele}
        row = ps.Series(rowDic)

        actualChrPos = annotationMerger.formatChrPosKey(row)

        self.assertEqual(expectedChrPosKey, actualChrPos)