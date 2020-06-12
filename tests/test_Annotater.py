#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
from unittest import TestCase
import IOutilities
import Annotater
import pandas as ps

class TestFilter(TestCase):

    def test_givenVCF_dropDuplicates(self):
        vcfWithDuplicates = ps.DataFrame({
            '#chrom' : ['12', '12', '12'],
            'pos' : [2500, 2500, 2500],
            "id" : [ '.', '.', '.'],
            "ref" : [ 'c', 'c', 'c'],
            "alt" : [ 't', 't', 't'],
            'qual' : [ '.', '.', '.'],
            'filter' : [ '.', '.', '.'],
            'info' : ['.', '.', '.'],
        })
        tmpFile = "/tmp/vcf.tsv"
        vcfWithDuplicates.to_csv(tmpFile, sep='\t', index=False)
        Annotater.dropDuplicates(tmpFile)
        actualDf = ps.read_csv(open(tmpFile, 'r'), sep='\t')

        self.assertEquals(len(actualDf.index), 1)

    def test_givenIncorrectchroFormat_ReturnCorrectFormat(self):

        givenChr = "chr1"
        expectedChr = "1"

        actualChr = IOutilities.formatChromo(givenChr)

        self.assertEquals(expectedChr, actualChr)

    def test_GivenImproperInsertionFormat_When_formatImproperInsertionIsCalled_Return_NsuffixedAlleles(self):

        givenRef = '-'
        givenAlt = 'A'

        expectedRef = 'A'
        expectedAl = 'AA'

        actualAlleles = Annotater.formatImproperInserions(givenRef,givenAlt)

        self.assertEquals(expectedRef,actualAlleles[0])
        self.assertEquals(expectedAl, actualAlleles[1])
