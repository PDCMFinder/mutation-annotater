#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
from unittest import TestCase
import IOutilities
import Annotater

class TestFilter(TestCase):


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
