import unittest
import filter
import pandas as ps


class TestFilter(unittest.TestCase):

    tmpLogLocation = "/tmp/"

    def test_GivenEmptyDF_Then_RowLensAreEqual(self):
        #Given
        expected = ps.DataFrame()
        #When
        actualRow = filter.run(expected, self.tmpLogLocation)
        #Then
        self.assertEqual(type(expected), type(actualRow))
        self.assertEqual(len(expected), len(actualRow))
        self.assertEqual(expected.get(0),actualRow.get(0))

    def test_GivenOneEMBLRowWithHeader_Then_ReturnRow(self):
        #Given
        data = ["4_1804915_A/G","4:1804915","G","ENSG00000068078","ENST00000260795","Transcript","3_prime_UTR_variant", \
               "1509","-","-","-","-","rs1466726466", \
               "IMPACT=MODIFIER;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]
        colNames = ["#Uploaded_variation","Location","Allele","Gene","Feature","Feature_type","Consequence",
                    "cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","Extra"]

        expectedRow = ps.DataFrame([data], columns=colNames)

        #When
        actualRow = filter.run(expectedRow, self.tmpLogLocation)

        #Then

        Row1EqualsRow0inActual = actualRow.equals(expectedRow[0:])

        self.assertTrue(Row1EqualsRow0inActual)

    def test_GivenOneCanonicalEMBLrowWithHeader_Then_returnRow(self):
        # Given
        data = ["4_1804915_A/G", "4:1804915", "G", "ENSG00000068078", "ENST00000260795", "Transcript",
                "3_prime_UTR_variant", \
                "1509", "-", "-", "-", "-", "rs1466726466", \
                "IMPACT=MODIFIER;STRAND=1;CANONICAL=YES;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]
        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation",
                    "Extra"]

        expectedRow = ps.DataFrame([data], columns=colNames)

        # When
        actualRow = filter.run(expectedRow, self.tmpLogLocation)

        # Then

        Row1EqualsRow0inActual = actualRow.equals(expectedRow[0:])
        self.assertTrue(Row1EqualsRow0inActual)

    def test_GivenOneNCBInonCanonical_Then_returnRow(self):
        data = ["4_1804915_A/G", "4:1804915", "G", "2261", "NR_148971.1", "Transcript",
                "3_prime_UTR_variant", \
                "1509", "-", "-", "-", "-", "rs1466726466", \
                "IMPACT=MODIFIER;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]
        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation",
                    "Extra"]
        expectedRow = ps.DataFrame([data], columns=colNames)

        # When
        actualRow = filter.run(expectedRow, self.tmpLogLocation)

        # Then
        Row1EqualsRow0inActual = actualRow.equals(expectedRow[0:])
        self.assertTrue(Row1EqualsRow0inActual)

    def test_GivenOneNCBICanonical_Then_returnRow(self):

        data = ["4_1804915_A/G", "4:1804915", "G", "2261", "NR_148971.1", "Transcript",
                "3_prime_UTR_variant", \
                "1509", "-", "-", "-", "-", "rs1466726466", \
                "IMPACT=MODIFIER;CANONICAL=YES;STRAND=1;BIOTYPE=lncRNA;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]
        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]
        expectedRow = ps.DataFrame([data], columns=colNames)
        # When
        actualRow = filter.run(expectedRow, self.tmpLogLocation)

        # Then
        Row1EqualsRow0inActual = actualRow.equals(expectedRow[0:])
        self.assertTrue(Row1EqualsRow0inActual)

    def test_GivenOneNCBIAndOneEMBLrow_Then_returnBothRows(self):
        #Given
        data = ["4_1804915_A/G", "4:1804915", "G", "ENSG00000068078", "ENST00000260795", "Transcript",
                "3_prime_UTR_variant", \
                "1509", "-", "-", "-", "-", "rs1466726466", \
                "IMPACT=MODIFIER;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]

        data1 = ["4_1804915_A/G", "4:1804915", "G", "2261", "NR_148971.1", "Transcript",
                "3_prime_UTR_variant", \
                "1509", "-", "-", "-", "-", "rs1466726466", \
                "IMPACT=MODIFIER;STRAND=1;BIOTYPE=lncRNA;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]
        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]
        EMBLdf = ps.DataFrame([data], columns=colNames)
        NCBIdf = ps.DataFrame([data1], columns=colNames)
        expectedRows = ps.concat([EMBLdf,NCBIdf])

        # When
        actualRows = filter.run(expectedRows, self.tmpLogLocation)

        # Then
        Row1EqualsRow0inActual = actualRows.equals(expectedRows)

        self.assertTrue(Row1EqualsRow0inActual)

    def test_GivenOneNCBIAndOneEMBLrowWithMissingData_Then_returnBothRows(self):
        #Given
        data = ["4_1804915_A/G", "4:1804915", "G", "ENSG00000068078", "ENST00000260795", "Transcript",
                "3_prime_UTR_variant", \
                "1509", "-", "-", "-", "-", "rs1466726466", \
                "IMPACT=MODIFIER;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]

        data1 = ["4_1804915_A/G", "4:1804915", "G", "2261", "NR_148971.1", "Transcript",
                "3_prime_UTR_variant", \
                "1509", "-", "-", "-", "-", "rs1466726466", \
                "IMPACT=MODIFIER;CANONICAL=YES;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]
        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]
        EMBLdf = ps.DataFrame([data], columns=colNames)
        NCBIdf = ps.DataFrame([data1], columns=colNames)
        expectedRows = ps.concat([EMBLdf,NCBIdf])

        # When
        actualRows = filter.run(expectedRows, self.tmpLogLocation)

        # Then
        Row1EqualsRow0inActual = actualRows.equals(expectedRows)

        self.assertTrue(Row1EqualsRow0inActual)



    def test_GivenOneNCBIAndEMBLCanonicals_Then_ReturnBothRows(self):
        # Given
        data = ["4_1804915_A/G", "4:1804915", "G", "ENSG00000068078", "ENST00000260795", "Transcript",
                "3_prime_UTR_variant", \
                "1509", "-", "-", "-", "-", "rs1466726466", \
                "IMPACT=MODIFIER;STRAND=1;CANONICAL=YES;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]

        data1 = ["4_1804915_A/G", "4:1804915", "G", "2261", "NR_148971.1", "Transcript",
                 "3_prime_UTR_variant", \
                 "1509", "-", "-", "-", "-", "rs1466726466", \
                 "IMPACT=MODIFIER;STRAND=1;CANONICAL=YES;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]
        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]
        EMBLdf = ps.DataFrame([data], columns=colNames)
        NCBIdf = ps.DataFrame([data1], columns=colNames)
        expectedRows = ps.concat([EMBLdf, NCBIdf])

        # When
        actualRows = filter.run(expectedRows, self.tmpLogLocation)

        # Then
        Row1EqualsRow0inActual = actualRows.equals(expectedRows)

        self.assertTrue(Row1EqualsRow0inActual)

    def test_GivenOneCanonicalAndOneNonCanonical_Then_ReturnBothrows(self):
        # Given
        data = ["4_1804915_A/G", "4:1804915", "G", "ENSG00000068078", "ENST00000260795", "Transcript",
                "3_prime_UTR_variant", \
                "1509", "-", "-", "-", "-", "rs1466726466", \
                "IMPACT=MODIFIER;STRAND=1;CANONICAL=YES;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]

        data1 = ["4_1804915_A/G", "4:1804915", "G", "2261", "NR_148971.1", "Transcript",
                 "3_prime_UTR_variant", \
                 "1509", "-", "-", "-", "-", "rs1466726466", \
                 "IMPACT=MODIFIER;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]
        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]
        EMBLdf = ps.DataFrame([data], columns=colNames)
        ncNCBIdf = ps.DataFrame([data1], columns=colNames)
        expectedRows = ps.concat([EMBLdf, ncNCBIdf])

        # When
        actualRows = filter.run(expectedRows, self.tmpLogLocation)

        # Then
        Row1EqualsRow0inActual = actualRows.equals(expectedRows)

        self.assertTrue(Row1EqualsRow0inActual)

    def test_GivenNonCanonAndCanonical_Then_ReturnBothrows(self):
        # Given
        data = ["4_1804915_A/G", "4:1804915", "G", "ENSG00000068078", "ENST00000260795", "Transcript",
                "3_prime_UTR_variant", \
                "1509", "-", "-", "-", "-", "rs1466726466", \
                "IMPACT=MODIFIER;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]

        data1 = ["4_1804915_A/G", "4:1804915", "G", "2261", "NR_148971.1", "Transcript",
                 "3_prime_UTR_variant", \
                 "1509", "-", "-", "-", "-", "rs1466726466", \
                 "IMPACT=MODIFIER;CANONICAL=YES;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]
        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]
        EMBLdf = ps.DataFrame([data], columns=colNames)
        ncNCBIdf = ps.DataFrame([data1], columns=colNames)
        expectedRows = ps.concat([EMBLdf, ncNCBIdf])

        # When
        actualRows = filter.run(expectedRows, self.tmpLogLocation)

        # Then
        Row1EqualsRow0inActual = actualRows.equals(expectedRows)

        self.assertTrue(Row1EqualsRow0inActual)

    def test_GivenMultipleNonCanonicalRowsInEMBLWithOneProtienCoding_ReturnProtienCoding(self):

        protienCodingRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000346085", "Transcript",
         "frameshift_variant", "2371-2372", "1445-1446", "482", "A/AX", "gcg/gcGg",
         "rs1554248082,TMP_ESP_6_157100260_157100259",
         "IMPACT=HIGH;STRAND=1;VARIANT_CLASS=insertion;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=1;APPRIS=A2;ENSP=ENSP00000344546;TREMBL=A0A3F2YNW7;UNIPARC=UPI000E6F3DF6;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=2/21;DOMAINS=Low_complexity_(Seg):seg,PANTHER:PTHR12656,PANTHER:PTHR12656:SF11;HGVSc=ENST00000346085.10:c.1451dup;HGVSp=ENSP00000344546.5:p.Phe485LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        miscRNARow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000350026", "Transcript",
         "upstream_gene_variant", "1196-1197", "1196-1197", "399", "A/AX", "gcg/gcGg",
         "rs1554248082,TMP_ESP_6_157100260_157100259",
         "IMPACT=HIGH;STRAND=1;VARIANT_CLASS=insertion;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=lncRNA;TSL=1;APPRIS=P3;CCDS=CCDS5251.2;ENSP=ENSP00000055163;SWISSPROT=Q8NFD5;UNIPARC=UPI000058E2EA;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=1/19;DOMAINS=PANTHER:PTHR12656,PANTHER:PTHR12656:SF11,Low_complexity_(Seg):seg;HGVSc=ENST00000350026.10:c.1202dup;HGVSp=ENSP00000055163.7:p.Phe402LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        lncRNArow =  ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000414678", "Transcript",
         "downstream_gene_variant", "-", "-", "-", "-", "-", "rs1554248082,TMP_ESP_6_157100260_157100259",
         "IMPACT=MODIFIER;DISTANCE=307;STRAND=1;FLAGS=cds_start_NF;VARIANT_CLASS=insertion;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=miRNA;TSL=5;ENSP=ENSP00000412835;TREMBL=H0Y7H8;UNIPARC=UPI0001D3BCFD;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]

        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]

        protienCodingDF = ps.DataFrame([protienCodingRow], columns=colNames)
        miscDF = ps.DataFrame([miscRNARow], columns=colNames)
        lncDF = ps.DataFrame([lncRNArow], columns=colNames)
        inputDF = ps.concat([protienCodingDF,miscDF,lncDF])

        # When
        actualRow = filter.run(inputDF, self.tmpLogLocation)

        actualRowEqualsExpected = protienCodingDF.equals(actualRow)

        self.assertTrue(actualRowEqualsExpected)

    def test_GivenMultipleCanonicalRowsInEMBL_ReturnCanonical(self):
        canonicalRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000346085",
                            "Transcript",
                            "frameshift_variant", "2371-2372", "1445-1446", "482", "A/AX", "gcg/gcGg",
                            "rs1554248082,TMP_ESP_6_157100260_157100259",
                            "IMPACT=HIGH;STRAND=1;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=1;APPRIS=A2;ENSP=ENSP00000344546;TREMBL=A0A3F2YNW7;UNIPARC=UPI000E6F3DF6;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=2/21;DOMAINS=Low_complexity_(Seg):seg,PANTHER:PTHR12656,PANTHER:PTHR12656:SF11;HGVSc=ENST00000346085.10:c.1451dup;HGVSp=ENSP00000344546.5:p.Phe485LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        miscRNARow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000350026",
                      "Transcript",
                      "upstream_gene_variant", "1196-1197", "1196-1197", "399", "A/AX", "gcg/gcGg",
                      "rs1554248082,TMP_ESP_6_157100260_157100259",
                      "IMPACT=HIGH;STRAND=1;VARIANT_CLASS=insertion;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=lncRNA;TSL=1;APPRIS=P3;CCDS=CCDS5251.2;ENSP=ENSP00000055163;SWISSPROT=Q8NFD5;UNIPARC=UPI000058E2EA;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=1/19;DOMAINS=PANTHER:PTHR12656,PANTHER:PTHR12656:SF11,Low_complexity_(Seg):seg;HGVSc=ENST00000350026.10:c.1202dup;HGVSp=ENSP00000055163.7:p.Phe402LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        lncRNArow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000414678",
                     "Transcript",
                     "downstream_gene_variant", "-", "-", "-", "-", "-", "rs1554248082,TMP_ESP_6_157100260_157100259",
                     "IMPACT=MODIFIER;DISTANCE=307;STRAND=1;FLAGS=cds_start_NF;VARIANT_CLASS=insertion;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=miRNA;TSL=5;ENSP=ENSP00000412835;TREMBL=H0Y7H8;UNIPARC=UPI0001D3BCFD;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]

        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]

        protienCodingDF = ps.DataFrame([canonicalRow], columns=colNames)
        miscDF = ps.DataFrame([miscRNARow], columns=colNames)
        lncDF = ps.DataFrame([lncRNArow], columns=colNames)
        inputDF = ps.concat([protienCodingDF, miscDF, lncDF])

        # When
        actualRow = filter.run(inputDF, self.tmpLogLocation)

        actualRowEqualsExpected = protienCodingDF.equals(actualRow)

        self.assertTrue(actualRowEqualsExpected)

    def test_GivenMultipleCanonicalRowsInEMBLandOneProteinCoding_ReturnProteinCoding(self):
        canonicalRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000346085",
                        "Transcript",
                        "frameshift_variant", "2371-2372", "1445-1446", "482", "A/AX", "gcg/gcGg",
                        "rs1554248082,TMP_ESP_6_157100260_157100259",
                        "IMPACT=HIGH;STRAND=1;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=1;APPRIS=A2;ENSP=ENSP00000344546;TREMBL=A0A3F2YNW7;UNIPARC=UPI000E6F3DF6;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=2/21;DOMAINS=Low_complexity_(Seg):seg,PANTHER:PTHR12656,PANTHER:PTHR12656:SF11;HGVSc=ENST00000346085.10:c.1451dup;HGVSp=ENSP00000344546.5:p.Phe485LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        miscRNARow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000350026",
                      "Transcript",
                      "upstream_gene_variant", "1196-1197", "1196-1197", "399", "A/AX", "gcg/gcGg",
                      "rs1554248082,TMP_ESP_6_157100260_157100259",
                      "IMPACT=HIGH;STRAND=1;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=lncRNA;TSL=1;APPRIS=P3;CCDS=CCDS5251.2;ENSP=ENSP00000055163;SWISSPROT=Q8NFD5;UNIPARC=UPI000058E2EA;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=1/19;DOMAINS=PANTHER:PTHR12656,PANTHER:PTHR12656:SF11,Low_complexity_(Seg):seg;HGVSc=ENST00000350026.10:c.1202dup;HGVSp=ENSP00000055163.7:p.Phe402LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        lncRNArow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000414678",
                     "Transcript",
                     "downstream_gene_variant", "-", "-", "-", "-", "-", "rs1554248082,TMP_ESP_6_157100260_157100259",
                     "IMPACT=MODIFIER;DISTANCE=307;STRAND=1;FLAGS=cds_start_NF;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=miRNA;TSL=5;ENSP=ENSP00000412835;TREMBL=H0Y7H8;UNIPARC=UPI0001D3BCFD;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]

        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]

        protienCodingDF = ps.DataFrame([canonicalRow], columns=colNames)
        miscDF = ps.DataFrame([miscRNARow], columns=colNames)
        lncDF = ps.DataFrame([lncRNArow], columns=colNames)
        inputDF = ps.concat([protienCodingDF, miscDF, lncDF])

        # When
        actualRow = filter.run(inputDF, self.tmpLogLocation)

        actualRowEqualsExpected = protienCodingDF.equals(actualRow)

        self.assertTrue(actualRowEqualsExpected)

    def test_GivenAllProtienCodingRows_ReturnHighestImpactRow(self):
        moderateImpactRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000346085",
                        "Transcript",
                        "frameshift_variant", "2371-2372", "1445-1446", "482", "A/AX", "gcg/gcGg",
                        "rs1554248082,TMP_ESP_6_157100260_157100259",
                        "IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=1;APPRIS=A2;ENSP=ENSP00000344546;TREMBL=A0A3F2YNW7;UNIPARC=UPI000E6F3DF6;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=2/21;DOMAINS=Low_complexity_(Seg):seg,PANTHER:PTHR12656,PANTHER:PTHR12656:SF11;HGVSc=ENST00000346085.10:c.1451dup;HGVSp=ENSP00000344546.5:p.Phe485LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        modifierImpactRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000350026",
                      "Transcript",
                      "upstream_gene_variant", "1196-1197", "1196-1197", "399", "A/AX", "gcg/gcGg",
                      "rs1554248082,TMP_ESP_6_157100260_157100259",
                      "IMPACT=MODIFIER;STRAND=1;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=1;APPRIS=P3;CCDS=CCDS5251.2;ENSP=ENSP00000055163;SWISSPROT=Q8NFD5;UNIPARC=UPI000058E2EA;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=1/19;DOMAINS=PANTHER:PTHR12656,PANTHER:PTHR12656:SF11,Low_complexity_(Seg):seg;HGVSc=ENST00000350026.10:c.1202dup;HGVSp=ENSP00000055163.7:p.Phe402LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        lowImpactRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000414678",
                     "Transcript",
                     "downstream_gene_variant", "-", "-", "-", "-", "-", "rs1554248082,TMP_ESP_6_157100260_157100259",
                     "IMPACT=LOW;DISTANCE=307;STRAND=1;FLAGS=cds_start_NF;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=5;ENSP=ENSP00000412835;TREMBL=H0Y7H8;UNIPARC=UPI0001D3BCFD;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]

        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]
        expectedRow = ps.DataFrame([moderateImpactRow], columns=colNames)
        modifierImpact = ps.DataFrame([modifierImpactRow], columns=colNames)
        lowImpact = ps.DataFrame([lowImpactRow], columns=colNames)
        inputDF = ps.concat([lowImpact,expectedRow,modifierImpact])

        # When
        actualRow = filter.run(inputDF, self.tmpLogLocation)

        ps.testing.assert_frame_equal(actualRow,expectedRow)

    def test_GivenTieinScoring_ReturnFirstAnnotation(self):
        moderateImpactRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000346085",
                             "Transcript",
                             "frameshift_variant", "2371-2372", "1445-1446", "482", "A/AX", "gcg/gcGg",
                             "rs1554248082,TMP_ESP_6_157100260_157100259",
                             "IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=1;APPRIS=A2;ENSP=ENSP00000344546;TREMBL=A0A3F2YNW7;UNIPARC=UPI000E6F3DF6;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=2/21;DOMAINS=Low_complexity_(Seg):seg,PANTHER:PTHR12656,PANTHER:PTHR12656:SF11;HGVSc=ENST00000346085.10:c.1451dup;HGVSp=ENSP00000344546.5:p.Phe485LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        secondRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000350026",
                             "Transcript",
                             "upstream_gene_variant", "1196-1197", "1196-1197", "399", "A/AX", "gcg/gcGg",
                             "rs1554248082,TMP_ESP_6_157100260_157100259",
                             "IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=1;APPRIS=P3;CCDS=CCDS5251.2;ENSP=ENSP00000055163;SWISSPROT=Q8NFD5;UNIPARC=UPI000058E2EA;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=1/19;DOMAINS=PANTHER:PTHR12656,PANTHER:PTHR12656:SF11,Low_complexity_(Seg):seg;HGVSc=ENST00000350026.10:c.1202dup;HGVSp=ENSP00000055163.7:p.Phe402LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        finalRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000414678",
                        "Transcript",
                        "downstream_gene_variant", "-", "-", "-", "-", "-",
                        "rs1554248082,TMP_ESP_6_157100260_157100259",
                        "IMPACT=MODERATE;DISTANCE=307;STRAND=1;FLAGS=cds_start_NF;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=5;ENSP=ENSP00000412835;TREMBL=H0Y7H8;UNIPARC=UPI0001D3BCFD;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]

        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]
        expectedRow = ps.DataFrame([moderateImpactRow], columns=colNames)
        secondDF = ps.DataFrame([secondRow], columns=colNames)
        thirdDF = ps.DataFrame([finalRow], columns=colNames)
        inputDF = ps.concat([expectedRow, secondDF, thirdDF])

        # When
        actualRow = filter.run(inputDF, self.tmpLogLocation)

        ps.testing.assert_frame_equal(actualRow, expectedRow)

    def test_GivenMultipleNCBItranscriptsAndOneEMBL_AllCanon_Then_ReturnNCBIwithMatchingSymbol(self):
        moderateImpactRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000346085",
                             "Transcript",
                             "frameshift_variant", "2371-2372", "1445-1446", "482", "A/AX", "gcg/gcGg",
                             "rs1554248082,TMP_ESP_6_157100260_157100259",
                             "IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=1;APPRIS=A2;ENSP=ENSP00000344546;TREMBL=A0A3F2YNW7;UNIPARC=UPI000E6F3DF6;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=2/21;DOMAINS=Low_complexity_(Seg):seg,PANTHER:PTHR12656,PANTHER:PTHR12656:SF11;HGVSc=ENST00000346085.10:c.1451dup;HGVSp=ENSP00000344546.5:p.Phe485LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        secondRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000350026",
                     "Transcript",
                     "upstream_gene_variant", "1196-1197", "1196-1197", "399", "A/AX", "gcg/gcGg",
                     "rs1554248082,TMP_ESP_6_157100260_157100259",
                     "IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=1;APPRIS=P3;CCDS=CCDS5251.2;ENSP=ENSP00000055163;SWISSPROT=Q8NFD5;UNIPARC=UPI000058E2EA;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=1/19;DOMAINS=PANTHER:PTHR12656,PANTHER:PTHR12656:SF11,Low_complexity_(Seg):seg;HGVSc=ENST00000350026.10:c.1202dup;HGVSp=ENSP00000055163.7:p.Phe402LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        finalRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000414678",
                    "Transcript",
                    "downstream_gene_variant", "-", "-", "-", "-", "-",
                    "rs1554248082,TMP_ESP_6_157100260_157100259",
                    "IMPACT=MODERATE;DISTANCE=307;STRAND=1;FLAGS=cds_start_NF;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=5;ENSP=ENSP00000412835;TREMBL=H0Y7H8;UNIPARC=UPI0001D3BCFD;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        NCBIrow = ["4_1804915_A/G", "4:1804915", "G", "2261", "NR_148971.1", "Transcript",
                 "3_prime_UTR_variant", \
                 "1509", "-", "-", "-", "-", "rs1466726466", \
                 "IMPACT=MODIFIER;CANONICAL=YES;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]


        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]

        expectedRow = ps.DataFrame([moderateImpactRow], columns=colNames)
        secondDF = ps.DataFrame([secondRow], columns=colNames)
        thirdDF = ps.DataFrame([finalRow], columns=colNames)
        secondExpectedRow = ps.DataFrame([NCBIrow], columns=colNames)
        inputDF = ps.concat([expectedRow, secondDF, thirdDF, secondExpectedRow])

        expectedDF = ps.concat([expectedRow,secondExpectedRow])
        expectedDF = expectedDF.reset_index(drop=True)

        # When
        actualDF = filter.run(inputDF, self.tmpLogLocation)

        ps.testing.assert_frame_equal(actualDF, expectedDF)

    def test_GivenMultipleNCBIAndEMBLAllCanonWithOnlyTwoMatchingSymbols_Then_ReturnNCBIwithMatchingSymbolIrrespectiveOfRow(self):
        moderateImpactRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000346085",
                             "Transcript",
                             "frameshift_variant", "2371-2372", "1445-1446", "482", "A/AX", "gcg/gcGg",
                             "rs1554248082,TMP_ESP_6_157100260_157100259",
                             "IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=1;APPRIS=A2;ENSP=ENSP00000344546;TREMBL=A0A3F2YNW7;UNIPARC=UPI000E6F3DF6;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=2/21;DOMAINS=Low_complexity_(Seg):seg,PANTHER:PTHR12656,PANTHER:PTHR12656:SF11;HGVSc=ENST00000346085.10:c.1451dup;HGVSp=ENSP00000344546.5:p.Phe485LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        secondRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000350026",
                     "Transcript",
                     "upstream_gene_variant", "1196-1197", "1196-1197", "399", "A/AX", "gcg/gcGg",
                     "rs1554248082,TMP_ESP_6_157100260_157100259",
                     "IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=1;APPRIS=P3;CCDS=CCDS5251.2;ENSP=ENSP00000055163;SWISSPROT=Q8NFD5;UNIPARC=UPI000058E2EA;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=1/19;DOMAINS=PANTHER:PTHR12656,PANTHER:PTHR12656:SF11,Low_complexity_(Seg):seg;HGVSc=ENST00000350026.10:c.1202dup;HGVSp=ENSP00000055163.7:p.Phe402LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        finalRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000414678",
                    "Transcript",
                    "downstream_gene_variant", "-", "-", "-", "-", "-",
                    "rs1554248082,TMP_ESP_6_157100260_157100259",
                    "IMPACT=MODERATE;DISTANCE=307;STRAND=1;FLAGS=cds_start_NF;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=5;ENSP=ENSP00000412835;TREMBL=H0Y7H8;UNIPARC=UPI0001D3BCFD;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        NCBIrow = ["6_156779126_-/G", "6:156779125-156779126", "G", "2261", "NR_148971.1", "Transcript",
                   "3_prime_UTR_variant", \
                   "1509", "-", "-", "-", "-", "rs1466726466", \
                   "IMPACT=MODIFIER;CANONICAL=YES;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC;BIOTYPE=protein_coding "]

        NCBIrow2 = ["6_156779126_-/G", "6:156779125-156779126", "G", "2261", "NR_148971.1", "Transcript",
                 "3_prime_UTR_variant", \
                 "1509", "-", "-", "-", "-", "rs1466726466", \
                 "IMPACT=MODIFIER;CANONICAL=YES;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=NONMATCHING;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC;BIOTYPE=protein_coding "]

        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]

        expectedRow = ps.DataFrame([moderateImpactRow], columns=colNames)
        secondDF = ps.DataFrame([secondRow], columns=colNames)
        thirdDF = ps.DataFrame([finalRow], columns=colNames)
        secondExpectedRow = ps.DataFrame([NCBIrow], columns=colNames)
        fourthRow = ps.DataFrame([NCBIrow2], columns=colNames)
        inputDF = ps.concat([expectedRow, secondDF, thirdDF, fourthRow, secondExpectedRow])

        expectedDF = ps.concat([expectedRow, secondExpectedRow]).reset_index(drop=True)

        # When
        actualDF = filter.run(inputDF, self.tmpLogLocation)

        ps.testing.assert_frame_equal(actualDF, expectedDF)

    def test_GivenMultipleNCBIAndEMBLAllcanonWithTwoMatchingBiotypesAndNoOtherMatchingInformation_Then_ReturnNCBIwithMatchingBiotypes(
            self):
        moderateImpactRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000346085",
                             "Transcript",
                             "frameshift_variant", "2371-2372", "1445-1446", "482", "A/AX", "gcg/gcGg",
                             "rs1554248082,TMP_ESP_6_157100260_157100259",
                             "IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=TEST1;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=MATCH;TSL=1;APPRIS=A2;ENSP=ENSP00000344546;TREMBL=A0A3F2YNW7;UNIPARC=UPI000E6F3DF6;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=2/21;DOMAINS=Low_complexity_(Seg):seg,PANTHER:PTHR12656,PANTHER:PTHR12656:SF11;HGVSc=ENST00000346085.10:c.1451dup;HGVSp=ENSP00000344546.5:p.Phe485LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        secondRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000350026",
                     "Transcript",
                     "upstream_gene_variant", "1196-1197", "1196-1197", "399", "A/AX", "gcg/gcGg",
                     "rs1554248082,TMP_ESP_6_157100260_157100259",
                     "IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=TEST2;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=NONMATCH1;TSL=1;APPRIS=P3;CCDS=CCDS5251.2;ENSP=ENSP00000055163;SWISSPROT=Q8NFD5;UNIPARC=UPI000058E2EA;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=1/19;DOMAINS=PANTHER:PTHR12656,PANTHER:PTHR12656:SF11,Low_complexity_(Seg):seg;HGVSc=ENST00000350026.10:c.1202dup;HGVSp=ENSP00000055163.7:p.Phe402LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        finalRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000414678",
                    "Transcript",
                    "downstream_gene_variant", "-", "-", "-", "-", "-",
                    "rs1554248082,TMP_ESP_6_157100260_157100259",
                    "IMPACT=MODERATE;DISTANCE=307;STRAND=1;FLAGS=cds_start_NF;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=TEST3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=NONMATCH2;TSL=5;ENSP=ENSP00000412835;TREMBL=H0Y7H8;UNIPARC=UPI0001D3BCFD;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        NCBIrow = ["6_156779126_-/G", "6:156779125-156779126", "G", "2261", "NR_148971.1", "Transcript",
                   "3_prime_UTR_variant", \
                   "1509", "-", "-", "-", "-", "rs1466726466", \
                   "IMPACT=MODIFIER;CANONICAL=YES;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=TEST4;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC;BIOTYPE=MATCH"]

        NCBIrow2 = ["6_156779126_-/G", "6:156779125-156779126", "G", "2261", "NR_148971.1", "Transcript",
                    "3_prime_UTR_variant", \
                    "1509", "-", "-", "-", "-", "rs1466726466", \
                    "IMPACT=MODIFIER;CANONICAL=YES;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=TEST5;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC;BIOTYPE=protein_coding"]

        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]

        expectedRow = ps.DataFrame([moderateImpactRow], columns=colNames)
        secondDF = ps.DataFrame([secondRow], columns=colNames)
        thirdDF = ps.DataFrame([finalRow], columns=colNames)
        secondExpectedRow = ps.DataFrame([NCBIrow], columns=colNames)
        fourthRow = ps.DataFrame([NCBIrow2], columns=colNames)
        inputDF = ps.concat([expectedRow, secondDF, thirdDF, fourthRow, secondExpectedRow])

        expectedDF = ps.concat([expectedRow, secondExpectedRow]).reset_index(drop=True)

        # When
        actualDF = filter.run(inputDF, self.tmpLogLocation)

        ps.testing.assert_frame_equal(actualDF, expectedDF)