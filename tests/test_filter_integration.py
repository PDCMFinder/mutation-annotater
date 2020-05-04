import unittest
import AnnotationFilter
import pandas as ps


class TestFilter(unittest.TestCase):

    tmpLogLocation = "/tmp/"
    tmpFileName = "filterTest"

    def test_GivenEmptyDF_Then_RowLensAreEqual(self):

        #Given
        expected = ps.DataFrame()
        #When
        actualRow = AnnotationFilter.run(expected, self.tmpFileName, self.tmpLogLocation)
        #Then
        self.assertEqual(type(expected), type(actualRow))
        self.assertEqual(len(expected), len(actualRow))
        self.assertEqual(expected.get(0), actualRow.get(0))

    def test_GivenOneEMBLRowWithHeader_Then_ReturnRow(self):
        data = ["4_1804915_A/G","4:1804915","G","ENSG00000068078","ENST00000260795","Transcript","3_prime_UTR_variant", \
               "1509","-","-","-","-","rs1466726466", \
               "IMPACT=MODIFIER;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]
        colNames = ["#Uploaded_variation","Location","Allele","Gene","Feature","Feature_type","Consequence",
                    "cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","Extra"]

        expectedRow = ps.DataFrame([data], columns=colNames)
        actualRow = AnnotationFilter.run(expectedRow, self.tmpFileName, self.tmpLogLocation)
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
        actualRow = AnnotationFilter.run(expectedRow, self.tmpFileName, self.tmpLogLocation)

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
        actualRow = AnnotationFilter.run(expectedRow, self.tmpFileName, self.tmpLogLocation)

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
        actualRow = AnnotationFilter.run(expectedRow, self.tmpFileName, self.tmpLogLocation)

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
        actualRows = AnnotationFilter.run(expectedRows, self.tmpFileName, self.tmpLogLocation)

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
        actualRows = AnnotationFilter.run(expectedRows, self.tmpFileName, self.tmpLogLocation)

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
        actualRows = AnnotationFilter.run(expectedRows, self.tmpFileName, self.tmpLogLocation)

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
        actualRows = AnnotationFilter.run(expectedRows, self.tmpFileName, self.tmpLogLocation)

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
        actualRows = AnnotationFilter.run(expectedRows, self.tmpFileName, self.tmpLogLocation)

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


        actualRow = AnnotationFilter.run(inputDF, self.tmpFileName, self.tmpLogLocation)

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


        actualRow = AnnotationFilter.run(inputDF, self.tmpFileName, self.tmpLogLocation)

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


        actualRow = AnnotationFilter.run(inputDF, self.tmpFileName, self.tmpLogLocation)

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

        actualRow = AnnotationFilter.run(inputDF, self.tmpFileName, self.tmpLogLocation)

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

        actualRow = AnnotationFilter.run(inputDF, self.tmpFileName, self.tmpLogLocation)

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

        actualDF = AnnotationFilter.run(inputDF, self.tmpFileName, self.tmpLogLocation)

        ps.testing.assert_frame_equal(actualDF, expectedDF)

    def test_GivenMultipleNCBIAndEMBLAllCanonWithOnlyTwoMatchingSymbols_Then_ReturnNCBIwithMatchingSymbolIrrespectiveOfRow(self):
        moderateImpactRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000346085",
                             "Transcript",
                             "frameshift_variant", "2371-2372", "1445-1446", "482", "A/AX", "gcg/gcGg",
                             "rs1554248082,TMP_ESP_6_157100260_157100259",
                             "IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=MATCH;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=1;APPRIS=A2;ENSP=ENSP00000344546;TREMBL=A0A3F2YNW7;UNIPARC=UPI000E6F3DF6;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=2/21;DOMAINS=Low_complexity_(Seg):seg,PANTHER:PTHR12656,PANTHER:PTHR12656:SF11;HGVSc=ENST00000346085.10:c.1451dup;HGVSp=ENSP00000344546.5:p.Phe485LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        secondRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000350026",
                     "Transcript",
                     "upstream_gene_variant", "1196-1197", "1196-1197", "399", "A/AX", "gcg/gcGg",
                     "rs1554248082,TMP_ESP_6_157100260_157100259",
                     "IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=NONMATCH;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=1;APPRIS=P3;CCDS=CCDS5251.2;ENSP=ENSP00000055163;SWISSPROT=Q8NFD5;UNIPARC=UPI000058E2EA;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;EXON=1/19;DOMAINS=PANTHER:PTHR12656,PANTHER:PTHR12656:SF11,Low_complexity_(Seg):seg;HGVSc=ENST00000350026.10:c.1202dup;HGVSp=ENSP00000055163.7:p.Phe402LeufsTer133;HGVS_OFFSET=6;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        finalRow = ["6_156779126_-/G", "6:156779125-156779126", "G", "ENSG00000049618", "ENST00000414678",
                    "Transcript",
                    "downstream_gene_variant", "-", "-", "-", "-", "-",
                    "rs1554248082,TMP_ESP_6_157100260_157100259",
                    "IMPACT=MODERATE;DISTANCE=307;STRAND=1;FLAGS=cds_start_NF;VARIANT_CLASS=insertion;CANONICAL=YES;SYMBOL=ARID1B;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:18040;BIOTYPE=protein_coding;TSL=5;ENSP=ENSP00000412835;TREMBL=H0Y7H8;UNIPARC=UPI0001D3BCFD;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-;GENE_PHENO=1;AA_AF=0.005795;EA_AF=0.006561;MAX_AF=0.006561;MAX_AF_POPS=EA,EA;CLIN_SIG=pathogenic;PHENO=1,0"]
        NCBIrow = ["6_156779126_-/G", "6:156779125-156779126", "G", "2261", "NR_148971.1", "Transcript",
                   "3_prime_UTR_variant", \
                   "1509", "-", "-", "-", "-", "rs1466726466", \
                   "IMPACT=MODIFIER;CANONICAL=YES;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=MATCH;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC;BIOTYPE=protein_coding "]

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
        actualDF = AnnotationFilter.run(inputDF, self.tmpFileName, self.tmpLogLocation)

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
                    "IMPACT=MODIFIER;CANONICAL=YES;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=TEST5;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC;BIOTYPE=MISSENSE"]

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

        actualDF = AnnotationFilter.run(inputDF, self.tmpFileName, self.tmpLogLocation)
        ps.testing.assert_frame_equal(actualDF, expectedDF)

    def test_RealLifeError_When_givenToFilter_DoNotReturnCloneBasedSymbol(self):

        actualFailingData = [ [ "X_83508838_C/T","X:83508838","T","ENSG00000279437","ENST00000625081","Transcript","non_coding_transcript_exon_variant","377","-","-","-","-","-","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=Z82170.1;SYMBOL_SOURCE=Clone_based_ensembl_gene;BIOTYPE=TEC;CANONICAL=YES;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;EXON=1/1;HGVSc=ENST00000625081.1:n.377G>A" ],
[ "X_83508838_C/T","X:83508838","T","ENSG00000196767","ENST00000644024","Transcript","missense_variant","549","514","172","H/Y","Cac/Tac","-","IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=POU3F4;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:9217;BIOTYPE=protein_coding;CANONICAL=YES;MANE=NM_000307.5;APPRIS=P1;CCDS=CCDS14450.1;ENSP=ENSP00000495996;TREMBL=A0A2R8Y739;UNIPARC=UPI0000131D8A;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;GENE_PHENO=1;SIFT=tolerated(0.17);PolyPhen=benign(0.216);EXON=1/1;DOMAINS=PIRSF:PIRSF002629,PANTHER:PTHR11636,PANTHER:PTHR11636:SF83;HGVSc=ENST00000644024.2:c.514C>T;HGVSp=ENSP00000495996.1:p.His172Tyr" ],
[ "X_83508838_C/T","X:83508838","T","5456","NM_000307.5","Transcript","missense_variant","549","514","172","H/Y","Cac/Tac","-","IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=POU3F4;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:9217;BIOTYPE=protein_coding;CANONICAL=YES;ENSP=NP_000298.3;SOURCE=RefSeq;GIVEN_REF=C;USED_REF=C;SIFT=tolerated(0.17);PolyPhen=benign(0.216);EXON=1/1;HGVSc=NM_000307.5:c.514C>T;HGVSp=NP_000298.3:p.His172Tyr" ],
[ "X_83508838_C/T","X:83508838","T","107985635","XR_001755986.1","Transcript","upstream_gene_variant","-","-","-","-","-","-","IMPACT=MODIFIER;DISTANCE=2022;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=LOC107985635;SYMBOL_SOURCE=EntrezGene;BIOTYPE=lncRNA;CANONICAL=YES;SOURCE=RefSeq;GIVEN_REF=C;USED_REF=C" ],
[ "X_83508838_C/T","X:83508838","T","-","ENSR00000909988","RegulatoryFeature","regulatory_region_variant","-","-","-","-","-","-","IMPACT=MODIFIER;VARIANT_CLASS=SNV;BIOTYPE=promoter" ],
[ "X_83508838_C/T","X:83508838","T","ENSG00000279437","ENST00000625081","Transcript","non_coding_transcript_exon_variant","377","-","-","-","-","-","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=Z82170.1;SYMBOL_SOURCE=Clone_based_ensembl_gene;BIOTYPE=TEC;CANONICAL=YES;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;EXON=1/1;HGVSc=ENST00000625081.1:n.377G>A" ],
[ "X_83508838_C/T","X:83508838","T","ENSG00000196767","ENST00000644024","Transcript","missense_variant","549","514","172","H/Y","Cac/Tac","-","IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=POU3F4;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:9217;BIOTYPE=protein_coding;CANONICAL=YES;MANE=NM_000307.5;APPRIS=P1;CCDS=CCDS14450.1;ENSP=ENSP00000495996;TREMBL=A0A2R8Y739;UNIPARC=UPI0000131D8A;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;GENE_PHENO=1;SIFT=tolerated(0.17);PolyPhen=benign(0.216);EXON=1/1;DOMAINS=PIRSF:PIRSF002629,PANTHER:PTHR11636,PANTHER:PTHR11636:SF83;HGVSc=ENST00000644024.2:c.514C>T;HGVSp=ENSP00000495996.1:p.His172Tyr" ],
[ "X_83508838_C/T","X:83508838","T","5456","NM_000307.5","Transcript","missense_variant","549","514","172","H/Y","Cac/Tac","-","IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=POU3F4;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:9217;BIOTYPE=protein_coding;CANONICAL=YES;ENSP=NP_000298.3;SOURCE=RefSeq;GIVEN_REF=C;USED_REF=C;SIFT=tolerated(0.17);PolyPhen=benign(0.216);EXON=1/1;HGVSc=NM_000307.5:c.514C>T;HGVSp=NP_000298.3:p.His172Tyr" ],
[ "X_83508838_C/T","X:83508838","T","107985635","XR_001755986.1","Transcript","upstream_gene_variant","-","-","-","-","-","-","IMPACT=MODIFIER;DISTANCE=2022;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=LOC107985635;SYMBOL_SOURCE=EntrezGene;BIOTYPE=lncRNA;CANONICAL=YES;SOURCE=RefSeq;GIVEN_REF=C;USED_REF=C" ],
[ "X_83508838_C/T","X:83508838","T","-","ENSR00000909988","RegulatoryFeature","regulatory_region_variant","-","-","-","-","-","-","IMPACT=MODIFIER;VARIANT_CLASS=SNV;BIOTYPE=promoter" ],
[ "X_83508838_C/T","X:83508838","T","ENSG00000279437","ENST00000625081","Transcript","non_coding_transcript_exon_variant","377","-","-","-","-","-","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=Z82170.1;SYMBOL_SOURCE=Clone_based_ensembl_gene;BIOTYPE=TEC;CANONICAL=YES;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;EXON=1/1;HGVSc=ENST00000625081.1:n.377G>A" ],
[ "X_83508838_C/T","X:83508838","T","ENSG00000196767","ENST00000644024","Transcript","missense_variant","549","514","172","H/Y","Cac/Tac","-","IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=POU3F4;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:9217;BIOTYPE=protein_coding;CANONICAL=YES;MANE=NM_000307.5;APPRIS=P1;CCDS=CCDS14450.1;ENSP=ENSP00000495996;TREMBL=A0A2R8Y739;UNIPARC=UPI0000131D8A;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;GENE_PHENO=1;SIFT=tolerated(0.17);PolyPhen=benign(0.216);EXON=1/1;DOMAINS=PIRSF:PIRSF002629,PANTHER:PTHR11636,PANTHER:PTHR11636:SF83;HGVSc=ENST00000644024.2:c.514C>T;HGVSp=ENSP00000495996.1:p.His172Tyr" ],
[ "X_83508838_C/T","X:83508838","T","5456","NM_000307.5","Transcript","missense_variant","549","514","172","H/Y","Cac/Tac","-","IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=POU3F4;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:9217;BIOTYPE=protein_coding;CANONICAL=YES;ENSP=NP_000298.3;SOURCE=RefSeq;GIVEN_REF=C;USED_REF=C;SIFT=tolerated(0.17);PolyPhen=benign(0.216);EXON=1/1;HGVSc=NM_000307.5:c.514C>T;HGVSp=NP_000298.3:p.His172Tyr" ],
[ "X_83508838_C/T","X:83508838","T","107985635","XR_001755986.1","Transcript","upstream_gene_variant","-","-","-","-","-","-","IMPACT=MODIFIER;DISTANCE=2022;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=LOC107985635;SYMBOL_SOURCE=EntrezGene;BIOTYPE=lncRNA;CANONICAL=YES;SOURCE=RefSeq;GIVEN_REF=C;USED_REF=C" ],
[ "X_83508838_C/T","X:83508838","T","-","ENSR00000909988","RegulatoryFeature","regulatory_region_variant","-","-","-","-","-","-","IMPACT=MODIFIER;VARIANT_CLASS=SNV;BIOTYPE=promoter" ],
[ "X_83508838_C/T","X:83508838","T","ENSG00000279437","ENST00000625081","Transcript","non_coding_transcript_exon_variant","377","-","-","-","-","-","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=Z82170.1;SYMBOL_SOURCE=Clone_based_ensembl_gene;BIOTYPE=TEC;CANONICAL=YES;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;EXON=1/1;HGVSc=ENST00000625081.1:n.377G>A" ],
[ "X_83508838_C/T","X:83508838","T","ENSG00000196767","ENST00000644024","Transcript","missense_variant","549","514","172","H/Y","Cac/Tac","-","IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=POU3F4;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:9217;BIOTYPE=protein_coding;CANONICAL=YES;MANE=NM_000307.5;APPRIS=P1;CCDS=CCDS14450.1;ENSP=ENSP00000495996;TREMBL=A0A2R8Y739;UNIPARC=UPI0000131D8A;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;GENE_PHENO=1;SIFT=tolerated(0.17);PolyPhen=benign(0.216);EXON=1/1;DOMAINS=PIRSF:PIRSF002629,PANTHER:PTHR11636,PANTHER:PTHR11636:SF83;HGVSc=ENST00000644024.2:c.514C>T;HGVSp=ENSP00000495996.1:p.His172Tyr" ],
[ "X_83508838_C/T","X:83508838","T","5456","NM_000307.5","Transcript","missense_variant","549","514","172","H/Y","Cac/Tac","-","IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=POU3F4;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:9217;BIOTYPE=protein_coding;CANONICAL=YES;ENSP=NP_000298.3;SOURCE=RefSeq;GIVEN_REF=C;USED_REF=C;SIFT=tolerated(0.17);PolyPhen=benign(0.216);EXON=1/1;HGVSc=NM_000307.5:c.514C>T;HGVSp=NP_000298.3:p.His172Tyr" ],
[ "X_83508838_C/T","X:83508838","T","107985635","XR_001755986.1","Transcript","upstream_gene_variant","-","-","-","-","-","-","IMPACT=MODIFIER;DISTANCE=2022;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=LOC107985635;SYMBOL_SOURCE=EntrezGene;BIOTYPE=lncRNA;CANONICAL=YES;SOURCE=RefSeq;GIVEN_REF=C;USED_REF=C" ],
[ "X_83508838_C/T","X:83508838","T","-","ENSR00000909988","RegulatoryFeature","regulatory_region_variant","-","-","-","-","-","-","IMPACT=MODIFIER;VARIANT_CLASS=SNV;BIOTYPE=promoter" ] ]


        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_var  iation",
                    "Extra"]

        testDF = ps.DataFrame(actualFailingData,columns=colNames)

        actualResults = AnnotationFilter.run(testDF, self.tmpFileName, self.tmpLogLocation)

        self.assertFalse("SOURCE=Clone_based_ensembl_gene" in actualResults.at[0,'Extra'])

    def test_bug02_When_returnPositiveTest(self):

        actualFailingAnnotations = [["1_156293306_A/G","1:156293306","G","ENSG00000163472","ENST00000295694","Transcript","downstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=863;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=TMEM79;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:28196;BIOTYPE=protein_coding;CANONICAL=YES;TSL=1;APPRIS=P2;CCDS=CCDS1138.1;ENSP=ENSP00000295694;SWISSPROT=Q9BSE2;UNIPARC=UPI000006F977;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000163472","ENST00000357501","Transcript","downstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=1358;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=TMEM79;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:28196;BIOTYPE=protein_coding;TSL=3;APPRIS=A2;ENSP=ENSP00000350100;TREMBL=Q5SZX3;UNIPARC=UPI0000470159;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000362007","Transcript","intron_variant","-","-","-","-","-rs10908494","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=protein_coding;CANONICAL=YES;MANE=NM_144580.3;TSL=1;APPRIS=P3;CCDS=CCDS1139.1;ENSP=ENSP00000354553;SWISSPROT=Q8WWB7;UNIPARC=UPI00000361F7;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;INTRON=5/5;HGVSc=ENST00000362007.6:c.1055+14T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000368264","Transcript","downstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=117;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=nonsense_mediated_decay;TSL=3;ENSP=ENSP00000357247;TREMBL=X6R6A3;UNIPARC=UPI000047015A;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000163472","ENST00000405535","Transcript","downstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=863;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=TMEM79;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:28196;BIOTYPE=protein_coding;MANE=NM_032323.3;TSL=1;APPRIS=P2;CCDS=CCDS1138.1;ENSP=ENSP00000384748;SWISSPROT=Q9BSE2;UNIPARC=UPI000006F977;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000163472","ENST00000456810","Transcript","downstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=1781;STRAND=1;FLAGS=cds_end_NF;VARIANT_CLASS=SNV;SYMBOL=TMEM79;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:28196;BIOTYPE=protein_coding;TSL=3;ENSP=ENSP00000394736;TREMBL=Q5SZW9;UNIPARC=UPI0000470158;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000461597","Transcript","intron_variant,NMD_transcript_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;STRAND=-1;FLAGS=cds_start_NF;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=nonsense_mediated_decay;TSL=3;ENSP=ENSP00000476260;TREMBL=V9GXZ9;UNIPARC=UPI0003B927A9;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;INTRON=1/2;HGVSc=ENST00000461597.1:c.62+14T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000163472","ENST00000463670","Transcript","downstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=863;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=TMEM79;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:28196;BIOTYPE=lncRNA;TSL=2;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000472870","Transcript","intron_variant,NMD_transcript_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=nonsense_mediated_decay;TSL=5;ENSP=ENSP00000477180;TREMBL=V9GYX8;UNIPARC=UPI0000D4CC7E;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;INTRON=2/3;HGVSc=ENST00000472870.5:c.379-186T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000476177","Transcript","downstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=67;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=nonsense_mediated_decay;TSL=5;ENSP=ENSP00000477253;TREMBL=V9GYZ8;UNIPARC=UPI0003B927DD;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000479084","Transcript","downstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=212;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=retained_intron;TSL=2;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000480968","Transcript","upstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=1800;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=lncRNA;TSL=5;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000481050","Transcript","intron_variant","-","-","-","-","-rs10908494","IMPACT=MODIFIER;STRAND=-1;FLAGS=cds_start_NF;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=protein_coding;TSL=3;ENSP=ENSP00000477187;TREMBL=V9GYY1;UNIPARC=UPI0003B928CF;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;INTRON=2/2;HGVSc=ENST00000481050.1:c.169+14T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000482579","Transcript","intron_variant,non_coding_transcript_variant-","-","-","-","-","rs10908494","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=lncRNA;TSL=3;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;INTRON=1/1;HGVSc=ENST00000482579.1:n.51-97T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000484214","Transcript","downstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=311;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=retained_intron;TSL=2;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000163472","ENST00000485135","Transcript","downstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=121;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=TMEM79;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:28196;BIOTYPE=lncRNA;TSL=5;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000163472","ENST00000495881","Transcript","downstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=863;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=TMEM79;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:28196;BIOTYPE=lncRNA;TSL=2;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000497831","Transcript","upstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=1756;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=lncRNA;TSL=3;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000497955","Transcript","downstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=78;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=retained_intron;TSL=2;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000612353","Transcript","intron_variant","-","-","-","-","-rs10908494","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=protein_coding;TSL=5;CCDS=CCDS72949.1;ENSP=ENSP00000483691;TREMBL=A0A087X0W3;UNIPARC=UPI00024EDF24;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;INTRON=4/4;HGVSc=ENST00000612353.4:c.799-97T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000614643","Transcript","intron_variant","-","-","-","-","-rs10908494","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=protein_coding;TSL=2;APPRIS=A2;CCDS=CCDS72947.1;ENSP=ENSP00000480936;SWISSPROT=Q8WWB7;UNIPARC=UPI00017A72C9;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;INTRON=5/5;HGVSc=ENST00000614643.4:c.812+14T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000622703","Transcript","intron_variant","-","-","-","-","-rs10908494","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=protein_coding;TSL=2;CCDS=CCDS72948.1;ENSP=ENSP00000479149;TREMBL=A0A087WV34;UNIPARC=UPI0000D4E1E1;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;INTRON=4/4;HGVSc=ENST00000622703.4:c.797+14T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","ENSG00000198715","ENST00000647767","Transcript","intron_variant,NMD_transcript_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:29436;BIOTYPE=nonsense_mediated_decay;CCDS=CCDS1139.1;ENSP=ENSP00000497576;SWISSPROT=Q8WWB7;UNIPARC=UPI00000361F7;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A;INTRON=5/6;HGVSc=ENST00000647767.1:c.1055+14T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","112770","NM_001256604.1","Transcript","intron_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:29436;BIOTYPE=protein_coding;ENSP=NP_001243533.1;SOURCE=RefSeq;GIVEN_REF=A;USED_REF=A;BAM_EDIT=OK;INTRON=5/5;HGVSc=NM_001256604.1:c.908+14T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","112770","NM_001256605.1","Transcript","intron_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:29436;BIOTYPE=protein_coding;ENSP=NP_001243534.1;SOURCE=RefSeq;GIVEN_REF=A;USED_REF=A;BAM_EDIT=OK;INTRON=4/4;HGVSc=NM_001256605.1:c.797+14T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","112770","NM_001256608.1","Transcript","intron_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:29436;BIOTYPE=protein_coding;ENSP=NP_001243537.1;SOURCE=RefSeq;GIVEN_REF=A;USED_REF=A;BAM_EDIT=OK;INTRON=4/4;HGVSc=NM_001256608.1:c.799-97T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","112770","NM_001256609.1","Transcript","intron_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:29436;BIOTYPE=protein_coding;ENSP=NP_001243538.1;SOURCE=RefSeq;GIVEN_REF=A;USED_REF=A;BAM_EDIT=OK;INTRON=5/5;HGVSc=NM_001256609.1:c.812+14T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","23381","NM_001323617.1","Transcript","upstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=1808;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=SMG5;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:24644;BIOTYPE=protein_coding;CANONICAL=YES;ENSP=NP_001310546.1;SOURCE=RefSeq;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","84283","NM_032323.3","Transcript","downstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=863;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=TMEM79;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:28196;BIOTYPE=protein_coding;ENSP=NP_115699.1;SOURCE=RefSeq;GIVEN_REF=A;USED_REF=A;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","112770","NM_144580.3","Transcript","intron_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:29436;BIOTYPE=protein_coding;ENSP=NP_653181.1;SOURCE=RefSeq;GIVEN_REF=A;USED_REF=A;INTRON=5/5;HGVSc=NM_144580.3:c.1055+14T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","84283","NR_026678.1","Transcript","downstream_gene_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;DISTANCE=863;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=TMEM79;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:28196;BIOTYPE=misc_RNA;CANONICAL=YES;SOURCE=RefSeq;GIVEN_REF=A;USED_REF=A;BAM_EDIT=OK;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","112770","XM_011509117.3","Transcript","intron_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:29436;BIOTYPE=protein_coding;ENSP=XP_011507419.1;SOURCE=RefSeq;GIVEN_REF=A;USED_REF=A;INTRON=4/4;HGVSc=XM_011509117.3:c.524+14T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"],
        ["1_156293306_A/G","1:156293306","G","112770","XR_426762.3","Transcript","intron_variant,non_coding_transcript_variant","-","-","-","-","-","rs10908494","IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=GLMP;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:29436;BIOTYPE=misc_RNA;CANONICAL=YES;SOURCE=RefSeq;GIVEN_REF=A;USED_REF=A;INTRON=5/5;HGVSc=XR_426762.3:n.936-97T>C;AF=0.2065;AFR_AF=0.0923;AMR_AF=0.2795;EAS_AF=0.2272;EUR_AF=0.2694;SAS_AF=0.2229;AA_AF=0.1232;EA_AF=0.2929;gnomAD_AF=0.2707;gnomAD_AFR_AF=0.1174;gnomAD_AMR_AF=0.2785;gnomAD_ASJ_AF=0.3119;gnomAD_EAS_AF=0.2511;gnomAD_FIN_AF=0.3001;gnomAD_NFE_AF=0.2928;gnomAD_OTH_AF=0.2865;gnomAD_SAS_AF=0.2359;MAX_AF=0.3119;MAX_AF_POPS=gnomAD_ASJ"]]

        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]

        testDF = ps.DataFrame(actualFailingAnnotations, columns=colNames)

        actualResults = AnnotationFilter.run(testDF, self.tmpFileName, self.tmpLogLocation)

        self.assertFalse("SOURCE=Clone_based_ensembl_gene" in actualResults.at[0, 'Extra'])

    def test_ifMultipleEMmblArePassedWithoutNcbiRows_when_givenToFilter_returnOnlyOneEnsembleMaybe(self):
        actualFailingData = [["X_83508838_C/T", "X:83508838", "T", "ENSG00000279437", "ENST00000625081", "Transcript",
                              "non_coding_transcript_exon_variant", "377", "-", "-", "-", "-", "-",
                              "IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=Z82170.1;SYMBOL_SOURCE=Clone_based_ensembl_gene;BIOTYPE=TEC;CANONICAL=YES;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;EXON=1/1;HGVSc=ENST00000625081.1:n.377G>A"],
                             ["X_83508838_C/T", "X:83508838", "T", "ENSG00000196767", "ENST00000644024", "Transcript",
                              "missense_variant", "549", "514", "172", "H/Y", "Cac/Tac", "-",
                              "IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=POU3F4;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:9217;BIOTYPE=protein_coding;CANONICAL=YES;MANE=NM_000307.5;APPRIS=P1;CCDS=CCDS14450.1;ENSP=ENSP00000495996;TREMBL=A0A2R8Y739;UNIPARC=UPI0000131D8A;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;GENE_PHENO=1;SIFT=tolerated(0.17);PolyPhen=benign(0.216);EXON=1/1;DOMAINS=PIRSF:PIRSF002629,PANTHER:PTHR11636,PANTHER:PTHR11636:SF83;HGVSc=ENST00000644024.2:c.514C>T;HGVSp=ENSP00000495996.1:p.His172Tyr"],
                             ["X_83508838_C/T", "X:83508838", "T", "-", "ENSR00000909988", "RegulatoryFeature",
                              "regulatory_region_variant", "-", "-", "-", "-", "-", "-",
                              "IMPACT=MODIFIER;VARIANT_CLASS=SNV;BIOTYPE=promoter"],
                             ["X_83508838_C/T", "X:83508838", "T", "ENSG00000279437", "ENST00000625081", "Transcript",
                              "non_coding_transcript_exon_variant", "377", "-", "-", "-", "-", "-",
                              "IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=Z82170.1;SYMBOL_SOURCE=Clone_based_ensembl_gene;BIOTYPE=TEC;CANONICAL=YES;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;EXON=1/1;HGVSc=ENST00000625081.1:n.377G>A"],
                             ["X_83508838_C/T", "X:83508838", "T", "ENSG00000196767", "ENST00000644024", "Transcript",
                              "missense_variant", "549", "514", "172", "H/Y", "Cac/Tac", "-",
                              "IMPACT=MODERATE;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=POU3F4;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:9217;BIOTYPE=protein_coding;CANONICAL=YES;MANE=NM_000307.5;APPRIS=P1;CCDS=CCDS14450.1;ENSP=ENSP00000495996;TREMBL=A0A2R8Y739;UNIPARC=UPI0000131D8A;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;GENE_PHENO=1;SIFT=tolerated(0.17);PolyPhen=benign(0.216);EXON=1/1;DOMAINS=PIRSF:PIRSF002629,PANTHER:PTHR11636,PANTHER:PTHR11636:SF83;HGVSc=ENST00000644024.2:c.514C>T;HGVSp=ENSP00000495996.1:p.His172Tyr"],
                             ["X_83508838_C/T", "X:83508838", "T", "-", "ENSR00000909988", "RegulatoryFeature",
                              "regulatory_region_variant", "-", "-", "-", "-", "-", "-",
                              "IMPACT=MODIFIER;VARIANT_CLASS=SNV;BIOTYPE=promoter"],
                             ["X_83508838_C/T", "X:83508838", "T", "ENSG00000279437", "ENST00000625081", "Transcript",
                              "non_coding_transcript_exon_variant", "377", "-", "-", "-", "-", "-",
                              "IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=Z82170.1;SYMBOL_SOURCE=Clone_based_ensembl_gene;BIOTYPE=TEC;CANONICAL=YES;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;EXON=1/1;HGVSc=ENST00000625081.1:n.377G>A"]]

        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]

        testDF = ps.DataFrame(actualFailingData,columns=colNames)

        actualResults = AnnotationFilter.run(testDF, self.tmpFileName, self.tmpLogLocation)

        self.assertEquals(len(actualResults), 1)