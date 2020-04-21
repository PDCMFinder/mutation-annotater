from unittest import TestCase
import AnnotationMerger
import pandas as ps


class TestFilter(TestCase):

    def test_GivenChrPosKeyWithoutMatchingKeys_WhenChrPosKeyIsused_Then_ReturnChrPosKey(self):
        chromo = '1'
        seq_start = '10000'
        refAllele = 'A'
        altAllele = 'T'

        expectedChrPosKey = "{0}:{1}".format(chromo, seq_start)

        rowDic = {'chromosome': chromo, 'seq_start_position': seq_start, 'ref_allele': refAllele,
                  'alt_allele': altAllele}
        row = ps.Series(rowDic)

        actualChrPos = AnnotationMerger.formatChrPosKey(row)

        self.assertEqual(expectedChrPosKey, actualChrPos)

    def test_GivenChrPosKeyWithInsertionWithMatchingFirstNucleotides_WhenChrPosKeyIsused_Then_ReturnShiftedAndHyphonatedResult(
            self):
        chromo = '1'
        seq_start = '10000'
        seq_shift = '10001'
        refAllele = 'A'
        altAllele = 'ATT'

        expectedChrPosKey = "{0}:{1}-{2}".format(chromo, seq_start, seq_shift)

        rowDic = {'chromosome': chromo, 'seq_start_position': seq_start, 'ref_allele': refAllele,
                  'alt_allele': altAllele}
        row = ps.Series(rowDic)

        actualChrPos = AnnotationMerger.formatChrPosKey(row)

        self.assertEqual(expectedChrPosKey, actualChrPos)

    def test_GivenTwoMatchingColumns_WhenRowBuilderIsCalled_Then_ReturnCorrectlyBuiltRow(self):
        emblGene = "ENSG00000068078"
        emblFeature = "ENST00000260795"
        ncbiGene = "9817"
        ncbiFeature = 'NM_203500.2'

        emblRow = ["4_1804915_A/G", "4:1804915", "G", emblGene, emblFeature, "Transcript",
                   "3_prime_UTR_variant", \
                   "1509", "-", "-", "-", "-", "rs1466726466", \
                   "IMPACT=MODIFIER;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]
        ncbiRow = {'Feature_type': 'Transcript', 'Extra': 'IMPACT=HIGH;STRAND=-1;VARIANT_'
                                                          'CLASS=deletion;SYMBOL=KEAP1;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:23177;BIOTYPE=protein_coding;CANONICAL=YES;ENSP=NP_987096.1;SOURCE=RefSeq;GIVEN_REF=T;USED_REF=T;EXON=3/6;HGVSc=NM_203500.2:c.1225del;HGVSp=NP_987096.1:p.Met409Ter',
                   '#Uploaded_variation': '19_10491677_T/-', 'cDNA_position': '1388', 'Feature': ncbiFeature,
                   'Codons': 'Atg/tg', 'Existing_variation': '-', 'Location': '19:10491677', 'CDS_position': '1225',
                   'Protein_position': '409', 'Consequence': 'frameshift_variant', 'Allele': '-', 'Gene': ncbiGene,
                   'Amino_acids': 'M/X'}
        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]
        tsvInputRow = {'alt_allele': '-', 'read_depth': '11572', 'ensembl_gene_id': '', 'rs_id_Variant': '',
                       'Platform': 'Targeted NGS', 'ucsc_gene_id': '',
                       'ref_allele': 'T', 'chromosome': '19', 'ncbi_gene_id': '', 'nucleotide_change': '1225delA',
                       'hgnc_symbol': 'KEAP1', 'consequence': 'stop gain', 'Allele_frequency': '96,47', 'Passage': '27',
                       'amino_acid_change': 'M409X', 'Model_ID': 'LCF16',
                       'datasource': 'CURIE-LC', 'ensembl_transcript_id': '', 'Sample_ID': 'LCF16p27:26/08/2016',
                       'seq_start_position': '10491677', 'genome_assembly': 'GRCh38', 'sample_origin': 'xenograft'}

        EMBLdf = ps.DataFrame([emblRow], columns=colNames)
        NCBIdf = ps.DataFrame([ncbiRow], columns=colNames)
        inputRows = ps.concat([EMBLdf, NCBIdf])

        actualRow = AnnotationMerger.buildFinalTemplate(inputRows, tsvInputRow)

        self.assertEquals(actualRow[20], ncbiGene)
        self.assertEquals(actualRow[21], ncbiFeature)
        self.assertEquals(actualRow[22], emblGene)
        self.assertEquals(actualRow[23], emblFeature)

    def test_GivenTwoEmblRow_When_buildFinalTemplateIsCalled_Return_BuiltRowWithoutNCBIdata(self):
        emblGene = "ENSG00000068078"
        emblFeature = "ENST00000260795"
        emblGene2 = "ENSG0000000999"
        emblFeature2 = 'ENST000000111'

        emblRow = ["4_1804915_A/G", "4:1804915", "G", emblGene, emblFeature, "Transcript",
                   "3_prime_UTR_variant",
                   "1509", "-", "-", "-", "-", "rs1466726466",
                   "IMPACT=MODIFIER;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]
        emblRow2 = ["4_1804915_A/G", "4:1804916", "G", emblGene2, emblFeature2, "Transcript",
                    "missense",
                    "1509", "-", "-", "-", "-", "rs1466726466",
                    "IMPACT=High;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=FGFR3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC "]

        colNames = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                    "Existing_variation",
                    "Extra"]
        tsvInputRow = {'alt_allele': '-', 'read_depth': '11572', 'ensembl_gene_id': '', 'rs_id_Variant': '',
                       'Platform': 'Targeted NGS', 'ucsc_gene_id': '',
                       'ref_allele': 'T', 'chromosome': '19', 'ncbi_gene_id': '', 'nucleotide_change': '1225delA',
                       'hgnc_symbol': 'KEAP1', 'consequence': 'stop gain', 'Allele_frequency': '96,47', 'Passage': '27',
                       'amino_acid_change': 'M409X', 'Model_ID': 'LCF16',
                       'datasource': 'CURIE-LC', 'ensembl_transcript_id': '', 'Sample_ID': 'LCF16p27:26/08/2016',
                       'seq_start_position': '10491677', 'genome_assembly': 'GRCh38', 'sample_origin': 'xenograft'}

        inputRows = ps.DataFrame([emblRow, emblRow2], columns=colNames)
        actualRow = AnnotationMerger.buildFinalTemplate(inputRows, tsvInputRow)

        self.assertEquals(actualRow[20], None)
        self.assertEquals(actualRow[21], None)
        self.assertEquals(actualRow[22], emblGene)
        self.assertEquals(actualRow[23], emblFeature)

    def test_blankDfFromFilter_When_buildFinalTemplateIsCalled_Return_blankList(self):
        tsvInputRow = {'alt_allele': '-', 'read_depth': '11572', 'ensembl_gene_id': '', 'rs_id_Variant': '',
                       'Platform': 'Targeted NGS', 'ucsc_gene_id': '',
                       'ref_allele': 'T', 'chromosome': '19', 'ncbi_gene_id': '', 'nucleotide_change': '1225delA',
                       'hgnc_symbol': 'KEAP1', 'consequence': 'stop gain', 'Allele_frequency': '96,47', 'Passage': '27',
                       'amino_acid_change': 'M409X', 'Model_ID': 'LCF16',
                       'datasource': 'CURIE-LC', 'ensembl_transcript_id': '', 'Sample_ID': 'LCF16p27:26/08/2016',
                       'seq_start_position': '10491677', 'genome_assembly': 'GRCh38', 'sample_origin': 'xenograft'}

        inputRows = ps.DataFrame()
        actualRow = AnnotationMerger.buildFinalTemplate(inputRows, tsvInputRow)
        self.assertEquals(len(actualRow), 0)

    def test_blankRowsFilter_When_buildFinalTemplateIsCalled_Return_blankList(self):
        blankrow1 = []
        blankrow2 = []
        tsvInputRow = {'alt_allele': '-', 'read_depth': '11572', 'ensembl_gene_id': '', 'rs_id_Variant': '',
                       'Platform': 'Targeted NGS', 'ucsc_gene_id': '',
                       'ref_allele': 'T', 'chromosome': '19', 'ncbi_gene_id': '', 'nucleotide_change': '1225delA',
                       'hgnc_symbol': 'KEAP1', 'consequence': 'stop gain', 'Allele_frequency': '96,47',
                       'Passage': '27',
                       'amino_acid_change': 'M409X', 'Model_ID': 'LCF16',
                       'datasource': 'CURIE-LC', 'ensembl_transcript_id': '', 'Sample_ID': 'LCF16p27:26/08/2016',
                       'seq_start_position': '10491677', 'genome_assembly': 'GRCh38', 'sample_origin': 'xenograft'}

        inputRows = ps.DataFrame([blankrow1, blankrow2])
        actualRow = AnnotationMerger.buildFinalTemplate(inputRows, tsvInputRow)
        self.assertEquals(len(actualRow), 0)

    def test_blankRowsWithBlankCellsFilter_When_buildFinalTemplateIsCalled_Return_blankList(self):
        blankrow1 = ["", "", "", "", "", "", "", "", "", "", "", "", "", "", ""]
        blankrow2 = ["", "", "", "", "", "", "", "", "", "", "", "", "", "", ""]
        tsvInputRow = {'alt_allele': '-', 'read_depth': '11572', 'ensembl_gene_id': '', 'rs_id_Variant': '',
                       'Platform': 'Targeted NGS', 'ucsc_gene_id': '',
                       'ref_allele': 'T', 'chromosome': '19', 'ncbi_gene_id': '', 'nucleotide_change': '1225delA',
                       'hgnc_symbol': 'KEAP1', 'consequence': 'stop gain', 'Allele_frequency': '96,47',
                       'Passage': '27',
                       'amino_acid_change': 'M409X', 'Model_ID': 'LCF16',
                       'datasource': 'CURIE-LC', 'ensembl_transcript_id': '', 'Sample_ID': 'LCF16p27:26/08/2016',
                       'seq_start_position': '10491677', 'genome_assembly': 'GRCh38', 'sample_origin': 'xenograft'}

        inputRows = ps.DataFrame([blankrow1, blankrow2])
        actualRow = AnnotationMerger.buildFinalTemplate(inputRows, tsvInputRow)
        self.assertEquals(len(actualRow), 0)
