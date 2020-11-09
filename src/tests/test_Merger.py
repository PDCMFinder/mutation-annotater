from unittest import TestCase
from src import AnnotationMerger
import pandas as ps


class TestFilter(TestCase):

    def test_GivenChrPosKeyWithInsertionWithMatchingFirstNucleotides_WhenChrPosKeyIsused_Then_ReturnShiftedAndHyphonatedResult(
            self):
        chromo = '1'
        seq_start = '10000'
        refAllele = 'A'
        altAllele = 'ATT'
        expectedAnnoId = "{0}_{1}_{2}_{3}".format(chromo, seq_start, refAllele, altAllele)

        rowDic = {'chromosome': chromo, 'seq_start_position': seq_start, 'ref_allele': refAllele,
                  'alt_allele': altAllele}
        row = ps.Series(rowDic)

        annoId = AnnotationMerger.createAnnotationKey(row)

        self.assertEqual(annoId,expectedAnnoId)

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

    def test_WhenDepracatedHeadersAreUse_Return_correspondingValue(self):
        deprecatedTsvInputRow = {'alt_allele': '-', 'read_depth': '11572', 'ensembl_gene_id': '', 'rs_id_Variant': '',
                       'Platform': 'Targeted NGS', 'ucsc_gene_id': '',
                       'ref_allele': 'T', 'chromosome': '19', 'ncbi_gene_id': '', 'nucleotide_change': '1225delA',
                       'hgnc_symbol': 'KEAP1', 'consequence': 'stop gain', 'Allele_frequency': '96,47', 'Passage': '27',
                       'amino_acid_change': 'M409X', 'Model_ID': 'LCF16',
                       'datasource': 'CURIE-LC', 'ensembl_transcript_id': '', 'Sample_ID': 'LCF16p27:26/08/2016',
                       'seq_start_position': '10491677', 'genome_assembly': 'GRCh38', 'sample_origin': 'xenograft'}
        tsvInputRow = {'alt_allele': '-', 'read_depth': '11572', 'ensembl_gene_id': '', 'rs_id_Variant': '',
                                 'platform': 'Targeted NGS', 'ucsc_gene_id': '',
                                 'ref_allele': 'T', 'chromosome': '19', 'ncbi_gene_id': '',
                                 'nucleotide_change': '1225delA',
                                 'symbol': 'KEAP1', 'consequence': 'stop gain', 'allele_frequency': '96,47',
                                 'passage': '27',
                                 'amino_acid_change': 'M409X', 'model_id': 'LCF16',
                                 'datasource': 'CURIE-LC', 'ensembl_transcript_id': '',
                                 'sample_ID': 'LCF16p27:26/08/2016',
                                 'seq_start_position': '10491677', 'genome_assembly': 'GRCh38',
                                 'sample_origin': 'xenograft'}

        actualDeprecatedHeader = AnnotationMerger.getEitherFromRow(deprecatedTsvInputRow, "Model_ID", "model_id")
        actualUpdatedHeader = AnnotationMerger.getEitherFromRow(tsvInputRow, "Model_ID", "model_id")

        self.assertEquals(actualDeprecatedHeader, 'LCF16')
        self.assertEquals(actualUpdatedHeader, 'LCF16')

    def test_Given_MultipleSNVs_WhenMatchingKeys_returnRowsWithMatchChrosomeAndRefAndAltAllele(self):
        expectedRef = "C"
        expectedAlt = "T"

        annoFileHeaders = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
          "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation",
          "Extra"]

        annoFile = [["12_25245350_C_A", "12:25245350", "A", "ENSG00000133703", "ENST00000311936", "Transcript",
          "missense_variant", "225", "35", "12", "G/V", "gGt/gTt",
          "rs121913529,COSM1135366,COSM1140133,COSM1140134,COSM12657,COSM49168,COSM520,COSM521,COSM522",
          "IMPACT=MODERATE;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=KRAS;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:6407;BIOTYPE=protein_coding;TSL=1;APPRIS=A1;CCDS=CCDS8702.1;ENSP=ENSP00000308495;SWISSPROT=P01116;TREMBL=A0A024RAV5;UNIPARC=UPI0000001252;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;GENE_PHENO=1;SIFT=deleterious(0);PolyPhen=probably_damaging(0.96);EXON=2/5;DOMAINS=Low_complexity_(Seg):seg,PROSITE_profiles:PS51421,CDD:cd04138,PANTHER:PTHR24070,PANTHER:PTHR24070:SF186,Pfam:PF00071,TIGRFAM:TIGR00231,Gene3D:3.40.50.300,SMART:SM00173,SMART:SM00174,SMART:SM00175,Superfamily:SSF52540,Prints:PR00449;HGVSc=ENST00000311936.8:c.35G>T;HGVSp=ENSP00000308495.3:p.Gly12Val;gnomAD_AF=0;gnomAD_AFR_AF=0;gnomAD_AMR_AF=0;gnomAD_ASJ_AF=0;gnomAD_EAS_AF=0;gnomAD_FIN_AF=0;gnomAD_NFE_AF=0;gnomAD_OTH_AF=0;gnomAD_SAS_AF=0;MAX_AF=0;MAX_AF_POPS=gnomAD_AFR,gnomAD_AMR,gnomAD_ASJ,gnomAD_EAS,gnomAD_FIN,gnomAD_NFE,gnomAD_OTH,gnomAD_SAS;CLIN_SIG=pathogenic,pathogenic/likely_pathogenic,likely_pathogenic;SOMATIC=0,1,1,1,1,1,1,1,1;PHENO=1,1,1,1,1,1,1,1,1;PUBMED=25157968,19029981,22499344,17332249,21079152,2278970,3122217,12460918,16434492,19075190,22407852,17384584,18794081,19047918,19255327,19773371,23406027,22683711,19018267,21975775,17704260,15696205,16361624,16618717,18316791,19114683,19679400,20921462,20921465,21228335,19794967,21398618,23182985,7773929,8439212,15842656,17910045,19358724,19881948,20609353,20805368,20949522,21169357,22025163,22235099,22282465,22897852,23014527,25044103,26372703,27872090,29525983,30463544"],
         ["12_25245350_C/_", "12:25245350", "A", "ENSG00000133703", "ENST00000556131", "Transcript",
          "missense_variant", "212", "35", "12", "G/V", "gGt/gTt",
          "rs121913529,COSM1135366,COSM1140133,COSM1140134,COSM12657,COSM49168,COSM520,COSM521,COSM522",
          "IMPACT=MODERATE;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=KRAS;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:6407;BIOTYPE=protein_coding;TSL=1;ENSP=ENSP00000451856;TREMBL=G3V4K2;UNIPARC=UPI00001FB7B5;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;GENE_PHENO=1;SIFT=deleterious(0);PolyPhen=possibly_damaging(0.884);EXON=2/3;DOMAINS=Low_complexity_(Seg):seg,PROSITE_profiles:PS51421,PANTHER:PTHR24070,PANTHER:PTHR24070:SF397,Pfam:PF00071,Gene3D:3.40.50.300,Superfamily:SSF52540,Prints:PR00449;HGVSc=ENST00000556131.1:c.35G>T;HGVSp=ENSP00000451856.1:p.Gly12Val;gnomAD_AF=0;gnomAD_AFR_AF=0;gnomAD_AMR_AF=0;gnomAD_ASJ_AF=0;gnomAD_EAS_AF=0;gnomAD_FIN_AF=0;gnomAD_NFE_AF=0;gnomAD_OTH_AF=0;gnomAD_SAS_AF=0;MAX_AF=0;MAX_AF_POPS=gnomAD_AFR,gnomAD_AMR,gnomAD_ASJ,gnomAD_EAS,gnomAD_FIN,gnomAD_NFE,gnomAD_OTH,gnomAD_SAS;CLIN_SIG=pathogenic,pathogenic/likely_pathogenic,likely_pathogenic;SOMATIC=0,1,1,1,1,1,1,1,1;PHENO=1,1,1,1,1,1,1,1,1;PUBMED=25157968,19029981,22499344,17332249,21079152,2278970,3122217,12460918,16434492,19075190,22407852,17384584,18794081,19047918,19255327,19773371,23406027,22683711,19018267,21975775,17704260,15696205,16361624,16618717,18316791,19114683,19679400,20921462,20921465,21228335,19794967,21398618,23182985,7773929,8439212,15842656,17910045,19358724,19881948,20609353,20805368,20949522,21169357,22025163,22235099,22282465,22897852,23014527,25044103,26372703,27872090,29525983,30463544"],
         ["12_25245350_C_A", "12:25245350", "A", "ENSG00000133703", "ENST00000557334", "Transcript",
          "missense_variant", "232", "35", "12", "G/V", "gGt/gTt",
          "rs121913529,COSM1135366,COSM1140133,COSM1140134,COSM12657,COSM49168,COSM520,COSM521,COSM522",
          "IMPACT=MODERATE;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=KRAS;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:6407;BIOTYPE=protein_coding;TSL=5;ENSP=ENSP00000452512;TREMBL=G3V5T7;UNIPARC=UPI00021CF477;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;GENE_PHENO=1;SIFT=deleterious_low_confidence(0);PolyPhen=probably_damaging(0.979);EXON=2/3;DOMAINS=Gene3D:3.40.50.300,Pfam:PF00071,Prints:PR00449,PROSITE_profiles:PS51421,PANTHER:PTHR24070,PANTHER:PTHR24070:SF186,SMART:SM00173,Superfamily:SSF52540,Low_complexity_(Seg):seg;HGVSc=ENST00000557334.5:c.35G>T;HGVSp=ENSP00000452512.1:p.Gly12Val;gnomAD_AF=0;gnomAD_AFR_AF=0;gnomAD_AMR_AF=0;gnomAD_ASJ_AF=0;gnomAD_EAS_AF=0;gnomAD_FIN_AF=0;gnomAD_NFE_AF=0;gnomAD_OTH_AF=0;gnomAD_SAS_AF=0;MAX_AF=0;MAX_AF_POPS=gnomAD_AFR,gnomAD_AMR,gnomAD_ASJ,gnomAD_EAS,gnomAD_FIN,gnomAD_NFE,gnomAD_OTH,gnomAD_SAS;CLIN_SIG=pathogenic,pathogenic/likely_pathogenic,likely_pathogenic;SOMATIC=0,1,1,1,1,1,1,1,1;PHENO=1,1,1,1,1,1,1,1,1;PUBMED=25157968,19029981,22499344,17332249,21079152,2278970,3122217,12460918,16434492,19075190,22407852,17384584,18794081,19047918,19255327,19773371,23406027,22683711,19018267,21975775,17704260,15696205,16361624,16618717,18316791,19114683,19679400,20921462,20921465,21228335,19794967,21398618,23182985,7773929,8439212,15842656,17910045,19358724,19881948,20609353,20805368,20949522,21169357,22025163,22235099,22282465,22897852,23014527,25044103,26372703,27872090,29525983,30463544"],
         ["12_25245350_C_T", "12:25245350", "T", "ENSG00000133703", "ENST00000311936", "Transcript",
          "missense_variant", "225", "35", "12", "G/D", "gGt/gAt",
          "rs121913529,COSM1135366,COSM1140133,COSM1140134,COSM12657,COSM49168,COSM520,COSM521,COSM522",
          "IMPACT=MODERATE;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=KRAS;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:6407;BIOTYPE=protein_coding;TSL=1;APPRIS=A1;CCDS=CCDS8702.1;ENSP=ENSP00000308495;SWISSPROT=P01116;TREMBL=A0A024RAV5;UNIPARC=UPI0000001252;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;GENE_PHENO=1;SIFT=deleterious(0);PolyPhen=benign(0.21);EXON=2/5;DOMAINS=Low_complexity_(Seg):seg,PROSITE_profiles:PS51421,CDD:cd04138,PANTHER:PTHR24070,PANTHER:PTHR24070:SF186,Pfam:PF00071,TIGRFAM:TIGR00231,Gene3D:3.40.50.300,SMART:SM00173,SMART:SM00174,SMART:SM00175,Superfamily:SSF52540,Prints:PR00449;HGVSc=ENST00000311936.8:c.35G>A;HGVSp=ENSP00000308495.3:p.Gly12Asp;gnomAD_AF=4.011e-06;gnomAD_AFR_AF=0;gnomAD_AMR_AF=0;gnomAD_ASJ_AF=0;gnomAD_EAS_AF=0;gnomAD_FIN_AF=0;gnomAD_NFE_AF=8.883e-06;gnomAD_OTH_AF=0;gnomAD_SAS_AF=0;MAX_AF=8.883e-06;MAX_AF_POPS=gnomAD_NFE;CLIN_SIG=pathogenic,likely_pathogenic;SOMATIC=0,1,1,1,1,1,1,1,1;PHENO=1,1,1,1,1,1,1,1,1;PUBMED=25157968,19029981,22499344,17332249,21079152,2278970,3122217,12460918,16434492,19075190,22407852,17384584,18794081,19047918,19255327,19773371,23406027,22683711,19018267,21975775,17704260,15696205,16361624,16618717,18316791,19114683,19679400,20921462,20921465,21228335,19794967,21398618,23182985,7773929,8439212,15842656,17910045,19358724,19881948,20609353,20805368,20949522,21169357,22025163,22235099,22282465,22897852,23014527,25044103,26372703,27872090,29525983,30463544"],
         ["12_25245350_C_T", "12:25245350", "T", "ENSG00000133703", "ENST00000556131", "Transcript",
          "missense_variant", "212", "35", "12", "G/D", "gGt/gAt",
          "rs121913529,COSM1135366,COSM1140133,COSM1140134,COSM12657,COSM49168,COSM520,COSM521,COSM522",
          "IMPACT=MODERATE;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=KRAS;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:6407;BIOTYPE=protein_coding;TSL=1;ENSP=ENSP00000451856;TREMBL=G3V4K2;UNIPARC=UPI00001FB7B5;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;GENE_PHENO=1;SIFT=deleterious(0);PolyPhen=benign(0.039);EXON=2/3;DOMAINS=Low_complexity_(Seg):seg,PROSITE_profiles:PS51421,PANTHER:PTHR24070,PANTHER:PTHR24070:SF397,Pfam:PF00071,Gene3D:3.40.50.300,Superfamily:SSF52540,Prints:PR00449;HGVSc=ENST00000556131.1:c.35G>A;HGVSp=ENSP00000451856.1:p.Gly12Asp;gnomAD_AF=4.011e-06;gnomAD_AFR_AF=0;gnomAD_AMR_AF=0;gnomAD_ASJ_AF=0;gnomAD_EAS_AF=0;gnomAD_FIN_AF=0;gnomAD_NFE_AF=8.883e-06;gnomAD_OTH_AF=0;gnomAD_SAS_AF=0;MAX_AF=8.883e-06;MAX_AF_POPS=gnomAD_NFE;CLIN_SIG=pathogenic,likely_pathogenic;SOMATIC=0,1,1,1,1,1,1,1,1;PHENO=1,1,1,1,1,1,1,1,1;PUBMED=25157968,19029981,22499344,17332249,21079152,2278970,3122217,12460918,16434492,19075190,22407852,17384584,18794081,19047918,19255327,19773371,23406027,22683711,19018267,21975775,17704260,15696205,16361624,16618717,18316791,19114683,19679400,20921462,20921465,21228335,19794967,21398618,23182985,7773929,8439212,15842656,17910045,19358724,19881948,20609353,20805368,20949522,21169357,22025163,22235099,22282465,22897852,23014527,25044103,26372703,27872090,29525983,30463544"],
         ["12_25245350_C_T", "12:25245350", "T", "ENSG00000133703", "ENST00000557334", "Transcript",
          "missense_variant", "232", "35", "12", "G/D", "gGt/gAt",
          "rs121913529,COSM1135366,COSM1140133,COSM1140134,COSM12657,COSM49168,COSM520,COSM521,COSM522",
          "IMPACT=MODERATE;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=KRAS;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:6407;BIOTYPE=protein_coding;TSL=5;ENSP=ENSP00000452512;TREMBL=G3V5T7;UNIPARC=UPI00021CF477;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;GENE_PHENO=1;SIFT=deleterious_low_confidence(0);PolyPhen=probably_damaging(0.969);EXON=2/3;DOMAINS=Gene3D:3.40.50.300,Pfam:PF00071,Prints:PR00449,PROSITE_profiles:PS51421,PANTHER:PTHR24070,PANTHER:PTHR24070:SF186,SMART:SM00173,Superfamily:SSF52540,Low_complexity_(Seg):seg;HGVSc=ENST00000557334.5:c.35G>A;HGVSp=ENSP00000452512.1:p.Gly12Asp;gnomAD_AF=4.011e-06;gnomAD_AFR_AF=0;gnomAD_AMR_AF=0;gnomAD_ASJ_AF=0;gnomAD_EAS_AF=0;gnomAD_FIN_AF=0;gnomAD_NFE_AF=8.883e-06;gnomAD_OTH_AF=0;gnomAD_SAS_AF=0;MAX_AF=8.883e-06;MAX_AF_POPS=gnomAD_NFE;CLIN_SIG=pathogenic,likely_pathogenic;SOMATIC=0,1,1,1,1,1,1,1,1;PHENO=1,1,1,1,1,1,1,1,1;PUBMED=25157968,19029981,22499344,17332249,21079152,2278970,3122217,12460918,16434492,19075190,22407852,17384584,18794081,19047918,19255327,19773371,23406027,22683711,19018267,21975775,17704260,15696205,16361624,16618717,18316791,19114683,19679400,20921462,20921465,21228335,19794967,21398618,23182985,7773929,8439212,15842656,17910045,19358724,19881948,20609353,20805368,20949522,21169357,22025163,22235099,22282465,22897852,23014527,25044103,26372703,27872090,29525983,30463544"],
         ["12_25245350_C_G", "12:25245350", "G", "ENSG00000133703", "ENST00000311936", "Transcript",
          "missense_variant", "225", "35", "12", "G/A", "gGt/gCt",
          "rs121913529,COSM1135366,COSM1140133,COSM1140134,COSM12657,COSM49168,COSM520,COSM521,COSM522",
          "IMPACT=MODERATE;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=KRAS;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:6407;BIOTYPE=protein_coding;TSL=1;APPRIS=A1;CCDS=CCDS8702.1;ENSP=ENSP00000308495;SWISSPROT=P01116;TREMBL=A0A024RAV5;UNIPARC=UPI0000001252;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;GENE_PHENO=1;SIFT=deleterious(0.02);PolyPhen=possibly_damaging(0.597);EXON=2/5;DOMAINS=Low_complexity_(Seg):seg,PROSITE_profiles:PS51421,CDD:cd04138,PANTHER:PTHR24070,PANTHER:PTHR24070:SF186,Pfam:PF00071,TIGRFAM:TIGR00231,Gene3D:3.40.50.300,SMART:SM00173,SMART:SM00174,SMART:SM00175,Superfamily:SSF52540,Prints:PR00449;HGVSc=ENST00000311936.8:c.35G>C;HGVSp=ENSP00000308495.3:p.Gly12Ala;gnomAD_AF=0;gnomAD_AFR_AF=0;gnomAD_AMR_AF=0;gnomAD_ASJ_AF=0;gnomAD_EAS_AF=0;gnomAD_FIN_AF=0;gnomAD_NFE_AF=0;gnomAD_OTH_AF=0;gnomAD_SAS_AF=0;MAX_AF=0;MAX_AF_POPS=gnomAD_AFR,gnomAD_AMR,gnomAD_ASJ,gnomAD_EAS,gnomAD_FIN,gnomAD_NFE,gnomAD_OTH,gnomAD_SAS;CLIN_SIG=not_provided,likely_pathogenic,pathogenic;SOMATIC=0,1,1,1,1,1,1,1,1;PHENO=1,1,1,1,1,1,1,1,1;PUBMED=25157968,19029981,22499344,17332249,21079152,2278970,3122217,12460918,16434492,19075190,22407852,17384584,18794081,19047918,19255327,19773371,23406027,22683711,19018267,21975775,17704260,15696205,16361624,16618717,18316791,19114683,19679400,20921462,20921465,21228335,19794967,21398618,23182985,7773929,8439212,15842656,17910045,19358724,19881948,20609353,20805368,20949522,21169357,22025163,22235099,22282465,22897852,23014527,25044103,26372703,27872090,29525983,30463544"],
         ["12_25245350_C_G", "12:25245350", "G", "ENSG00000133703", "ENST00000556131", "Transcript",
          "missense_variant", "212", "35", "12", "G/A", "gGt/gCt",
          "rs121913529,COSM1135366,COSM1140133,COSM1140134,COSM12657,COSM49168,COSM520,COSM521,COSM522",
          "IMPACT=MODERATE;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=KRAS;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:6407;BIOTYPE=protein_coding;TSL=1;ENSP=ENSP00000451856;TREMBL=G3V4K2;UNIPARC=UPI00001FB7B5;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;GENE_PHENO=1;SIFT=deleterious(0.04);PolyPhen=possibly_damaging(0.614);EXON=2/3;DOMAINS=Low_complexity_(Seg):seg,PROSITE_profiles:PS51421,PANTHER:PTHR24070,PANTHER:PTHR24070:SF397,Pfam:PF00071,Gene3D:3.40.50.300,Superfamily:SSF52540,Prints:PR00449;HGVSc=ENST00000556131.1:c.35G>C;HGVSp=ENSP00000451856.1:p.Gly12Ala;gnomAD_AF=0;gnomAD_AFR_AF=0;gnomAD_AMR_AF=0;gnomAD_ASJ_AF=0;gnomAD_EAS_AF=0;gnomAD_FIN_AF=0;gnomAD_NFE_AF=0;gnomAD_OTH_AF=0;gnomAD_SAS_AF=0;MAX_AF=0;MAX_AF_POPS=gnomAD_AFR,gnomAD_AMR,gnomAD_ASJ,gnomAD_EAS,gnomAD_FIN,gnomAD_NFE,gnomAD_OTH,gnomAD_SAS;CLIN_SIG=not_provided,likely_pathogenic,pathogenic;SOMATIC=0,1,1,1,1,1,1,1,1;PHENO=1,1,1,1,1,1,1,1,1;PUBMED=25157968,19029981,22499344,17332249,21079152,2278970,3122217,12460918,16434492,19075190,22407852,17384584,18794081,19047918,19255327,19773371,23406027,22683711,19018267,21975775,17704260,15696205,16361624,16618717,18316791,19114683,19679400,20921462,20921465,21228335,19794967,21398618,23182985,7773929,8439212,15842656,17910045,19358724,19881948,20609353,20805368,20949522,21169357,22025163,22235099,22282465,22897852,23014527,25044103,26372703,27872090,29525983,30463544"],
         ["12_25245350_C_G", "12:25245350", "G", "ENSG00000133703", "ENST00000557334", "Transcript",
          "missense_variant", "232", "35", "12", "G/A", "gGt/gCt",
          "rs121913529,COSM1135366,COSM1140133,COSM1140134,COSM12657,COSM49168,COSM520,COSM521,COSM522",
          "IMPACT=MODERATE;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=KRAS;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:6407;BIOTYPE=protein_coding;TSL=5;ENSP=ENSP00000452512;TREMBL=G3V5T7;UNIPARC=UPI00021CF477;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C;GENE_PHENO=1;SIFT=tolerated_low_confidence(0.13);PolyPhen=possibly_damaging(0.899);EXON=2/3;DOMAINS=Gene3D:3.40.50.300,Pfam:PF00071,Prints:PR00449,PROSITE_profiles:PS51421,PANTHER:PTHR24070,PANTHER:PTHR24070:SF186,SMART:SM00173,Superfamily:SSF52540,Low_complexity_(Seg):seg;HGVSc=ENST00000557334.5:c.35G>C;HGVSp=ENSP00000452512.1:p.Gly12Ala;gnomAD_AF=0;gnomAD_AFR_AF=0;gnomAD_AMR_AF=0;gnomAD_ASJ_AF=0;gnomAD_EAS_AF=0;gnomAD_FIN_AF=0;gnomAD_NFE_AF=0;gnomAD_OTH_AF=0;gnomAD_SAS_AF=0;MAX_AF=0;MAX_AF_POPS=gnomAD_AFR,gnomAD_AMR,gnomAD_ASJ,gnomAD_EAS,gnomAD_FIN,gnomAD_NFE,gnomAD_OTH,gnomAD_SAS;CLIN_SIG=not_provided,likely_pathogenic,pathogenic;SOMATIC=0,1,1,1,1,1,1,1,1;PHENO=1,1,1,1,1,1,1,1,1;PUBMED=25157968,19029981,22499344,17332249,21079152,2278970,3122217,12460918,16434492,19075190,22407852,17384584,18794081,19047918,19255327,19773371,23406027,22683711,19018267,21975775,17704260,15696205,16361624,16618717,18316791,19114683,19679400,20921462,20921465,21228335,19794967,21398618,23182985,7773929,8439212,15842656,17910045,19358724,19881948,20609353,20805368,20949522,21169357,22025163,22235099,22282465,22897852,23014527,25044103,26372703,27872090,29525983,30463544"]]


        headers = ["datasource","Model_ID","Sample_ID","sample_origin","Passage","host_strain_name","hgnc_symbol","amino_acid_change","nucleotide_change","consequence","read_depth","Allele_frequency",
        "chromosome","seq_start_position","ref_allele","alt_allele","ucsc_gene_id","ncbi_gene_id","ensembl_gene_id","ensembl_transcript_id","rs_id_Variant","genome_assembly","Platform"]
        data = ["IRCC-CRC","CRC1082PR","CRC1082PRX0A01001TUMD05000","xenograft","1","NOD.Cg-Prkdc<sup>scid</sup> Il2rg<sup>tm1Wjl</sup>/SzJ","KRAS","G12D","","missense","","0.96","12","25245350",
         expectedRef,expectedAlt,"","","","","rs121913529","GRCh38","TargetedNGS_MUT"]

        testAnnotations = ps.DataFrame(annoFile, columns=annoFileHeaders)
        row = ps.DataFrame([data], columns=headers).to_dict("index").pop(0)

        actualMatches = AnnotationMerger.returnMatchingRows(row, testAnnotations)
        actualAllele = actualMatches['Allele'].unique()
        self.assertEquals(len(actualAllele), 1)

def test_givenNonMatchingFirstNucleotideAlleles_WhenformatAlleles_ThenReturnInput(self):
    basicRefAllele = "A"
    basicAltAllele = "T"
    row1 = {"ref_allele": basicRefAllele, "alt_allele": basicAltAllele}

    refAllele = "ATC"
    altAllele = "CTT"
    row2 = {"alt_allele" : altAllele, "ref_allele" : refAllele}

    actualRef1, actualAlt1 = AnnotationMerger.formatAlleles(row1)
    actualRef2, actualAlt2 = AnnotationMerger.formatAlleles(row2)

    self.assertEquals(actualRef1, basicRefAllele)
    self.assertEquals(actualRef2, refAllele)
    self.assertEquals(actualAlt1, basicAltAllele)
    self.assertEquals(actualAlt2, altAllele)

def test_givenNonMatchingFirstNucleotideAlleles_WhenformatAlleles_ThenReturnInput(self):
    basicRefAllele = "A"
    basicAltAllele = "T"
    row1 = {"ref_allele": basicRefAllele, "alt_allele": basicAltAllele}

    refAllele = "ATC"
    altAllele = "CTT"
    row2 = {"alt_allele" : altAllele, "ref_allele" : refAllele}

    actualRef1, actualAlt1 = AnnotationMerger.formatAlleles(row1)
    actualRef2, actualAlt2 = AnnotationMerger.formatAlleles(row2)

    self.assertEquals(actualRef1, basicRefAllele)
    self.assertEquals(actualRef2, refAllele)
    self.assertEquals(actualAlt1, basicAltAllele)
    self.assertEquals(actualAlt2, altAllele)


def test_givenMatchingFirstNucleotideAlleles_WhenformatAlleles_ThenReturnFirstDroppedAlleles(self):
    refAllele = "NC"
    altAllele = "NT"
    row = {"alt_allele" : altAllele, "ref_allele" : refAllele}
    refAllele2 = "NATC"
    altAllele2 = "NCTT"
    row2 = {"alt_allele": altAllele2, "ref_allele": refAllele2}

    expectedRefAllele1 = "C"
    expectedAltAllele1 = "T"
    expectedRefAllele2 = "ATC"
    expectedAltAllele2 = "CTT"

    actualRef1, actualAlt1 = AnnotationMerger.formatAlleles(row)
    actualRef2, actualAlt2 = AnnotationMerger.formatAlleles(row2)

    self.assertEquals(actualRef1, expectedRefAllele1)
    self.assertEquals(actualRef2, expectedRefAllele2)
    self.assertEquals(actualAlt1, expectedAltAllele1)
    self.assertEquals(actualAlt2, expectedAltAllele2)

def test_buildFinalAnnotation_realDroppedAnnotations_Debug(self):


    AnnotationMerger.buildFinalTemplate()