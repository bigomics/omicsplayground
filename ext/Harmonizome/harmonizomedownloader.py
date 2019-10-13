#!/usr/local/bin/python3

"""---------------------------------------------------------------------------
A utility script for downloading data from the Harmonizome, with the ability
to configure which datasets and which download types from which to download.
Note that all content decompressed is roughly 30GB. The default is to not
decompress the files on download.
Dependencies:
- requests is an HTTP library with an easy-to-use API:
  http://docs.python-requests.org/en/latest/
---------------------------------------------------------------------------"""


import requests
import os
import zlib


def _download_file(response, filename):
    with open(filename, 'wb') as f:
        for chunk in response.iter_content(chunk_size=1024):
            f.write(chunk)


def _download_and_decompress_file(response, filename):
    decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
    filename = filename[:-3]
    with open(filename, 'w+') as f:
        while True:
            chunk = response.raw.read(1024)
            if not chunk:
                break
            string = decompressor.decompress(chunk)
            f.write(string)


def download_datasets(selected_datasets, selected_downloads, decompress=False):
    for dataset, path in selected_datasets:
        if not os.path.exists(dataset):
            os.mkdir(dataset)

        for downloadable in selected_downloads:
            url = 'https://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/%s/%s' %\
                  (path, downloadable)
            response = requests.get(url, stream=True)
            filename = '%s/%s' % (dataset, downloadable)

            # Not every dataset has all downloadables.
            if response.status_code != 200:
                continue

            if decompress and 'txt.gz' in filename:
                _download_and_decompress_file(response, filename)
            else:
                _download_file(response, filename)

        print('%s downloaded.' % dataset)


if __name__ == '__main__':
    # Uncomment a dataset to download it.
    download_datasets([
        # ('Achilles Cell Line Gene Essentiality Profiles', 'achilles'),
        # ('Allen Brain Atlas Adult Human Brain Tissue Gene Expression Profiles', 'brainatlasadulthuman'),
        # ('Allen Brain Atlas Adult Mouse Brain Tissue Gene Expression Profiles', 'brainatlasadultmouse'),
        # ('Allen Brain Atlas Developing Human Brain Tissue Gene Expression Profiles by Microarray', 'brainatlasdevelopmentalhumanmicroarray'),
        # ('Allen Brain Atlas Developing Human Brain Tissue Gene Expression Profiles by RNA-seq', 'brainatlasdevelopmentalhumanrnaseq'),
        # ('Allen Brain Atlas Prenatal Human Brain Tissue Gene Expression Profiles', 'brainatlasprenatalhuman'),
        # ('BIND Biomolecular Interactions', 'bind'),
        # ('BioGPS Cell Line Gene Expression Profiles', 'biogpsnci60'),
        # ('BioGPS Human Cell Type and Tissue Gene Expression Profiles', 'biogpshuman'),
        # ('BioGPS Mouse Cell Type and Tissue Gene Expression Profiles', 'biogpsmouse'),
        # ('BioGRID Protein-Protein Interactions', 'biogrid'),
        # ('Biocarta Pathways', 'biocarta'),
        # ('CCLE Cell Line Gene CNV Profiles', 'cclecnv'),
        # ('CCLE Cell Line Gene Expression Profiles', 'cclemrna'),
        # ('CCLE Cell Line Gene Mutation Profiles', 'cclemut'),
        # ('CHEA Transcription Factor Binding Site Profiles', 'chea'),
        # ('CHEA Transcription Factor Targets', 'cheappi'),
        # ('CMAP Signatures of Differentially Expressed Genes for Small Molecules', 'cmap'),
        # ('COMPARTMENTS Curated Protein Localization Evidence Scores', 'jensencompartmentcurated'),
        # ('COMPARTMENTS Experimental Protein Localization Evidence Scores', 'jensencompartmentexpts'),
        # ('COMPARTMENTS Text-mining Protein Localization Evidence Scores', 'jensencompartmenttextmining'),
        # ('CORUM Protein Complexes', 'corum'),
        # ('COSMIC Cell Line Gene CNV Profiles', 'cosmiccnv'),
        # ('COSMIC Cell Line Gene Mutation Profiles', 'cosmicmut'),
        # ('CTD Gene-Chemical Interactions', 'ctdchemical'),
        # ('CTD Gene-Disease Associations', 'ctddisease'),
        # ('ClinVar SNP-Phenotype Associations', 'clinvar'),
        # ('Combined Pathways Pathways', 'combinedpathways'),
        # ('dbGAP Gene-Trait Associations', 'dbgap'),
        # ('DEPOD Substrates of Phosphatases', 'depod'),
        # ('DIP Protein-Protein Interactions', 'dip'),
        # ('DISEASES Curated Gene-Disease Assocation Evidence Scores', 'jensendiseasecurated'),
        # ('DISEASES Experimental Gene-Disease Assocation Evidence Scores', 'jensendiseaseexpts'),
        # ('DISEASES Text-mining Gene-Disease Assocation Evidence Scores', 'jensendiseasetextmining'),
        # ('DrugBank Drug Targets', 'drugbank'),
        # ('ENCODE Histone Modification Site Profiles', 'encodehm'),
        # ('ENCODE Transcription Factor Binding Site Profiles', 'encodetf'),
        # ('ENCODE Transcription Factor Targets', 'encodetfppi'),
        # ('ESCAPE Omics Signatures of Genes and Proteins for Stem Cells', 'escape'),
        # ('GAD Gene-Disease Associations', 'gad'),
        # ('GAD High Level Gene-Disease Associations', 'gadhighlevel'),
        # ('GDSC Cell Line Gene Expression Profiles', 'gdsc'),
        # ('GEO Signatures of Differentially Expressed Genes for Diseases', 'geodisease'),
        # ('GEO Signatures of Differentially Expressed Genes for Gene Perturbations', 'geogene'),
        # ('GEO Signatures of Differentially Expressed Genes for Kinase Perturbations', 'geokinase'),
        # ('GEO Signatures of Differentially Expressed Genes for Small Molecules', 'geochemical'),
        # ('GEO Signatures of Differentially Expressed Genes for Transcription Factor Perturbations', 'geotf'),
        # ('GEO Signatures of Differentially Expressed Genes for Viral Infections', 'geovirus'),
        # ('GO Biological Process Annotations', 'gobp'),
        # ('GO Cellular Component Annotations', 'gocc'),
        # ('GO Molecular Function Annotations', 'gomf'),
        # ('GTEx Tissue Gene Expression Profiles', 'gtextissue'),
        # ('GTEx Tissue Sample Gene Expression Profiles', 'gtexsample'),
        # ('GTEx eQTL', 'gtexeqtl'),
        # ('GWAS Catalog SNP-Phenotype Associations', 'gwascatalog'),
        # ('GWASdb SNP-Disease Associations', 'gwasdbdisease'),
        # ('GWASdb SNP-Phenotype Associations', 'gwasdbphenotype'),
        # ('GeneRIF Biological Term Annotations', 'generif'),
        # ('GeneSigDB Published Gene Signatures', 'genesigdb'),
        # ('Graph of Medicine EHR Text-mining Clinical Term Annotations', 'graphofmedicine'),
        # ('Guide to Pharmacology Chemical Ligands of Receptors', 'guidetopharmchemical'),
        # ('Guide to Pharmacology Protein Ligands of Receptors', 'guidetopharmprotein'),
        # ('HMDB Metabolites of Enzymes', 'hmdb'),
        # ('HPA Cell Line Gene Expression Profiles', 'hpacelllines'),
        # ('HPA Tissue Gene Expression Profiles', 'hpatissuesmrna'),
        # ('HPA Tissue Protein Expression Profiles', 'hpatissuesprotein'),
        # ('HPA Tissue Sample Gene Expression Profiles', 'hpasamples'),
        # ('HPM Cell Type and Tissue Protein Expression Profiles', 'hpm'),
        # ('HPO Gene-Disease Associations', 'hpo'),
        # ('HPRD Protein-Protein Interactions', 'hprd'),
        # ('Heiser et al., PNAS, 2011 Cell Line Gene Expression Profiles', 'heiser'),
        # ('HuGE Navigator Gene-Phenotype Associations', 'hugenavigator'),
        # ('Hub Proteins Protein-Protein Interactions', 'hubs'),
        # ('HumanCyc Biomolecular Interactions', 'humancycppi'),
        # ('HumanCyc Pathways', 'humancyc'),
        # ('IntAct Biomolecular Interactions', 'intact'),
        # ('InterPro Predicted Protein Domain Annotations', 'interpro'),
        # ('JASPAR Predicted Transcription Factor Targets', 'jasparpwm'),
        # ('KEA Substrates of Kinases', 'kea'),
        # ('KEGG Biomolecular Interactions', 'keggppi'),
        # ('KEGG Pathways', 'kegg'),
        # ('Kinativ Kinase Inhibitor Bioactivity Profiles', 'kinativ'),
        # ('KinomeScan Kinase Inhibitor Targets', 'kinomescan'),
        # ('Klijn et al., Nat. Biotechnol., 2015 Cell Line Gene CNV Profiles', 'klijncnv'),
        # ('Klijn et al., Nat. Biotechnol., 2015 Cell Line Gene Expression Profiles', 'klijnmrna'),
        # ('Klijn et al., Nat. Biotechnol., 2015 Cell Line Gene Mutation Profiles', 'klijnmut'),
        # ('LINCS L1000 CMAP Signatures of Differentially Expressed Genes for Gene Knockdowns', 'lincscmapgene'),
        # ('LINCS L1000 CMAP Signatures of Differentially Expressed Genes for Small Molecules', 'lincscmapchemical'),
        # ('LOCATE Curated Protein Localization Annotations', 'locate'),
        # ('LOCATE Predicted Protein Localization Annotations', 'locatepredicted'),
        # ('MPO Gene-Phenotype Associations', 'mgimpo'),
        # ('MSigDB Cancer Gene Co-expression Modules', 'msigdbcomp'),
        # ('MSigDB Signatures of Differentially Expressed Genes for Cancer Gene Perturbations', 'msigdbonc'),
        # ('MiRTarBase microRNA Targets', 'mirtarbase'),
        # ('MotifMap Predicted Transcription Factor Targets', 'motifmap'),
        # ('NURSA Protein Complexes', 'nursa'),
        # ('NURSA Protein-Protein Interactions', 'nursappi'),
        # ('OMIM Gene-Disease Associations', 'omim'),
        # ('PANTHER Biomolecular Interactions', 'pantherppi'),
        # ('PANTHER Pathways', 'panther'),
        # ('PID Biomolecular Interactions', 'pidppi'),
        # ('PID Pathways', 'pid'),
        # ('Pathway Commons Protein-Protein Interactions', 'pc'),
        # ('PhosphoSitePlus Phosphosite-Disease Associations', 'phosphositeplusdisease'),
        # ('PhosphoSitePlus Substrates of Kinases', 'phosphositeplus'),
        # ('Phosphosite Textmining Biological Term Annotations', 'phosphositetextmining'),
        # ('ProteomicsDB Cell Type and Tissue Protein Expression Profiles', 'proteomicsdb'),
        # ('Reactome Biomolecular Interactions', 'reactomeppi'),
        # ('Reactome Pathways', 'reactome'),
        # ('Recon X Predicted Biomolecular Interactions', 'reconx'),
        # ('Roadmap Epigenomics Cell and Tissue DNA Accessibility Profiles', 'epigenomicsdnaaccessibility'),
        # ('Roadmap Epigenomics Cell and Tissue DNA Methylation Profiles', 'epigenomicsdnamethylation'),
        # ('Roadmap Epigenomics Cell and Tissue Gene Expression Profiles', 'epigenomicsmrna'),
        # ('Roadmap Epigenomics Histone Modification Site Profiles', 'epigenomicshm'),
        # ('SILAC Phosphoproteomics Signatures of Differentially Phosphorylated Proteins for Drugs', 'silacdrug'),
        # ('SILAC Phosphoproteomics Signatures of Differentially Phosphorylated Proteins for Gene Perturbations', 'silacgene'),
        # ('SILAC Phosphoproteomics Signatures of Differentially Phosphorylated Proteins for Protein Ligands', 'silacligand'),
        # ('SNPedia SNP-Phenotype Associations', 'snpedia'),
        # ('TCGA Signatures of Differentially Expressed Genes for Tumors', 'tcga'),
        # ('TISSUES Curated Tissue Protein Expression Evidence Scores', 'jensentissuecurated'),
        # ('TISSUES Experimental Tissue Protein Expression Evidence Scores', 'jensentissueexpts'),
        # ('TISSUES Text-mining Tissue Protein Expression Evidence Scores', 'jensentissuetextmining'),
        # ('TRANSFAC Curated Transcription Factor Targets', 'transfac'),
        # ('TRANSFAC Predicted Transcription Factor Targets', 'transfacpwm'),
        # ('TargetScan Predicted Conserved microRNA Targets', 'targetscan'),
        # ('TargetScan Predicted Nonconserved microRNA Targets', 'targetscannonconserved'),
        # ('Virus MINT Protein-Viral Protein Interactions', 'virusmintppi'),
        # ('Virus MINT Protein-Virus Interactions', 'virusmint'),
        # ('Wikipathways Pathways', 'wikipathways'),
    ], [
        # 'gene_attribute_matrix.txt.gz',
        # 'gene_attribute_edges.txt.gz',
        # 'gene_set_library_crisp.txt.gz',
        # 'gene_set_library_up_crisp.txt.gz',
        # 'gene_set_library_dn_crisp.txt.gz',
        # 'attribute_set_library_crisp.txt.gz',
        # 'attribute_set_library_up_crisp.txt.gz',
        # 'attribute_set_library_dn_crisp.txt.gz',
        # 'gene_similarity_matrix_cosine.txt.gz',
        # 'attribute_similarity_matrix_cosine.txt.gz',
        # 'gene_list_terms.txt.gz',
        # 'attribute_list_entries.txt.gz',
        # 'processing_script.m'
    ])