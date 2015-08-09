"""Template script to analyze the changes in the genome of of a sample or multiple samples from the .gd file from breseq"""


# define search paths manually
import sys
# dependency dirs
sys.path.append('C:/Users/dmccloskey-sbrg/Documents/GitHub/sequencing_analysis')
sys.path.append('C:/Users/dmccloskey-sbrg/Documents/GitHub/io_utilities')
sys.path.append('C:/Users/dmccloskey-sbrg/Documents/GitHub/sequencing_utilities')
sys.path.append('C:/Users/dmccloskey-sbrg/Documents/GitHub/calculate_utilities')

from sequencing_analysis.genome_diff import genome_diff
from sequencing_analysis.mutations_lineage import mutations_lineage
from sequencing_analysis.mutations_endpoints import mutations_endpoints
from sequencing_analysis.mutations_heatmap import mutations_heatmap
from sequencing_analysis.gff_coverage import gff_coverage

#import and parse .gd files, filter the mutations, and export the filtered mutations tables to .csv
filenames_I = ['C:/Users/dmccloskey-sbrg/Documents/GitHub/sbaas/sbaas/data/tests/analysis_resequencing/breseq/ALEevo04-tpiA/output.gd',
             'C:/Users/dmccloskey-sbrg/Documents/GitHub/sbaas/sbaas/data/tests/analysis_resequencing/breseq/Evo04tpiAEvo01EP/output.gd',
             'C:/Users/dmccloskey-sbrg/Documents/GitHub/sbaas/sbaas/data/tests/analysis_resequencing/breseq/Evo04tpiAEvo02EP/output.gd',
             'C:/Users/dmccloskey-sbrg/Documents/GitHub/sbaas/sbaas/data/tests/analysis_resequencing/breseq/Evo04tpiAEvo03EP/output.gd',
             'C:/Users/dmccloskey-sbrg/Documents/GitHub/sbaas/sbaas/data/tests/analysis_resequencing/breseq/Evo04tpiAEvo04EP/output.gd',
             ];
samplenames = ['Evo04tpiA',
               'Evo04tpiAEvo01EP',
               'Evo04tpiAEvo02EP',
               'Evo04tpiAEvo03EP',
               'Evo04tpiAEvo04EP',
               ];
gd = genome_diff();
for cnt,filename in enumerate(filenames_I):
    gd.import_gd(filename,experiment_id='ALEsKOs01',sample_name=samplenames[cnt]);
    #use default settings to filter the data
    gd.filter_mutations_population(); 
    #annotate the genome using genbank and MG1655
    gd.annotate_mutations(table_I = 'mutationsFiltered',
                          ref_genome_I = 'C:/Users/dmccloskey-sbrg/Documents/GitHub/sbaas/sbaas/data/U00096.2.gb',
                          ref_I = 'genbank',
                          geneReference_I = 'C:/Users/dmccloskey-sbrg/Documents/GitHub/sbaas_workspace/sbaas_workspace/workspace_data/_input/150527_MG1655_geneReference.csv',
                          biologicalmaterial_id_I = 'MG1655')
    #export the filtered and annotated data for later use
    filename_O = samplenames[cnt] + ".csv"
    gd.export_mutationsAnnotated(filename_O);
    gd.clear_data();

#analyze a single lineage
mutationslineage = mutations_lineage(lineage_name_I = 'tpiAEvo01',
                                     lineage_sample_names_I = {0:"Evo04tpiA",1:"Evo04tpiAEvo01EP"});
mutationslineage.import_mutations(['Evo04tpiA.csv','Evo04tpiAEvo01EP.csv']);

#visualize all mutations prior to filtering out mutations that do not persist to the end-point:
mutationslineage.annotate_mutations(table_I = 'mutations',
                          ref_genome_I = 'C:/Users/dmccloskey-sbrg/Documents/GitHub/sbaas/sbaas/data/U00096.2.gb',
                          ref_I = 'genbank',
                          geneReference_I = 'C:/Users/dmccloskey-sbrg/Documents/GitHub/sbaas_workspace/sbaas_workspace/workspace_data/_input/150527_MG1655_geneReference.csv',
                          biologicalmaterial_id_I = 'MG1655')
mutationslineage.export_lineage_js();

#filter out mutations that do not persist to the end-point:
mutationslineage.analyze_lineage_population();
mutationslineage.export_mutationsLineage('lineage.csv');
mutationslineage.annotate_mutations(table_I = 'mutationsLineage',
                          ref_genome_I = 'C:/Users/dmccloskey-sbrg/Documents/GitHub/sbaas/sbaas/data/U00096.2.gb',
                          ref_I = 'genbank',
                          geneReference_I = 'C:/Users/dmccloskey-sbrg/Documents/GitHub/sbaas_workspace/sbaas_workspace/workspace_data/_input/150527_MG1655_geneReference.csv',
                          biologicalmaterial_id_I = 'MG1655')
mutationslineage.export_mutationsAnnotated('lineage_annotated.csv');
mutationslineage.export_lineage_js();

#analyze endpoint replicates
mutationsendpoints = mutations_endpoints(endpoint_name_I = 'tpiAEP',
                                         endpoint_sample_names_I = ['Evo04tpiAEvo01EP','Evo04tpiAEvo02EP','Evo04tpiAEvo03EP','Evo04tpiAEvo04EP']);
mutationsendpoints.import_mutations(['Evo04tpiAEvo01EP.csv','Evo04tpiAEvo02EP.csv','Evo04tpiAEvo03EP.csv','Evo04tpiAEvo04EP.csv']);
mutationsendpoints.analyze_endpoints_population();
mutationsendpoints.export_mutationsEndpoints('endpoints.csv');
mutationsendpoints.annotate_mutations(table_I = 'mutationsEndpoints',
                          ref_genome_I = 'C:/Users/dmccloskey-sbrg/Documents/GitHub/sbaas/sbaas/data/U00096.2.gb',
                          ref_I = 'genbank',
                          geneReference_I = 'C:/Users/dmccloskey-sbrg/Documents/GitHub/sbaas_workspace/sbaas_workspace/workspace_data/_input/150527_MG1655_geneReference.csv',
                          biologicalmaterial_id_I = 'MG1655')
mutationsendpoints.export_mutationsAnnotated('endpoints_annotated.csv');

#generate a heatmap for the strains
heatmap = mutations_heatmap(sample_names_I = samplenames);
heatmap.import_mutations(['Evo04tpiA.csv','Evo04tpiAEvo01EP.csv','Evo04tpiAEvo02EP.csv','Evo04tpiAEvo03EP.csv','Evo04tpiAEvo04EP.csv']);
heatmap.make_heatmap(mutation_id_exclusion_list=['MOB_insA-/-uspC_1977510',
                    'SNP_ylbE_547694',
                    'SNP_yifN_3957957',
                    'DEL_corA_3999668',
                    'MOB_tdk_1292255',
                    'SNP_rpoB_4182566',
                    'INS__4294403',
                    'DEL_pyrE-/-rph_3813882',
                    'SNP_wcaA_2130811'])
heatmap.export_heatmap_js();