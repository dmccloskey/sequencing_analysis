"""Template script to analyze the coverage of a sample or multiple samples from the .gff file"""

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

#analyze the coverage for a particular strain
gffcoverage = gff_coverage();
gffcoverage.extract_coverage_fromGff(gff_file = '//proline/Users/dmccloskey/Resequencing_DNA/Evo04ptsHIcrrEvo04EP/Evo04ptsHIcrrEvo04EP/data/Evo04ptsHIcrrEvo04EP_reference.gff',
            strand_start = 0,strand_stop = 4640000,
            scale_factor = False,downsample_factor = 2000,
            experiment_id_I = 'ALEsKOs01',
            sample_name_I = 'Evo04ptsHIcrrEvo04EP');
# calculate the coverage statistics
gffcoverage.calculate_coverageStats_fromGff(gff_file = '//proline/Users/dmccloskey/Resequencing_DNA/Evo04ptsHIcrrEvo04EP/Evo04ptsHIcrrEvo04EP/data/Evo04ptsHIcrrEvo04EP_reference.gff',
            strand_start = 0,strand_stop = 4640000,
            scale_factor = False,downsample_factor = 0,
            experiment_id_I = 'ALEsKOs01',
            sample_name_I = 'Evo04ptsHIcrrEvo04EP')
gffcoverage.export_coverageStats('Evo04ptsHIcrrEvo04EP_coverage.csv');
gffcoverage.export_coverage_js();
# find amplifications
gffcoverage.findAndCalculate_amplificationStats_fromGff(gff_file = '//proline/Users/dmccloskey/Resequencing_DNA/Evo04ptsHIcrrEvo04EP/Evo04ptsHIcrrEvo04EP/data/Evo04ptsHIcrrEvo04EP_reference.gff',
            strand_start = 0,strand_stop = 4640000,
            scale_factor = True,downsample_factor = 200,
            reads_min=1.25,reads_max=4.0, indices_min=5000,consecutive_tol=50,
            experiment_id_I = 'ALEsKOs01',
            sample_name_I = 'Evo04ptsHIcrrEvo04EP');
gffcoverage.export_amplificationStats('Evo04ptsHIcrrEvo04EP_amplificationStats.csv');
gffcoverage.annotate_amplifications(ref_genome_I = 'C:/Users/dmccloskey-sbrg/Documents/GitHub/sbaas/sbaas/data/U00096.2.gb',
                          ref_I = 'genbank',
                          geneReference_I = 'C:/Users/dmccloskey-sbrg/Documents/GitHub/sbaas_workspace/sbaas_workspace/workspace_data/_input/150527_MG1655_geneReference.csv',
                          biologicalmaterial_id_I = 'MG1655')
gffcoverage.export_amplificationAnnotations('Evo04ptsHIcrrEvo04EP_amplificationAnnotations.csv');
gffcoverage.export_amplifications_js();