from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData
from .genome_diff import genome_diff

class mutations(genome_diff):
    #def __init__(self,mutations_I = []):
    #    if mutations_I:
    #        self.mutations=mutations_I;
    #    else:
    #        self.mutations = [];
    
    #def format_mutationData(self,mutationData_I):
    #    """converts '{}' to {}"""
    #    for mutationData in mutationData_I:
    #        if type(mutationData['mutation_data'])==type('string'):
    #            mutationData['mutation_data'] = eval(mutationData['mutation_data']);
    #    return mutationData;

    def import_mutations(self,filenames_I=[]):
        """import mutations for multiple experiments/samples"""
        io = base_importData();
        for filename in filenames_I:
            io.read_csv(filename);
            self.mutations.extend(self.format_mutationData(io.data));
            io.clear_data();

    def _get_mutationsBySampleName(self,data_I,sample_name_I):
        """return mutations that match the sample_name"""
        mutations_O = [];
        for d in data_I:
            if d['sample_name']==sample_name_I:
                mutations_O.append(d);
        return mutations_O;

    def _get_mutationsByExperimentIDAndSampleName(self,data_I,experiment_id_I,sample_name_I):
        """return mutations that match the experiment_id and sample_name"""
        mutations_O = [];
        for d in data_I:
            if d['experiment_id']==experiment_id_I and d['sample_name']==sample_name_I:
                mutations_O.append(d);
        return mutations_O;