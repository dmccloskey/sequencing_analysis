from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData
from .genome_diff import genome_diff
from .genome_diff_mutations import mutations

class mutations_lineage(mutations):
    def __init__(self,lineage_name_I = None,lineage_sample_names_I={}, mutations_I = [],mutationsLineage_I=[]):
        '''
        INPUT:
        mutations_I
        lineage_name_I = name of the lineage
        lineage_sample_names_I = dict of key/value pairs specifying the intermediate # and the sample_name
                                e.g. {0:sample_name,1:sample_name,2:sample_name,...,n:sample_name}
                                     where n is the end-point strain
        '''
        if mutations_I:
            self.mutations=mutations_I;
        else:
            self.mutations = [];
        if lineage_name_I:
            self.lineage_name=lineage_name_I;
        else:
            self.lineage_name = None;
        if lineage_sample_names_I:
            self.lineage_sample_names=lineage_sample_names_I;
        else:
            self.lineage_sample_names = {};
        if mutationsLineage_I:
            self.mutationsLineage=mutationsLineage_I;
        else:
            self.mutationsLineage=[];

    def analyze_lineage_population(self,lineage_name_I=None,lineage_sample_names_I={}):
        '''Analyze a lineage to identify the following:
        1. conserved mutations
        2. changes in frequency of mutations
        3. hitch-hiker mutations
        Input:
           lineage_name = string,
           lineage_sample_names = {0:sample_name,1:sample_name,2:sample_name,...,n:sample_name}
                                     n is the end-point strain
        Output:
        '''
        if lineage_name_I:
            self.lineage_name = lineage_name_I;
        if lineage_sample_names_I:
            self.lineage_sample_names = lineage_sample_names_I;
        data_O = [];
        mutations = self.mutations;
        strain_lineage = {self.lineage_name:self.lineage_sample_names};
        for lineage_name,strain in strain_lineage.items():
            print('analyzing lineage ' + lineage_name);
            lineage = list(strain.keys());
            end_point = max(lineage)
            # query end data:
            end_mutations = [];
            end_mutations = self._get_mutationsBySampleName(mutations,strain[end_point]);
            intermediates = [i for i in lineage if i!=end_point];
            intermediate_mutations = [];
            for intermediate in intermediates:
                print('analyzing intermediate ' + str(intermediate));
                # query intermediate data:
                intermediate_mutations = [];
                intermediate_mutations = self._get_mutationsBySampleName(mutations,strain[intermediate]);
                # extract mutation lineages
                data_O.extend(self._extract_mutationsLineage(lineage_name,end_mutations,intermediate_mutations,intermediate,end_point));
        # record the data
        self.mutationsLineage = data_O;

    def _extract_mutationsLineage(self,lineage_name,end_mutations,intermediate_mutations,intermediate,end_point):
        """extract out mutation lineages based on the end-point mutation"""
        data_O = [];
        for end_cnt,end_mutation in enumerate(end_mutations):
            print('end mutation type/position ' + end_mutation['mutation_data']['type'] + '/' + str(end_mutation['mutation_data']['position']));
            for inter_cnt,intermediate_mutation in enumerate(intermediate_mutations):
                print('intermediate mutation type/position ' + intermediate_mutation['mutation_data']['type'] + '/' + str(intermediate_mutation['mutation_data']['position']));
                if intermediate == 0 and inter_cnt == 0:
                    #copy end_point data (only once per strain lineage)
                    data_tmp = {};
                    data_tmp['experiment_id'] = end_mutation['experiment_id'];
                    data_tmp['sample_name'] = end_mutation['sample_name'];
                    data_tmp['intermediate'] = end_point;
                    frequency = 1.0;
                    if 'frequency' in end_mutation['mutation_data']: frequency = end_mutation['mutation_data']['frequency'];
                    data_tmp['mutation_frequency'] = frequency
                    data_tmp['mutation_position'] = end_mutation['mutation_data']['position']
                    data_tmp['mutation_type'] = end_mutation['mutation_data']['type']
                    data_tmp['lineage_name'] = lineage_name;
                    data_tmp['mutation_data'] = end_mutation['mutation_data'];
                    data_O.append(data_tmp);
                # find the mutation in the intermediates
                # filter by mutation type-specific criteria
                match = {};
                if end_mutation['mutation_data']['type'] == 'SNP':
                    if end_mutation['mutation_data']['type']==intermediate_mutation['mutation_data']['type'] and \
                        end_mutation['mutation_data']['position']==intermediate_mutation['mutation_data']['position'] and \
                        end_mutation['mutation_data']['new_seq']==intermediate_mutation['mutation_data']['new_seq']:
                        match = intermediate_mutation;
                elif end_mutation['mutation_data']['type'] == 'SUB':
                    if end_mutation['mutation_data']['type']==intermediate_mutation['mutation_data']['type'] and \
                        end_mutation['mutation_data']['position']==intermediate_mutation['mutation_data']['position'] and \
                        end_mutation['mutation_data']['size']==intermediate_mutation['mutation_data']['size'] and \
                        end_mutation['mutation_data']['new_seq']==intermediate_mutation['mutation_data']['new_seq']:
                        match = intermediate_mutation;
                elif end_mutation['mutation_data']['type'] == 'DEL':
                    if end_mutation['mutation_data']['type']==intermediate_mutation['mutation_data']['type'] and \
                        end_mutation['mutation_data']['position']==intermediate_mutation['mutation_data']['position'] and \
                        end_mutation['mutation_data']['size']==intermediate_mutation['mutation_data']['size']:
                        match = intermediate_mutation;
                elif end_mutation['mutation_data']['type'] == 'INS':
                    if end_mutation['mutation_data']['type']==intermediate_mutation['mutation_data']['type'] and \
                        end_mutation['mutation_data']['position']==intermediate_mutation['mutation_data']['position'] and \
                        end_mutation['mutation_data']['new_seq']==intermediate_mutation['mutation_data']['new_seq']:
                        match = intermediate_mutation;
                elif end_mutation['mutation_data']['type'] == 'MOB':
                    if end_mutation['mutation_data']['type']==intermediate_mutation['mutation_data']['type'] and \
                        end_mutation['mutation_data']['repeat_name']==intermediate_mutation['mutation_data']['repeat_name'] and \
                        end_mutation['mutation_data']['strand']==intermediate_mutation['mutation_data']['strand'] and \
                        end_mutation['mutation_data']['duplication_size']==intermediate_mutation['mutation_data']['duplication_size']:
                        match = intermediate_mutation;
                elif end_mutation['mutation_data']['type'] == 'AMP':
                    if end_mutation['mutation_data']['type']==intermediate_mutation['mutation_data']['type'] and \
                        end_mutation['mutation_data']['position']==intermediate_mutation['mutation_data']['position'] and \
                        end_mutation['mutation_data']['size']==intermediate_mutation['mutation_data']['size'] and \
                        end_mutation['mutation_data']['new_copy_number']==intermediate_mutation['mutation_data']['new_copy_number']:
                        match = intermediate_mutation;
                elif end_mutation['mutation_data']['type'] == 'CON':
                    if end_mutation['mutation_data']['type']==intermediate_mutation['mutation_data']['type'] and \
                        end_mutation['mutation_data']['position']==intermediate_mutation['mutation_data']['position'] and \
                        end_mutation['mutation_data']['size']==intermediate_mutation['mutation_data']['size'] and \
                        end_mutation['mutation_data']['region']==intermediate_mutation['mutation_data']['region']:
                        match = intermediate_mutation;
                elif end_mutation['mutation_data']['type'] == 'INV':
                    if end_mutation['mutation_data']['type']==intermediate_mutation['mutation_data']['type'] and \
                        end_mutation['mutation_data']['position']==intermediate_mutation['mutation_data']['position'] and \
                        end_mutation['mutation_data']['size']==intermediate_mutation['mutation_data']['size']:
                        match = intermediate_mutation;
                else:
                    print('unknown mutation type');
                if match:
                    data_tmp = {};
                    data_tmp['experiment_id'] = match['experiment_id'];
                    data_tmp['sample_name'] = match['sample_name'];
                    data_tmp['intermediate'] = intermediate;
                    frequency = 1.0;
                    if 'frequency' in match['mutation_data']: frequency = match['mutation_data']['frequency'];
                    data_tmp['mutation_frequency'] = frequency
                    data_tmp['mutation_position'] = match['mutation_data']['position']
                    data_tmp['mutation_type'] = match['mutation_data']['type']
                    data_tmp['lineage_name'] = lineage_name;
                    data_tmp['mutation_data'] = match['mutation_data'];
                    data_O.append(data_tmp);
        return data_O;
    
    def export_mutationsLineage(self,filename_O):
        """export mutationsLineage"""
        io = base_exportData(self.mutationsLineage);
        io.write_dict2csv(filename_O);

    def clear_data(self):
        del self.mutations[:];
        self.lineage_name = None;
        self.lineage_sample_names = {};
        del self.mutationsLineage[:];

    def import_mutationsLineage(self,filename_I):
        """import mutationsLineage"""
        io = base_importData();
        io.read_csv(filename_I);
        self.mutationsLineage = io.data;