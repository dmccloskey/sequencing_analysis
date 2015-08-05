from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData
from .genome_diff import genome_diff
from .genome_diff_mutations import mutations

class mutations_endpoints(mutations):
    def __init__(self,endpoint_name_I = None, endpoint_sample_names_I=[], mutations_I = [],mutationsEndpoints_I=[]):
        '''
        INPUT:
        endpoint_name_I = name of the endpoint grouping
        endpoint_name_I = list of strains
        '''
        if mutations_I:
            self.mutations=mutations_I;
        else:
            self.mutations = [];
        if endpoint_name_I:
            self.endpoint_name=endpoint_name_I;
        else:
            self.endpoint_name = None;
        if endpoint_sample_names_I:
            self.endpoint_sample_names=endpoint_sample_names_I;
        else:
            self.endpoint_sample_names = [];
        if mutationsEndpoints_I:
            self.mutationsEndpoints=mutationsEndpoints_I;
        else:
            self.mutationsEndpoints=[];
    
    def analyze_endpoints_population(self,endpoint_name_I=None,endpoint_sample_names_I=[]):
        '''Analyze endpoints to identify the following:
        1. conserved mutations among replicates
        2. unique mutations among replicates
        Input:
           endpoint_name = string 
           endpoint_sample_names = [sample_name_1,sample_name_2,sample_name_3,...]
        #Output:
        '''
        if endpoint_name_I:
            self.endpoint_name = endpoint_name_I;
        if endpoint_sample_names_I:
            self.endpoint_sample_names = endpoint_sample_names_I;
        print('Executing analyzeEndpoints_population...')
        endpoint_name = self.endpoint_name;
        strains = self.endpoint_sample_names;
        mutations = self.mutations;
        data_O = [];

        print('analyzing endpoint ' + endpoint_name);
        analyzed_strain1 = []; # strain1s that have been analyzed
        matched_mutations = {};
        for strain1 in strains:
            # query strain 1 data:
            strain1_mutations = [];
            strain1_mutations = self._get_mutationsBySampleName(mutations,strain1);
            analyzed_strain1.append(strain1);
            analyzed_strain1_mutations = []; # mutations from strain 1 that have been analyzed
            analyzed_strain2_mutations_all = []; # all mutations from strain 2 that have been analyzed
            strain2_cnt = 0;
            for strain2 in strains:
                if strain2 == strain1: continue; # do not compare the same strain to itself
                print('comparing ' + strain1 + ' to ' + strain2);
                # query strain 1 data:
                strain2_mutations = [];
                strain2_mutations = self._get_mutationsBySampleName(mutations,strain2);
                analyzed_strain2_mutations = []; # mutations from strain 2 that have been analyzed
                # extract common mutations
                analyzed_strain1_mutations_tmp = [];
                analyzed_strain2_mutations_tmp = [];
                matched_mutations_tmp = {}
                data_tmp = [];
                matched_mutations_tmp,\
                    analyzed_strain1_mutations_tmp,\
                    analyzed_strain2_mutations_tmp,\
                    data_tmp= self._extract_commonMutations(matched_mutations,\
                        analyzed_strain1_mutations,\
                        analyzed_strain2_mutations,\
                        strain1_mutations,strain2_mutations,
                        strain1,strain2_cnt,endpoint_name);
                
                analyzed_strain1_mutations.extend(analyzed_strain1_mutations_tmp)
                analyzed_strain2_mutations.extend(analyzed_strain2_mutations_tmp)
                analyzed_strain2_mutations_all.append(analyzed_strain2_mutations);
                matched_mutations.update(matched_mutations_tmp);
                data_O.extend(data_tmp);
                #update the strain2 counter
                strain2_cnt+=1
            # extract unique mutations
            data_tmp = [];
            data_tmp = self._extract_uniqueMutations(analyzed_strain1_mutations,analyzed_strain2_mutations_all,strain1_mutations,endpoint_name);
            data_O.extend(data_tmp);
        # record the data
        self.mutationsEndpoints = data_O;

    def _extract_uniqueMutations(self,analyzed_strain1_mutations,analyzed_strain2_mutations_all,strain1_mutations,endpoint_name):
        '''Extract out unique mutations'''
        data_O = [];
        for analyzed_strain1_mutation in analyzed_strain1_mutations:
            isUnique_bool = True;
            isConserved_cnt = 0;
            for analyzed_strain2_mutations_cnt,analyzed_strain2_mutations in enumerate(analyzed_strain2_mutations_all):
                for analyzed_strain2_mutation in analyzed_strain2_mutations:
                    if analyzed_strain1_mutation == analyzed_strain2_mutation:
                        isUnique_bool = False;
                        isConserved_cnt += 1;
            if isUnique_bool:
                for strain1_mutation_cnt,strain1_mutation in enumerate(strain1_mutations):
                    if (strain1_mutation['mutation_data']['type'],strain1_mutation['mutation_data']['position'])==analyzed_strain1_mutation:
                        data_tmp = {};
                        data_tmp['experiment_id'] = strain1_mutation['experiment_id'];
                        data_tmp['sample_name'] = strain1_mutation['sample_name'];
                        frequency = 1.0;
                        if 'frequency' in strain1_mutation['mutation_data']: frequency = strain1_mutation['mutation_data']['frequency'];
                        data_tmp['mutation_frequency'] = frequency
                        data_tmp['mutation_position'] = strain1_mutation['mutation_data']['position']
                        data_tmp['mutation_type'] = strain1_mutation['mutation_data']['type']
                        data_tmp['endpoint_name'] = endpoint_name;
                        data_tmp['mutation_data'] = strain1_mutation['mutation_data'];
                        data_tmp['isUnique'] = True;
                        data_O.append(data_tmp);
        return data_O;
    def _extract_commonMutations(self,matched_mutations,
                        analyzed_strain1_mutations,
                        analyzed_strain2_mutations,
                        strain1_mutations,strain2_mutations,
                        strain1,
                        strain2_cnt,
                        endpoint_name):
        '''extract out mutations that are in common between strains'''
        data_O = [];
        for strain1_mutation_cnt,strain1_mutation in enumerate(strain1_mutations):
            print('strain1 mutation type/position ' + strain1_mutation['mutation_data']['type'] + '/' + str(strain1_mutation['mutation_data']['position']));
            if strain2_cnt == 0: # record strain 1 mutations only once for all strain 2 mutations
                analyzed_strain1_mutations.append((strain1_mutation['mutation_data']['type'],strain1_mutation['mutation_data']['position']));
            for strain2_mutation_cnt,strain2_mutation in enumerate(strain2_mutations):
                print('strain2 mutation type/position ' + strain2_mutation['mutation_data']['type'] + '/' + str(strain2_mutation['mutation_data']['position']));
                if strain2_mutation_cnt == 0 and \
                    (strain1,strain1_mutation['mutation_data']['type'],strain1_mutation['mutation_data']['position']) not in matched_mutations:
                    matched_mutations[(strain1,strain1_mutation['mutation_data']['type'],strain1_mutation['mutation_data']['position'])] = 0;
                # find the mutations that are common to strain1 and strain2
                # filter by mutation type-specific criteria
                match = {};
                if strain1_mutation['mutation_data']['type'] == 'SNP':
                    if strain1_mutation['mutation_data']['type']==strain2_mutation['mutation_data']['type'] and \
                        strain1_mutation['mutation_data']['position']==strain2_mutation['mutation_data']['position'] and \
                        strain1_mutation['mutation_data']['new_seq']==strain2_mutation['mutation_data']['new_seq']:
                        match = strain1_mutation;
                elif strain1_mutation['mutation_data']['type'] == 'SUB':
                    if strain1_mutation['mutation_data']['type']==strain2_mutation['mutation_data']['type'] and \
                        strain1_mutation['mutation_data']['position']==strain2_mutation['mutation_data']['position'] and \
                        strain1_mutation['mutation_data']['size']==strain2_mutation['mutation_data']['size'] and \
                        strain1_mutation['mutation_data']['new_seq']==strain2_mutation['mutation_data']['new_seq']:
                        match = strain1_mutation;
                elif strain1_mutation['mutation_data']['type'] == 'DEL':
                    if strain1_mutation['mutation_data']['type']==strain2_mutation['mutation_data']['type'] and \
                        strain1_mutation['mutation_data']['position']==strain2_mutation['mutation_data']['position'] and \
                        strain1_mutation['mutation_data']['size']==strain2_mutation['mutation_data']['size']:
                        match = strain1_mutation;
                elif strain1_mutation['mutation_data']['type'] == 'INS':
                    if strain1_mutation['mutation_data']['type']==strain2_mutation['mutation_data']['type'] and \
                        strain1_mutation['mutation_data']['position']==strain2_mutation['mutation_data']['position'] and \
                        strain1_mutation['mutation_data']['new_seq']==strain2_mutation['mutation_data']['new_seq']:
                        match = strain1_mutation;
                elif strain1_mutation['mutation_data']['type'] == 'MOB':
                    if strain1_mutation['mutation_data']['type']==strain2_mutation['mutation_data']['type'] and \
                        strain1_mutation['mutation_data']['position']==strain2_mutation['mutation_data']['position'] and \
                        strain1_mutation['mutation_data']['repeat_name']==strain2_mutation['mutation_data']['repeat_name'] and \
                        strain1_mutation['mutation_data']['strand']==strain2_mutation['mutation_data']['strand'] and \
                        strain1_mutation['mutation_data']['duplication_size']==strain2_mutation['mutation_data']['duplication_size']:
                        match = strain1_mutation;
                elif strain1_mutation['mutation_data']['type'] == 'AMP':
                    if strain1_mutation['mutation_data']['type']==strain2_mutation['mutation_data']['type'] and \
                        strain1_mutation['mutation_data']['position']==strain2_mutation['mutation_data']['position'] and \
                        strain1_mutation['mutation_data']['size']==strain2_mutation['mutation_data']['size'] and \
                        strain1_mutation['mutation_data']['new_copy_number']==strain2_mutation['mutation_data']['new_copy_number']:
                        match = strain1_mutation;
                elif strain1_mutation['mutation_data']['type'] == 'CON':
                    if strain1_mutation['mutation_data']['type']==strain2_mutation['mutation_data']['type'] and \
                        strain1_mutation['mutation_data']['position']==strain2_mutation['mutation_data']['position'] and \
                        strain1_mutation['mutation_data']['size']==strain2_mutation['mutation_data']['size'] and \
                        strain1_mutation['mutation_data']['region']==strain2_mutation['mutation_data']['region']:
                        match = strain1_mutation;
                elif strain1_mutation['mutation_data']['type'] == 'INV':
                    if strain1_mutation['mutation_data']['type']==strain2_mutation['mutation_data']['type'] and \
                        strain1_mutation['mutation_data']['position']==strain2_mutation['mutation_data']['position'] and \
                        strain1_mutation['mutation_data']['size']==strain2_mutation['mutation_data']['size']:
                        match = strain1_mutation;
                else:
                    print('unknown mutation type');
                if match and \
                        matched_mutations[(strain1,strain1_mutation['mutation_data']['type'],strain1_mutation['mutation_data']['position'])] == 0:
                    # check that the mutation combination and pairs of strains have not already been analyzed
                    data_tmp = {};
                    data_tmp['experiment_id'] = match['experiment_id'];
                    data_tmp['sample_name'] = match['sample_name'];
                    frequency = 1.0;
                    if 'frequency' in match['mutation_data']: frequency = match['mutation_data']['frequency'];
                    data_tmp['mutation_frequency'] = frequency
                    data_tmp['mutation_position'] = match['mutation_data']['position']
                    data_tmp['mutation_type'] = match['mutation_data']['type']
                    data_tmp['endpoint_name'] = endpoint_name;
                    data_tmp['mutation_data'] = match['mutation_data'];
                    data_tmp['isUnique'] = False;
                    data_O.append(data_tmp);
                    matched_mutations[(strain1,strain1_mutation['mutation_data']['type'],strain1_mutation['mutation_data']['position'])] += 1;
                if strain1_mutation_cnt == 0: # record strain 2 mutations only once for all strain 1 mutations
                    analyzed_strain2_mutations.append((strain2_mutation['mutation_data']['type'],strain2_mutation['mutation_data']['position']));
        return matched_mutations,analyzed_strain1_mutations,analyzed_strain2_mutations,data_O;

    def clear_data(self):
        del self.mutations[:];
        self.endpoint_name = None;
        del self.endpoint_sample_names[:];
        del self.mutationsEndpoints[:];
    
    def export_mutationsEndpoints(self,filename_O):
        """export mutationsEndpoints"""
        io = base_exportData(self.mutationsEndpoints);
        io.write_dict2csv(filename_O);

    def import_mutationsEndpoints(self,filename_I):
        """import mutationsEndpoints"""
        io = base_importData();
        io.read_csv(filename_I);
        self.mutationsEndpoints = io.data;