from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData
from .genome_diff import genome_diff
from .genome_diff_mutations import mutations
import json

class mutations_lineage(mutations):
    def __init__(self,lineage_name_I = None,lineage_sample_names_I={},mutations_I = [],mutationsLineage_I=[]):
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
        self.mutationsLineage = self.format_mutationData(io.data);

    def export_lineage_js(self,mutation_id_exclusion_list=[],data_dir_I="tmp"):
        '''export lineage to js file
        Default: will export the mutationsAnnotated table if mutationsLineage have been annotated'''

        intermediates = list(self.lineage_sample_names.keys());
        sample_names = list(self.lineage_sample_names.values());

        # get the lineage information
        lineage_data = [];
        if self.mutationsAnnotated:
            lineage_data = self.mutationsAnnotated;
        else:
            lineage_data = self.mutationsLineage;
        # generate a unique mutation_id if it does not exist
        if not 'mutation_id' in lineage_data[0]:
            for d in lineage_data:
                d['mutation_id'] = self._make_mutationID(d['mutation_genes'],d['mutation_type'],d['mutation_position']);
        mutation_ids = [x['mutation_id'] for x in lineage_data];
        mutation_ids_screened = [x for x in mutation_ids if x not in mutation_id_exclusion_list];
        mutation_ids_unique = list(set(mutation_ids_screened));
        # get mutation information for all unique mutations
        mutation_ids_uniqueInfo = [];
        for mutation_id in mutation_ids_unique:
            for mutation in lineage_data:
                if mutation_id == mutation['mutation_id']:
                    tmp = {};
                    tmp['mutation_id']=mutation['mutation_id'];
                    tmp['mutation_frequency']=mutation['mutation_frequency'];
                    if mutation['mutation_genes']:
                        tmp['mutation_genes']=";".join([x for x in mutation['mutation_genes'] if x is not None]);
                    else: tmp['mutation_genes']=mutation['mutation_genes'];
                    if mutation['mutation_position']:
                        tmp['mutation_position']=mutation['mutation_position'];
                    else: tmp['mutation_position']=mutation['mutation_position'];
                    if mutation['mutation_annotations']:
                        tmp['mutation_annotations']=";".join([x for x in mutation['mutation_annotations'] if x is not None]);
                    else: tmp['mutation_annotations']=mutation['mutation_annotations'];
                    if mutation['mutation_locations']:
                        tmp['mutation_locations']=";".join([x for x in mutation['mutation_locations'] if x is not None]);
                    else: tmp['mutation_locations']=mutation['mutation_locations'];
                    if mutation['mutation_links']:
                        tmp['mutation_links']=";".join([x for x in mutation['mutation_links'] if x is not None]);
                    else: tmp['mutation_links']=mutation['mutation_links'];
                    tmp['mutation_type']=mutation['mutation_type'];
                    tmp['used_']=True;
                    tmp['comment_']=None;
                    mutation_ids_uniqueInfo.append(tmp);          
        data_O = [];
        # add in 0.0 frequency for mutations that are not found
        for sample_name_cnt,sample_name in enumerate(sample_names):
            for mutation_id in mutation_ids_uniqueInfo:
                tmp = {};
                tmp_fitted = {};
                tmp['mutation_id']=mutation_id['mutation_id']
                tmp['intermediate']=intermediates[sample_name_cnt]
                #tmp['experiment_id']=experiment_ids[sample_name_cnt]
                tmp['sample_name']=sample_name
                tmp['mutation_frequency']=0.0;  
                tmp['mutation_genes']=mutation_id['mutation_genes'];
                tmp['mutation_position']=mutation_id['mutation_position'];
                tmp['mutation_annotations']=mutation_id['mutation_annotations'];
                tmp['mutation_locations']=mutation_id['mutation_locations'];
                tmp['mutation_links']=mutation_id['mutation_links'];
                tmp['mutation_type']=mutation_id['mutation_type'];
                tmp['used_']=mutation_id['used_'];
                for mutation in lineage_data:
                    if sample_name == mutation['sample_name'] and mutation_id['mutation_id'] == mutation['mutation_id']:
                        tmp['mutation_frequency']=mutation['mutation_frequency'];
                        break;
                data_O.append(tmp);
        # dump chart parameters to a js files
        data1_keys = [
                    #'experiment_id',
                    #'lineage_name',
                    'sample_name',
                    'mutation_id',
                    #'mutation_frequency',
                    'mutation_type',
                    'mutation_position',
                    #'mutation_data',
                    #'mutation_annotations',
                    'mutation_genes',
                    #'mutation_links',
                    'mutation_locations'
                    ];
        data1_nestkeys = ['mutation_id'];
        data1_keymap = {'xdata':'intermediate',
                        'ydata':'mutation_frequency',
                        'serieslabel':'mutation_id',
                        'featureslabel':''};
        # make the data object
        dataobject_O = [{"data":data_O,"datakeys":data1_keys,"datanestkeys":data1_nestkeys},{"data":data_O,"datakeys":data1_keys,"datanestkeys":data1_nestkeys}];
        # make the tile parameter objects
        formtileparameters_O = {'tileheader':'Filter menu','tiletype':'html','tileid':"filtermenu1",'rowid':"row1",'colid':"col1",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-4"};
        formparameters_O = {'htmlid':'filtermenuform1',"htmltype":'form_01',"formsubmitbuttonidtext":{'id':'submit1','text':'submit'},"formresetbuttonidtext":{'id':'reset1','text':'reset'},"formupdatebuttonidtext":{'id':'update1','text':'update'}};
        formtileparameters_O.update(formparameters_O);
        svgparameters_O = {"svgtype":'scatterlineplot2d_01',"svgkeymap":[data1_keymap,data1_keymap],
                            'svgid':'svg1',
                            "svgmargin":{ 'top': 50, 'right': 150, 'bottom': 50, 'left': 50 },
                            "svgwidth":500,"svgheight":350,
                            "svgx1axislabel":"intermediate","svgy1axislabel":"frequency",
    						'svgformtileid':'filtermenu1','svgresetbuttonid':'reset1','svgsubmitbuttonid':'submit1'};
        svgtileparameters_O = {'tileheader':'Population mutation frequency','tiletype':'svg','tileid':"tile2",'rowid':"row1",'colid':"col2",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-8"};
        svgtileparameters_O.update(svgparameters_O);
        tableparameters_O = {"tabletype":'responsivetable_01',
                    'tableid':'table1',
                    "tablefilters":None,
                    "tableclass":"table  table-condensed table-hover",
    			    'tableformtileid':'filtermenu1','tableresetbuttonid':'reset1','tablesubmitbuttonid':'submit1'};
        tabletileparameters_O = {'tileheader':'Population mutation frequency','tiletype':'table','tileid':"tile3",'rowid':"row2",'colid':"col1",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-12"};
        tabletileparameters_O.update(tableparameters_O);
        parametersobject_O = [formtileparameters_O,svgtileparameters_O,tabletileparameters_O];
        tile2datamap_O = {"filtermenu1":[0],"tile2":[0,1],"tile3":[0]};
        # dump the data to a json file
        data_str = 'var ' + 'data' + ' = ' + json.dumps(dataobject_O) + ';';
        parameters_str = 'var ' + 'parameters' + ' = ' + json.dumps(parametersobject_O) + ';';
        tile2datamap_str = 'var ' + 'tile2datamap' + ' = ' + json.dumps(tile2datamap_O) + ';';
        if data_dir_I=='tmp':
            filename_str = 'ddt_data.js'
        elif data_dir_I=='data_json':
            data_json_O = data_str + '\n' + parameters_str + '\n' + tile2datamap_str;
            return data_json_O;
        with open(filename_str,'w') as file:
            file.write(data_str);
            file.write(parameters_str);
            file.write(tile2datamap_str);