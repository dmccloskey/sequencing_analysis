from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData
from .genome_diff_mutations import mutations
from calculate_utilities.base_calculate import base_calculate
import numpy
import json

class mutations_heatmap(mutations):
    def __init__(self,heatmap_I=[],dendrogram_col_I=[],dendrogram_row_I=[],mutations_I=[],sample_names_I=[]):
        '''
        INPUT:
        mutations_I
        sample_names_I
        heatmap_I = heatmap data
        dendrogram_I = dendrogram data
        '''
        if mutations_I:
            self.mutations=mutations_I;
        else:
            self.mutations = [];
        if sample_names_I:
            self.sample_names=sample_names_I;
        else:
            self.sample_names = [];
        if heatmap_I:
            self.heatmap=heatmap_I;
        else:
            self.heatmap = [];
        if dendrogram_col_I:
            self.dendrogram_col=dendrogram_col_I;
        else:
            self.dendrogram_col = [];
        if dendrogram_row_I:
            self.dendrogram_row=dendrogram_row_I;
        else:
            self.dendrogram_row = [];
    def make_heatmap(self, mutation_id_exclusion_list=[],max_position=4000000,
                row_pdist_metric_I='euclidean',row_linkage_method_I='complete',
                col_pdist_metric_I='euclidean',col_linkage_method_I='complete'):
        '''Execute hierarchical cluster on row and column data'''

        print('executing heatmap...');
        calculate = base_calculate();

        # partition into variables:
        mutation_data = self.mutations;
        sample_names = self.sample_names;
        mutation_data_O = [];
        mutation_ids_all = [];
        for end_cnt,mutation in enumerate(mutation_data):
            if int(mutation['mutation_position']) > max_position: #ignore positions great than 4000000
                continue;
            # mutation id
            mutation_id = '';
            mutation_id = self._make_mutationID(mutation['mutation_genes'],mutation['mutation_type'],int(mutation['mutation_position']))
            tmp = {};
            tmp.update(mutation);
            tmp.update({'mutation_id':mutation_id});
            mutation_data_O.append(tmp);
            mutation_ids_all.append(mutation_id);
        mutation_ids_all_unique = list(set(mutation_ids_all));
        mutation_ids = [x for x in mutation_ids_all_unique if not x in mutation_id_exclusion_list];
        # generate the frequency matrix data structure (mutation x intermediate)
        data_O = numpy.zeros((len(sample_names),len(mutation_ids)));
        samples=[];
        # order 2: groups each sample by mutation (intermediate x mutation)
        for sample_name_cnt,sample_name in enumerate(sample_names): #all samples for intermediate j / mutation i
            samples.append(sample_name); # corresponding label from hierarchical clustering
            for mutation_cnt,mutation in enumerate(mutation_ids): #all mutations i for intermediate j
                for row in mutation_data_O:
                    if row['mutation_id'] == mutation and row['sample_name'] == sample_name:
                        data_O[sample_name_cnt,mutation_cnt] = row['mutation_frequency'];
        # generate the clustering for the heatmap
        heatmap_O = [];
        dendrogram_col_O = {};
        dendrogram_row_O = {};
        heatmap_O,dendrogram_col_O,dendrogram_row_O = calculate.heatmap(data_O,samples,mutation_ids,
                row_pdist_metric_I=row_pdist_metric_I,row_linkage_method_I=row_linkage_method_I,
                col_pdist_metric_I=col_pdist_metric_I,col_linkage_method_I=col_linkage_method_I);
        # record the data
        self.heatmap = heatmap_O;
        self.dendrogram_col = dendrogram_col_O;
        self.dendrogram_row = dendrogram_row_O;

    def _make_mutationID(self,mutation_genes,mutation_type,mutation_position):
        '''return a unique mutation id string'''
        mutation_genes_str = '';
        for gene in mutation_genes:
            mutation_genes_str = mutation_genes_str + gene + '-/-'
        mutation_genes_str = mutation_genes_str[:-3];
        mutation_id = mutation_type + '_' + mutation_genes_str + '_' + str(mutation_position);
        return mutation_id;

    def clear_data(self):
        del self.mutations[:];
        del self.heatmap[:];
        del self.dendrogram[:];
        
    def export_heatmap_js(self,data_dir_I="tmp"):
        """export heatmap to js file"""

        #get the heatmap data for the analysis
        data_O = self.heatmap;
        # dump chart parameters to a js files
        data1_keys = [
            'analysis_id',
                      'row_label','col_label','row_index','col_index','row_leaves','col_leaves',
                'col_pdist_metric','row_pdist_metric','col_linkage_method','row_linkage_method',
                'value_units']
        data1_nestkeys = ['analysis_id'];
        data1_keymap = {'xdata':'row_leaves','ydata':'col_leaves','zdata':'value',
                'rowslabel':'row_label','columnslabel':'col_label',
                'rowsindex':'row_index','columnsindex':'col_index',
                'rowsleaves':'row_leaves','columnsleaves':'col_leaves'};
        # make the data object
        dataobject_O = [{"data":data_O,"datakeys":data1_keys,"datanestkeys":data1_nestkeys}];
        # make the tile parameter objects
        svgparameters_O = {"svgtype":'heatmap2d_01',"svgkeymap":[data1_keymap],
                            'svgid':'svg1',
                             'svgcellsize':18,'svgmargin':{ 'top': 200, 'right': 50, 'bottom': 100, 'left': 200 },
                            'svgcolorscale':'quantile',
                            'svgcolorcategory':'heatmap10',
                            'svgcolordomain':[0,1],
                            'svgcolordatalabel':'value',
                            'svgdatalisttileid':'tile1'};
        svgtileparameters_O = {'tileheader':'heatmap','tiletype':'svg','tileid':"tile2",'rowid':"row2",'colid':"col1",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-12"};
        svgtileparameters_O.update(svgparameters_O);
        formtileparameters_O = {'tileheader':'filter menu','tiletype':'html','tileid':"tile1",'rowid':"row1",'colid':"col1",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-4"
            
            };
        formparameters_O = {'htmlid':'datalist1','htmltype':'datalist_01','datalist': [{'value':'hclust','text':'by cluster'},
                            {'value':'probecontrast','text':'by row and column'},
                            {'value':'probe','text':'by row'},
                            {'value':'contrast','text':'by column'},
                            {'value':'custom','text':'by value'}]}
        formtileparameters_O.update(formparameters_O);
        parametersobject_O = [formtileparameters_O,svgtileparameters_O];
        tile2datamap_O = {"tile1":[0],"tile2":[0]};
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