from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData
from .genes_fpkm_tracking import genes_fpkm_tracking
import copy

class gene_exp_diff(genes_fpkm_tracking):
    def __init__(self,geneExpDiff_I=None,genesFpkmTracking_I=None):
        if geneExpDiff_I:
            self.geneExpDiff = geneExpDiff_I;
        else:
            self.geneExpDiff = [];
        if genesFpkmTracking_I:
            self.genesFpkmTracking = genesFpkmTracking_I;
        else:
            self.genesFpkmTracking = [];

    def import_geneExpDiff(self,filename_I, experiment_id_1_I=None,experiment_id_2_I=None,sample_name_abbreviation_1_I=None,sample_name_abbreviation_2_I=None):
        """import geneExpDiff
        INPUT:
        filename_I = input filename
        
        OPTIONAL INPUT:
        the following are optional for analyzing a single sample,
        but required when analyzing multiple samples

        experiment_id_1_I = string, name of the experiment that generated the samples
        experiment_id_I = string, name of the experiment that generated the samples
        sample_name_abbreviation_1_I = string, name of the sample
        sample_name_abbreviation_2_I = string, name of the sample
        """
        io = base_importData();
        io.read_tab(filename_I);
        geneExpDiff = self.format_geneExpDiff(io.data);
        for d in geneExpDiff:
            d['experiment_id_1'] = experiment_id_1_I;
            d['experiment_id_2'] = experiment_id_2_I;
            d['sample_name_abbreviation_1'] = sample_name_abbreviation_1_I;
            d['sample_name_abbreviation_2'] = sample_name_abbreviation_2_I;
        self.geneExpDiff = geneExpDiff;

    def import_genesFpkmTracking(self,filename_I,analysis_id_I = None,experiment_id_1_I=None,experiment_id_2_I=None,sample_name_abbreviation_1_I=None,sample_name_abbreviation_2_I=None):
        """import genes.fpkm_tracking from cuffdiff
        INPUT:
        filename_I = input filename
        
        OPTIONAL INPUT:
        the following are optional for analyzing a single sample,
        but required when analyzing multiple samples

        experiment_id_1_I = string, name of the experiment that generated the samples
        experiment_id_I = string, name of the experiment that generated the samples
        sample_name_abbreviation_1_I = string, name of the sample
        sample_name_abbreviation_2_I = string, name of the sample
        """
        io = base_importData();
        io.read_tab(filename_I);
        genesFpkmTracking = self.format_geneExpDiff(io.data);
        genesFpkmTracking_O = [];
        for d in genesFpkmTracking:
            tmp = copy.copy(d);
            #1
            tmp['experiment_id'] = experiment_id_1_I;
            tmp['analysis_id'] = analysis_id_I;
            tmp['FPKM'] = d[sample_name_abbreviation_1_I+'_FPKM'];
            tmp['FPKM_conf_lo'] = d[sample_name_abbreviation_1_I+'_conf_lo'];
            tmp['FPKM_conf_hi'] = d[sample_name_abbreviation_1_I+'_conf_hi'];
            tmp['FPKM_status'] = d[sample_name_abbreviation_1_I+'_status'];
            tmp['sample_name_abbreviation'] = sample_name_abbreviation_1_I;
            tmp['used_'] = True;
            tmp['comment_'] = None;
            genesFpkmTracking_O.append(tmp);
            #2
            tmp = copy.copy(d);
            tmp['analysis_id'] = analysis_id_I;
            tmp['experiment_id'] = experiment_id_2_I;
            tmp['FPKM'] = d[sample_name_abbreviation_2_I+'_FPKM'];
            tmp['FPKM_conf_lo'] = d[sample_name_abbreviation_2_I+'_conf_lo'];
            tmp['FPKM_conf_hi'] = d[sample_name_abbreviation_2_I+'_conf_hi'];
            tmp['FPKM_status'] = d[sample_name_abbreviation_2_I+'_status'];
            tmp['sample_name_abbreviation'] = sample_name_abbreviation_2_I;
            tmp['used_'] = True;
            tmp['comment_'] = None;
            genesFpkmTracking_O.append(tmp);
        self.genesFpkmTracking = self.format_genesFpkmTracking(genesFpkmTracking_O);

    def export_geneExpDiff(self,filename_O):
        """export geneExpDiff"""
        io = base_exportData(self.geneExpDiff);
        io.write_dict2csv(filename_O);

    def format_nanAndInf(self,str_I):
        """Converts +/-nan to None and +/-inf to None"""
        data_O = None;
        if str_I == '-nan':
            data_O = None;
        elif str_I == 'nan':
            data_O = None;
        elif str_I == '-inf':
            data_O = None;
        elif str_I == 'inf':
            data_O = None;
        return data_O;
    
    def format_geneExpDiff(self,fpkmTracking_I):
        """formats raw string input into their appropriate values"""
        for fpkmTracking in fpkmTracking_I:
            if 'log2(fold_change)' in fpkmTracking and type(fpkmTracking['log2(fold_change)'])==type('string'):
                if 'nan' in fpkmTracking['log2(fold_change)'] or 'inf' in fpkmTracking['log2(fold_change)']:
                    fpkmTracking['log2(fold_change)'] = self.format_nanAndInf(fpkmTracking['log2(fold_change)']);
                else:
                    fpkmTracking['log2(fold_change)'] = eval(fpkmTracking['log2(fold_change)']);
            if 'test_stat' in fpkmTracking and type(fpkmTracking['test_stat'])==type('string'):
                if 'nan' in fpkmTracking['test_stat'] or 'inf' in fpkmTracking['test_stat']:
                    fpkmTracking['test_stat'] = self.format_nanAndInf(fpkmTracking['test_stat']);
                else:
                    fpkmTracking['test_stat'] = eval(fpkmTracking['test_stat']);
            if 'value_1' in fpkmTracking and type(fpkmTracking['value_1'])==type('string'):
                fpkmTracking['value_1'] = eval(fpkmTracking['value_1']);
            if 'value_2' in fpkmTracking and type(fpkmTracking['value_2'])==type('string'):
                fpkmTracking['value_2'] = eval(fpkmTracking['value_2']);
            if 'p_value' in fpkmTracking and type(fpkmTracking['p_value'])==type('string'):
                fpkmTracking['p_value'] = eval(fpkmTracking['p_value']);
            if 'q_value' in fpkmTracking and type(fpkmTracking['q_value'])==type('string'):
                fpkmTracking['q_value'] = eval(fpkmTracking['q_value']);
        return fpkmTracking_I;
        
    def export_geneExpDiff_js(self,data_dir_I='tmp'):
        '''Export data for a volcano plot'''
        
        #get the data for the analysis
        data_O = [];
        data_O = self.geneExpDiff;
        # transform p_value to -log10(p_value)
        for d in data_O:
            if not '-log10(p_value)' in d:
                d['-log10(p_value)']=-log(d['p_value'],10.0);
        # make the data parameters
        data1_keys = ['experiment_id_1','experiment_id_2','sample_name_abbreviation_1','sample_name_abbreviation_2','gene','log2(fold_change)','-log10(p_value)','significant'
                    ];
        data1_nestkeys = ['experiment_id_1'];
        data1_keymap = {'ydata':'-log10(p_value)',
                        'xdata':'log2(fold_change)',
                        'serieslabel':'',
                        'featureslabel':'gene'};
        # make the data object
        dataobject_O = [{"data":data_O,"datakeys":data1_keys,"datanestkeys":data1_nestkeys}];
        # make the tile parameter objects
        formtileparameters_O = {'tileheader':'Filter menu','tiletype':'html','tileid':"filtermenu1",'rowid':"row1",'colid':"col1",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-4"};
        formparameters_O = {'htmlid':'filtermenuform1',"htmltype":'form_01',"formsubmitbuttonidtext":{'id':'submit1','text':'submit'},"formresetbuttonidtext":{'id':'reset1','text':'reset'},"formupdatebuttonidtext":{'id':'update1','text':'update'}};
        formtileparameters_O.update(formparameters_O);
        svgparameters_O = {"svgtype":'volcanoplot2d_01',"svgkeymap":[data1_keymap],
                            'svgid':'svg1',
                            "svgmargin":{ 'top': 50, 'right': 50, 'bottom': 50, 'left': 50 },
                            "svgwidth":500,"svgheight":350,
                            "svgx1axislabel":'Fold Change [log2(FC)]',"svgy1axislabel":'Probability [-log10(P)]',
    						'svgformtileid':'filtermenu1','svgresetbuttonid':'reset1','svgsubmitbuttonid':'submit1'};
        svgtileparameters_O = {'tileheader':'Volcano plot','tiletype':'svg','tileid':"tile2",'rowid':"row1",'colid':"col2",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-8"};
        svgtileparameters_O.update(svgparameters_O);
        tableparameters_O = {"tabletype":'responsivetable_01',
                    'tableid':'table1',
                    "tablefilters":None,
                    "tableclass":"table  table-condensed table-hover",
    			    'tableformtileid':'filtermenu1','tableresetbuttonid':'reset1','tablesubmitbuttonid':'submit1'};
        tabletileparameters_O = {'tileheader':'pairWiseTest','tiletype':'table','tileid':"tile3",'rowid':"row2",'colid':"col1",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-12"};
        tabletileparameters_O.update(tableparameters_O);
        parametersobject_O = [formtileparameters_O,svgtileparameters_O,tabletileparameters_O];
        tile2datamap_O = {"filtermenu1":[0],"tile2":[0],"tile3":[0]};
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