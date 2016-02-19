from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData
from .genes_fpkm_tracking import genes_fpkm_tracking
from python_statistics.calculate_interface import calculate_interface

class fpkms(genes_fpkm_tracking):
    def __init__(self,genesFpkmTracking_I=None,genesFpkmTrackingStats_I=[]):
        if genesFpkmTracking_I:
            self.genesFpkmTracking = genesFpkmTracking_I;
        else:
            self.genesFpkmTracking = [];
        if genesFpkmTrackingStats_I:
            self.genesFpkmTrackingStats = genesFpkmTrackingStats_I;
        else:
            self.genesFpkmTrackingStats = [];

    def import_genesFpkmTracking(self,filename_I,sample_name_abbreviation_I):
        """import geneExpDiff
        INPUT:
        filename_I = list of input filenames (.csv output from genes_fpkm_tracking.export_genesFpkmTracking;
        
        OPTIONAL INPUT:
        the following are optional for analyzing a single sample,
        but required when analyzing multiple samples

        sample_name_abbreviation_I = list of strings, sample name abbreviation indicating which sample_names/filenames belong to which replicate group
                                    (default = '')
        """
        io = base_importData();
        for cnt,filename in enumerate(filename_I):
            io.read_csv(filename);
            genesFpkmTracking = self.format_genesFpkmTracking(io.data);
            for d in genesFpkmTracking:
                if sample_name_abbreviation_I:
                    d['sample_name_abbreviation'] = sample_name_abbreviation_I[cnt];
                else:
                    d['sample_name_abbreviation'] = '';
            self.genesFpkmTracking.extend(genesFpkmTracking);

    def export_genesFpkmTracking(self,filename_O):
        """export genesFpkmTracking"""
        io = base_exportData(self.genesFpkmTracking);
        io.write_dict2csv(filename_O);

    def calculate_genesFpkmTrackingStats(self,
                experiment_id_I = None,
                sample_name_I = None,):
        """calculate statistics of replicate samples from genesFpkmTracking

        INPUT:
        OPTION INPUT:
        experiment_id_I = limiter for the experiment_id
        sample_name_I = limiter for the sample_name
        
        """
        
        data_O=[];
        stats_O=[];
        experiment_id = experiment_id_I;
        sn = sample_name_I;
        genesFpkmTracking = self.genesFpkmTracking;
        calculate = calculate_interface();
        # get the uniqueSampleNameAbbreviations
        sna_unique = self._get_uniqueSampleNameAbbreviations();
        for sna in sna_unique:
            data_tmp = [];
            data_tmp = self._get_rowsBySampleNameAbbreviation(sna);
            # calculate using scipy
            data_ave_O, data_var_O, data_lb_O, data_ub_O = None, None, None, None;
            data_ave_O, data_var_O, data_lb_O, data_ub_O = calculate.calculate_ave_var(data_tmp,confidence_I = 0.95);
            # calculate the interquartile range
            min_O, max_O, median_O, iq_1_O, iq_3_O = None, None, None, None, None;
            min_O, max_O, median_O, iq_1_O, iq_3_O=calculate.calculate_interquartiles(data_tmp);

    def _get_uniqueSampleNameAbbreviations(self,data_I,experiment_id_I=None,sample_name_I=[]):
        """return unique sample_name_abbreviations
       OPTION: that match the experiment_id and sample_name"""
        sna_all = [];
        for d in data_I:
            if experiment_id_I and sample_name_I and d['experiment_id']==experiment_id_I and d['sample_name']==sample_name_I:
                sna_all.append(d['sample_name_abbreviation']);
            else:
                sna_all.append(d['sample_name_abbreviation']);
        sna_unique = [];
        sna_unique = list(set(sna_all));
        sna_unique.sort();
        return sna_unique;

    def _get_rowsBySampleNameAbbreviation(self,sample_name_abbreviation_I,data_I,experiment_id_I=None,sample_name_I=[]):
        """return unique sample_name_abbreviations
        INPUT:
        sample_name_abbreviation_I = string;       
        OPTION: that match the experiment_id and sample_name"""
        data_O = [];
        for d in data_I:
            if experiment_id_I and sample_name_I and d['experiment_id']==experiment_id_I and d['sample_name']==sample_name_I:
                if d['sample_name_abbreviation'] == sample_name_abbreviation_I:
                    data_O.append(d);
            elif d['sample_name_abbreviation'] == sample_name_abbreviation_I:
                data_O.append(d);
        return data_O;

    def export_genesFpkmTracking_js(self,data_dir_I='tmp'):
        '''Export data for a box and whiskers plot'''

        data_O = [];
        data_O = self.genesFpkmTracking
        # dump chart parameters to a js files
        data1_keys = ['experiment_id','sample_name','gene_short_name'
                    ];
        data1_nestkeys = ['gene_short_name'];
        data1_keymap = {'xdata':'gene_short_name',
                        'ydata':'FPKM',
                        'ydatalb':'FPKM_conf_lo',
                        'ydataub':'FPKM_conf_hi',
                        'serieslabel':'sample_name',
                        'featureslabel':'gene_short_name'};
        # make the data object
        dataobject_O = [{"data":data_O,"datakeys":data1_keys,"datanestkeys":data1_nestkeys}];
        # make the tile parameter objects
        formtileparameters_O = {'tileheader':'Filter menu','tiletype':'html','tileid':"filtermenu1",'rowid':"row1",'colid':"col1",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-4"};
        formparameters_O = {'htmlid':'filtermenuform1',"htmltype":'form_01',"formsubmitbuttonidtext":{'id':'submit1','text':'submit'},"formresetbuttonidtext":{'id':'reset1','text':'reset'},"formupdatebuttonidtext":{'id':'update1','text':'update'}};
        formtileparameters_O.update(formparameters_O);
        svgparameters_O = {"svgtype":'boxandwhiskersplot2d_01',"svgkeymap":[data1_keymap],
                            'svgid':'svg1',
                            "svgmargin":{ 'top': 50, 'right': 150, 'bottom': 50, 'left': 50 },
                            "svgwidth":500,"svgheight":350,
                            "svgx1axislabel":"gene","svgy1axislabel":"FPKM",
    						'svgformtileid':'filtermenu1','svgresetbuttonid':'reset1','svgsubmitbuttonid':'submit1'};
        svgtileparameters_O = {'tileheader':'Custom box and whiskers plot','tiletype':'svg','tileid':"tile2",'rowid':"row1",'colid':"col2",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-8"};
        svgtileparameters_O.update(svgparameters_O);
        tableparameters_O = {"tabletype":'responsivetable_01',
                    'tableid':'table1',
                    "tablefilters":None,
                    "tableclass":"table  table-condensed table-hover",
    			    'tableformtileid':'filtermenu1','tableresetbuttonid':'reset1','tablesubmitbuttonid':'submit1'};
        tabletileparameters_O = {'tileheader':'FPKM','tiletype':'table','tileid':"tile3",'rowid':"row2",'colid':"col1",
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
            