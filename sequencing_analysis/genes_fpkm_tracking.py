from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData

class genes_fpkm_tracking():
    def __init__(self,genesFpkmTracking_I=None):
        if genesFpkmTracking_I:
            self.genesFpkmTracking = genesFpkmTracking_I;
        else:
            self.genesFpkmTracking = [];

    def import_genesFpkmTracking(self,filename_I,experiment_id_I=None,sample_name_I=None):
        """import geneExpDiff
        INPUT:
        filename_I = input filename
        
        OPTIONAL INPUT:
        the following are optional for analyzing a single sample,
        but required when analyzing multiple samples

        experiment_id_I = string, name of the experiment that generated the sample
        sample_name_I = string, name of the sample
        """
        io = base_importData();
        io.read_csv(filename_I);
        genesFpkmTracking = self.format_genesFpkmTracking(io.data);
        for d in genesFpkmTracking:
            d['experiment_id'] = experiment_id_I;
            d['sample_name'] = sample_name_I;
        self.genesFpkmTracking = genesFpkmTracking;

    def export_genesFpkmTracking(self,filename_O):
        """export genesFpkmTracking"""
        io = base_exportData(self.genesFpkmTracking);
        io.write_dict2csv(filename_O);
    
    def format_genesFpkmTracking(self,fpkmTracking_I):
        """formats raw string input into their appropriate values"""
        for fpkmTracking in fpkmTracking_I:
            if 'FPKM' in fpkmTracking and type(fpkmTracking['FPKM'])==type('string'):
                fpkmTracking['FPKM'] = eval(fpkmTracking['FPKM']);
            if 'FPKM_conf_lo' in fpkmTracking and type(fpkmTracking['FPKM_conf_lo'])==type('string'):
                fpkmTracking['FPKM_conf_lo'] = eval(fpkmTracking['FPKM_conf_lo']);
            if 'FPKM_conf_hi' in fpkmTracking and type(fpkmTracking['FPKM_conf_hi'])==type('string'):
                fpkmTracking['FPKM_conf_hi'] = eval(fpkmTracking['FPKM_conf_hi']);
        return fpkmTracking_I;

    def export_genesFpkmTracking_js(self,data_dir_I='tmp'):
        '''Export data for a box and whiskers plot'''

        data_O = [];
        data_O = self.genesFpkmTracking
        # dump chart parameters to a js files
        data1_keys = ['analysis_id','experiment_id','sample_name','gene'
                    ];
        data1_nestkeys = ['gene'];
        data1_keymap = {'xdata':'gene',
                        'ydata':'FPKM',
                        'ydatalb':'FPKM_conf_lo',
                        'ydataub':'FPKM_conf_hi',
                        'serieslabel':'sample_name',
                        'featureslabel':'gene'};
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