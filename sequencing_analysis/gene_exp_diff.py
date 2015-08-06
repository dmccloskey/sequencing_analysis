from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData

class gene_exp_diff():
    def __init__(self,geneExpDiff_I=None):
        if geneExpDiff_I:
            self.geneExpDiff = geneExpDiff_I;
        else:
            self.geneExpDiff = [];

    def import_geneExpDiff(self,filename_I,experiment_id_1_I=None,experiment_id_2_I=None,sample_name_abbreviation_1_I=None,sample_abbreviation_2_I=None):
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
        io.read_csv(filename_I);
        geneExpDiff = self.format_geneExpDiff(io.data);
        for d in geneExpDiff:
            d['experiment_id_1'] = experiment_id_1_I;
            d['experiment_id'] = experiment_id_2_I;
            d['sample_name_abbreviation_1'] = sample_name_abbreviation_1_I;
            d['sample_name_abbreviation_2'] = sample_name_abbreviation_2_I;
        self.geneExpDiff = geneExpDiff;

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