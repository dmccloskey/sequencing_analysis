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