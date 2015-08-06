from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData
from .genes_fpkm_tracking import genes_fpkm_tracking
from calculate_utilities.base_calculate import base_calculate

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
        filename_I = list of input filename
        
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
        calculate = base_calculate();
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
            