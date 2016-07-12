from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData
from .genes_fpkm_tracking import genes_fpkm_tracking
import copy

class genes_countFPKMattr_table(genes_fpkm_tracking):
    '''Helper class to parse the output from cuffnorm
    http://cole-trapnell-lab.github.io/cufflinks/cuffnorm/#cuffnorm-output
    '''

    def __init__(self,countTable_I=None,fpkmTable_I=None,attrTable_I=None):
        if countTable_I:
            self.countTable = countTable_I;
        else:
            self.countTable = [];
        if fpkmTable_I:
            self.fpkmTable = fpkmTable_I;
        else:
            self.fpkmTable = [];
        if attrTable_I:
            self.attrTable = attrTable_I;
        else:
            self.attrTable = [];

    def import_countTable(self,filename_I):
        """import count table
        INPUT:
        filename_I = input filename
        """
        countTable = self.import_countOrFPKMTable(filename_I);
        self.countTable = countTable;

    def import_fpkmTable(self,filename_I):
        """import FPKM table
        INPUT:
        filename_I = input filename
        """
        fpkmTable = self.import_countOrFPKMTable(filename_I);
        self.fpkmTable = fpkmTable;

    def import_attrTable(self,filename_I):
        """import attr table
        INPUT:
        filename_I = input filename
        """
        #import and format the data
        io = base_importData();
        io.read_tab(filename_I);
        attrTable = self.format_genesFpkmTracking(io.data);
        self.attrTable = attrTable;
        
    def import_countOrFPKMTable(
        self,filename_I):
        """import count or FPKM table
        INPUT:
        filename_I = input filename
        """
        #import and format the data
        io = base_importData();
        io.read_tab(filename_I);
        countOrFPKMTable = self.format_countOrFPKMTable(io.data);
        return countOrFPKMTable;

    def reformat_attrTable(
        self):
        """reformat attr tables into a dictionary
        for rapid alignment of attr table with tracking_id

        """
        #format into a dictionary of rows for quick aligning with the tracking_id
        if self.attrTable: attrTable = self.attrTable[:];
        else: attrTable = [];

        attrTable_dict = {};
        for row in attrTable:
            attrTable_dict[row['tracking_id']] = row;
        return attrTable_dict;

    def reformat_countTable(
        self,analysis_id_I=None,sna2experimentID_I=None,
        sna2sns_I=None):
        """reformat count table into a flattened table of sample_names/values

        OPTIONAL INPUT:
        the following are optional for analyzing a single sample,
        but required when analyzing multiple samples

        analysis_id_I = string, name of the grouping of samples
        sna2experimentID_I = dict, mapping of cuffnorm sample label to desired experiment_id
        sna2sns_I = dict, mapping of cuffnorm sample label to desired sample name
        """
        if self.countTable: countTable = self.countTable[:];
        else: countTable = [];

        countTable_flat = self.reformat_countOrFPKMTable(
            countOrFPKMTable_I=countTable,
            analysis_id_I=analysis_id_I,
            sna2experimentID_I=sna2experimentID_I,
            sna2sns_I=sna2sns_I,
            count_or_FPKM = 'count');
        return countTable_flat;

    def reformat_fpkmTable(
        self,analysis_id_I=None,sna2experimentID_I=None,
        sna2sns_I=None):
        """reformat fpkm table into flattened table of sample_names/values

        OPTIONAL INPUT:
        the following are optional for analyzing a single sample,
        but required when analyzing multiple samples

        analysis_id_I = string, name of the grouping of samples
        sna2experimentID_I = dict, mapping of cuffnorm sample label to desired experiment_id
        sna2sns_I = dict, mapping of cuffnorm sample label to desired sample name
        """
        if self.fpkmTable: fpkmTable = self.fpkmTable[:];
        else: fpkmTable = [];

        fpkmTable_flat = self.reformat_countOrFPKMTable(
            countOrFPKMTable_I=fpkmTable,
            analysis_id_I=analysis_id_I,
            sna2experimentID_I=sna2experimentID_I,
            sna2sns_I=sna2sns_I,
            count_or_FPKM = 'fpkm');
        return fpkmTable_flat;

    def reformat_countOrFPKMTable(
        self,
        countOrFPKMTable_I=None,
        analysis_id_I=None,
        sna2experimentID_I=None,
        sna2sns_I=None,
        count_or_FPKM = 'count'):
        """reformat count or FPKM tables into flattened table of sample_names/values
        for rapid alignment of attr table with tracking_id

        OPTIONAL INPUT:
        the following are optional for analyzing a single sample,
        but required when analyzing multiple samples

        analysis_id_I = string, name of the grouping of samples
        sna2experimentID_I = dict, mapping of cuffnorm sample label to desired experiment_id
        sna2sns_I = dict, mapping of cuffnorm sample label to desired sample name
        count_or_FPKM = string, normalized 'count' or raw 'FPKM'
        """
        #format into a dictionary of rows for quick aligning with the tracking_id
        countOrFPKMTable_flat = [];
        for row in countOrFPKMTable_I:
            for k,v in row.items():
                if k=='tracking_id':continue;
                tmp = {};
                tmp['analysis_id'] = analysis_id_I;
                tmp['tracking_id'] = row['tracking_id'];

                sample_name_lst = k.split('_');
                sample_name_base = '_'.join(sample_name_lst[:-1]);
                sample_name_rep = eval(sample_name_lst[-1]);
                if sna2experimentID_I: 
                    experiment_id = sna2experimentID_I[sample_name_base];
                else:
                    experiment_id=None;
                tmp['experiment_id'] = experiment_id;
                if sna2sns_I: 
                    sample_name = sna2sns_I[sample_name_base][sample_name_rep];
                else:
                    sample_name=k;
                tmp['sample_name'] = sample_name;

                tmp['value'] = v;
                tmp['value_units'] = count_or_FPKM;
                tmp['used_'] = True;
                tmp['comment_'] = None;
                countOrFPKMTable_flat.append(tmp);
        return countOrFPKMTable_flat;

    def alignAndReformat_countFPKMattrTables(
        self,analysis_id_I=None,sna2experimentID_I=None,
        sna2sns_I=None):
        """reformat count or FPKM tables into a dictionary
        for rapid alignment of attr table with tracking_id

        OPTIONAL INPUT:
        the following are optional for analyzing a single sample,
        but required when analyzing multiple samples

        analysis_id_I = string, name of the grouping of samples
        sna2experimentID_I = dict, mapping of cuffnorm sample label to desired experiment_id
        sna2sns_I = dict, mapping of cuffnorm sample label to desired sample name
        """
        #reformat
        countTable_flat = self.reformat_countTable(
            analysis_id_I=analysis_id_I,
            sna2experimentID_I=sna2experimentID_I,
            sna2sns_I=sna2sns_I,);
        fpkmTable_flat = self.reformat_fpkmTable(
            analysis_id_I=analysis_id_I,
            sna2experimentID_I=sna2experimentID_I,
            sna2sns_I=sna2sns_I,);
        attrTable_dict = self.reformat_attrTable();
        #align
        countAndFpkmTable_aligned = [];
        for row in countTable_flat[:]:
            row.update(attrTable_dict[row['tracking_id']]);
            countAndFpkmTable_aligned.append(row);
        for row in fpkmTable_flat[:]:
            row.update(attrTable_dict[row['tracking_id']]);
            countAndFpkmTable_aligned.append(row);
        return countAndFpkmTable_aligned;
    
    def format_countOrFPKMTable(self,fpkmTracking_I):
        """formats raw string input into their appropriate values"""
        for fpkmTracking in fpkmTracking_I:
            for k,v in fpkmTracking.items():
                if k=='tracking_id' and type(fpkmTracking['tracking_id'])==type('string'):
                    pass;
                elif k!='tracking_id' and type(fpkmTracking[k])==type('string'):
                    fpkmTracking[k] = eval(v);
        return fpkmTracking_I;
        