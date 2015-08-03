from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData
from .general_feature_format import general_feature_format
from calculate_utilities.base_calculate import base_calculate
from Bio import SeqIO
from Bio import Entrez

class gff_coverage(general_feature_format):
    def __init__(self,gff_file_I,plus_I,minus_I,plus_high_regions_I, minus_high_regions_I,
                 coverage_I,coverageStats_I,
                 amplifications_I,amplificationStats_I,amplificationAnnotations):
        if gff_file_I:
            self.gff_file = gff_file_I;
        else:
            self.gff_file = None;
        if plus_I:
            self.plus = plus_I;
        else:
            self.plus = None;
        if minus_I:
            self.minus = minus_I;
        else:
            self.minus = None;
        if plus_high_regions_I:
            self.plus_high_regions = plus_high_regions_I;
        else:
            self.plus_high_regions = None;
        if minus_high_regions_I:
            self.minus_high_regions = minus_high_regions_I;
        else:
            self.minus_high_regions = None;
        
        if coverage_I:
            self.coverage = coverage_I;
        else:
            self.coverage = [];
        if coverageStats_I:
            self.coverageStats = coverageStats_I;
        else:
            self.coverageStats = [];
        
        if amplifications_I:
            self.amplifications = amplifications_I;
        else:
            self.amplifications = [];
        if amplificationStats_I:
            self.amplificationStats = amplificationStats_I;
        else:
            self.amplificationStats = [];

    def find_amplifications_fromGff(self,gff_file,
                strand_start, strand_stop,
                experiment_id_I = None,
                sample_name_I = None,
                scale_factor=True, downsample_factor=0,
                reads_min=1.5,reads_max=5.0,
                indices_min=200,consecutive_tol=10):
        """find amplifications from the gff file
        INPUT:
        strand_start = index of the start position
        strand_stop = index of the stop position
        scale_factor = boolean, if true, reads will be normalized to have 100 max
        downsample_factor = integer, factor to downsample the points to
        reads_min = minimum number of reads to identify an amplification
        reads_max = maximum number of reads to identify an amplification
        indices_min : minimum number of points of a high coverage region
        consecutive_tol: maximum number of consecutive points that do not meet the coverage_min/max criteria that can be included a high coverage region

        OPTION INPUT:
        experiment_id_I = tag for the experiment from which the sample came
        sample_name_I = tag for the sample name
        
        """
        data_O=[];
        experiment_id = experiment_id_I;
        sn = sample_name_I;
        # get the data_dir
        self.set_gffFile(gff_file);
        # extract the strands
        self.extract_strandsFromGff(strand_start, strand_stop, scale=scale_factor, downsample=downsample_factor)
        # find high coverage regions
        plus_high_region_indices,minus_high_region_indices = self.find_highCoverageRegions(coverage_min=reads_min,coverage_max=reads_max,points_min=indices_min,consecutive_tol=consecutive_tol);
        # record high coverage regions
        # + strand
        iter = 0;
        for index,reads in self.plus_high_regions.iteritems():
            if index > plus_high_region_indices[iter]['stop']:
                iter+=1;
            data_O.append({
            #'analysis_id':analysis_id,
            'experiment_id':experiment_id,
            'sample_name':sn,
            'genome_chromosome':1, #default
            'genome_strand':'+',
            'genome_index':int(index),
            'strand_start':strand_start,
            'strand_stop':strand_stop,
            'reads':float(reads),
            'reads_min':reads_min,
            'reads_max':reads_max,
            'indices_min':indices_min,
            'consecutive_tol':consecutive_tol,
            'scale_factor':scale_factor,
            'downsample_factor':downsample_factor,
            'amplification_start':int(plus_high_region_indices[iter]['start']),
            'amplification_stop':int(plus_high_region_indices[iter]['stop']),
            'used_':True,
            'comment_':None
                });
        # - strand
        iter = 0;
        for index,reads in self.minus_high_regions.iteritems():
            if index > minus_high_region_indices[iter]['stop']:
                iter+=1;
            data_O.append({
            #'analysis_id':analysis_id,
            'experiment_id':experiment_id,
            'sample_name':sn,
            'genome_chromosome':1, #default
            'genome_strand':'-',
            'genome_index':int(index),
            'strand_start':strand_start,
            'strand_stop':strand_stop,
            'reads':float(reads),
            'reads_min':reads_min,
            'reads_max':reads_max,
            'indices_min':indices_min,
            'consecutive_tol':consecutive_tol,
            'scale_factor':scale_factor,
            'downsample_factor':downsample_factor,
            'amplification_start':int(minus_high_region_indices[iter]['start']),
            'amplification_stop':int(minus_high_region_indices[iter]['stop']),
            'used_':True,
            'comment_':None
                });
        self.amplifications = data_O;

    def _get_strandsByStrand(self,data_I,strand_I):
        """return all data for the given strand"""
        data_O = [];
        for d in data_I:
            if data_I['genome_strand'] == strand_I:
                data_O.append(d);
        return data_O;

    def _get_strandsByChromosomeAndStrand(self,data_I,chromosome_I,strand_I):
        """return all data for the given chromosome and strand"""
        data_O = [];
        for d in data_I:
            if data_I['genome_chromosome'] == chromosome_I and data_I['genome_strand'] == strand_I:
                data_O.append(d);
        return data_O;

    def import_coverage(self,filename):
        '''import coverage from csv'''
        io = base_importData();
        io.read_csv(filename);
        self.coverage=io.data;

    def import_coverageStats(self,filename):
        '''import coverageStats from csv'''
        io = base_importData();
        io.read_csv(filename);
        self.coverageStats=io.data;

    def import_amplifications(self,filename):
        '''import amplifications from csv'''
        io = base_importData();
        io.read_csv(filename);
        self.amplifications=io.data;

    def import_amplificationStats(self,filename):
        '''import amplificationStats from csv'''
        io = base_importData();
        io.read_csv(filename);
        self.amplificationStats=io.data;

    def import_amplificationAnnotations(self,filename):
        '''import amplificationAnnotations from csv'''
        io = base_importData();
        io.read_csv(filename);
        self.amplificationAnnotations=io.data;

    def export_coverage(self,filename_O):
        """export coverage"""
        io = base_exportData(self.coverage);
        io.write_dict2csv(filename_O);

    def export_coverageStats(self,filename_O):
        """export coverageStats"""
        io = base_exportData(self.coverageStats);
        io.write_dict2csv(filename_O);

    def export_amplifications(self,filename_O):
        """export amplifications"""
        io = base_exportData(self.amplifications);
        io.write_dict2csv(filename_O);

    def export_amplificationStats(self,filename_O):
        """export amplificationStats"""
        io = base_exportData(self.amplificationStats);
        io.write_dict2csv(filename_O);

    def export_amplificationAnnotations(self,filename_O):
        """export amplificationAnnotations"""
        io = base_exportData(self.amplificationAnnotations);
        io.write_dict2csv(filename_O);

    def clear_data(self):
        self.gff_file = None;
        self.minus = None;
        self.plus = None;
        self.plus_high_regions = None;
        self.minus_high_regions = None;
        del self.coverage[:];
        del self.coverageStats[:];
        del self.amplifications[:];
        del self.amplificationStats[:];
        del self.amplificationAnnotations[:];

    def findAndCalculate_amplificationStats_fromGff(self,gff_file,
                strand_start, strand_stop,
                experiment_id_I = None,
                sample_name_I = None,
                scale_factor=True, downsample_factor=0,
                reads_min=1.5,reads_max=5.0,
                indices_min=200,consecutive_tol=10):
        """find amplifications from the gff file and calculate their statistics

        INPUT:
        strand_start = index of the start position
        strand_stop = index of the stop position
        scale_factor = boolean, if true, reads will be normalized to have 100 max
        downsample_factor = integer, factor to downsample the points to
        reads_min = minimum number of reads to identify an amplification
        reads_max = maximum number of reads to identify an amplification
        indices_min : minimum number of points of a high coverage region
        consecutive_tol: maximum number of consecutive points that do not meet the coverage_min/max criteria that can be included a high coverage region

        OPTION INPUT:
        experiment_id_I = tag for the experiment from which the sample came
        sample_name_I = tag for the sample name
        
        """
        data_O=[];
        stats_O=[];
        experiment_id = experiment_id_I;
        sn = sample_name_I;
        calculate = base_calculate();
        # get the data_dir
        self.set_gffFile(gff_file);
        # extract the strands
        self.extract_strandsFromGff(strand_start, strand_stop, scale=scale_factor, downsample=downsample_factor)
        # find high coverage regions
        plus_high_region_indices,minus_high_region_indices = self.find_highCoverageRegions(coverage_min=reads_min,coverage_max=reads_max,points_min=indices_min,consecutive_tol=consecutive_tol);
        
        # record the means for later use
        plus_mean,minus_mean = plus.mean(),minus.mean();
        plus_min,minus_min = plus.min(),minus.min();
        plus_max,minus_max = plus.max(),minus.max();
        # find high coverage regions
        plus_high_region_indices,minus_high_region_indices,plus_high_regions, minus_high_regions = find_highCoverageRegions(plus,minus,coverage_min=reads_min,coverage_max=reads_max,points_min=indices_min,consecutive_tol=consecutive_tol)
        # calculate stats on the high coverage regions
        # + strand
        for row_cnt,row in enumerate(plus_high_region_indices):
            plus_region = plus_high_regions[(plus_high_regions.index>=row['start']) & (plus_high_regions.index<=row['stop'])]
            # calculate using scipy
            data_ave_O, data_var_O, data_lb_O, data_ub_O = None, None, None, None;
            data_ave_O, data_var_O, data_lb_O, data_ub_O = calculate.calculate_ave_var(plus_region.values,confidence_I = 0.95);
            # calculate the interquartile range
            min_O, max_O, median_O, iq_1_O, iq_3_O = None, None, None, None, None;
            min_O, max_O, median_O, iq_1_O, iq_3_O=calculate.calculate_interquartiles(plus_region.values);
            # record data
            stats_O.append({
                #'analysis_id':analysis_id,
                'experiment_id':experiment_id,
                'sample_name':sn,
                'genome_chromosome':1,
                'genome_strand':'plus',
                'strand_start':strand_start,
                'strand_stop':strand_stop,
                'reads_min':min_O,
                'reads_max':max_O,
                'reads_lb':data_lb_O,
                'reads_ub':data_ub_O,
                'reads_iq1':iq_1_O,
                'reads_iq3':iq_3_O,
                'reads_median':median_O,
                'reads_mean':data_ave_O,
                'reads_var':data_var_O,
                'reads_n':len(plus_region.values),
                'amplification_start':int(row['start']),
                'amplification_stop':int(row['stop']),
                'used_':True,
                'comment_':None
                })
            # downsample
            collapse_factor = None;
            if downsample_factor > 1:
                collapse_factor = int((row['stop'] - row['start']) / downsample_factor)
            if collapse_factor and collapse_factor > 1:
                plus_region = plus_region.groupby(lambda x: x // collapse_factor).mean()
                plus_region.index *= collapse_factor
            # add mean to index before and after the amplification start and stop, respectively (for visualization)
            if downsample_factor > 1 and row_cnt==0:
                #plus_region[strand_start]=plus_mean;
                #plus_region[strand_stop]=plus_mean;
                data_O.append({
                    #'analysis_id':analysis_id,
                    'experiment_id':experiment_id,
                    'sample_name':sn,
                    'genome_chromosome':1, #default
                    'genome_strand':'plus_mean',
                    #'genome_index':int(strand_start),
                    'genome_index':int(row['start']-1),
                    'strand_start':strand_start,
                    'strand_stop':strand_stop,
                    'reads':plus_mean,
                    'reads_min':reads_min,
                    'reads_max':reads_max,
                    'indices_min':indices_min,
                    'consecutive_tol':consecutive_tol,
                    'scale_factor':scale_factor,
                    'downsample_factor':downsample_factor,
                    'amplification_start':strand_start,
                    'amplification_stop':strand_stop,
                    'used_':True,
                    'comment_':'mean reads of the plus strand'
                    });
            if downsample_factor > 1 and row_cnt==len(plus_high_region_indices)-1:
                data_O.append({
                    #'analysis_id':analysis_id,
                    'experiment_id':experiment_id,
                    'sample_name':sn,
                    'genome_chromosome':1, #default
                    'genome_strand':'plus_mean',
                    #'genome_index':int(strand_stop),
                    'genome_index':int(row['stop']+1),
                    'strand_start':strand_start,
                    'strand_stop':strand_stop,
                    'reads':plus_mean,
                    'reads_min':reads_min,
                    'reads_max':reads_max,
                    'indices_min':indices_min,
                    'consecutive_tol':consecutive_tol,
                    'scale_factor':scale_factor,
                    'downsample_factor':downsample_factor,
                    'amplification_start':strand_start,
                    'amplification_stop':strand_stop,
                    'used_':True,
                    'comment_':'mean reads of the plus strand'
                    });
            ## add zeros to strand start and stop, respectively (for visualization)
            #if downsample_factor > 1:
            #    plus_region[row['start']-1]=plus_mean;
            #    plus_region[row['stop']+1]=plus_mean;
            # record high coverage regions
            for index,reads in plus_region.iteritems():
                data_O.append({
                    #'analysis_id':analysis_id,
                    'experiment_id':experiment_id,
                    'sample_name':sn,
                    'genome_chromosome':1, #default
                    'genome_strand':'plus',
                    'genome_index':int(index),
                    'strand_start':strand_start,
                    'strand_stop':strand_stop,
                    'reads':float(reads),
                    'reads_min':reads_min,
                    'reads_max':reads_max,
                    'indices_min':indices_min,
                    'consecutive_tol':consecutive_tol,
                    'scale_factor':scale_factor,
                    'downsample_factor':downsample_factor,
                    'amplification_start':int(row['start']),
                    'amplification_stop':int(row['stop']),
                    'used_':True,
                    'comment_':None
                });
        # - strand
        for row_cnt,row in enumerate(minus_high_region_indices):
            minus_region = minus_high_regions[(minus_high_regions.index>=row['start']) & (minus_high_regions.index<=row['stop'])]
            # calculate using scipy
            data_ave_O, data_var_O, data_lb_O, data_ub_O = None, None, None, None;
            data_ave_O, data_var_O, data_lb_O, data_ub_O = calculate.calculate_ave_var(minus_region.values,confidence_I = 0.95);
            # calculate the interquartile range
            min_O, max_O, median_O, iq_1_O, iq_3_O = None, None, None, None, None;
            min_O, max_O, median_O, iq_1_O, iq_3_O=calculate.calculate_interquartiles(minus_region.values);
            # record data
            stats_O.append({
                #'analysis_id':analysis_id,
                'experiment_id':experiment_id,
                'sample_name':sn,
                'genome_chromosome':1,
                'genome_strand':'minus',
                'strand_start':strand_start,
                'strand_stop':strand_stop,
                'reads_min':min_O,
                'reads_max':max_O,
                'reads_lb':data_lb_O,
                'reads_ub':data_ub_O,
                'reads_iq1':iq_1_O,
                'reads_iq3':iq_3_O,
                'reads_median':median_O,
                'reads_mean':data_ave_O,
                'reads_var':data_var_O,
                'reads_n':len(minus_region.values),
                'amplification_start':int(row['start']),
                'amplification_stop':int(row['stop']),
                'used_':True,
                'comment_':None
                })
            # downsample
            collapse_factor = None;
            if downsample_factor > 1:
                collapse_factor = int((row['stop'] - row['start']) / downsample_factor)
            if collapse_factor and collapse_factor > 1:
                minus_region = minus_region.groupby(lambda x: x // collapse_factor).mean()
                minus_region.index *= collapse_factor
            # add mean to index before and after the amplification start and stop, respectively (for visualization)
            if downsample_factor > 1 and row_cnt==0:
                #minus_region[strand_start]=minus_mean;
                #minus_region[strand_stop]=minus_mean;
                data_O.append({
                    #'analysis_id':analysis_id,
                    'experiment_id':experiment_id,
                    'sample_name':sn,
                    'genome_chromosome':1, #default
                    'genome_strand':'minus_mean',
                    #'genome_index':int(strand_start),
                    'genome_index':int(row['start']-1),
                    'strand_start':strand_start,
                    'strand_stop':strand_stop,
                    'reads':minus_mean,
                    'reads_min':reads_min,
                    'reads_max':reads_max,
                    'indices_min':indices_min,
                    'consecutive_tol':consecutive_tol,
                    'scale_factor':scale_factor,
                    'downsample_factor':downsample_factor,
                    'amplification_start':strand_start,
                    'amplification_stop':strand_stop,
                    'used_':True,
                    'comment_':'mean reads of the minus strand'
                    });
            if downsample_factor > 1 and row_cnt==len(minus_high_region_indices)-1:
                data_O.append({
                    #'analysis_id':analysis_id,
                    'experiment_id':experiment_id,
                    'sample_name':sn,
                    'genome_chromosome':1, #default
                    'genome_strand':'minus_mean',
                    #'genome_index':int(strand_stop),
                    'genome_index':int(row['stop']+1),
                    'strand_start':strand_start,
                    'strand_stop':strand_stop,
                    'reads':minus_mean,
                    'reads_min':reads_min,
                    'reads_max':reads_max,
                    'indices_min':indices_min,
                    'consecutive_tol':consecutive_tol,
                    'scale_factor':scale_factor,
                    'downsample_factor':downsample_factor,
                    'amplification_start':strand_start,
                    'amplification_stop':strand_stop,
                    'used_':True,
                    'comment_':'mean reads of the minus strand'
                    });
            ## add zeros to strand start and stop, respectively (for visualization)
            #if downsample_factor > 1:
            #    minus_region[row['start']-1]=minus_mean;
            #    minus_region[row['stop']+1]=minus_mean;
            # record high coverage regions
            for index,reads in minus_region.iteritems():
                data_O.append({
                    #'analysis_id':analysis_id,
                    'experiment_id':experiment_id,
                    'sample_name':sn,
                    'genome_chromosome':1, #default
                    'genome_strand':'minus',
                    'genome_index':int(index),
                    'strand_start':strand_start,
                    'strand_stop':strand_stop,
                    'reads':float(reads),
                    'reads_min':reads_min,
                    'reads_max':reads_max,
                    'indices_min':indices_min,
                    'consecutive_tol':consecutive_tol,
                    'scale_factor':scale_factor,
                    'downsample_factor':downsample_factor,
                    'amplification_start':int(row['start']),
                    'amplification_stop':int(row['stop']),
                    'used_':True,
                    'comment_':None});
        #record the data
        self.amplifications = data_O;
        self.amplificationStats = stats_O;

    def annotate_amplifications(self,ref_genome_I='U00096.2.gb',
                           ref_I = 'genbank',geneReference_I=None,biologicalmaterial_id_I='MG1655'):
        """annotate amplificaitons from reference
        ref_genome_I = reference genome to use for the annotation
        ref_I = reference database
        geneReference_I = filename for the gene reference table
        biologicalmaterial_id_I = biologicalmatereial_id for the geneReference to use for the annotation (required to generate an ecocyc link)
        """
        genomeannotation = genome_annotation(ref_genome_I,ref_I,geneReference_I);
        data_O = [];


