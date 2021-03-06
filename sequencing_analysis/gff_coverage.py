from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData
from .general_feature_format import general_feature_format
from python_statistics.calculate_interface import calculate_interface
from .genome_annotations import genome_annotations
import json

class gff_coverage(general_feature_format):
    def __init__(self,gff_file_I = None,plus_I = None,minus_I = None,
                 plus_high_regions_I = None, minus_high_regions_I = None,
                 coverage_I = None,coverageStats_I = None,
                 amplifications_I = None,amplificationStats_I = None,
                 amplificationAnnotations_I = None):
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
        if amplificationAnnotations_I:
            self.amplificationAnnotations = amplificationAnnotations_I;
        else:
            self.amplificationAnnotations = [];

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
            if d['genome_strand'] == strand_I:
                data_O.append(d);
        return data_O;

    def _get_strandsByChromosomeAndStrand(self,data_I,chromosome_I,strand_I):
        """return all data for the given chromosome and strand"""
        data_O = [];
        for d in data_I:
            if d['genome_chromosome'] == chromosome_I and d['genome_strand'] == strand_I:
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
        calculate = calculate_interface();
        # get the data_dir
        self.set_gffFile(gff_file);
        # extract the strands
        self.extract_strandsFromGff(strand_start, strand_stop, scale=scale_factor, downsample=0)
        # find high coverage regions
        plus_high_region_indices,minus_high_region_indices = self.find_highCoverageRegions(coverage_min=reads_min,coverage_max=reads_max,points_min=indices_min,consecutive_tol=consecutive_tol);
        
        # record the means for later use
        plus_mean,minus_mean = self.plus.mean(),self.minus.mean();
        plus_min,minus_min = self.plus.min(),self.minus.min();
        plus_max,minus_max = self.plus.max(),self.minus.max();
        # calculate stats on the high coverage regions
        # + strand
        for row_cnt,row in enumerate(plus_high_region_indices):
            plus_region = self.plus_high_regions[(self.plus_high_regions.index>=row['start']) & (self.plus_high_regions.index<=row['stop'])]
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
            minus_region = self.minus_high_regions[(self.minus_high_regions.index>=row['start']) & (self.minus_high_regions.index<=row['stop'])]
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
        genomeannotation = genome_annotations(ref_genome_I,ref_I,geneReference_I);
        
        data_O = [];
        # get amplification regions
        amplificationStats = self.amplificationStats;
        # annotate each region
        for row in amplificationStats:
            # annotate each mutation based on the position
            annotations = [];
            annotations = genomeannotation._find_genesInRegion(row['amplification_start'],row['amplification_start'])
            for annotation in annotations:
                # record the data
                tmp = {
                    'experiment_id':row['experiment_id'],
                    'sample_name':row['sample_name'],
                    'genome_chromosome':row['genome_chromosome'],
                    'genome_strand':row['genome_strand'],
                    'strand_start':row['strand_start'],
                    'strand_stop':row['strand_stop'],
                    'amplification_start':row['amplification_start'],
                    'amplification_stop':row['amplification_stop'],
                    'used_':True,
                    'comment_':None};
                tmp['feature_genes'] = annotation['gene']
                tmp['feature_locations'] = annotation['location']
                tmp['feature_annotations'] = annotation['product']
                tmp['feature_start'] = annotation['start'];
                tmp['feature_stop'] = annotation['stop'];
                tmp['feature_types'] = annotation['type']
                # generate a link to ecogene for the genes
                tmp['feature_links'] = [];
                for bnumber in annotation['locus_tag']:
                    if bnumber:
                        ecogenes = [];
                        ecogenes = genomeannotation._get_ecogenesByBiologicalmaterialIDAndOrderedLocusName(biologicalmaterial_id_I,bnumber);
                        if ecogenes:
                            ecogene = ecogenes[0];
                            ecogene_link = genomeannotation._generate_httplink2gene_ecogene(ecogene['ecogene_accession_number']);
                            tmp['feature_links'].append(ecogene_link)
                        else: print('no ecogene_accession_number found for ordered_locus_location ' + bnumber);
                data_O.append(tmp);
        self.amplificationAnnotations = data_O;

    def _get_chromosomes(self,data_I,experiment_id_I=None,sample_name_I=None):
        """return all chromosomes"""
        data_O = [];
        for d in data_I:
            if experiment_id_I and sample_name_I:
                if d['experiment_id'] == experiment_id_I and d['sample_name'] == sample_name_I:
                    data_O.append(d['genome_chromosome']);
            else:
                data_O.append(d['genome_chromosome']);
        return data_O;

    def _get_strandsByChromosome(self,data_I,chromosome_I,experiment_id_I=None,sample_name_I=None):
        """return strands for the given chromosome"""
        data_O = [];
        for d in data_I:
            if experiment_id_I and sample_name_I:
                if d['experiment_id'] == experiment_id_I and d['sample_name'] == sample_name_I and d['genome_chromosome'] == chromosome_I:
                    data_O.append(d['genome_strand']);
            elif d['genome_chromosome'] == chromosome_I:
                data_O.append(d['genome_strand']);
        return data_O;

    def _get_startAndStopsByChromosomeAndStrand(self,data_I,chromosome_I,strand_I,experiment_id_I=None,sample_name_I=None):
        """return strand start and stop positions for the given chromosome and strand"""
        genomic_starts,genomic_stops = [],[]
        for d in data_I:
            if experiment_id_I and sample_name_I:
                if d['experiment_id'] == experiment_id_I and d['sample_name'] == sample_name_I and d['genome_chromosome'] == chromosome_I and d['genome_strand'] == strand_I:
                    genomic_starts.append(d['strand_start']);
                    genomic_stops.append(d['strand_stop']);
            elif d['genome_chromosome'] == chromosome_I and d['genome_strand'] == strand_I:
                genomic_starts.append(d['strand_start']);
                genomic_stops.append(d['strand_stop']);
        return genomic_starts,genomic_stops;

    def _get_amplificationRegionsByChromosomeAndStrand(self,data_I,chromosome_I,strand_I,experiment_id_I=None,sample_name_I=None):
        """return strand start and stop positions for the given chromosome and strand"""
        genomic_starts,genomic_stops = [],[]
        for d in data_I:
            if experiment_id_I and sample_name_I:
                if d['experiment_id'] == experiment_id_I and d['sample_name'] == sample_name_I and d['genome_chromosome'] == chromosome_I and d['genome_strand'] == strand_I:
                    genomic_starts.append(d['amplification_start']);
                    genomic_stops.append(d['amplification_stop']);
            elif d['genome_chromosome'] == chromosome_I and d['genome_strand'] == strand_I:
                genomic_starts.append(d['amplification_start']);
                genomic_stops.append(d['amplification_stop']);
        return genomic_starts,genomic_stops;

    def _get_amplificationRegions(self,data_I,experiment_id_I=None,sample_name_I=None):
        """return strand start and stop positions"""
        genomic_starts,genomic_stops = [],[]
        for d in data_I:
            if experiment_id_I and sample_name_I:
                if d['experiment_id'] == experiment_id_I and d['sample_name'] == sample_name_I:
                    genomic_starts.append(d['amplification_start']);
                    genomic_stops.append(d['amplification_stop']);
            else:
                genomic_starts.append(d['amplification_start']);
                genomic_stops.append(d['amplification_stop']);
        return genomic_starts,genomic_stops;

    def extract_coverage_fromGff(self,gff_file, 
         strand_start,strand_stop,scale_factor=True,downsample_factor=2000,
         experiment_id_I=None, sample_name_I=None):
        """extract coverage (genome position and reads) from .gff
        INPUT:
        strand_start = index of the start position
        strand_stop = index of the stop position
        scale_factor = boolean, if true, reads will be normalized to have 100 max
        downsample_factor = integer, factor to downsample the points to
     
        OPTION INPUT:
        experiment_id_I = tag for the experiment from which the sample came
        sample_name_I = tag for the sample name
        
        """
        self.set_gffFile(gff_file);
        filename = self.gff_file;
        experiment_id = experiment_id_I;
        sample_name = sample_name_I;
        # parse the gff file into pandas dataframes
        self.extract_strandsFromGff(strand_start, strand_stop, scale=scale_factor, downsample=downsample_factor)
        # split into seperate data structures based on the destined table add
        coverage_data = [];
        if not self.plus.empty:
            for index,reads in self.plus.iteritems():
                coverage_data.append({
                                #'analysis_id':analysis_id,
                                'experiment_id':experiment_id,
                                'sample_name':sample_name,
                                'data_dir':filename,
                                'genome_chromosome':1, #default
                                'genome_strand':'plus',
                                'genome_index':int(index),
                                'strand_start':strand_start,
                                'strand_stop':strand_stop,
                                'reads':float(reads),
                                'scale_factor':scale_factor,
                                'downsample_factor':downsample_factor,
                                'used_':True,
                                'comment_':None});
        if not self.minus.empty:
            for index,reads in self.minus.iteritems():
                coverage_data.append({
                                #'analysis_id':analysis_id,
                                'experiment_id':experiment_id,
                                'sample_name':sample_name,
                                'data_dir':filename,
                                'genome_chromosome':1, #default
                                'genome_strand':'minus',
                                'genome_index':int(index),
                                'strand_start':strand_start,
                                'strand_stop':strand_stop,
                                'reads':float(reads),
                                'scale_factor':scale_factor,
                                'downsample_factor':downsample_factor,
                                'used_':True,
                                'comment_':None});
        # add data to the database:
        self.coverage = coverage_data;

    def calculate_coverageStats_fromGff(self,gff_file, 
         strand_start,strand_stop,scale_factor=True,downsample_factor=2000,
         experiment_id_I=None, sample_name_I=None):
        """extract coverage (genome position and reads) from .gff
        INPUT:
        strand_start = index of the start position
        strand_stop = index of the stop position
        scale_factor = boolean, if true, reads will be normalized to have 100 max
        downsample_factor = integer, factor to downsample the points to
     
        OPTION INPUT:
        experiment_id_I = tag for the experiment from which the sample came
        sample_name_I = tag for the sample name
        
        """
        calculate = calculate_interface();

        self.set_gffFile(gff_file);
        filename = self.gff_file;
        experiment_id = experiment_id_I;
        sn = sample_name_I;
        # parse the gff file into pandas dataframes
        self.extract_strandsFromGff(strand_start, strand_stop, scale=scale_factor, downsample=downsample_factor)
        # split into seperate data structures based on the destined table add
        coverageStats_data = [];
        # plus strand
        # calculate using scipy
        data_ave_O, data_var_O, data_lb_O, data_ub_O = None, None, None, None;
        data_ave_O, data_var_O, data_lb_O, data_ub_O = calculate.calculate_ave_var(self.plus.values,confidence_I = 0.95);
        # calculate the interquartile range
        min_O, max_O, median_O, iq_1_O, iq_3_O = None, None, None, None, None;
        min_O, max_O, median_O, iq_1_O, iq_3_O=calculate.calculate_interquartiles(self.plus.values);
        # record data
        coverageStats_data.append({
            #'analysis_id':analysis_id,
            'experiment_id':experiment_id,
            'sample_name':sn,
            'genome_chromosome':1,
            'genome_strand':'plus',
            'strand_start':strand_start,
            'strand_stop':strand_stop,
            'reads_min':int(min_O),
            'reads_max':int(max_O),
            'reads_lb':data_lb_O,
            'reads_ub':data_ub_O,
            'reads_iq1':iq_1_O,
            'reads_iq3':iq_3_O,
            'reads_median':median_O,
            'reads_mean':data_ave_O,
            'reads_var':data_var_O,
            'reads_n':len(self.plus.values),
            'used_':True,
            'comment_':None});
        # minus strand
        # calculate using scipy
        data_ave_O, data_var_O, data_lb_O, data_ub_O = None, None, None, None;
        data_ave_O, data_var_O, data_lb_O, data_ub_O = calculate.calculate_ave_var(self.minus.values,confidence_I = 0.95);
        # calculate the interquartile range
        min_O, max_O, median_O, iq_1_O, iq_3_O = None, None, None, None, None;
        min_O, max_O, median_O, iq_1_O, iq_3_O=calculate.calculate_interquartiles(self.minus.values);
        # record data
        coverageStats_data.append({
            #'analysis_id':analysis_id,
            'experiment_id':experiment_id,
            'sample_name':sn,
            'genome_chromosome':1,
            'genome_strand':'minus',
            'strand_start':strand_start,
            'strand_stop':strand_stop,
            'reads_min':int(min_O),
            'reads_max':int(max_O),
            'reads_lb':data_lb_O,
            'reads_ub':data_ub_O,
            'reads_iq1':iq_1_O,
            'reads_iq3':iq_3_O,
            'reads_median':median_O,
            'reads_mean':data_ave_O,
            'reads_var':data_var_O,
            'reads_n':len(self.minus.values),
            'used_':True,
            'comment_':None});
        # record the data
        self.coverageStats = coverageStats_data;
        
    def export_amplifications_js(self,data_dir_I="tmp"):
        """export amplifications and statistics to js file"""

        #get the data for the analysis
        data1_O = [];
        data2_O = [];
        data3_O = [];
        data1_O = self._make_sampleNameStrand(self.amplifications);
        data2_O = self._make_sampleNameStrand(self.amplificationStats);
        data3_O = self._make_sampleNameStrand(self.amplificationAnnotations);
        # dump chart parameters to a js files
        data1_keys = ['experiment_id',
                    'sample_name',
                    'genome_chromosome',
                    'genome_strand',
                    'amplification_start',
                    'amplification_stop',
                    'sample_name_strand',
                    ]
        data1_nestkeys = [
                        #'sample_name',
                        'genome_strand'
                        ];
        data1_keymap = {'xdata':'genome_index',
                        'ydata':'reads',
                        'serieslabel':'sample_name_strand',#custom for vis
                        #'serieslabel':'genome_strand',
                        'featureslabel':'reads'};
        data2_keys = ['experiment_id',
                'sample_name',
                'genome_chromosome',
                'genome_strand',
                #'reads_min',
                #'reads_max',
                #'reads_lb',
                #'reads_ub',
                #'reads_iq1',
                #'reads_iq3',
                #'reads_median',
                #'reads_mean',
                #'reads_var',
                #'reads_n',
                'amplification_start',
                'amplification_stop',
                    ]
        data2_nestkeys = ['sample_name'];
        data2_keymap = {'xdata':'genome_index',
                        'ydata':'reads',
                        'serieslabel':'genome_strand',
                        'featureslabel':'reads'};
        data3_keys = ['experiment_id',
                'sample_name',
                'genome_chromosome',
                'genome_strand',
                'feature_annotations',
                'feature_genes',
                'feature_locations',
                'feature_links',
                'feature_start',
                'feature_stop',
                'feature_types',
                'amplification_start',
                'amplification_stop',
                    ]
        data3_nestkeys = ['sample_name'];
        data3_keymap = {'xdata':'genome_index',
                        'ydata':'reads',
                        'serieslabel':'genome_strand',
                        'featureslabel':'reads'};
        # make the data object
        dataobject_O = [{"data":data1_O,"datakeys":data1_keys,"datanestkeys":data1_nestkeys},
                        {"data":data2_O,"datakeys":data2_keys,"datanestkeys":data2_nestkeys},
                        {"data":data3_O,"datakeys":data3_keys,"datanestkeys":data3_nestkeys}
                        ];
        # make the tile parameter objects
        # linked set #1
        formtileparameters_O = {'tileheader':'Filter menu','tiletype':'html','tileid':"filtermenu1",'rowid':"row1",'colid':"col1",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-4"};
        formparameters_O = {'htmlid':'filtermenuform1',"htmltype":'form_01',"formsubmitbuttonidtext":{'id':'submit1','text':'submit'},"formresetbuttonidtext":{'id':'reset1','text':'reset'},"formupdatebuttonidtext":{'id':'update1','text':'update'}};
        formtileparameters_O.update(formparameters_O);
        svgparameters_O = {"svgtype":'scatterplot2d_01',"svgkeymap":[data1_keymap,data1_keymap],
                            'svgid':'svg1',
                            "svgmargin":{ 'top': 50, 'right': 150, 'bottom': 50, 'left': 50 },
                            "svgwidth":500,"svgheight":350,
                            "svgx1axislabel":"index","svgy1axislabel":"reads",
    						'svgformtileid':'filtermenu1','svgresetbuttonid':'reset1','svgsubmitbuttonid':'submit1',
                            "svgx1axistickformat":".2e",
                            "svgx1axisticktextattr":{"transform":"matrix(0,1,-1,0,16,6)",
                                                     #"transform":'rotate(90)',"transform":'translate(0,10)'
                                                     },
                            "svgx1axisticktextstyle":{"text-anchor":"start"}
                            };
        svgtileparameters_O = {'tileheader':'Amplifications','tiletype':'svg','tileid':"tile2",'rowid':"row1",'colid':"col2",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-8"};
        svgtileparameters_O.update(svgparameters_O);
        # linked set #2
        formtileparameters2_O = {'tileheader':'Filter menu','tiletype':'html','tileid':"filtermenu2",'rowid':"row2",'colid':"col1",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-4"};
        formparameters2_O = {'htmlid':'filtermenuform2',"htmltype":'form_01',"formsubmitbuttonidtext":{'id':'submit2','text':'submit'},"formresetbuttonidtext":{'id':'reset2','text':'reset'},"formupdatebuttonidtext":{'id':'update2','text':'update'}};
        formtileparameters2_O.update(formparameters2_O);
        tableparameters_O = {"tabletype":'responsivetable_01',
                    'tableid':'table1',
                    "tablefilters":None,
                    "tableclass":"table  table-condensed table-hover",
    			    'tableformtileid':'filtermenu2','tableresetbuttonid':'reset2','tablesubmitbuttonid':'submit2'};
        tabletileparameters_O = {'tileheader':'Amplification statistics','tiletype':'table','tileid':"tile3",'rowid':"row2",'colid':"col2",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-12"};
        tabletileparameters_O.update(tableparameters_O);
        # linked set #3
        formtileparameters3_O = {'tileheader':'Filter menu','tiletype':'html','tileid':"filtermenu3",'rowid':"row3",'colid':"col1",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-4"};
        formparameters3_O = {'htmlid':'filtermenuform3',"htmltype":'form_01',"formsubmitbuttonidtext":{'id':'submit3','text':'submit'},"formresetbuttonidtext":{'id':'reset3','text':'reset'},"formupdatebuttonidtext":{'id':'update3','text':'update'}};
        formtileparameters3_O.update(formparameters3_O);
        tableparameters2_O = {"tabletype":'responsivetable_01',
                    'tableid':'table2',
                    "tablefilters":None,
                    "tableclass":"table  table-condensed table-hover",
    			    'tableformtileid':'filtermenu3','tableresetbuttonid':'reset3','tablesubmitbuttonid':'submit3'};
        tabletileparameters2_O = {'tileheader':'Amplification annotations','tiletype':'table','tileid':"tile4",'rowid':"row3",'colid':"col2",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-12"};
        tabletileparameters2_O.update(tableparameters2_O);
        parametersobject_O = [formtileparameters_O,svgtileparameters_O,formtileparameters2_O,tabletileparameters_O,formtileparameters3_O,tabletileparameters2_O];
        tile2datamap_O = {"filtermenu1":[0],"tile2":[0,0],"tile3":[1],"tile4":[2],"filtermenu2":[1],"filtermenu3":[2]};
        filtermenuobject_O = [{"filtermenuid":"filtermenu1","filtermenuhtmlid":"filtermenuform1",
                "filtermenusubmitbuttonid":"submit1","filtermenuresetbuttonid":"reset1",
                "filtermenuupdatebuttonid":"update1"},{"filtermenuid":"filtermenu2","filtermenuhtmlid":"filtermenuform2",
                "filtermenusubmitbuttonid":"submit2","filtermenuresetbuttonid":"reset2",
                "filtermenuupdatebuttonid":"update2"},{"filtermenuid":"filtermenu3","filtermenuhtmlid":"filtermenuform3",
                "filtermenusubmitbuttonid":"submit3","filtermenuresetbuttonid":"reset3",
                "filtermenuupdatebuttonid":"update3"}];
        # dump the data to a json file
        data_str = 'var ' + 'data' + ' = ' + json.dumps(dataobject_O) + ';';
        parameters_str = 'var ' + 'parameters' + ' = ' + json.dumps(parametersobject_O) + ';';
        tile2datamap_str = 'var ' + 'tile2datamap' + ' = ' + json.dumps(tile2datamap_O) + ';';
        filtermenu_str = 'var ' + 'filtermenu' + ' = ' + json.dumps(filtermenuobject_O) + ';';
        if data_dir_I=='tmp':
            filename_str = 'ddt_data.js'
        elif data_dir_I=='data_json':
            data_json_O = data_str + '\n' + parameters_str + '\n' + tile2datamap_str + '\n' + filtermenu_str;
            return data_json_O;
        with open(filename_str,'w') as file:
            file.write(data_str);
            file.write(parameters_str);
            file.write(tile2datamap_str);
            file.write(filtermenu_str);

    def export_coverage_js(self,data_dir_I="tmp"):
        """exportcoverage data to js file"""

        #get the data for the analysis
        data1_O = [];
        data2_O = [];
        data1_O = self._make_sampleNameStrand(self.coverage);
        data2_O = self._make_sampleNameStrand(self.coverageStats);
        # dump chart parameters to a js files
        data1_keys = ['experiment_id',
                    'sample_name',
                    'genome_chromosome',
                    'genome_strand',
                    'sample_name_strand'
                    ]
        data1_nestkeys = [
                        #'sample_name',
                        'genome_strand'
                        ];
        data1_keymap = {'xdata':'genome_index',
                        'ydata':'reads',
                        'serieslabel':'sample_name_strand',#custom for vis
                        'featureslabel':'reads'};
        data2_keys = ['experiment_id',
                'sample_name',
                'genome_chromosome',
                'genome_strand',
                #'strand_start',
                #'strand_stop',
                #'reads_min',
                #'reads_max',
                #'reads_lb',
                #'reads_ub',
                #'reads_iq1',
                #'reads_iq3',
                #'reads_median',
                #'reads_mean',
                #'reads_var',
                #'reads_n',
                'amplification_start',
                'amplification_stop',
                'used_',
                'comment_'
                    ]
        data2_nestkeys = ['sample_name'];
        data2_keymap = {'xdata':'genome_index',
                        'ydata':'reads',
                        'serieslabel':'genome_strand',
                        'featureslabel':'reads'};
        # make the data object
        dataobject_O = [{"data":data1_O,"datakeys":data1_keys,"datanestkeys":data1_nestkeys},{"data":data2_O,"datakeys":data2_keys,"datanestkeys":data2_nestkeys}];
        # make the tile parameter objects
        formtileparameters_O = {'tileheader':'Filter menu','tiletype':'html','tileid':"filtermenu1",'rowid':"row1",'colid':"col1",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-4"};
        formparameters_O = {'htmlid':'filtermenuform1',"htmltype":'form_01',"formsubmitbuttonidtext":{'id':'submit1','text':'submit'},"formresetbuttonidtext":{'id':'reset1','text':'reset'},"formupdatebuttonidtext":{'id':'update1','text':'update'}};
        formtileparameters_O.update(formparameters_O);
        svgparameters_O = {"svgtype":'scatterplot2d_01',"svgkeymap":[data1_keymap,data1_keymap],
                            'svgid':'svg1',
                            "svgmargin":{ 'top': 50, 'right': 150, 'bottom': 50, 'left': 50 },
                            "svgwidth":500,"svgheight":350,
                            "svgx1axislabel":"index","svgy1axislabel":"reads",
    						'svgformtileid':'filtermenu1','svgresetbuttonid':'reset1','svgsubmitbuttonid':'submit1',
                            "svgx1axistickformat":".2e",
                            "svgx1axisticktextattr":{"transform":"matrix(0,1,-1,0,16,6)",
                                                     #"transform":'rotate(90)',"transform":'translate(0,10)'
                                                     },
                            "svgx1axisticktextstyle":{"text-anchor":"start"}
                            };
        svgtileparameters_O = {'tileheader':'Resequencing coverage','tiletype':'svg','tileid':"tile2",'rowid':"row1",'colid':"col2",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-8"};
        svgtileparameters_O.update(svgparameters_O);
        tableparameters_O = {"tabletype":'responsivetable_01',
                    'tableid':'table1',
                    "tablefilters":None,
                    "tableclass":"table  table-condensed table-hover",
    			    'tableformtileid':'filtermenu1','tableresetbuttonid':'reset1','tablesubmitbuttonid':'submit1'};
        tabletileparameters_O = {'tileheader':'Resequencing coverage statistics','tiletype':'table','tileid':"tile3",'rowid':"row2",'colid':"col1",
            'tileclass':"panel panel-default",'rowclass':"row",'colclass':"col-sm-12"};
        tabletileparameters_O.update(tableparameters_O);
        parametersobject_O = [formtileparameters_O,svgtileparameters_O,tabletileparameters_O];
        tile2datamap_O = {"filtermenu1":[0],"tile2":[0,0],"tile3":[1]};
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

    def _make_sampleNameStrand(self,coverage_I):
        """generate a unique sample/strand name for visualization"""
        for d in coverage_I:
            d['sample_name_strand']="_".join([d['sample_name'],d['genome_strand']]);
        return coverage_I;


