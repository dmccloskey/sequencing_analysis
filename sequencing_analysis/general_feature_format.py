from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData
from sequencing_utilities.makegff import write_samfile_to_gff
import pandas, numpy
from matplotlib.pyplot import subplot, fill_between, xlabel, xlim, \
    ylim, setp, savefig, figure

class general_feature_format():
    #.gff
    # NOTE: multiple chromosomes not yet supported
    def __init__(self,gff_file_I,plus_I,minus_I,plus_high_regions_I, minus_high_regions_I):
        
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

    def convert_bam2gff(self,bam_file):
        '''convert a bam file to a gff file'''

        filename = filename.replace('.bam','.gff');
        extract_strandsFromGff(bam_file,filename,separate_strand=False);
        print(bam_file + " converted to " + filename);
        # record the .gff filename
        self.gff_file = filename;

    def set_gffFile(self,gff_file):
        '''set the name of the gff file'''
        self.gff_file = gff_file;

    def extract_strandsFromGff(self, left, right, scale=True, downsample=0):
        """convert a gff file to a table of positions and reads

        Input:
        gff_file: gff file to read
        left: left position to start analysis
        right: right position to end analysis
        scale: reads will be normalized to have 100 max
        downsample: the number of positions to downsample to

        Output:
        plus: table [index,reads] for the plus strand
        minus: table [index,reads] for the minus strand
        """
        gff_file = self.gff_file;
        # sometimes the first line is a comment which pandas can't handle
        skiprows = 0
        with open(gff_file, "r") as infile:
            if infile.read(1) == "#":
                skiprows = 1
        table = pandas.read_table(gff_file, header=None,
            usecols=[0, 2, 3, 4, 5, 6], comment="#", skiprows=skiprows,
            names=["chromosome", "name", "leftpos", "rightpos", "reads", "strand"])
        table = table[(table.rightpos >= left) & (table.leftpos <= right)]
        # TODO - detect if chromsome_plus and chromosome_minus
        if len(table.chromosome.unique()) > 1:
            raise Exception("multiple chromosomes not supported")
        if (table.leftpos == table.rightpos).all():  # each line is one point
            table = table[["leftpos", "reads", "strand"]]
        table_plus = table[table.strand == "+"].set_index("leftpos")
        table_minus = table[table.strand == "-"].set_index("leftpos")
        # fill missing values with 0
        filler = pandas.Series([range(left, right + 1)], [range(left, right + 1)])
        table_plus["filler"] = 1
        table_minus["filler"] = 1
        table_plus.fillna(0)
        table_minus.fillna(0)
        # extract only the series we need
        plus = table_plus.reads
        minus = table_minus.reads.abs()  # in case stored negative
        if scale:
            plus *= 100. / plus.max()
            minus *= 100. / minus.max()
        # downsample
        collapse_factor = None;
        if downsample > 1:
            collapse_factor = int((right - left) / downsample)
        if collapse_factor and collapse_factor > 1:
            plus = plus.groupby(lambda x: x // collapse_factor).mean()
            plus.index *= collapse_factor
            minus = minus.groupby(lambda x: x // collapse_factor).mean()
            minus.index *= collapse_factor
        #return plus,minus;
        self.plus = plus;
        self.minus = minus;

    def find_highCoverageRegions(self, coverage_min=1.5,coverage_max=5.0,
                points_min=200,consecutive_tol=10):
        '''Find regions of high coverage

        Input:
        plus: plus strand
        minus: minus strand
        coverage_min: minimum coverage of a high coverage region
        coverage_max: maximum coverage of a high coverage region
        points_min: minimum number of points of a high coverage region
        consecutive_tol: maximum number of consecutive points that do not meet the coverage_min/max criteria that can be included a high coverage region

        Output:
        '''

        plus = self.plus;
        minus = self.minus;

        # find indices that are within the min/max coverage
        plus_high = plus[(plus>=coverage_min*float(plus.mean())) & (plus<=coverage_max*float(plus.mean()))];
        minus_high = minus[(minus>=coverage_min*float(minus.mean())) & (minus<=coverage_max*float(minus.mean()))];
        # find indices that meet the minimum number of consecutive points and consecutive point count tolerance limit
        plus_high_regions = plus_high*0;
        minus_high_regions = minus_high*0;
        plus_high_region_summary = [];
        minus_high_region_summary = [];
        cnt = 0;
        for index,value in plus_high.iteritems():
            # initialize all tmp variables and iterators
            if cnt==0: 
                high_regions_tmp = plus_high*0;
                previous = index;
                consecutive_points = 0;
                start=index;
                cnt += 1;
                continue;
            # check for regions of high coverage that
            if index-previous<=consecutive_tol: # 1 satisfy the consecutive tol
                # case for all points but the last point
                # update current high region and iterate counters
                if cnt<len(plus_high)-1:
                    high_regions_tmp[index]=value;
                    consecutive_points += 1;
                    previous = index;
                # case for the last point
                # update all high regions (points_min satisfied)
                elif consecutive_points >= points_min: # 2 satisfy the consecutive points
                    plus_high_regions=record_highCoverageRegions(high_regions_tmp,plus_high_regions)
                    plus_high_region_summary.append({'start':start,'stop':index,'mean':high_regions_tmp[high_regions_tmp>0].mean()});    
            else:
                # update all high regions and reset (points_min satisfied)
                if consecutive_points >= points_min:
                    plus_high_regions=record_highCoverageRegions(high_regions_tmp,plus_high_regions)
                    plus_high_region_summary.append({'start':start,'stop':index,'mean':high_regions_tmp[high_regions_tmp>0].mean()});
                    previous = index;
                    consecutive_points = 0;
                    start=index;
                    high_regions_tmp = plus_high*0;
                    # reset variables and counters
                else:
                    high_regions_tmp = plus_high*0;
                    previous = index;
                    consecutive_points = 0;
                    start=index;
            # update the counter
            cnt += 1;
        cnt = 0;
        for index,value in minus_high.iteritems():
            # initialize all tmp variables and iterators
            if cnt==0: 
                high_regions_tmp = minus_high*0;
                previous = index;
                consecutive_points = 0;
                start=index;
                cnt += 1;
                continue;
            # check for regions of high coverage that
            if index-previous<=consecutive_tol: # 1 satisfy the consecutive tol
                # case for all points but the last point
                # update current high region and iterate counters
                if cnt<len(minus_high)-1:
                    high_regions_tmp[index]=value;
                    consecutive_points += 1;
                    previous = index;
                # case for the last point
                # update all high regions (points_min satisfied)
                elif consecutive_points >= points_min: # 2 satisfy the consecutive points
                    minus_high_regions=self.record_highCoverageRegions(high_regions_tmp,minus_high_regions)
                    minus_high_region_summary.append({'start':start,'stop':index,'mean':high_regions_tmp[high_regions_tmp>0].mean()});    
            else:
                # update all high regions and reset (points_min satisfied)
                if consecutive_points >= points_min:
                    minus_high_regions=self.record_highCoverageRegions(high_regions_tmp,minus_high_regions)
                    minus_high_region_summary.append({'start':start,'stop':index,'mean':high_regions_tmp[high_regions_tmp>0].mean()});
                    previous = index;
                    consecutive_points = 0;
                    start=index;
                    high_regions_tmp = minus_high*0;
                    # reset variables and counters
                else:
                    high_regions_tmp = minus_high*0;
                    previous = index;
                    consecutive_points = 0;
                    start=index;
            # update the counter
            cnt += 1;
        # filter out all 0 values
        plus_high_regions = plus_high_regions[plus_high_regions>0];
        minus_high_regions = minus_high_regions[minus_high_regions>0];
        # record the data
        self.plus_high_regions = plus_high_regions;
        self.minus_high_regions = minus_high_regions;

        #return plus_high_region_summary,minus_high_region_summary,plus_high_regions, minus_high_regions
        return plus_high_region_summary,minus_high_region_summary

    def record_highCoverageRegions(self,high_region,high_regions):
        '''Record high coverage region index and reads to master list of all high regions indices and reads'''
        high_regions += high_region;
        return high_regions;
    
    def plot_strands(plus,minus,left,right,scale=True, output=None,
                  coverage_max=5.0, downsample=2000):
        '''Plot strand coverage using matplotlib

        Input:
        plus: table [index,reads] for the plus strand
        minus: table [index,reads] for the minus strand
        left: left position to start analysis
        right: right position to end analysis
        scale: reads will be normalized to have 100 max
        output: output file
        coverage_max: maximum read coverage to plot
        downsample: the number of positions to downsample to

        Output:
        plus: table [index,reads] for the plus strand
        minus: table [index,reads] for the minus strand'''

        plus = self.plus;
        minus = self.minus;
        ## fill missing values with 0
        #filler = pandas.Series([range(left, right + 1)])
        #plus = filler + plus;
        #minus = filler + minus;
        # scale
        if scale:
            plus *= 100. / plus.max()
            minus *= 100. / minus.max()
        # downsample
        collapse_factor = None;
        if downsample > 1:
            collapse_factor = int((right - left) / downsample)
        if collapse_factor and collapse_factor > 1:
            plus = plus.groupby(lambda x: x // collapse_factor).mean()
            plus.index *= collapse_factor
            minus = minus.groupby(lambda x: x // collapse_factor).mean()
            minus.index *= collapse_factor;
        # plot the figure:
        fig=figure()
        ax = subplot(1,1,1)
        fill_between(minus.index.values, 0, minus.values, color="orange", alpha=0.75)
        fill_between(plus.index.values, 0, plus.values, color="blue", alpha=0.5)
        xlim(left, right)
        if scale:
            ylim(0, 100)
        else:
            ylim(ymin=0)
        # hide extra text and labels
        ax.grid(axis="x")
        ax.ticklabel_format(axis="x", style="sci", scilimits=(1, 3))
        if output:
            savefig(output, transparent=True)

    def clear_data(self):
        self.gff_file = None;
        self.minus = None;
        self.plus = None;
        self.plus_high_regions = None;
        self.minus_high_regions = None;